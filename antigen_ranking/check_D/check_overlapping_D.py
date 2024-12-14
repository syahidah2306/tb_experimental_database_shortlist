#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd

# Step 1: Load the data
fasta_file_path = './mycobacterium_tuberculosis_H37Rv_proteins_mycodb_D.fasta'
#epitope_data = pd.read_excel('./Data/epitope_score_20240819_7criteria,inc_uniq_mtb, bce,-locus.xlsx')
epitope_data = pd.read_excel( "./epitope_score_20240819_7criteria,inc_uniq_mtb, bce,-locus.xlsx",engine='openpyxl')


#get_ipython().system('head $fasta_file_path')
epitope_data


# In[ ]:


from Bio import SeqIO, SearchIO
from Bio.Blast.Applications import NcbiblastpCommandline
import tempfile
import os

# Function to perform a local BLAST alignment using blastp
def blast_sequence(tce_sequence, protein_sequence):
    # Create temporary files for the query and subject sequences
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as query_file:
        query_file.write(f">Query\n{tce_sequence}")
        query_file.close()
    
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as subject_file:
        subject_file.write(f">Subject\n{protein_sequence}")
        subject_file.close()
    
    # Create a temporary file for BLAST output
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as output_file:
        output_file.close()

    # Define the BLAST command with parameters adjusted for short sequences
    blastp_cline = NcbiblastpCommandline(
        num_threads=4, 
        query=query_file.name, 
        subject=subject_file.name, 
        outfmt=5, # XML format for easier parsing
        out=output_file.name,
        word_size=2, # Adjusted for short sequences
        gapopen=11,
        gapextend=1,
        evalue=0.01, # A bit strict e-value to reduce false positive
        matrix="BLOSUM62",
        max_hsps=1, # One alignment per pair
        max_target_seqs=1, # One target sequence (subject)
        soft_masking="false",
        ungapped=False
    )

    # Run the BLAST command
    stdout, stderr = blastp_cline()

    # Parse the BLAST output to extract alignments
    alignments = []
    blast_qresult = SearchIO.read(output_file.name, "blast-xml")
    
    for hit in blast_qresult:
        for hsp in hit.hsps:
            # Calculate mismatches manually by comparing aligned sequences
            mismatches = sum(1 for a, b in zip(hsp.query.seq, hsp.hit.seq) if a != b)
            
            # Check if the full length of the query is matched and allow up to 3 mismatches
            if hsp.aln_span >= len(tce_sequence) and mismatches <= 3:
                alignments.append((hsp.hit_start, hsp.hit_end))
    
    # Clean up temporary files
    os.remove(query_file.name)
    os.remove(subject_file.name)
    os.remove(output_file.name)
    
    return alignments

# Function to print the result for each protein sequence
def print_result(header, lowercase_sequence, matched_tces, total_score, total_matched_len, total_lowercase_len, total_epitope_score):
    print(f"Header: {header}")
    print(f"Lowercased Sequence: {lowercase_sequence}")
    print(f"Matched TCEs: {matched_tces}")
    print(f"Total Score: {total_score}")
    print(f"Total Matched Length: {total_matched_len}")
    print(f"Total Lowercase Length: {total_lowercase_len}")
    print(f"Total Epitope Score: {total_epitope_score}")
    print("\n" + "-"*80 + "\n")


# In[ ]:


from Bio import SeqIO
from Bio.Seq import Seq

# Store results for tabulation later
results = []

# Step 2: Process each protein sequence
for record in SeqIO.parse(fasta_file_path, "fasta"):
    protein_sequence = str(record.seq)
    header = record.description
    
    total_score = 0
    total_matched_len = 0
    lowercase_sequence = protein_sequence
    
    matched_tces = []
    matched_positions = []
    
    for _, row in epitope_data.iterrows():
        tce_sequence = row['TCE_sequence']
        score = row['score']
        
        alignments = blast_sequence(tce_sequence, protein_sequence)
        
        if alignments:
            matched_tces.append(tce_sequence)
            total_score += score
            total_matched_len += len(tce_sequence)
            
            for start, end in alignments:
                matched_positions.append((start, end))
    
    # Convert matched positions to lowercase in the protein sequence
    lowercase_sequence = list(lowercase_sequence)
    for start, end in matched_positions:
        lowercase_sequence[start:end] = [x.lower() for x in lowercase_sequence[start:end]]
    
    lowercase_sequence = "".join(lowercase_sequence)
    
    # Calculate the total lowercase length in the protein sequence
    total_lowercase_len = sum(1 for c in lowercase_sequence if c.islower())
    
    # Calculate the total epitope score
    if total_matched_len > 0:
        total_epitope_score = (total_score * total_lowercase_len) / total_matched_len
    else:
        total_epitope_score = 0
    
    # Append the result to the list
    results.append({
        'Header': header,
        'Lowercased_Sequence': lowercase_sequence,
        'Matched_TCEs': matched_tces,
        'Total_Score': total_score,
        'Total_Matched_Length': total_matched_len,
        'Total_Lowercase_Length': total_lowercase_len,
        'Total_Epitope_Score': total_epitope_score
    })
    
    # Step 7: Print the result for the current protein
    print_result(header, lowercase_sequence, matched_tces, total_score, total_matched_len, total_lowercase_len, total_epitope_score)


# In[ ]:


# Step 8: Tabulate and sort the results by total_epitope_score
results_df = pd.DataFrame(results)
sorted_results_df = results_df.sort_values(by='Total_Epitope_Score', ascending=False)

# Final Print: Display the tabulated results
print("\nFinal Tabulated Results (Sorted by Total Epitope Score):\n")
sorted_results_df[['Header', 'Total_Epitope_Score', 'Total_Score', 'Total_Matched_Length', 'Total_Lowercase_Length']]

# Export the DataFrame to an Excel file
sorted_results_df.to_excel('./output_antigen.xlsx', index=False,engine='openpyxl')


# In[ ]:




