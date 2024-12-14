''''
Script 4: This script to enable the prediction of the MHCI supertypes and allele using IEDB consesus algorithm but with initial
algorithm to regulate the prediction length. This script is the continuity of script 3 which conclude the promiscous criteria for the T and B epitopes.
ADDITION: This script is full length script for the prediction after combined with other previous scripts.
'''
import os
import re
import funcprog as fx
from collections import OrderedDict
import pickle
import shutil  # module to handle files that os can't
import time
import statistics
from itertools import chain

report_stat = 1 # 0 = print Error only; 1 = print with minimal info; 2 = print up to variable content; 3 = print with iteration
cmd_stat = 1# 0 = cmd do not run, 1 = cmd run

mhc_all_swc = True # false for batches, but if want to run all tce, this is true
mhc_1x_swc = True
mhc_2x_swc = True
mhc_mutual_swc = True
mhc1_debug_swc = True #default = False, only to debug set to True (which will override other switches)
mhc2_debug_swc = True  #default = False, only to debug set to True (which will override other switches)
debug_stat = 0 # 0 = debug do not run, 1 = debug run
rm_stat = 0 # to remove the database created # 0 = remove do not run, 1 = remove run

HEG_swc1 = True
HEG_swc2 = True  #tce blast with mtb seq from NCBI db
popcov_swc = True
tce_pickle_swc = True
IR_swc = True  #contains def "tce_all_dos" , need to be always open even for debug
hyp_swc = True
sau_swc = True
grd_smry_tble_swc = True

# input file for script 1
mtb_myco_table = "Mycobacterium_tuberculosis_H37Rv_txt_v4_myco_table.txt" #prev file = ID_Mycobacterium_tuberculosis_H37Rv_txt_v3.txt"
mtb_ncbi_table = "mtb_h37rv_protein_table.txt" #prev file = ID_GCF_000195955.2_ASM19595v2_feature_table20190527.txt" # this file also used in script 3
mice = "combined_mtb_expressed_in_mice.txt" #prev file name= ID_mice_temporal_expreesion_profile_sorted.txt
human = "combined_mtb_expressed_in_human.txt" #prev file name= combined_genes_notredundant.txt. This current file was added with heg2_genes.txt
latency = "combined_mtb_expressed_in_lab.txt" #prevlatency_genes.txt"
mtb_ncbi = "GCF_000195955.2_ASM19595v2_protein.faa.out" #this file contains 3906 proteins (mtb ori NCBI)
mtb_mycodb = "Mycobacterium_tuberculosis_H37Rv_proteins_mycodb.out2" #"Mycobacterium_tuberculosis_H37Rv_proteins_v3.fasta"  # mtb database originated from mycobrowser, 4098 proteins.

# input file for script 3a
mtb_heg_rsid = "mtb_heg_rsid.faa" # contained 1180 HEG protein
heg_blastp_tce = "blastp_tce_heg_op.txt" #blast output
heg_blastp_bce = "blastp_bce_heg_op.txt" #blast output
tce_all_ep_file = "tce_h37rv_pos_allmhc.faa" #TCE downloaded from IEDB (in FAA file format), this file also used in script 4  # prev_file = export_table__TCell_out.faa
database_mtb_np = 'db/' + mtb_heg_rsid # make a new directory for the database so that easy to delete after used # contained 900+ HEG protein

#input file for script 3b - the process was done in SAU section
mtb_heg_rvid = "mtb_heg_rvid.faa" # contained 2500+ HEG protein
database_mtb_rv = 'db/' + mtb_heg_rvid
#heg_rvid_blastp_tce = "blastp_tce_heg_rvid_op.txt" #blast output
# input file for script
#ort_ip_file = "orthologous_15.faa" #prev_file = matched_h37rv_cds.faa.out # o/p of 173 mtb strains as input for blastp with tce
#hyp_tce_result = "blastp_ort_tce_op.txt"  # blastp file produced internally

tce_mhc1_ep = "" # prev_file =  epitope_table_Tcell_mhc1.faa # this file changed according subset
tce_mhc2_ep = ""  # prev_file =  epitope_table_Tcell_mhc2.faa # this file changed according subset
if mhc_all_swc:
    tce_mhc1_ep = "tce_h37rv_pos_allmhc.faa" # all mtb mhc1 epitopes = "epitope_mhc1_20201117.faa"
    tce_mhc2_ep = "tce_h37rv_pos_allmhc.faa" # all mtb mhc2 epitopes = "epitope_mhc2_20201117.faa"
if mhc_1x_swc:
    tce_mhc1_ep = "tce_h37rv_pos_mhc1x.faa" # all mtb mhc1x epitopes = epitope_mhc1x_20201117_H-I.faa
    tce_mhc2_ep = ""
if mhc_2x_swc:
    tce_mhc1_ep = ""
    tce_mhc2_ep = "tce_h37rv_pos_mhc2x.faa" # all mtb mhc2x epitopes = epitope_mhc2x_20201117_D.faa
if mhc_mutual_swc:
    tce_mhc1_ep = "tce_h37rv_pos_mhcnonspec_mutual.faa"  # all mtb mutual epitopes = epitope_mhc1-2_mutual_20201117_A-G.faa
    tce_mhc2_ep = "tce_h37rv_pos_mhcnonspec_mutual.fa"

bce_ep = "epitope_bce_20201117.faa" #othervise, bce_pos_mtb_h37Rv =  bce_h37rv_pos.faa # # prev_file = export_table__BCell_out.faa
# mhc1_netmhc_alleles = "consensus_mhc1.txt"
# mhc1_netmhc_alleles = "mhc_i_netmhccons_alleles_full.txt" #mhc1_recommended = netmhcpan,
mhc1_iedbrec_alleles = "mhc1_alleles_more1percent_sorted.txt" # alleles consensus+ IEDB suggested = "human_alleles_total_mhc1_IEDB_rec_sorted.txt"
# MHC2 alleles used was in the text, but can be review in file named = mhc2_alleles_in_used.txt
# mhc2_cnss = "MHCII_IEDB_recommeded_alleles_2_5.txt"
mhc2_iedbrec_alleles = "human_alleles_total_mhc2_IEDB_rec_sorted.txt"
# mhc2_IEDBrec_netmhcpan = "MHCII_IEDB_recommeded_alleles.txt"
blastp_bce_tce = "blastp_bce_op.txt"  # this blastp file used in IR section
# processing required variable
hs = "human_protein_GRCh38.faa"
biom = "human_microbiome_project_sequence_1.faa"
database_mtb_ncbi = 'db/' + mtb_ncbi
database_tce = 'db/' + tce_all_ep_file
database_mtb_ncbi = 'db/' + mtb_ncbi
database_mtb_myco = 'db/' + mtb_mycodb
database_hs = 'db/' + hs
database_biom = 'db/' + biom

### xxxxxxxxxxxxxxxxxxxxxxxxxxxxx SCRIPT 1: CONVERSION OF HIGHLY EXPRESSED GENES xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx###
if HEG_swc1:
    if report_stat >= 0: print("\n\n############## SCRIPT 1: CONVERSION OF HIGHLY EXPRESSED GENES ################\n")

    # read input file for mtb_mycobrowser
    fx.file_val(mtb_myco_table)
    mtb_ctr = 0
    mtb_myco_table = open(mtb_myco_table).readlines()
    for i in mtb_myco_table:
        if i.startswith("Rv"):
            mtb_ctr += 1
    if report_stat == 1: print("mtb_myco_table = ", mtb_myco_table)
    if report_stat == 1: print("mtb_ctr = ", mtb_ctr)

    #to split the string in the file and put in dictionary
    split_cntr = 0
    split_mtb_myco = {}
    for i in range(len(mtb_myco_table)):
        temp = mtb_myco_table[i].split()
        split_mtb_myco[temp[0].strip()] = temp[1] #to add into dictionary
        split_cntr += 1
    if report_stat == 1: print("split_mtb_myco", split_mtb_myco, "\n")
    if report_stat == 1: print("split_cntr = ", split_cntr)

    #read input file Mtb_NCBI
    fx.file_val(mtb_ncbi_table)
    NCBI_ctr = 0
    mtb_NCBI = open(mtb_ncbi_table).readlines()
    for i in mtb_NCBI:
        if i.startswith("Rv"):
            mtb_ctr += 1
    if report_stat == 1: print("mtb_NCBI = ", mtb_NCBI)
    if report_stat == 1: print("mtb_ctr = ", NCBI_ctr)

    #to split the string in the file and put in dictionary
    sNCBI_cntr = 0
    split_mtb_NCBI = {}
    for i in range(len(mtb_NCBI)):
        temp = mtb_NCBI[i].split()
        split_mtb_NCBI[temp[0]] = temp[1].strip()
        sNCBI_cntr += 1
        #print("split_mtb_NCBI[", temp[0], "] = ", temp[1])
    if report_stat == 1: print("split_mtb_NCBI", split_mtb_NCBI, "\n")
    if report_stat == 1: print("split_NCBI_cntr = ", sNCBI_cntr)

    ##TODO: to check the code for Myco db dictionary
    fx.file_val(mtb_mycodb)
    mtb_seq_mycodb_list = open(mtb_mycodb).readlines()
    mtb_seq_mycodb_dict = {}
    for rv_line in range(len(mtb_seq_mycodb_list)):
        if re.match(r"^>(Rv\w+)", mtb_seq_mycodb_list[rv_line]):
            rv_id = re.findall(r"^>(Rv\w+)", mtb_seq_mycodb_list[rv_line])[0]
        if not re.match(r"^>(Rv\w+)", mtb_seq_mycodb_list[rv_line]):
            rv_seq = mtb_seq_mycodb_list[rv_line].rstrip()
            mtb_seq_mycodb_dict[rv_id] = rv_seq
    if report_stat == 1:
        print("mtb_seq_mycodb_dict:", len(mtb_seq_mycodb_dict))
        print("mtb_seq_mycodb_dict:", mtb_seq_mycodb_dict.keys())

    #read input file mice
    fx.file_val(mice)
    mice = open(mice).readlines()
    mice_list = []
    for column in mice:
        col_mice = column.split()[0] #to cover \t, \r and \n
        mice_str = col_mice.rstrip()  #rstrip is to remove \n trail
        mice_list.append(mice_str)
    mtb_gene_mice_rvid = list(set(mice_list))
    if report_stat == 1: print("mice list unique  = ", len(mtb_gene_mice_rvid), mtb_gene_mice_rvid)

    mtb_gene_mice_rsid = []
    missing_mice = []
    missing_mice_myco = []
    with open("mice_genes_NCBI.txt", 'w') as mice_genes_NCBI:
        for i in mtb_gene_mice_rvid:
            try:
                translate_m = split_mtb_NCBI[i]
                mtb_gene_mice_rsid.append(translate_m)
                mice_genes_NCBI.write('%s\n' % translate_m)
            except KeyError:
                missing_mice.append(i)

            try:
                translate_m = mtb_seq_mycodb_dict[i]
                #TODO:append translate_m (conversion from mycobrowser db) into a variable

            except KeyError:
                missing_mice_myco.append(i)
    mice_genes_NCBI.close()
    if report_stat == 1: print("translated mice NCBI :", len(mtb_gene_mice_rsid), mtb_gene_mice_rsid)
    if report_stat == 1: print("missing mice NCBI:", len(missing_mice), missing_mice)
    if report_stat == 1: print("missing mice myco:", len(missing_mice_myco), missing_mice_myco)


    #read input file human (convert from list to set)
    fx.file_val(human)
    human = open(human).readlines()
    human_list = []
    for line in human:
        col0 = line.split()[0] #to cover \t, \r and \n
        human_str = col0.rstrip()
        human_list.append(human_str)
    if report_stat == 1: print("human_list = ", len(human_list), human_list)

    mtb_gene_human_rvid = []  #convert from common name to Rv
    missing_humanRv = []
    missing_human_myco = []
    with open("human_genes_rsid.txt", 'w') as human_genes_rsid:
        for i in human_list:
            try:
                translate_Rv = split_mtb_myco[i]
                mtb_gene_human_rvid.append(translate_Rv)
                human_genes_rsid.write('%s\n' % translate_Rv)
            except KeyError:
                if re.match(r"^Rv.*", i):
                    mtb_gene_human_rvid.append(i) #to add RV genes into the list regardless using Rv name
                else:
                    missing_humanRv.append(i) #missing because it detect common name to rv name not otherwise  #just to debug list, not to remove
        human_genes_rsid.close()
    mtb_gene_human_rvid = list(set(mtb_gene_human_rvid))
    missing_humanRv_uniq = list(set(missing_humanRv))
    if report_stat == 1: print("human_Rv_myco = ", len(mtb_gene_human_rvid), mtb_gene_human_rvid)
    if report_stat == 1: print("human list unique myco :", len(mtb_gene_human_rvid), mtb_gene_human_rvid)
    if report_stat == 1: print("missing humanRv uniq myco:", len(missing_humanRv_uniq), missing_humanRv_uniq)

    # convert from mtb genes in human from Rv to NCBI name
    mtb_gene_human_rsid = []
    missing_humanNCBI = []
    with open("human_genes_NCBI.txt", 'w') as human_genes_NCBI:
        for i in mtb_gene_human_rvid:
            try:
                translate2 = split_mtb_NCBI[i]
                mtb_gene_human_rsid.append(translate2)
                human_genes_NCBI.write('%s\n' % translate2)
            except KeyError:
                missing_humanNCBI.append(i)
    human_genes_NCBI.close()
    if report_stat == 1: print("human_NCBI = ", len(mtb_gene_human_rsid), mtb_gene_human_rsid)
    if report_stat == 1: print("missing key human NCBI :", len(missing_humanNCBI), missing_humanNCBI)

    file = open("NCBI_human","w")
    for i in mtb_gene_human_rsid:
        file.write(i)
        file.write('\n')
    file.close()

    #read input file lab list
    fx.file_val(latency)
    lab_list = []
    lab_f = open(latency).readlines()
    for i in lab_f:
        colL = i.split()[0]  #to cover \t, \r and \n
        lab_str = colL.rstrip()
        lab_list.append(lab_str)
    if report_stat == 1: print("lab_list :", len(lab_list), lab_list)

    lab_Rv_list = []  #to standardise common name into Rv name
    missing_lab_Rv = []
    with open("latency_genes_rsid.txt", 'w') as latency_genes_rsid:
        for i in lab_list:
            try:
                translate = split_mtb_myco[i]
                lab_Rv_list.append(translate)
                latency_genes_rsid.write('%s\n' % translate)
            except KeyError:
                if re.match(r"^Rv.*", i):
                    lab_Rv_list.append(i)  # to add RV genes into the list regardless using Rv name
                else:
                    missing_lab_Rv.append(i)  #just to debug the #of RV list not converted
        mtb_gene_lab_rvid = list(set(lab_Rv_list))  # make list to set (unique0 list)
    latency_genes_rsid.close()

    if report_stat == 1: print("lab_Rv_list_myco:", len(lab_Rv_list), lab_Rv_list)
    if report_stat == 1: print("lab_rv_list_uniq_myco :", len(mtb_gene_lab_rvid), mtb_gene_lab_rvid)
    if report_stat == 1: print("missing_lab_Rv_myco :", len(missing_lab_Rv), missing_lab_Rv)

    mtb_gene_lab_rsid = []
    missing_L = []
    with open("latency_genes_NCBI.txt", 'w') as latency_genes_NCBI:
        for i in mtb_gene_lab_rvid:
            try:
                translate_L = split_mtb_NCBI[i]
                mtb_gene_lab_rsid.append(translate_L)
                latency_genes_NCBI.write('%s\n' % translate_L)
            except KeyError:
                missing_L.append(i)
    latency_genes_NCBI.close()
    if report_stat == 1: print("Translated latency NCBI :", len(mtb_gene_lab_rsid), mtb_gene_lab_rsid)
    if report_stat == 1: print("missing latency NCBI :", len(missing_L), missing_L)

    #making variable to put all Rv_id
    mtb_heg_rvid_set = set()  # this list is unique
    missing_mtb_heg_rvid_set = set()
    for rv_id in mtb_seq_mycodb_dict.keys():
        if rv_id in mtb_gene_mice_rvid or rv_id in mtb_gene_human_rvid or rv_id in mtb_gene_lab_rvid:
            mtb_heg_rvid_set.add(rv_id)
        else:
            missing_mtb_heg_rvid_set.add(rv_id)
    mtb_heg_rvid_list = list(mtb_heg_rvid_set)
    if report_stat == 1: print("mtb_heg_rvid_list:", len(mtb_heg_rvid_list), mtb_heg_rvid_list)
    if report_stat == 1: print("missing_mtb_heg_rvid_set:", len(list(missing_mtb_heg_rvid_set)), list(missing_mtb_heg_rvid_set))


    # making a file of Rvid genes with its sequence
    with open("mtb_heg_rvid.faa",'w') as mtb_heg_rvid_op_file:  #mtb_full_dict(variable)- header(NP id does not have >) but in file have >
        mtb_heg_rvid_op_dict = OrderedDict()
        for rv_id in mtb_heg_rvid_list:
            if rv_id in mtb_seq_mycodb_dict:
                rv_seq = mtb_seq_mycodb_dict[rv_id]
                #if NP_list_str in mtb_full_dict:
                mtb_heg_rvid_op_dict[rv_id] = rv_seq
                mtb_heg_rvid_op_file.write('>%s\n%s\n' % (rv_id, rv_seq))
    mtb_heg_rvid_op_file.close()
    if report_stat == 1: print("mtb_heg_rvid_op:", mtb_heg_rvid_op_dict)




    #make an NP dictionary with its sequence to be filtered with HEG list
    fx.file_val(mtb_ncbi)
    Mf = open(mtb_ncbi, "r")
    mtb_list = Mf.read().splitlines()
    mtb_full_dict = OrderedDict()
    for line in range(len(mtb_list)):
        if re.match(r"^>(\S+) ", mtb_list[line]):
            NP_id_str = re.findall(r"^>(\S+) ", mtb_list[line])[0]
        if not re.match(r"^>(\S+) ", mtb_list[line]):
            NP_seq_str = mtb_list[line]
            mtb_full_dict[NP_id_str] = NP_seq_str
    if report_stat == 0: print("mtb_full_dict:", mtb_full_dict)

    # to match with heg list and shortlisted the NP list (with NP seq)
    mtb_heg_list = set()  # this list is unique
    for NP_id_str in mtb_full_dict.keys():
        if NP_id_str in mtb_gene_mice_rsid:
            mtb_heg_list.add(NP_id_str)
        if NP_id_str in mtb_gene_human_rsid:
            mtb_heg_list.add(NP_id_str)
        if NP_id_str in mtb_gene_lab_rsid:
            mtb_heg_list.add(NP_id_str)
    mtb_heg_list_uniq = list(set(mtb_heg_list))
    if report_stat == 1: print("mtb_heg_list NCBI:", len(mtb_heg_list), mtb_heg_list)

    #to make a file for mtb_heg_dict
    with open("mtb_heg_rsid.faa",'w') as mtb_heg_op_file:  #mtb_full_dict(variable)- header(NP id does not have >) but in file have >
        mtb_heg_dict = OrderedDict()
        for NP_list_str in mtb_heg_list:
            if NP_list_str in mtb_full_dict:
                NP_seq = mtb_full_dict[NP_list_str]
                #if NP_list_str in mtb_full_dict:
                mtb_heg_dict[NP_list_str] = NP_seq
                mtb_heg_op_file.write('>%s\n%s\n' % (NP_list_str, NP_seq))
    mtb_heg_op_file.close()
    if report_stat == 1: print("mtb_heg_dict NCBI:", len(mtb_heg_dict), mtb_heg_dict)


    '''
    cross_check = []
    missing_cross_check = []
    for i in range(len(missing_mice)):
        for j in range(len(missing_human2)):
            for k in range(len(missing_L)):
                if missing_mice[i] == missing_human2[j] and missing_mice[i] == missing_L[k]:
                    cross_check.append(missing_mice[i])
    if report_stat == 1: print("cross_check:", len(cross_check), cross_check)
    if report_stat == 1: print("missing_cross_check:", len(missing_cross_check), missing_cross_check)'''


### xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx SCRIPT 2&3b: HEG FOR OVERLAPPING TCE AND BCE WITH MTB ORIGINATED FROM NCBI xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx###
if HEG_swc2:
    if report_stat >= 0: print(
        "\n\n############## SCRIPT 2&3b: HEG FOR OVERLAPPING TCE AND BCE WITH MTB ORIGINATED FROM NCBI  ####################\n")
    '''SCRIPT 3:This script is the modified script 2, which there is some strategy amendment from script 2.
    This script is used to process the blastp result of
    1) epitope_mhc_all_20201117.faa vs mtb_heg_rsid.faa prev_file: NCBI_filtered_genes_matchseq.faa (filtered from NCBI_filtered_genes.out.txt)
    2) epitope_bce_20201117.faa vs mtb_heg_rsid.faa prev_file: NCBI_filtered_genes_matchseq.faa (filtered from NCBI_filtered_genes.out.txt) 
    3) epitope_mhc_all_20201117.faa vs mtb_heg_rvid_file.faa 
    '''

    startscript3_time = time.time()
    if report_stat >= 0: print("\nIdentify the reference protein (H37Rv) sequences that contain epitope region.\n")

    # make NCBI database from fasta format file
    #cmd_script3_1 = "makeblastdb -in " + mtb_heg_rsid + " -out " + database_mtb_np + " -parse_seqids -dbtype prot"
    cmd_script3_1 = "makeblastdb -in " + mtb_heg_rvid + " -out " + database_mtb_rv + " -parse_seqids -dbtype prot"
    if cmd_stat == 1: os.system(cmd_script3_1)

    # run blast mtb_seq (HEG-mycobrowser) as query, tce_ep as subject
    cmd_script3_2 = "blastp -num_threads 6 -task blastp-short -query " + tce_all_ep_file + " -db " + database_mtb_rv + " -out blastp_tce_heg_op.txt"  #
    if cmd_stat == 1: os.system(cmd_script3_2)
    print(cmd_script3_2)
    # run blast mtb_seq (HEG-mycobrowser) as query, bce_ep as subject
    cmd_script3_3 = "blastp -num_threads 6 -task blastp-short -query " + bce_ep + " -db " + database_mtb_rv + " -out blastp_bce_heg_op.txt"
    if cmd_stat == 1: os.system(cmd_script3_3)
    print(cmd_script3_3)


    # processing T cell vs HEG-NCBI blastp o/p file
    fx.file_val(heg_blastp_tce)
    blastp_tce_f = open(heg_blastp_tce, "r")
    linesT = blastp_tce_f.read().splitlines()
    total_lineT = len(linesT)
    blastp_tce_f.close()
    blastlinenumsT = []  # list of line where the query start
    i = 0
    for i in range(total_lineT):
        if re.match(r"^Query= (\S+)", linesT[i]):
            blastlinenumsT.append(i)
    total_queryT = len(blastlinenumsT)
    if report_stat == 1: print("Total queries = ", total_queryT)
    blastlinenumsT.append(total_lineT)  # append EOL into blastlinenums
    if report_stat == 1: print("blastlinenums T = ", blastlinenumsT, "\n")
    # get the lines where query start
    if report_stat >= 0: print("Processing T cell blast result.")

    # process query by queryT
    qu_T = [''] * total_queryT  # list of result according to query heg_blastp_tce
    for i in range(total_queryT):
        temp = []
        for l in range(blastlinenumsT[i], blastlinenumsT[i + 1]):
            temp.append(linesT[l])
        if report_stat == 2: print("Total lines for query heg_blastp_tce #", i + 1, " are ", len(temp))
        qu_T[i] = temp
        del temp

    '''###only for debug to check for index in list 'qu_T'
    for i in range(len(TCE_blasted_f)):
        if report_stat == 1: print("\n******Query TCE_blasted_f " + str(i+1) + "******\n")
        for j in range(len(qu_T[i])):
            lineT=qu_T[i][j]
            if report_stat == 2: print("qu_T[" , i,"][" , j, "] " , lineT)
    if report_stat == 1: print('\n\n\n')'''




    # processing B cell vs HEG-NCBI blastp o/p file
    fx.file_val(heg_blastp_bce)
    blastp_bce_f = open(heg_blastp_bce, "r")
    linesB = blastp_bce_f.read().splitlines()
    total_lineB = len(linesB)
    blastp_bce_f.close()
    blastlinenumsB = []  # list of line where the query start
    i = 0
    for i in range(total_lineB):
        if re.match(r"^Query= (\S+)", linesB[i]):
            blastlinenumsB.append(i)
    total_queryB = len(blastlinenumsB)
    if report_stat == 1: print("Total queries = ", total_queryB)
    blastlinenumsB.append(total_lineB)  # append EOL into blastlinenums
    if report_stat == 1: print("blastlinenums B = ", blastlinenumsB, "\n")
    # get the lines where query start
    if report_stat >= 0: print("Processing B cell blast result.")

    # process query by queryB
    qu_B = [''] * total_queryB  # list of result according to query Bcell
    for i in range(total_queryB):
        temp = []
        for l in range(blastlinenumsB[i], blastlinenumsB[i + 1]):
            temp.append(linesB[l])
        if report_stat == 1: print("Total lines for query Bcell #", i + 1, " are ", len(temp))
        qu_B[i] = temp
        del temp

    '''
    ###only for debug to check for index in list 'qu_B'
    for i in range(len(Bcell)):
        if report_stat == 1: print("\n******Query Bcell " + str(i+1) + "******\n")
        for j in range(len(qu_B[i])):
            lineB=qu_B[i][j]
            if report_stat == 1: print("qu_B[" , i,"][" , j, "] " , lineB)
    if report_stat == 1: print('\n\n\n')'''

    # Declaring mtbfiles(mtb_smry_dolol) as a dictionary (summary dictionary)(for debuging)
    # Structure: mtb_smry_dolol = {prtn_id:[seq_normal, seq_coded,[TCE_id, qu_seq,sub_seq, st,sp], [BCE_id, qu_seq, sub_seq, st,sp]]}
    mtb_smry_dolol = {}
    fx.file_val(mtb_heg_rvid)
    Mf = open(mtb_heg_rvid, "r")
    mtb_heg = Mf.read().splitlines()
    Mf.close()
    prtn_id = ""
    for i in range(len(mtb_heg)):
        if re.match(r"^>(\w+)", mtb_heg[i]):
            prtn_id = re.findall(r"^>(\w+)", mtb_heg[i])[0]
        if not re.match(r"^>(\w+)", mtb_heg[i]):
            seq_ori = mtb_heg[i].rstrip()
            seq_coded = ""
            tce_id_lol = []
            bce_id_lol = []
            mtb_smry_dolol[prtn_id] = [seq_ori, seq_coded, tce_id_lol, bce_id_lol]


    # Declaring TCE table (tce_smry_dolol) as a dictionary
    # Structure: mtb_smry_dolol = {Np_id:[seq_normal, seq_coded,[TCE_id, qu_seq,sub_seq, st,sp], [BCE_id, qu_seq, sub_seq, st,sp]]}
    # Structure: tce_smry_dolol = {TCE_id:[tce_qu_seq,[[prtn_id, tce_loc]], [[prtn_id, bce_base_ovlp]] }
    # None= bases overlapped/total tce length
    tce_smry_dolol = {}
    fx.file_val(tce_all_ep_file)
    Tf = open(tce_all_ep_file, "r")
    tce_all_ep = Tf.read().splitlines()
    Tf.close()
    TCE_id = ""
    for i in range(len(tce_all_ep)):
        if re.match(r"^>(\S+)", tce_all_ep[i]):
            TCE_id = re.findall(r"^>(\S+)", tce_all_ep[i])[0]
        if not re.match(r"^>(\S+)", tce_all_ep[i]):
            tce_qu_seq = tce_all_ep[i].rstrip()
            tce_loc_lol = []
            bce_base_lol = []
            tce_smry_dolol[TCE_id] = [tce_qu_seq, tce_loc_lol, bce_base_lol]

    '''###only for debug to check for index in list 'mtb_smry_dolol'
    ctr = 0
    for i in mtb_smry_dolol:
        if report_stat == 1: print("\n******Query", ctr, "******\n")
        for j in range(len(mtb_smry_dolol[i])):
            #line = mtb_smry_dolol[i][j]
    #       length = len(mtb_smry_dolol[i])
    #        if report_stat == 1: print("test mtb_smry_dolol:", len(mtb_smry_dolol[i]))
    #        for k in range(len(mtb_smry_dolol[i][j])):
            if j == 0:
                line = mtb_smry_dolol[i][j]
                if report_stat == 2: print("mtb_smry_dolol[", i,"][" , j, "]:\n", line)
        ctr+=1
    if report_stat == 1: print('\n\n\n')

    ###only for debug to check for index in list 'tce_smry_dolol'
    ctr = 0
    for i in tce_smry_dolol:
        if report_stat == 1: print("\n******Query", ctr, "******\n")
        for j in range(len(tce_smry_dolol[i])):
            #line = tce_smry_dolol[i][j]
    #       length = len(tce_smry_dolol[i])
    #        if report_stat == 1: print("test tce_smry_dolol:", len(tce_smry_dolol[i]))
    #        for k in range(len(tce_smry_dolol[i][j])):
            if j == 0:
                line = tce_smry_dolol[i][j]
                if report_stat == 2: print("tce_smry_dolol[", i,"][" , j, "]:\n", line)
        ctr+=1
    if report_stat == 1: print('\n\n\n')'''


    # to insert info into mtbseq dictionary
    def process(total_query, mtb_smry_dolol, qu_p, b_or_t, tce_smry_dolol={}):
        if report_stat == 1: print("linesT:", total_query, "mtb_smry_dolol:", mtb_smry_dolol)
        # process each query line by line
        i = 0
        while i < total_query:
            j = 0
            while j < (len(qu_p[i])):
                line = qu_p[i][j]  # this should be the first line of a  query
                while not re.match(r"^Query=\s*(\S+)", line):
                    j += 1
                    line = qu_p[i][j]
                q_id = re.findall(r"^Query=\s*(\S+)", line)[0]
                if report_stat == 1: print("\nProcessing Query ID = ", q_id)
                j += 1
                line = qu_p[i][j]
                while not re.match(r"^Length=(\d+)", line):
                    j += 1
                    line = qu_p[i][j]
                q_len = re.findall(r"^Length=(\d+)", line)[0]
                j += 1
                line = qu_p[i][j]
                while not (re.match(r"^Sequences producing", line) or re.match(r"(.+)No hits found(.+)", line)):
                    j += 1
                    line = qu_p[i][j]

                if re.match(r"^Sequences producing", line):
                    while not re.match(r"^>\s*(\S+)", line):
                        if report_stat == 2: print("i:", i, "j:,", j, "line:", line)
                        j += 1
                        line = qu_p[i][j]
                    while re.match(r"^>\s*(\S+)", line):
                        hit_id = re.findall(r"^>\s*(\S+)", line)[0]
                        j += 1
                        line = qu_p[i][j]
                        while not re.match(r"^ Identities = ((\d+)/(\d+)\s\((\d+)\%\))", line):
                            j += 1
                            line = qu_p[i][j]
                        x = re.findall(r"^ Identities = (\d+)/(\d+)\s\((\d+)\%\)", line)[0]
                        mismatch = int(x[1]) - int(x[0])
                        hit_subject_mlen = x[1]

                        # accepted
                        if int(hit_subject_mlen) >= int(q_len) and int(mismatch) <= 3:
                            if hit_id not in mtb_smry_dolol:
                                if report_stat == 2: print("hit id cannot be found: ", hit_id)
                            j += 1
                            line = qu_p[i][j]
                            while not re.match(r"^Query\s+\d+\s+(\S+)\s+\d+", line):
                                j += 1
                                line = qu_p[i][j]
                                qu_seq = re.findall(r"^Query\s+\d+\s+(\S+)\s+\d+", line)[0]
                            while not re.match(r"^Sbjct\s+(\d+)\s+\S+\s+(\d+)", line):
                                j += 1
                                line = qu_p[i][j]
                            if report_stat == 2: print("line:", line)
                            if report_stat == 2: print("list:", re.search(r"^Sbjct\s+(\d+)\s+\S+\s+(\d+)", line))
                            start = re.search(r"^Sbjct\s+(\d+)\s+(\S+)\s+(\d+)", line).group(1)
                            sub_seq = re.search(r"^Sbjct\s+(\d+)\s+(\S+)\s+(\d+)", line).group(2)
                            end = re.search(r"^Sbjct\s+(\d+)\s+(\S+)\s+(\d+)", line).group(3)
                            if report_stat == 1: print("start:", start)
                            if report_stat == 1: print("end:", end)
                            if b_or_t == 1:
                                mtb_smry_dolol[hit_id][3].append([q_id, qu_seq, sub_seq, start, end])  # to insert bce info

                            if b_or_t == 3:
                                # if q_id == "189575" and hit_id == "NP_217478.1":
                                # print("")
                                mtb_smry_dolol[hit_id][2].append([q_id, qu_seq, sub_seq, start, end])  # to insert tce info
                                # This structure is to accommodate for repetition of the same hit_id in a same q_id
                                tce_loc = str(start) + "-" + str(end)
                                tce_smry_dolol[q_id][1].append([hit_id, tce_loc])

                            # discarded query
                            # if int(hit_subject_mlen) < int(q_len) or int(mismatch) > 3:
                            if report_stat == 2: print("")

                        while not (re.match(r"^>\s*(\S+)", line) or re.match(r"^Lambda", line)):
                            j += 1
                            line = qu_p[i][j]

                # no hit
                if re.match(r"(.+)No hits found(.+)", line):
                    if report_stat == 2:   print("")
                    if report_stat == 2:  print("No hit found for Query ID=", q_id)

                while not j == (len(qu_p[i])):
                    j += 1

            i += 1


    process(total_queryB, mtb_smry_dolol, qu_B, 1)
    process(total_queryT, mtb_smry_dolol, qu_T, 3, tce_smry_dolol)

    if report_stat == 2: print("mtb_smry_dolol:", mtb_smry_dolol)
    if report_stat == 2: print("mtb_smry_dolol:BCE_id_lol:", mtb_smry_dolol[prtn_id][3])
    if report_stat == 2: print("mtb_smry_dolol:TCE_id_lol:", mtb_smry_dolol[prtn_id][2])
    if report_stat == 1: print("Done processing blast result of all queries.")
    if report_stat == 1: print("Proceed with the processing the H37Rv antigen that contain epitopes.")


    # to split B and T epitopes from its subject sequence
    for prtn_id in mtb_smry_dolol:
        # if i != 'NP_214954.1': # debuging
        # continue
        EpB_loc = []  # "start.end",..
        EpT_loc = []  # "start.end",..
        if len(mtb_smry_dolol[prtn_id][3]) > 0:
            for bce_list in mtb_smry_dolol[prtn_id][3]:
                EpB_loc.append(bce_list[3] + "-" + bce_list[4])
            if report_stat == 2: print("name:", prtn_id, " EpB_loc:", EpB_loc)

        if len(mtb_smry_dolol[prtn_id][2]) > 0:
            for tce_list in mtb_smry_dolol[prtn_id][2]:
                EpT_loc.append(tce_list[3] + "-" + tce_list[4])
            if report_stat == 2: print("name:", prtn_id, " EpT_loc:", EpT_loc)

        if prtn_id == "NP_217478.1":
            print("")

        # to mark the NP that have the overlapping B epitopes
        seq_digit_list = [0] * len(mtb_smry_dolol[prtn_id][0])
        for j in range(len(EpB_loc)):
            EpB_loc_strt = int(EpB_loc[j].split("-")[0]) - 1  # index strt 0
            EpB_loc_end = int(
                EpB_loc[j].split("-")[1])  # don't need to minus 1 coz in list slicing, end already minus 1
            seq_digit_list[EpB_loc_strt:EpB_loc_end] = [1] * (EpB_loc_end - EpB_loc_strt)  # end- strt x enough 1

        # to mark the NP that have the overlapping T epitopes and combined with B epitopes
        for j in range(len(EpT_loc)):
            EpT_loc_strt = int(EpT_loc[j].split("-")[0]) - 1  # index strt 0
            EpT_loc_end = int(
                EpT_loc[j].split("-")[1])  # don't need to minus 1 coz in list slicing, end already minus 1
            for k in range(EpT_loc_strt, EpT_loc_end):
                try:
                    if seq_digit_list[k] == 0 or seq_digit_list[k] == 1:
                        seq_digit_list[k] = seq_digit_list[k] + 2
                except:
                    if report_stat == 2: print("prtn_id:", prtn_id, "k:", k, "lengt:", len(seq_digit_list), "EpT loc:",
                                               EpT_loc)

        if report_stat == 2: print("name:", prtn_id, "Seq_digit:", seq_digit_list)
        mtb_smry_dolol[prtn_id][1] = seq_digit_list

    # to print table for NCBI_summary dictionary
    tbl_ctr = 0
    if report_stat == 1: print("NCBI_smry_table")
    if report_stat == 1: print(
        "Table counter no \t NCBI ID \t seq_ori \t seq_digit \t Tce_dict_id \t Tce_dict_qu_seq \t Tce_dict_sub_seq \t Tce_dict_start \t Tce_dict_stop \t Bce_dict_id \t Bce_dict_qu_seq \t Bce_dict_sub_seq \t Bce_dict_start \t Bce_dict_stop")
    for prtn_id in mtb_smry_dolol:
        tbl_ctr += 1
        seq_ori = mtb_smry_dolol[prtn_id][0]
        seq_digit = mtb_smry_dolol[prtn_id][1]
        Tce_list_ctr = len(mtb_smry_dolol[prtn_id][2])
        Bce_list_ctr = len(mtb_smry_dolol[prtn_id][3])

        if Tce_list_ctr == 0:
            continue

        Tce_dict_string = "{"
        for tce_list in mtb_smry_dolol[prtn_id][2]:
            q_id, qu_seq, sub_seq, start, stop = tce_list
            Tce_dict_string = Tce_dict_string + q_id + ":" + "[" + qu_seq + "," + sub_seq + "," + start + "," + stop + "]\t"
        Tce_dict_string = Tce_dict_string + "}"

        Bce_dict_string = "{"
        if Bce_list_ctr == 0:
            Bce_dict_string = "{no BCE}"

        for bce_list in mtb_smry_dolol[prtn_id][3]:
            q_id, qu_seq, sub_seq, start, stop = bce_list
            Bce_dict_string = Bce_dict_string + q_id + ":" + "[" + qu_seq + "," + sub_seq + "," + start + "," + stop + "]\t"
        Bce_dict_string = Bce_dict_string + "}"

        if len(seq_ori) != len(seq_digit):
            if report_stat == 1: print("Original sequence and coded sequence is not match")

        if report_stat == 1: print(tbl_ctr, prtn_id, "\t", seq_ori, "\t", seq_digit, "\t", Tce_dict_string, "\t",
                                   Bce_dict_string)

    # Structure: mtb_smry_dolol = {Np_id:[seq_normal, seq_coded,[TCE_id, qu_seq,sub_seq, st,sp], [BCE_id, qu_seq, sub_seq, st,sp]]}
    # Structure: tce_smry_dolol = {TCE_id:[tce_qu_seq,[[prtn_id, tce_loc]], [[prtn_id, bce_base_ovlp]] }
    # counting bce_ovlp within each tce
    for tce_id in tce_smry_dolol:
        tce_loc_lol = tce_smry_dolol[tce_id][1]

        if len(tce_loc_lol) > 0:
            # if tce_id == "189575" and tce_loc_lol[0][0] == "NP_217478.1":
            # print("")
            for tce_loc_list in tce_smry_dolol[tce_id][1]:
                prtn_id = tce_loc_list[0]

                # debugger: check if tce has more than 1 occurences in the same protein
                #       if len(tce_loc_list) > 1:
                #          if report_stat == 2: print("check this out")

                tce_loc_strt = int(tce_loc_list[1].split("-")[0]) - 1  # -1 to offset for 0 base index
                tce_loc_end = int(
                    tce_loc_list[1].split("-")[1])  # no need -1 bcz slicing only take b4 this location
                tce_code_list = mtb_smry_dolol[prtn_id][1][tce_loc_strt:tce_loc_end]
                if tce_loc_end == len(
                        mtb_smry_dolol[prtn_id][0]):  # to handle cases where tce_loc_end is equal tce_seq_end
                    tce_code_list = mtb_smry_dolol[prtn_id][1][tce_loc_strt:]
                if report_stat == 1: print("tce_id=", tce_id, ":", prtn_id, "coded:", mtb_smry_dolol[prtn_id][1])
                if sum(tce_code_list) == len(
                        tce_code_list) * 2:  # if the sum of tce base code == tce len * 2, the is no bce & tce overlapped
                    tce_loc_list.append("No BCE Overlap")
                    # tce_loc_strt + 1 to reset the 0 index, thus start count index with 1
                    if report_stat == 1: print("tce_id=", tce_id, ":", prtn_id, "tce_loc_start:", tce_loc_strt + 1,
                                               "tce_loc_end:", tce_loc_end, "tce_code_list:", tce_code_list,
                                               "=> No BCE Overlapped")
                elif sum(tce_code_list) >= len(tce_code_list) * 2:
                    tce_bce_ovlp_ratio = str(sum(tce_code_list) - len(tce_code_list) * 2) + "/" + str(
                        len(tce_code_list))
                    tce_loc_list.append(tce_bce_ovlp_ratio)
                    # tce_loc_strt + 1 to reset the 0 index, thus start count index with 1
                    if report_stat == 1: print("tce_id=", tce_id, ":", prtn_id, "tce_loc_start:", tce_loc_strt + 1,
                                               "tce_loc_end:", tce_loc_end, "tce_code_list:", tce_code_list,
                                               "tce_bce_ovlp_ratio:", tce_bce_ovlp_ratio)
            if report_stat == 1: print("tce_id=", tce_id, ":", prtn_id, "tce_loc_list", tce_loc_list)
        else:
            if report_stat == 1: print("tce_id=" + tce_id + ", no prtn_id, tce_loc_list: []")

    # Structure: mtb_smry_dolol = {Np_id:[seq_normal, seq_coded,[TCE_id, qu_seq,sub_seq, st,sp], [BCE_id, qu_seq, sub_seq, st,sp]}
    # Structure: tce_smry_dolol = {TCE_id:[tce_qu_seq,[[prtn_id, tce_loc]], [[prtn_id, bce_base_ovlp]] }
    # to print table for TCE_summary dictionary
    tbl_ctr = 0
    if report_stat == 1: print("TCE_smry_table")
    if report_stat == 0: print("Table counter no \t TCE ID \t tce_qu_seq \t prtn_id \t tce_count_NP \t bce_base_ovlp")
    for tce_id in tce_smry_dolol:
        tce_qu_seq = tce_smry_dolol[tce_id][0]
        tce_loc_lol = tce_smry_dolol[tce_id][1]
        tbl_ctr += 1
        tce_np_ctr_dict = {}
        tce_np_bceovlp_dict = {}

        for tce_loc_list in tce_loc_lol:
            prtn_id = tce_loc_list[0]
            bce_ovlp = tce_loc_list[2]
            if prtn_id in tce_np_ctr_dict:
                tce_np_ctr_dict[prtn_id] += 1
                tce_np_bceovlp_dict[prtn_id] += ", " + bce_ovlp
            else:
                tce_np_ctr_dict[prtn_id] = 1
                tce_np_bceovlp_dict[prtn_id] = bce_ovlp

        NP_id_list = []
        tce_count_NP_list = []
        bce_base_list = []
        for prtn_id, ctr in tce_np_ctr_dict.items():
            NP_id_list.append(prtn_id)
            tce_count_NP_list.append(ctr)
            bce_base_list.append(tce_np_bceovlp_dict[prtn_id])

        if report_stat == 1: print(tbl_ctr, "\t", tce_id, "\t", tce_qu_seq, "\t", NP_id_list, "\t",
                                   tce_count_NP_list, "\t", bce_base_list)

    if cmd_stat == 1:
        filehandler_scpt3_1 = open('mtb_smry_dolol.obj','wb')  # to save memory to file so do not need to run cmd everytime
        pickle.dump(mtb_smry_dolol, filehandler_scpt3_1)
        filehandler_scpt3_1.close()

    if cmd_stat == 1:
        filehandler_scpt3_2 = open('tce_smry_dolol.obj', 'wb')  # to save memory to file so do not need to run cmd everytime
        pickle.dump(tce_smry_dolol, filehandler_scpt3_2)
        filehandler_scpt3_2.close()
    if report_stat == 1: print("--- script3 processing in:%s seconds ---" % (time.time() - startscript3_time))

filehandler_scpt3_1 = open('mtb_smry_dolol.obj', 'rb')  # comment this first and uncomment it after run cmd
mtb_smry_dolol = pickle.load(filehandler_scpt3_1)
filehandler_scpt3_1.close()

filehandler_scpt3_2 = open('tce_smry_dolol.obj', 'rb')  # comment this first and uncomment it after run cmd
tce_smry_dolol = pickle.load(filehandler_scpt3_2)
filehandler_scpt3_2.close()


### xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx SCRIPT 4a : MHC1 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx###
if mhc_all_swc or mhc_1x_swc or mhc_mutual_swc:
    if report_stat >= 0: print("\n\n############## SCRIPT 4a: MHC1 PREDICTION ####################\n")
    # processing tce_mhc1_ep file
    fx.file_val(tce_mhc1_ep)
    tce_mhc1_f = open(tce_mhc1_ep, "r")
    tce_mhc1 = tce_mhc1_f.read().splitlines()
    # tce_gbl_dict = {8:{tce_id:seq}, 9:{tce_id:seq, tce_id:seq}}
    tce_mhc1_gbl_dod = {'7': {}, '8': {}, '9': {}, '10': {}, '11': {}}  # gbl = group by length

    tce_id = ""
    tce_mhc1_dos = {}
    for i in range(len(tce_mhc1)):
        if re.match(r"^>(\S+)", tce_mhc1[i]):
            tce_id = re.findall(r"^>(\S+)", tce_mhc1[i])[0]
        if not re.match(r"^>(\S+)", tce_mhc1[i]):
            seq = tce_mhc1[i].rstrip()
            len_seq = len(seq)
            if len_seq >= 11:
                tce_seq_dict = tce_mhc1_gbl_dod['11']  #assign memory location, not value
                tce_seq_dict[tce_id] = seq
            else:
                tce_seq_dict = tce_mhc1_gbl_dod[str(len_seq)]
                tce_seq_dict[tce_id] = seq
            tce_mhc1_dos[tce_id] = seq
    if report_stat==3: print("tce_mhc1_gbl_dod", tce_mhc1_gbl_dod)

    #to debug the tce_gbl_dict
    for tce_len in tce_mhc1_gbl_dod:
        if debug_stat == 1: print("tce_gbl_dict_no:", tce_len)
        for seq in tce_mhc1_gbl_dod[tce_len]:
            if debug_stat == 1: print(tce_len, tce_mhc1_gbl_dod[tce_len][seq])

    mhc1_alle_gbl_dict = OrderedDict()  # to initiallise mhc1_alle_gbl_dict
    for length in range(8, 12): #12 because len needed is 11
        mhc1_alle_gbl_dict[str(length)] = []

    ''' #ONLY FOR NETMHCCONS METHOD 
    # processing mhc1_cnns --> change to IEDB_recommended file  # mhc1_alle_gbl_dict= {8:[mhc, mhc, mhc], 9:[mhc, mhc, mhc]}
    fx.file_val(mhc1_netmhc_alleles)
    mhc1_netmhc_alleles = open(mhc1_netmhc_alleles, "r")
    mhc1_netmhc_alleles_file = mhc1_netmhc_alleles.read().splitlines()
    
    for i in range(len(mhc1_netmhc_alleles_file)):
        mhc, mhc_len = mhc1_netmhc_alleles_file[i].split(",")
        mhc1_alle_gbl_dict[mhc_len].append(mhc)
    if report_stat == 1: print("mhc1_alle_gbl_dict", mhc1_alle_gbl_dict)

    # to debug the mhc1_alle_gbl_dict
    for mhc_len in mhc1_alle_gbl_dict:
        if debug_stat == 1: print("mhc1_gbl_dict_no:", mhc_len)
        for seq in range(len(mhc1_alle_gbl_dict[mhc_len])):
            if debug_stat == 1: print(mhc_len, mhc1_alle_gbl_dict[mhc_len][seq])
    '''
    fx.file_val(mhc1_iedbrec_alleles)
    mhc1_iedbrec_alleles_file = open(mhc1_iedbrec_alleles, "r")
    mhc1_iedbrec_alleles_list = mhc1_iedbrec_alleles_file.read().splitlines()
    startmhc1_time = time.time()

    # this dictionary is for the conversion of seq_no using tce file (IEDB) to convert the seq id in the mhc1_output file.
    # seqno_to_tceid_gbl_dod = {'8': {seq_no->tce_id}, '9': {}, '10': {}, '11': {}}
    seqno_to_tceid_gbl_dod = {'8': {}, '9': {}, '10': {}, '11': {}}  #minimum allele length for tce to run in predictor

    # to produce mhc1 parameter, one length by one length
    for mhc_len_str in mhc1_alle_gbl_dict.keys(): # MHC1 HLA length have 8-11
        # to produce file from tce_mhc1_gbl_dod (mhc(one) to tce(many)), one lenght have many tce_id and seq
        with open("tce_gbl_file.faa", 'w') as tce_file:
            seq_no_ctr = 0
            for tce_lngth_str, tce_seq_dict in tce_mhc1_gbl_dod.items():
                if int(tce_lngth_str) >= int(mhc_len_str):
                    for tce_id, tce_seq in tce_seq_dict.items():
                        if report_stat == 2: print("tce_lngth_str:", tce_lngth_str, " tce:", tce_id, " tce_seq:", tce_seq, "\n")
                        seq_no_ctr += 1
                        seqno_to_tceid_gbl_dod[mhc_len_str][seq_no_ctr] = tce_id
                        tce_file.write('>%s\n%s\n' % (tce_id, tce_seq))
        tce_file.close()

        # bcoz allele>20k cannot run this command in bash in single run, decide to loop this 30x coz in total there is 2915 alleles for each length
        # in mhc_i_netmhccons_alleles_full.txt file
        if not os.path.exists('mhc1_list'):
            os.makedirs('mhc1_list')  # if mhc_1x_swc or mhc_all_swc or mhc_mutual_swc: as a switch to run in linux

        ''' #now we decide to change the list of allele from 2915 to 90+, thus no need to use index anymore
        for idx in range(1, 31):
            mhc_list_curr = mhc1_alle_gbl_dict[mhc_len_str]
            if idx < 30:
                mhc_list_curr = mhc_list_curr[(idx - 1) * 100:(idx * 100)]
            else:
                mhc_list_curr = mhc_list_curr[(idx - 1) * 100:] # to handle the last length because it only have until 2915
            param1_str = ",".join([str(x) for x in mhc_list_curr])
            param2_str = (str(mhc_len_str) + ',') * len(mhc_list_curr)
            param2_str = param2_str[:-1]  # It slices the string to omit the last "," character
            if report_stat == 2: print("param1:", param1_str)
            if report_stat == 2: print("param2:", param2_str)

            fx.file_val("tce_gbl_file.faa") # to enable auto run in server
            output_file = "mhc1_output" + mhc_len_str + ".txt"
            raw_op_file = "mhc1_raw_output" + mhc_len_str + "_idx" + str(idx) + ".txt"
            cmd = "../mhc_i_3.1.1/src/predict_binding.py netmhccons " + param1_str + " " + param2_str + " tce_gbl_file.faa > mhc1_list/" + raw_op_file
            if report_stat == 1: print("cmd", cmd)
            os.system(cmd)
        
        # this step is to filter useful columns with percentile rank <=1%
        ttl_raw_op_file = "mhc1_raw_output" + mhc_len_str + "_idx*.txt"
        cmd2 = "perl -F\"\\t\" -lane 'print \"$F[0]\\t$F[1]\\t$F[2]\\t$F[3]\\t$F[4]\\t$F[5]\\t$F[6]\\t$F[7]\" if $F[7] <= 1 && ($F[7]=~/^\\d/)' mhc1_list/" + ttl_raw_op_file + " | grep -v allele > mhc1_list/" + output_file
        
        # output_file has redundant contents due to it was extracted from 30 files. Check if the redundancy affect logic error below
        # need to remove the raw file after filtered  to a new file to save space
        cmd3 = "rm mhc1_list/" + ttl_raw_op_file
        
        '''

        param1_str = ",".join([str(x) for x in mhc1_iedbrec_alleles_list])
        param2_str = (str(mhc_len_str) + ',') * len(mhc1_iedbrec_alleles_list)
        param2_str = param2_str[:-1]  # It slices the string to omit the last "," character

        if report_stat == 2: print("param1:", param1_str)
        if report_stat == 2: print("param2:", param2_str)

        fx.file_val("tce_gbl_file.faa")  # to enable auto run in server
        output_file = "mhc1_output" + mhc_len_str + ".txt"
        raw_op_file = "mhc1_raw_output" + mhc_len_str + ".txt"
        cmd = "../mhc_i_3.1.1/src/predict_binding.py IEDB_recommended " + param1_str + " " + param2_str + " tce_gbl_file.faa > mhc1_list/" + raw_op_file
        if report_stat == 1: print("cmd", cmd)
        if cmd_stat == 1:os.system(cmd)

        # this step is to filter useful columns with percentile rank <=1%
        cmd2 = "perl -F\"\\t\" -lane 'print \"$F[0]\\t$F[1]\\t$F[2]\\t$F[3]\\t$F[4]\\t$F[5]\\t$F[6]\\t$F[7]\\t$F[9]\" if $F[9] <= 1 && ($F[9]=~/^\\d/)' mhc1_list/" + raw_op_file + " | grep -v allele > mhc1_list/" + output_file
        print("cmd2: ", cmd2)
        if cmd_stat == 1:os.system(cmd2)
        # output_file has redundant contents due to it was extracted from 30 files. Check if the redundancy affect logic error below
        # need to remove the raw file after filtered  to a new file to save space
        cmd3 = "rm mhc1_list/" + raw_op_file
        if cmd_stat == 1:os.system(cmd3)

    filehandler1 = open('seqno_to_tceid_gbl_dod.obj', 'wb')  # to save memory to file so do not need to run cmd everytime
    pickle.dump(seqno_to_tceid_gbl_dod, filehandler1)
    filehandler1.close()

# open when want to debug
if mhc1_debug_swc:
    filehandler1 = open('seqno_to_tceid_gbl_dod.obj', 'rb') #comment this first and uncomment it after run cmd
    seqno_to_tceid_gbl_dod = pickle.load(filehandler1)
    filehandler1.close()

# to process data file after the MHC1 prediction
# to open mhc_output file
# to create a dictionary for mhc1 output
# do the function to make sure the result filtered kept in the dictionary
# tce_mhc_op_dol = {tce_id:[mhc1_info_list[], mhc1_info_list[], mhc1_info_list[],.. ]}
# mhc1_debug_info_list consist of allele, consensus_percentile_rank and ANN_ic_50
tce_mhc_op_dol = {}  # mhcii result will also be kept here  ### pls open this even to close mhc1 loop
if mhc_all_swc or mhc_1x_swc or mhc_mutual_swc or mhc1_debug_swc:
    for mhc_gbl_int in seqno_to_tceid_gbl_dod.keys():
        output_file_name = "mhc1_list/mhc1_output" + str(mhc_gbl_int) + ".txt"  #"mhc1_list/mhc1_output"
        output_f = open(output_file_name, "r")
        output_mhc1 = output_f.read().splitlines()
        for line_num in range(len(output_mhc1)):
            if line_num == 0:  # skip header line
                continue
            else:
                line_cols = output_mhc1[line_num].split("\t")
                if report_stat == 1: print("line column", line_cols)
                seq_no = line_cols[1]
                print("seq_no: ", seq_no)
                tce_id = seqno_to_tceid_gbl_dod[mhc_gbl_int][int(seq_no)] # use seq_no because it has the same usage as key[ctr] above and its interchangeble
                allele = line_cols[0]
                tce_strt = line_cols[2]
                tce_end = line_cols[3]
                tce_lngth = line_cols[4]
                allele_seq = line_cols[5]
                pred_prcntile_rank = line_cols[7]
                # ANN_ic_50 = line_cols[7]
                mhc_debug_info_list = [allele, tce_strt, tce_end, tce_lngth, allele_seq, pred_prcntile_rank]
                if tce_id in tce_mhc_op_dol.keys():
                    tce_mhc_op_dol[tce_id][0].add(allele)
                    tce_mhc_op_dol[tce_id][1].append(mhc_debug_info_list)
                else:
                    tce_mhc_op_dol[tce_id] = [set(), []]
                    tce_mhc_op_dol[tce_id][0].add(allele)
                    tce_mhc_op_dol[tce_id][1].append(mhc_debug_info_list)
    if report_stat == 1: print("--- mhc1 processing in:%s seconds ---" % (time.time() - startmhc1_time))

### xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx SCRIPT 4b : MHC2 PREDICTION xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx###
if mhc_all_swc or mhc_2x_swc or mhc_mutual_swc:
    if report_stat >= 0: print("\n\n############## SCRIPT 4b: MHC2 PREDICTION ####################\n")
    # to open file for mhc2 allele and iterate to put in string
    # processing tce_mhc2_ep file
    fx.file_val(tce_mhc2_ep)
    tce_mhc2_f = open(tce_mhc2_ep, "r")
    tce_mhc2 = tce_mhc2_f.read().splitlines()
    # tce_gbl_dict = {15:{tce_id:seq}, 16:{tce_id:seq, tce_id:seq}}  # gbl = group by length
    tce_mhc2_gbl_dict = OrderedDict()
    for length in range(15, 31): #31 because len needed is 30
        tce_mhc2_gbl_dict[str(length)] = {}

    tce_id = ""
    tce_mhc2_dos = {}
    tce_mhc2_skip_dict = {}  # to handle the filtered tce id len <15 in printing loop
    for i in range(len(tce_mhc2)):
        if re.match(r"^>(\S+)", tce_mhc2[i]):
            tce_id = re.findall(r"^>(\S+)", tce_mhc2[i])[0]
        if not re.match(r"^>(\S+)", tce_mhc2[i]):
            seq = tce_mhc2[i].rstrip()
            len_seq = len(seq)
            if len_seq < 15:  # 15 bcz of the predictor limitation
                tce_mhc2_skip_dict[tce_id] = ""
                continue
            if len_seq >= 20:  # 20 is the max length for practicality & optimization
                tce_seq_dict = tce_mhc2_gbl_dict['20']
                tce_seq_dict[tce_id] = seq
            else:
                tce_seq_dict = tce_mhc2_gbl_dict[str(len_seq)]
                tce_seq_dict[tce_id] = seq
            tce_mhc2_dos[tce_id] = seq

    filehandler_mhc_skip = open('tce_mhc2_skip_dict.obj','wb')  # to save memory to file so do not need to run cmd everytime
    pickle.dump(tce_mhc2_skip_dict, filehandler_mhc_skip)
    filehandler_mhc_skip.close()

    if report_stat == 1: print("tce_mhc2_gbl_dict", tce_mhc2_gbl_dict)
    startmhc2_time = time.time()

    '''
    all_mhc2_cnss_str = ("HLA-DPA1*01/DPB1*04:01,HLA-DRB1*01:01,HLA-DRB1*04:04,HLA-DRB1*08:04,HLA-DRB1*11:21," 
    "HLA-DPA1*01:03/DPB1*02:01,HLA-DRB1*01:02,HLA-DRB1*04:05,HLA-DRB1*08:06,HLA-DRB1*11:28,"  
    "HLA-DPA1*02:01/DPB1*01:01,HLA-DRB1*03:01,HLA-DRB1*04:08,HLA-DRB1*08:13,HLA-DRB1*13:01,HLA-DRB1*15:01," 
    "HLA-DPA1*02:01/DPB1*05:01,HLA-DRB1*03:05,HLA-DRB1*04:10,HLA-DRB1*08:17,HLA-DRB1*13:02,HLA-DRB1*15:02,"  
    "HLA-DPA1*03:01/DPB1*04:02,HLA-DRB1*03:06,HLA-DRB1*04:21,HLA-DRB1*11:01,HLA-DRB1*13:04,HLA-DRB1*15:06,"  
    "HLA-DQA1*01:01/DQB1*05:01,HLA-DRB1*03:07,HLA-DRB1*04:23,HLA-DRB1*11:02,HLA-DRB1*13:05,HLA-DRB3*01:01,"
    "HLA-DQA1*01:02/DQB1*06:02,HLA-DRB1*03:08,HLA-DRB1*04:26,HLA-DRB1*11:04,HLA-DRB1*13:07,HLA-DRB4*01:01,"
    "HLA-DQA1*03:01/DQB1*03:02,HLA-DRB1*03:09,HLA-DRB1*07:01,HLA-DRB1*11:06,HLA-DRB1*13:11,HLA-DRB5*01:01,"
    "HLA-DQA1*04:01/DQB1*04:02,HLA-DRB1*07:03,HLA-DRB1*11:07,HLA-DRB1*13:21,HLA-DRB5*01:05,"
    "HLA-DQA1*05:01/DQB1*02:01,HLA-DRB1*04:01,HLA-DRB1*08:01,HLA-DRB1*11:14,HLA-DRB1*13:22,"
    "HLA-DQA1*05:01/DQB1*03:01,HLA-DRB1*04:02,HLA-DRB1*08:02,HLA-DRB1*11:20,HLA-DRB1*13:23,HLA-DRB1*13:27,HLA-DRB1*13:28")

    all_mhc2_cnss_list = all_mhc2_cnss_str.split(",")
    '''

    fx.file_val(mhc2_iedbrec_alleles)
    all_mhc2_cnss_f = open(mhc2_iedbrec_alleles, "r")
    all_mhc2_cnss_list = all_mhc2_cnss_f.read().splitlines()
    all_mhc2_cnss_str = ",".join([str(x) for x in all_mhc2_cnss_list])

    # this dictionary is for the conversion of seq_no using tce file (IEDB) to convert the seq id in the mhc1_output file.
    # seqno_to_tceid_gbl_dod = {seq_no->tce_id
    seqno_to_tceid_mhc2_gbl_dict = {'15': {}, '16': {},'17': {}, '18': {}, '19': {}, '20': {}}  # gbl = group}

    # to group the tce into groups of 15-to-20 and above for MHC2 processing
    # to produce a file for tce_lngth to run mhc2
    for tce_lngth_str, tce_seq_dict in tce_mhc2_gbl_dict.items():
        with open("tce_mhc2_gbl_file.faa", 'w') as tce_mhc2_file:
            seq_no_ctr = 0
            if int(tce_lngth_str) >= 15:
                for tce_id, tce_seq in tce_seq_dict.items():
                    seq_no_ctr += 1 # in o/p file predict_binding, seq_no start with 1
                    seqno_to_tceid_mhc2_gbl_dict[tce_lngth_str][seq_no_ctr] = tce_id
                    tce_mhc2_file.write('>%s\n%s\n' % (tce_id, tce_seq))
                    if report_stat == 2: print("tce_mhc2_file:", tce_mhc2_file)
        tce_mhc2_file.close()

        fx.file_val("tce_mhc2_gbl_file.faa") # to enable auto run in server
        output_file = "mhc2_output" + tce_lngth_str + ".txt"
        if int(tce_lngth_str) == 15:
            cmd = "../mhc_ii_3.1.2/mhc_II_binding.py IEDB_recommended" + " " + all_mhc2_cnss_str + " " + "tce_mhc2_gbl_file.faa 15 > " + output_file
        else:
            cmd = "../mhc_ii_3.1.2/mhc_II_binding.py IEDB_recommended" + " " + all_mhc2_cnss_str + " " + "tce_mhc2_gbl_file.faa 15-" + tce_lngth_str + " > " + output_file
        if report_stat == 1: print("cmd", cmd)
        if cmd_stat == 1:os.system(cmd)
        cmd2 = "sort " + output_file + " | uniq | perl -F\"\\t\" -lane 'print \"$F[0]\\t$F[1]\\t$F[2]\\t$F[3]\\t$F[4]\\t$F[5]\\t$F[6]\\t$F[7]\\t$F[13]\"' > " + "mhc2_output" + tce_lngth_str + "_filtered.txt"
        if cmd_stat == 1:os.system(cmd2)
        cmd3 = "mv " + "mhc2_output" + tce_lngth_str + "_filtered.txt " + output_file
        if cmd_stat == 1:os.system(cmd3)

    filehandler2 = open('seqno_to_tceid_mhc2_gbl_dict.obj', 'wb')  # to save memory to file so do not need to run cmd everytime
    pickle.dump(seqno_to_tceid_mhc2_gbl_dict, filehandler2)
    filehandler2.close()

if mhc2_debug_swc:
    filehandler2 = open('seqno_to_tceid_mhc2_gbl_dict.obj', 'rb') #comment this first and uncomment it after run cmd
    seqno_to_tceid_mhc2_gbl_dict = pickle.load(filehandler2)
    filehandler2.close()

# to process data file after the MHC2 prediction
# to create a dictionary for mhc2 output
# do the function to make sure the result filtered kept in the dictionary
if mhc_all_swc or mhc_2x_swc or mhc_mutual_swc or mhc2_debug_swc:
    for mhc_gbl_int in seqno_to_tceid_mhc2_gbl_dict.keys():
        output_file_name = "mhc2_output" + mhc_gbl_int + ".txt"
        output_f = open(output_file_name, "r")
        output_mhc2 = output_f.read().splitlines()
        for i in range(len(output_mhc2)):
            if i == 0:  # skip header line
                continue
            else:
                line_cols = output_mhc2[i].split("\t")
                # matched_list = re.findall(r"^(\S+)\t(\d+)\t(\d+\.\d+)\t(\d+\.\d+)", output_mhc2[i])
                # mhc2_output_dict = {tce_id:mhc2_info_list}
                # mhc2_info_list consist of allele, consensus_percentile_rank and SMM_ic_50
                seq_no = line_cols[1]
                tce_id = seqno_to_tceid_mhc2_gbl_dict[mhc_gbl_int][int(seq_no)] #use seq_no because it has the same usage as key[ctr] above and its interchangeble
                allele = line_cols[0]
                tce_strt = line_cols[2]
                tce_end = line_cols[3]
                tce_lngth = line_cols[4]
                allele_seq = line_cols[6]
                pred_prcntile_rank = line_cols[7]
                # SMM_ic_50 = line_cols[8]
                mhc_debug_info_list = [allele, tce_strt, tce_end, tce_lngth, allele_seq, pred_prcntile_rank]
                if float(pred_prcntile_rank) <= 10:
                    if tce_id in tce_mhc_op_dol.keys():
                        tce_mhc_op_dol[tce_id][0].add(allele) # will be index 1 or above cz list have ordering
                        tce_mhc_op_dol[tce_id][1].append(mhc_debug_info_list)
                    else:
                        tce_mhc_op_dol[tce_id] = [set(), []]
                        tce_mhc_op_dol[tce_id][0].add(allele) # will be index 0 cz list have ordering
                        tce_mhc_op_dol[tce_id][1].append(mhc_debug_info_list)

    if report_stat == 2: print("CHECK PICKLE! tce_mhc_op_dol:", tce_mhc_op_dol)
    if report_stat == 1: print("--- mhc2 processing in:%s seconds ---" % (time.time() - startmhc2_time))

filehandler3 = open('tce_mhc_op_dol.obj', 'wb')  # to save memory to file so do not need to run cmd everytime
pickle.dump(tce_mhc_op_dol, filehandler3)
filehandler3.close()

if mhc1_debug_swc or mhc2_debug_swc:
    filehandler3 = open('tce_mhc_op_dol.obj', 'rb') #comment this first and uncomment it after run cmd
    tce_mhc_op_dol = pickle.load(filehandler3)
    filehandler3.close()


### xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx SCRIPT 4c : POPULATION COVERAGE xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ###
if popcov_swc:
    if report_stat >= 0: print("\n\n############## SCRIPT 4c : POPULATION COVERAGE ####################\n")

    # tce_mhc_op_dol = {tce_id:[mhc1_info_list[], mhc1_info_list[], mhc1_info_list[],.. ]}
    # mhc1_debug_info_list consist mhc_debug_info_list = [allele, tce_strt, tce_end, tce_lngth, allele_seq, pred_prcntile_rank]
    if tce_pickle_swc:
        tce_mhc_op_dol = {}
        tce_mhc_op_dol_pickle_list = os.listdir('tce_mhc_op_dol_pickle_list')  # make a list of mhc_pred_pickle files
        for idx in range(len(tce_mhc_op_dol_pickle_list)):
            tce_mhc_op_dol_file = 'tce_mhc_op_dol_pickle_list/' + tce_mhc_op_dol_pickle_list[idx]
            filehandler_scpt3 = open(tce_mhc_op_dol_file, 'rb')  # comment this first and uncomment it after run cmd
            tce_mhc_op_dol_subset = pickle.load(filehandler_scpt3)
            filehandler_scpt3.close()
            for tce_id, mhc_debug_info_lol in tce_mhc_op_dol_subset.items():
                #if tce_id == "29256":
                    #print("pickle file name:", tce_mhc_op_dol_pickle_list[idx], "tce_id:", tce_id, "mhc debug info:", mhc_debug_info_list)
                tce_mhc_op_dol[tce_id] = mhc_debug_info_lol #TCE in all subset file is not overlap to each other

    # remove redundant mhc allele   #set return the list into tuple thus, have to change data type to list again by list(....)
    # In popcov perspective, we are more interested to look into collective alleles so that we will get higher % for all populations
    startpc_time = time.time()
    popcov_output_dod = {} # popcov_output_dod = {tce_id: {Alleles:set(), world:%, Asian:%, Southeast Asia:%, Malaysia:%}}
    for tce_id, allele_n_info in tce_mhc_op_dol.items():
        allele_set = allele_n_info[0]
        new_allele_set = set()
        for allele in allele_set:
            dbl_allele = allele.split('/') #if no "/" present, return a list with original strings
            if len(dbl_allele) == 2: #because the "/" invoving 2 alleles
                #print("pop_cov_op:dbl allele", dbl_allele)
                new_allele_set.add(dbl_allele[0])
                new_allele_set.add("HLA-" + dbl_allele[1])
            else:
                new_allele_set.add(allele)
            #print("new_allele_set:", new_allele_set)
        allele_set = new_allele_set
        mhc_debug_info_list = allele_n_info[1]
        popcov_output_dod[tce_id] = {}
        popcov_output_dod[tce_id]['Alleles'] = allele_set
        popcov_output_dod[tce_id]['World'] = None
        popcov_output_dod[tce_id]['Asian'] = None
        popcov_output_dod[tce_id]['Southeast Asia'] = None
        popcov_output_dod[tce_id]['Malaysia'] = None
        if report_stat == 1: print("popcov_output_dict_dodod:", popcov_output_dod)

    # to process the input file for population coverage prediction
    for tce_id in popcov_output_dod.keys():
        allele_set = popcov_output_dod[tce_id]['Alleles']
        with open("ip_population_cov.txt", 'w') as ip_pop_cov_file:
            ip_pop_cov_file.write('%s    %s\n' % (tce_id, ','.join(str(s) for s in allele_set)))
        ip_pop_cov_file.close()
        #for allele, info_dict in alleles_dict.items():
        cmd = "python ../population_coverage/calculate_population_coverage.py -p World,Asian,'Southeast Asia',Malaysia -c combined -f ip_population_cov.txt > op_population_cov.txt"
        if cmd_stat == 1: os.system(cmd)

        # processing output file
        try:
            pc_f = open("op_population_cov.txt", "r")
            lines = pc_f.read().splitlines()
            pc_f.close()

        except IOError:
            if report_stat == 1: print('op_population_cov.txt not found')

        total_line = len(lines)
        for i in range(total_line):
            if re.match(r"^World\t(\d+.\d+%)\t", lines[i]):
                percentage = round(float(re.findall(r"^World\t(\d+.\d+)%\t", lines[i])[0])/100,4)  #"{:.4f}".format = to convert into 0.0000 in string
                popcov_output_dod[tce_id]['World'] = str(percentage)
            if re.match(r"^Asian\t(\d+.\d+%)\t", lines[i]):
                percentage = round(float(re.findall(r"^Asian\t(\d+.\d+)%\t", lines[i])[0])/100, 4)
                popcov_output_dod[tce_id]['Asian'] = str(percentage)
            if re.match(r"^Southeast Asia\t(\d+.\d+%)\t", lines[i]):
                percentage = round(float(re.findall(r"^Southeast Asia\t(\d+.\d+)%\t", lines[i])[0])/100, 4)
                popcov_output_dod[tce_id]['Southeast Asia'] = str(percentage)
            if re.match(r"^Malaysia\t(\d+.\d+%)\t", lines[i]):
                percentage = round(float(re.findall(r"^Malaysia\t(\d+.\d+)%\t", lines[i])[0])/100, 4)
                popcov_output_dod[tce_id]['Malaysia'] = str(percentage)
        if report_stat == 1: print("popcov_output_dod[" + tce_id + "]" + str(popcov_output_dod[tce_id]))

    if cmd_stat == 1:
        filehandler4 = open('popcov_output_dod.obj', 'wb')  # to save memory to file so do not need to run cmd everytime
        pickle.dump(popcov_output_dod, filehandler4)
        filehandler4.close()

filehandler4 = open('popcov_output_dod.obj', 'rb') #comment this first and uncomment it after run cmd
popcov_output_dod = pickle.load(filehandler4)
filehandler4.close()


### xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx SCRIPT 4d: PREDICTING IMMUNOGENIC REGION USING BCE AND TCE OVERLAPPING  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ###
if IR_swc:
    if report_stat >= 0: print("\n\n############## SCRIPT 4d: PREDICTING IMMUNOGENIC REGION USING BCE AND TCE OVERLAPPING ##############\n")
    # to find the overlapping between BCE and TCE (without considering in the TCE proteins)
    startir_time = time.time()

    # make a dictionary for tce_id and tce_seq
    fx.file_val(tce_all_ep_file)
    tce_all_f = open(tce_all_ep_file, "r")
    tce_all_lines = tce_all_f.read().splitlines()
    tce_all_f.close()
    tce_id_str = ""
    tce_seq_str = ""
    tce_all_dos = OrderedDict()  # to make a dictionary for tce_seq for population coverage input file
    for line in range(len(tce_all_lines)):
        if re.match(r"^>(\d+)", tce_all_lines[line]):
            tce_id_str = re.findall(r"^>(\d+)", tce_all_lines[line])[0]
        if not re.match(r"^>(\d+)", tce_all_lines[line]):
            tce_seq_str = tce_all_lines[line].rstrip()
            tce_all_dos[tce_id_str] = tce_seq_str
    if report_stat == 1: print("tce_all_dos:", tce_all_dos)

    # run blast bce_ep as query, tce_ep as subject because the percentage is calculated based on query
    # blastp -query export_table__BCell_out.faa -subject export_table__TCell_out.faa -evalue 1e-2 -out blastp_bce_op.txt

    cmd_bce_blast = "blastp -task blastp-short -query " + bce_ep + " -subject " + tce_all_ep_file + " -evalue 1e-2 -out blastp_bce_e2.txt"
    if cmd_stat == 1: os.system(cmd_bce_blast)

    # processing output file
    try:
        blastp_tce_f = open("blastp_bce_e2.txt", "r")
        lines = blastp_tce_f.read().splitlines()
        blastp_tce_f.close()

    except IOError:
        if report_stat == 1: print('blastp_bce_e2.txt not found')


    blast_bce_f = open("blastp_bce_e2.txt", "r")
    # processing output file
    if report_stat == 1: print("Processing blastp_bce_tce result.")
    lines = blast_bce_f.read().splitlines()
    total_line = len(lines)
    # get the lines where query start
    blastlinenums = []  # list of line where the query start
    i = 0
    for i in range(total_line):
        if re.match(r"^Query= (\S+)", lines[i]):
            blastlinenums.append(i)
    total_query = len(blastlinenums)
    if report_stat == 1: print("Total queries = " , total_query)
    blastlinenums.append(total_line)  # append EOL into blastlinenums
    if report_stat == 1: print("blastlinenums = ",blastlinenums, "\n")

    # process query by query
    qu_p = [''] * total_query  # list of result according to query
    for i in range(total_query):
        temp = []
        for l in range(blastlinenums[i], blastlinenums[i + 1]):
            temp.append(lines[l])
            if report_stat == 1: print("Total lines for query #", i+1 , " are " , len(temp))
        qu_p[i] = temp
        del temp


    ###only for debug to check for index in list 'qu_p'
    for i in range(total_query):
        if report_stat == 1: print("\n******Query " + str(i+1) + "******\n")
        for j in range(len(qu_p[i])):
            line=qu_p[i][j]
            if report_stat == 1 : print("qu_p[" , i,"][" , j, "] " , line)
    if report_stat == 1: print('\n\n\n')
    ## comment until this line

    # to make a holder dictionary for ir data
    # ir_output_dolol[tce_id_str] = [bce_list]  # bce_list = [bce_id, tce_ovlp_bce, tce_prct_ovlp]
    # tce_ovlp_bce = sbjt_start + "\t" + tce_ovlp_bce[:int(sbjt_start)-1] + lc(alligned_bce_seq) + tce_ovlp_bce[int(sbjt_end):] + "\t" + sbjt_end
    ir_output_dolol = {}
    for tce_id_str in tce_all_dos.keys():
        ir_output_dolol[tce_id_str] = []

    # process each query line by line
    for i in range(total_query):
        j = 0
        while j < (len(qu_p[i])):
            line = qu_p[i][j]  # this should be the first line of a  query
            while not re.match(r"^Query= (\S+)", line):
                j += 1
                line = qu_p[i][j]
            bce_id = re.findall(r"^Query= (\S+)", line)[0]
            if report_stat == 2: print("\nProcessing Query ID = " , bce_id)
            j += 1
            line = qu_p[i][j]
            while not re.match(r"^Length=(\d+)", line):
                j += 1
                line = qu_p[i][j]
            bce_len = re.findall(r"^Length=(\d+)", line)[0]
            j += 1
            line = qu_p[i][j]
            while not (re.match(r"^Sequences producing", line) or re.match(r"(.+)No hits found(.+)", line)):
                j += 1
                line = qu_p[i][j]

            if re.match(r"^Sequences producing", line):
                while not re.match(r"^> (\S+)", line):
                    j += 1
                    line = qu_p[i][j]
                while re.match(r"^> (\S+)", line):
                    tce_id = re.findall(r"^> (\S+)", line)[0]
                    j += 1
                    line = qu_p[i][j]
                    tce_length = re.findall(r"^Length=(\d+)", line)[0]
                    while not re.match(r"^ Identities = ((\d+)/(\d+)\s\((\d+)\%\))", line):
                        j += 1
                        line = qu_p[i][j]
                    x = re.findall(r"^ Identities = (\d+)/(\d+)\s\((\d+)\%\)", line)[0]
                    mismatch = int(x[1]) - int(x[0])
                    hit_subject_mlen = x[1]

                    # accepted
                    if int(hit_subject_mlen) >= 8: #andint(mismatch) <= 3:
                        j += 1
                        line = qu_p[i][j]
                        while not re.match(r"^Query \s\d+\s+(\S+)\s+\d+", line):
                            j += 1
                            line = qu_p[i][j]
                        alligned_bce_seq = re.findall(r"^Query \s\d+\s+(\S+)\s+\d+", line)[0]
                        alligned_bce_seq = alligned_bce_seq.lower()
                        j += 1
                        line = qu_p[i][j]
                        while not re.match(r"^Sbjct \s(\d+)\s+(\S+)\s+(\d+)", line):
                            j += 1
                            line = qu_p[i][j]
                        x = re.findall(r"^Sbjct \s(\d+)\s+(\S+)\s+(\d+)", line)[0]
                        sbjt_start = x[0]
                        sbjt_end = x[2]
                        tce_ovlp_bce = tce_all_dos[tce_id] # tce_seq ori
                        tce_ovlp_bce = tce_ovlp_bce[:int(sbjt_start)-1] + alligned_bce_seq + tce_ovlp_bce[int(sbjt_end):]

                        tce_prct_ovlp = round(float(hit_subject_mlen) / float(tce_length) * 100, 2)
                        #if bce_id == "43332" and tce_id == "153":
                            #print("tce_prct_ovlp:" + str(tce_prct_ovlp) + " hit_subject_mlen:" + hit_subject_mlen + " tce_length:" + tce_length)
                        #    exit()
                        bce_list = [bce_id, tce_ovlp_bce, tce_prct_ovlp]
                        ir_output_dolol[tce_id].append(bce_list)
                        if report_stat == 1: print("bce id:", bce_id, "alligned_bce_seq:", alligned_bce_seq, "tce id:", tce_id, "tce_ori:", tce_all_dos[tce_id], "bce_tce_ovlp:", tce_ovlp_bce, str(bce_list))
                        #if int(bce_id) == 429190 and int(tce_id) == 58225:
                            #print("bce id:", bce_id, "alligned_bce_seq:", alligned_bce_seq, "tce id:", tce_id, "tce_ori:", tce_all_dos[tce_id], "bce_tce_ovlp:", tce_ovlp_bce, str(bce_list))
                            #val = input("press enter")

                    # discarded query
                    if int(hit_subject_mlen) < 8 or int(mismatch) > 3:
                        if report_stat == 1: print("hit subject len <8 or mismatched >3")

                    while not (re.match(r"^> (\S+)", line) or re.match(r"^Lambda", line)):
                        j += 1
                        line = qu_p[i][j]

            # no hit
            if re.match(r"(.+)No hits found(.+)", line):
                if report_stat == 2: print("No hit found for Query ID=", bce_id)


            while not j == (len(qu_p[i])):
                j += 1

    if report_stat == 1: print("Done processing blast result of all queries.")
    if report_stat == 1: print("Proceed with the processing the H37Rv antigen that contain epitopes.")
    blast_bce_f.close()

    if cmd_stat == 1:
        filehandler5 = open('ir_output_dolol.obj', 'wb')  # to save memory to file so do not need to run cmd everytime
        pickle.dump(ir_output_dolol, filehandler5)
        filehandler5.close()
    if report_stat == 1: print("--- ir processing in:%s seconds ---" % (time.time() - startir_time))

filehandler5 = open('ir_output_dolol.obj', 'rb') #comment this first and uncomment it after run cmd
ir_output_dolol = pickle.load(filehandler5)
filehandler5.close()



### xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx SCRIPT 4e: HYPERCONSERVATION OF MTB PROTEINS xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ###
if hyp_swc:
    if report_stat >= 0: print("\n\n############## SCRIPT 4e: HYPERCONSERVATION OF MTB PROTEINS ##############\n")
    # This section is to process the imported result of hyperconservation MTB strains (172 complete sequence strains)
    # The result of hyperconservation cds processing (antigen.fasta, 173 strains) the blast o/p file (matched_h37rv_cds.faa) was compared with mtb H37Rv cds reference strains (GCF_000195955.2_ASM19595v2_cds_from_genomic.fna.out) to recheck the length matched
    # The cds file was converted to faa file to run blastp

    starthyp_time = time.time()

    strain_protein_list = os.listdir('273_mtb_strain_protein')  # make a list of mhc_pred_pickle files

    # populate tce_strainprotien_dol using strain_fn loop, blast using tce_all_ep_file, need to check the blast output pattern matching
    tce_strainprotein_dol = {}  # tce_strainprotein_dol[tce] = [#count of matched epitope] list follow strain_fn ordering
    for tce in tce_all_dos.keys():
        tce_strainprotein_dol[tce] = [0] * len(strain_protein_list)

    for idx in range(len(strain_protein_list)):
        if not os.path.exists('db_strain'):
            os.makedirs('db_strain')
        strain_fn = '273_mtb_strain_protein/' + strain_protein_list[idx]
        strain_dir = 'db_strain/' + strain_protein_list[idx][:-4]  # 4 for (.faa) extension file
        cmd1 = "makeblastdb -in " + strain_fn + " -out " + strain_dir + " -parse_seqids -dbtype prot"
        if cmd_stat == 1: os.system(cmd1)

        cmd_strain_dir_blast = "blastp -task blastp-short -query " + tce_all_ep_file + " -db " + strain_dir + " -out blastp_" + strain_protein_list[idx][:-4] + ".txt"
        if cmd_stat == 1: os.system(cmd_strain_dir_blast)
        if rm_stat == 1: shutil.rmtree('db_strain')  # to delete database after use

        # processing T cell blastp with hyperconservation cds o/p file
        strain_dir_blastp_op = "blastp_" + strain_protein_list[idx][:-4] + ".txt"
        fx.file_val(strain_dir_blastp_op)
        strain_dir_op_f = open(strain_dir_blastp_op, "r")
        lines = strain_dir_op_f.read().splitlines()
        total_line = len(lines)
        strain_dir_op_f.close()
        os.remove(strain_dir_blastp_op)  # to remove the blastp file
        blastlinenums = []  # list of line where the query start
        i = 0
        for i in range(total_line):
            if re.match(r"^Query= (\S+)", lines[i]):
                blastlinenums.append(i)
        total_query = len(blastlinenums)
        if report_stat == 1: print("Total queries = ", total_query)
        blastlinenums.append(total_line)  # append EOL into blastlinenums
        if report_stat == 1: print("blastlinenums = ", blastlinenums, "\n")
        # get the lines where query start
        if report_stat >= 0: print("Processing strain_dir blastp result.")

        # process query by query
        qu_p = [''] * total_query  # list of result according to query blastp_ort_tce
        for i in range(total_query):
            temp = []
            for l in range(blastlinenums[i], blastlinenums[i + 1]):
                temp.append(lines[l])
            if report_stat == 2: print("Total lines for query blastp_ort_tce #", i + 1, " are ", len(temp))
            qu_p[i] = temp
            del temp

        i = 0
        while i < total_query:
            j = 0
            while j < (len(qu_p[i])):
                line = qu_p[i][j]  # this should be the first line of a  query
                while not re.match(r"^Query= (\S+)", line):
                    j += 1
                    line = qu_p[i][j]
                q_id = re.findall(r"^Query= (\S+)", line)[0]
                if report_stat == 1: print("\nProcessing Query ID = ", q_id)
                j += 1
                line = qu_p[i][j]
                while not re.match(r"^Length=(\d+)", line):
                    j += 1
                    line = qu_p[i][j]
                q_len = re.findall(r"^Length=(\d+)", line)[0]
                j += 1
                line = qu_p[i][j]
                while not (re.match(r"^Sequences producing", line) or re.match(r"(.+)No hits found(.+)", line)):
                    j += 1
                    line = qu_p[i][j]

                if re.match(r"^Sequences producing", line):
                    while not re.match(r"^>(\S+)", line):
                        if report_stat == 2: print("i:", i, "j:,", j, "line:", line)
                        j += 1
                        line = qu_p[i][j]
                    while re.match(r"^>(\S+)", line):
                        hit_id = re.findall(r"^>(\S+)", line)[0]
                        j += 1
                        line = qu_p[i][j]
                        while not re.match(r"^ Identities = ((\d+)/(\d+)\s\((\d+)\%\))", line):
                            j += 1
                            line = qu_p[i][j]
                        x = re.findall(r"^ Identities = (\d+)/(\d+)\s\((\d+)\%\)", line)[0]
                        mismatch = int(x[1]) - int(x[0])
                        hit_subject_mlen = x[1]

                        # accepted
                        if int(hit_subject_mlen) >= int(q_len) and int(mismatch) <= 3:
                            j += 1
                            line = qu_p[i][j]
                            while not re.match(r"^Query\s+\d+\s+(\S+)\s+\d+", line):
                                j += 1
                                line = qu_p[i][j]
                                qu_seq = re.findall(r"^Query\s+\d+\s+(\S+)\s+\d+", line)[0]
                            while not re.match(r"^Sbjct\s+(\d+)\s+\S+\s+(\d+)", line):
                                j += 1
                                line = qu_p[i][j]
                            if report_stat == 2: print("line:", line)
                            if report_stat == 2: print("list:", re.search(r"^Sbjct\s+(\d+)\s+\S+\s+(\d+)", line))
                            start = re.search(r"^Sbjct\s+(\d+)\s+(\S+)\s+(\d+)", line).group(1)
                            sub_seq = re.search(r"^Sbjct\s+(\d+)\s+(\S+)\s+(\d+)", line).group(2)
                            end = re.search(r"^Sbjct\s+(\d+)\s+(\S+)\s+(\d+)", line).group(3)
                            if report_stat == 1: print("start:", start)
                            if report_stat == 1: print("end:", end)

                            # if q_id == "189575" and hit_id == "NP_217478.1":
                            # print("")
                            tce_strainprotein_dol[q_id][idx] += 1  # to insert tce info
                            print("idx=" + str(idx), "tce_strainprotein_dol[" + q_id + "]=",
                                  tce_strainprotein_dol[q_id])

                        # discarded query
                        # if int(hit_subject_mlen) < int(q_len) or int(mismatch) > 3:
                        #    if report_stat == 2: print("")
                        while not (re.match(r"^>(\S+)", line) or re.match(r"^Lambda", line)):  # to skip line
                            j += 1
                            line = qu_p[i][j]

                # no hit
                if re.match(r"(.+)No hits found(.+)", line):
                    if report_stat == 2: print("")
                    if report_stat == 2: print("No hit found for Query ID=", q_id)

                while not j == (len(qu_p[i])):
                    j += 1

            i += 1

    if cmd_stat == 1:
        filehandler6 = open('tce_strainprotein_dol.obj', 'wb')  # to save memory to file so do not need to run cmd everytime
        pickle.dump(tce_strainprotein_dol, filehandler6)
        filehandler6.close()

    if report_stat == 1: print("--- ort processing in:%s seconds ---" % (time.time() - starthyp_time))

filehandler6 = open('tce_strainprotein_dol.obj', 'rb')  # comment this first and uncomment it after run cmd
tce_strainprotein_dol = pickle.load(filehandler6)
filehandler6.close()

### xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx SCRIPT 4f: MICROBIOTA AND HUMAN PROTEIN (SELF ANTIGENIC UNIVERSAL-SAU) xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ###
if sau_swc:
    if report_stat >= 0: print("\n\n############## SCRIPT 4f: MICROBIOTA AND HUMAN PROTEIN (SELF ANTIGENIC UNIVERSAL-SAU) ##############\n")
    # This section is for blast the epitopes with human proteins and human's microbiome protein.
    # sau = self antigenic universal, hs = "homo sapiens", biom = microbiota
    # sau_output_dodol ={'tce_id':{'mtb_tce_prtn':(mtb_prtn_ids), 'col_uniq_mtb':none,'hs_tce_prtn':(hs_prtn_ids), 'col_mtb_hs':none,
    # 					            'biom_tce_prtn':(biome_prtn_ids),'col_mtb_biom':none,'col_mtb_hs_biom':none}}

    # Initialise sau dictionary
    # we will use the same dict from above for tce i.e. tce_all_dos[tce_id_str] = tce_seq_str

    startsau_time = time.time()
    sau_output_dodol = OrderedDict()
    for tce_id_str, tce_seq_str in tce_all_dos.items():
        mtb_prtn_ids = []
        mtb_tce_rv_prtn = []
        hs_prtn_ids = []
        biome_prtn_ids = []
        sau_output_dodol[tce_id_str] = {'mtb_tce_prtn' : mtb_prtn_ids, 'col_uniq_mtb': None, 'mtb_tce_rv_prtn': mtb_tce_rv_prtn, 'hs_tce_prtn': hs_prtn_ids, 'col_mtb_hs': None, 'biom_tce_prtn': biome_prtn_ids, 'col_mtb_biom': None, 'col_mtb_hs_biom': None}
        if report_stat == 1: print("sau_output_dodol:[", tce_id_str, "]", sau_output_dodol[tce_id_str])

    # make database from fasta format file
    cmd1 = "makeblastdb -in " + mtb_ncbi + " -out " + database_mtb_ncbi + " -parse_seqids -dbtype prot"
    if cmd_stat == 1: os.system(cmd1)

    # make mycobrowser database from fasta format file
    cmd2 = "makeblastdb -in " + mtb_mycodb + " -out " + database_mtb_rv + " -parse_seqids -dbtype prot"
    if cmd_stat == 1: os.system(cmd2)

    cmd3 = "makeblastdb -in " + hs + " -out " + database_hs + " -parse_seqids -dbtype prot"
    if cmd_stat == 1: os.system(cmd3)

    cmd4 = "makeblastdb -in " + biom + " -out " + database_biom + " -parse_seqids -dbtype prot"
    if cmd_stat == 1: os.system(cmd4)
    print(cmd4)
    # makeblastdb -in human_microbiome_project_sequence_1.faa -out db/human_microbiome_project_sequence_1.faa -parse_seqids -dbtype prot


    def prss_blastp_opf(tce_id_str, tce_seq_str, type, database, sau_output_dodol):
        if report_stat == 1: print("Processing blast result.")

        f = open("tce_single.txt", "w")
        f.write('>%s\n%s\n' % (tce_id_str, tce_seq_str))
        f.close()
        cmd4 = "blastp -num_threads 6 -task blastp-short -query tce_single.txt -db " + database + " -out blastp_tce_op.txt"
        # run blast
        if cmd_stat == 1: os.system(cmd4)
        print("cmd4:", cmd4)


        # processing output file
        try:
            blastp_tce_f = open("blastp_tce_op.txt", "r")
            lines = blastp_tce_f.read().splitlines()
            blastp_tce_f.close()

        except IOError:
            if report_stat == 1: print('file not found')

        total_line = len(lines)
        # get the lines where query start
        blastlinenums = []  # list of line where the query start
        i = 0

        for i in range(total_line):
            if re.match(r"^Query=\s*([a-zA-Z0-9_.]+)", lines[i]):
                blastlinenums.append(i)
        total_query = len(blastlinenums)
        if report_stat == 2: print("Total queries = " , total_query)
        blastlinenums.append(total_line)  # append EOL into blastlinenums
        if report_stat == 2: print("blastlinenums = ",blastlinenums, "\n")

        # process query by query
        qu_p = [''] * total_query  # list of result according to query
        for i in range(total_query):
            temp = []
            for l in range(blastlinenums[i], blastlinenums[i + 1]):
                temp.append(lines[l])
                if report_stat == 2: print("Total lines for query #", i+1 , " are " , len(temp))
            qu_p[i] = temp
            del temp

        ###only for debug to check for index in list 'qu_p'
        for i in range(total_query):
            if debug_stat == 1: print("\n******Query " + str(i + 1) + "******\n")
            for j in range(len(qu_p[i])):
                line = qu_p[i][j]
                if debug_stat == 1: print("qu_p[", i, "][", j, "] ", line)
        if debug_stat == 1: print('\n\n\n')
        ## comment until this line

        if type == "biom":
            print("")
        accepted_hitid = []  # list of hit_id that is not discard
        # process each query line by line
        for i in range(total_query):
            j = 0
            while j < (len(qu_p[i])):
                line = qu_p[i][j]  # this should be the first line of a  query
                while not re.match(r"^Query=\s*([a-zA-Z0-9_.]+)", line):
                    j += 1
                    line = qu_p[i][j]
                bce_id = re.findall(r"^Query=\s*([a-zA-Z0-9_.]+)", line)[0]
                if report_stat == 2: print("\nProcessing Query ID = " ,bce_id)
                j += 1
                line = qu_p[i][j]
                while not re.match(r"^Length=(\d+)", line):
                    j += 1
                    line = qu_p[i][j]
                bce_len = re.findall(r"^Length=(\d+)", line)[0]
                j += 1
                line = qu_p[i][j]
                while not (re.match(r"^Sequences producing", line) or re.match(r"(.+)No hits found(.+)", line)):
                    j += 1
                    line = qu_p[i][j]


                if re.match(r"^Sequences producing", line):
                    while not re.match(r"^>\s*([a-zA-Z0-9_.]+)", line):
                        j += 1
                        line = qu_p[i][j]

                    while re.match(r"^>\s*([a-zA-Z0-9_.]+)", line):
                        hit_id = re.findall(r"^>\s*([a-zA-Z0-9_.]+)", line)[0]
                        j += 1
                        line = qu_p[i][j]
                        while not re.match(r"^ Identities = ((\d+)/(\d+)\s\((\d+)\%\))", line):
                            j += 1
                            line = qu_p[i][j]
                        x = re.findall(r"^ Identities = (\d+)/(\d+)\s\((\d+)\%\)", line)[0] #x == list of list
                        mismatch = int(x[1]) - int(x[0])
                        hit_subject_len = x[1]

                        # accepted
                        if int(hit_subject_len) >= int(bce_len) and int(mismatch) <= 3:
                            accepted_hitid.append(hit_id)
                            j += 1
                            line = qu_p[i][j]
                            while not re.match(r"^Query \s\d+\s+(\S+)\s+\d+", line):
                                j += 1
                                line = qu_p[i][j]
                            alligned_bce_seq = re.findall(r"^Query \s\d+\s+(\S+)\s+\d+", line)[0]
                            while not re.match(r"^Sbjct \s(\d+)\s+(\S+)\s+(\d+)", line):
                                j += 1
                                line = qu_p[i][j]
                            x = re.findall(r"^Sbjct \s(\d+)\s+(\S+)\s+(\d+)", line)[0]
                            # print("Done ")

                        # discarded query
                        if int(hit_subject_len) < int(bce_len) or int(mismatch) > 3:
                            if report_stat == 2: print("")
                        while not (re.match(r"^>\s*([a-zA-Z0-9_.]+)", line) or re.match(r"^Lambda", line)):
                            j += 1
                            line = qu_p[i][j]

                # no hit
                if re.match(r"(.+)No hits found(.+)", line):
                    if report_stat == 1: print("No hit found for Query ID=" + bce_id)

                while not j == (len(qu_p[i])):
                    j += 1

        if debug_stat == 1: print("only for debug to check for accepted hitid (list)")
        for i in range(len(accepted_hitid)):
            if debug_stat == 1: print("accepted_hitid[" , i,"][" , "] " ,)
        if debug_stat == 1: print("END DEBUGGING")

        # in this part, the accepted hit_id will be use as the value in the corresponding keys
        accepted_hitid = list(set(accepted_hitid))  # remove redundant hit_id   #set return the list into tuple thus, have to change data type to list again by list(....)
        accepted_hitid.sort()
        if type == "mtb_np":
            sau_output_dodol[tce_id_str]['mtb_tce_prtn'] = accepted_hitid
        if type == "mtb_rv":
            sau_output_dodol[tce_id_str]['mtb_tce_rv_prtn'] = accepted_hitid
        if type == "hs":
            sau_output_dodol[tce_id_str]['hs_tce_prtn'] = accepted_hitid
        if type == "biom":
            sau_output_dodol[tce_id_str]['biom_tce_prtn'] = accepted_hitid


    if report_stat == 0: print("SAU table: " )
    for tce_id_str, tce_seq_str in tce_all_dos.items():

        prss_blastp_opf(tce_id_str, tce_seq_str, "mtb_np", database_mtb_ncbi, sau_output_dodol)
        prss_blastp_opf(tce_id_str, tce_seq_str, "mtb_rv", database_mtb_rv, sau_output_dodol)
        prss_blastp_opf(tce_id_str, tce_seq_str, "hs", database_hs, sau_output_dodol)
        prss_blastp_opf(tce_id_str, tce_seq_str, "biom", database_biom, sau_output_dodol)
        found_in_mtb = found_in_hs = found_in_biom = False

        if len(sau_output_dodol[tce_id_str]['mtb_tce_prtn']) > 0:  # this condition for checking the logic for column below
            found_in_mtb = True

        if len(sau_output_dodol[tce_id_str]['hs_tce_prtn']) > 0:
            found_in_hs = True

        if len(sau_output_dodol[tce_id_str]['biom_tce_prtn']) > 0:
            found_in_biom = True

        if found_in_mtb == True and found_in_hs == False and found_in_biom == False: # 1st condition = only found in mtb, tce is unique to mtb only
            sau_output_dodol[tce_id_str]['col_uniq_mtb'] = "Yes"
            sau_output_dodol[tce_id_str]['col_mtb_hs'] = "No"
            sau_output_dodol[tce_id_str]['col_mtb_biom'] = "No"
            sau_output_dodol[tce_id_str]['col_mtb_hs_biom'] = "No"

        if found_in_mtb == True and found_in_hs == True and found_in_biom == False:  # 2nd condition = only found in mtb and hs , tce no longer unique to mtb
            sau_output_dodol[tce_id_str]['col_uniq_mtb'] = "No"
            sau_output_dodol[tce_id_str]['col_mtb_hs'] = "Yes"
            sau_output_dodol[tce_id_str]['col_mtb_biom'] = "No"
            sau_output_dodol[tce_id_str]['col_mtb_hs_biom'] = "No"

        if found_in_mtb == True and found_in_hs == False and found_in_biom == True: #3rd condition = only found in mtb and biom, tce no longer unique to mtb
            sau_output_dodol[tce_id_str]['col_uniq_mtb'] = "No"
            sau_output_dodol[tce_id_str]['col_mtb_hs'] = "No"
            sau_output_dodol[tce_id_str]['col_mtb_biom'] = "Yes"
            sau_output_dodol[tce_id_str]['col_mtb_hs_biom'] = "No"

        if found_in_mtb == True and found_in_hs == True and found_in_biom == True: #4th condition = found in mtb, hs and biom, tce no longer unique to mtb
            sau_output_dodol[tce_id_str]['col_uniq_mtb'] = "No"
            sau_output_dodol[tce_id_str]['col_mtb_hs'] = "No"
            sau_output_dodol[tce_id_str]['col_mtb_biom'] = "No"
            sau_output_dodol[tce_id_str]['col_mtb_hs_biom'] = "Yes"

        if found_in_mtb == False and found_in_hs == False and found_in_biom == False: #5th condition, tce is not found in mtb, or hs or biom.
            sau_output_dodol[tce_id_str]['col_uniq_mtb'] = None
            sau_output_dodol[tce_id_str]['col_mtb_hs'] = None
            sau_output_dodol[tce_id_str]['col_mtb_biom'] = None
            sau_output_dodol[tce_id_str]['col_mtb_hs_biom'] = None

        if report_stat ==1: print("sau_output_dodol:[" + tce_id_str + "]: " + str(sau_output_dodol[tce_id_str]))
        if report_stat ==1: print("tce_id:" + tce_id_str + "\t mtb_prtn_id:" + str(sau_output_dodol[tce_id_str]['mtb_tce_prtn']) +
                                  "\t mtb_prtn_rv_id:" + str(sau_output_dodol[tce_id_str]['mtb_tce_rv_prtn']) + "\t hs_tce_prtn:" +
                                  str(sau_output_dodol[tce_id_str]['hs_tce_prtn'])  + "\t biom_tce_prtn:" + str(sau_output_dodol[tce_id_str]['biom_tce_prtn']) +
                                  "\t col_uniq_mtb:" + str(sau_output_dodol[tce_id_str]['col_uniq_mtb']) + "\t col_mtb_hs:" + str(sau_output_dodol[tce_id_str]['col_mtb_hs']) +
                                  "\t col_mtb_biom:" + str(sau_output_dodol[tce_id_str]['col_mtb_biom']) + "\t col_mtb_hs_biom:" + str(sau_output_dodol[tce_id_str]['col_mtb_hs_biom']))
    if report_stat == 0: print("Done processing blast result of all tce_id for database mtb_ncbi, mtb_myco_table, hs and biom with all conditions.")

    if cmd_stat == 1:
        filehandler7 = open('sau_output_dodol.obj', 'wb')  # to save memory to file so do not need to run cmd everytime
        pickle.dump(sau_output_dodol, filehandler7)
        filehandler7.close()

    if report_stat == 1: print("--- sau processing in:%s seconds ---" % (time.time() - startsau_time))
    if rm_stat == 1: shutil.rmtree(database_mtb_ncbi)  # delete the database after used to save space
    if rm_stat == 1: shutil.rmtree(database_mtb_rv)
    if rm_stat == 1: shutil.rmtree(database_hs)
    if rm_stat == 1: shutil.rmtree(database_biom)

filehandler7 = open('sau_output_dodol.obj', 'rb') #comment this first and uncomment it after run cmd
sau_output_dodol = pickle.load(filehandler7)
filehandler7.close()

### xxxxxxxxxxxx SCRIPT 4g : PRINTING THE COMBINATION OF MHC1, MHC2, POPULATION COVERAGE, IR, ORTHOLOGOUS GENES AND SAU xxxxxxxxxxxx #
if grd_smry_tble_swc:

    if report_stat >= 0: print("\n\n############ SCRIPT 4g: PRINTING THE COMBINATION OF MHC1, MHC2, POPULATION COVERAGE, IR, ORTHOLOGOUS GENES AND SAU ##############\n")
    # to produce 2 table i.e. final summary table and debug table, based on supertypes and allele counted

    if report_stat == 0: print("summary table:")
    with open("grand_summary_table.txt", 'w') as tce_output_file_handler:
        sau_header = "Mtb_NP_Prtn_id\tMtb_NP_Prtn_count\tMtb_Rv_Prtn_id\tMtb_Rv_Prtn_count\tHs_Protein_id\tHS_Protein_count\tBiome_Protein_id\tBiom_Protein_count\tUnique_to_Mtb\tFound_in_Mtb-Hs\tFound_in_Mtb-Biom\tFound_in_Mtb-Hs-Biom\t"
        mhc_header = "Locus\tLocus_count\tAlleles\tAlleles_count\t"
        popcov_header = "World\tAsian\tSoutheast_Asia\tMalaysia\t"
        ir_header = "IR_matched_sequence\tIR_BCE_TCE_%matched\tIR_BCE_counter\t"
        hyp_header = "Hyp_strain_cov\tHyp_min,max,avg,sd\t"
        heg_header = "HEG_Rv_id\tHEG_Rv_loc\tHEG_Rv_count\tHEG_species/condition\tHEG_expr_score\t"
        file_header = "TCE_id\tTCE_sequence\t" + sau_header + mhc_header + popcov_header + ir_header + hyp_header + heg_header
        #\tTCE_matched_HEG\tTCE_location_HEG\tTCE_location_HEG_count\tTCE_matched_orthologous
        tce_output_file_handler.write('%s\n' % (file_header))
        print(file_header)

        master_tce_list = [] # bcoz need to run subset
        if mhc_all_swc:
            master_tce_list = tce_all_dos.keys()
        elif mhc_1x_swc or mhc_mutual_swc:
            master_tce_list = tce_mhc1_dos.keys()
        elif mhc_2x_swc:
            master_tce_list = tce_mhc2_dos.keys()
        else:
            master_tce_list = tce_all_dos.keys()

        if tce_pickle_swc:
            # to skip the tce_id list len <15 in mhc2
            filehandler_mhc_skip = open('tce_mhc2_skip_dict.obj', 'rb')  # comment this first and uncomment it after run cmd
            tce_mhc2_skip_dict = pickle.load(filehandler_mhc_skip)
            filehandler_mhc_skip.close()

        for tce_id in master_tce_list: # replace with tce_mhc1_dos (case #1-10) / tce_mhc2_dos (case #11) for selective printing
            # not all tce_id in tce_mhc_op_dol is present in tce_all_dos as we discard some tce_id due to unmeet criteria eg length <8
            # sau_output_dodol ={'tce_id':{'mtb_tce_prtn':(mtb_prtn_ids), 'col_uniq_mtb':none,'hs_tce_prtn':(hs_prtn_ids), 'col_mtb_hs':none,
            # 					            'biom_tce_prtn':(biome_prtn_ids),'col_mtb_biom':none,'col_mtb_hs_biom':none}}
            if tce_pickle_swc:
                if tce_id in tce_mhc2_skip_dict.keys():
                    continue

            sau_output_str = ""
            if tce_id in sau_output_dodol.keys():
                tce_id_sau = tce_id # need to reassign as tce_id_sau does not have value for "if" loop

                if sau_output_dodol[tce_id_sau]['col_uniq_mtb'] == None: sau_output_dodol[tce_id_sau]['col_uniq_mtb'] = "Nil"
                if sau_output_dodol[tce_id_sau]['col_mtb_hs'] == None: sau_output_dodol[tce_id_sau]['col_mtb_hs'] = "Nil"
                if sau_output_dodol[tce_id_sau]['col_mtb_biom'] == None: sau_output_dodol[tce_id_sau]['col_mtb_biom'] = "Nil"
                if sau_output_dodol[tce_id_sau]['col_mtb_hs_biom'] == None: sau_output_dodol[tce_id_sau]['col_mtb_hs_biom'] = "Nil"

                print_list = ["", "", "", "", "", "", "", "", "", "", "", ""]
                print_list[8] = sau_output_dodol[tce_id_sau]['col_uniq_mtb']
                print_list[9] = sau_output_dodol[tce_id_sau]['col_mtb_hs']
                print_list[10] = sau_output_dodol[tce_id_sau]['col_mtb_biom']
                print_list[11] = sau_output_dodol[tce_id_sau]['col_mtb_hs_biom']

                if len(sau_output_dodol[tce_id_sau]['mtb_tce_prtn']) > 0: # for condition case #1
                    print_list[0] = ','.join(str(mtb_tce) for mtb_tce in sau_output_dodol[tce_id_sau]['mtb_tce_prtn'])
                    print_list[1] = str(len(sau_output_dodol[tce_id_sau]['mtb_tce_prtn']))

                    if len(sau_output_dodol[tce_id_sau]['hs_tce_prtn']) > 0: # for condition case #2 and 4
                        print_list[4] = ','.join(str(hs_tce) for hs_tce in sau_output_dodol[tce_id_sau]['hs_tce_prtn'])
                        print_list[5] = str(len(sau_output_dodol[tce_id_sau]['hs_tce_prtn']))
                    else:
                        print_list[4] = "Nil"
                        print_list[5] = "0"

                    if len(sau_output_dodol[tce_id_sau]['biom_tce_prtn']) > 0:  # for condition case #3  and 4
                        print_list[6] = ','.join(str(biom_tce) for biom_tce in sau_output_dodol[tce_id_sau]['biom_tce_prtn'])
                        print_list[7] = str(len(sau_output_dodol[tce_id_sau]['biom_tce_prtn']))
                    else:
                        print_list[6] = "Nil"
                        print_list[7] = "0"
                else:  # for condition case #5
                    print_list[0:8] = ["Nil", "0", "Nil", "0", "Nil", "0", "Nil", "0"]

                if len(sau_output_dodol[tce_id_sau]['mtb_tce_rv_prtn']) > 0: # for condition case #1
                    print_list[2] = ','.join(str(mtb_tce) for mtb_tce in sau_output_dodol[tce_id_sau]['mtb_tce_rv_prtn'])
                    print_list[3] = str(len(sau_output_dodol[tce_id_sau]['mtb_tce_rv_prtn']))
                else:
                    print_list[2:4] = ["Nil", "0"]

                sau_output_str = '\t'.join(str(sau_list) for sau_list in print_list)
                print("sau_info:", sau_output_str)

            else:
                sau_output_str = "Nil\t0\tNil\t0\tNil\t0\tNil\t0\tNil\tNil\tNil\tNil"

            mhc_op_str = ""
            if tce_id in tce_mhc_op_dol.keys():
                mhc_info_lol = tce_mhc_op_dol[tce_id]
                allele_set = mhc_info_lol[0]
                debug_info_lol = mhc_info_lol[1]
                sum_table_list = [set(), 0, allele_set]  # [(locus), total_allele, [list_of_allele_info]]
                for allele in allele_set:
                    sum_table_list[1] += 1
                    if (re.findall(r"^HLA-A(\S+)", allele)):
                        sum_table_list[0].add("A")
                    elif (re.findall(r"^HLA-B(\S+)", allele)):
                        sum_table_list[0].add("B")
                    elif (re.findall(r"^HLA-C(\S+)", allele)):
                        sum_table_list[0].add("C")
                    elif (re.findall(r"^HLA-E(\S+)", allele)):
                        sum_table_list[0].add("E")
                    elif (re.findall(r"^HLA-DP(\S+)", allele)):
                        sum_table_list[0].add("DP")
                    elif (re.findall(r"^HLA-DQ(\S+)", allele)):
                        sum_table_list[0].add("DQ")
                    elif (re.findall(r"^HLA-DR(\S+)", allele)):
                        sum_table_list[0].add("DR")
                sum_table_list_locus = sorted(sum_table_list[0])
                sum_table_list_allele = sorted(sum_table_list[2])

                mhc_op_str = ','.join(str(s) for s in sum_table_list_locus) + "\t" + str(len(sum_table_list_locus)) + "\t" + ','.join(str(s) for s in sum_table_list_allele) + "\t" + str(sum_table_list[1])
                if report_stat == 1: print("tce_id\tlocus\ttotal_locus\ttotal_allele\t[list_of_allele_info]")

                if len(debug_info_lol) > 0:
                    for debug_info in debug_info_lol:
                        allele, tce_strt, tce_end, tce_lngth, allele_seq, pred_prcntile_rank = debug_info
                        debug_info_str = allele + "\t" + tce_strt + "\t" + tce_end + "\t" + tce_lngth + "\t" + allele_seq + "\t" + pred_prcntile_rank
                else:
                    debug_info_str = "" + allele + "\tNil\tNil\tNil\tNil\tNil"
                if report_stat == 2: print("mhc_debug_info_str:", debug_info_str)
            else:
                mhc_op_str = "Nil\t0\tNil\t0"  # for tce_id with no predict binding result

            if report_stat == 1: print("mhc_op_str:", mhc_op_str)

            # popcov_output_dod = {tce_id: {Alleles:set(), world:%, Asian:%, Southeast Asia:%, Malaysia:%}}
            popcov_info_str = ""
            if tce_id in popcov_output_dod.keys():
                # allele_set = popcov_output_dod[tce_id]["Alleles"] #do not need as it already printed in mhc section
                world = popcov_output_dod[tce_id]["World"]
                if world == None: world = "0.0000"
                asian = popcov_output_dod[tce_id]["Asian"]
                if asian == None: asian = "0.0000"
                southeast_asia = popcov_output_dod[tce_id]["Southeast Asia"]
                if southeast_asia == None: southeast_asia = "0.0000"
                malaysia = popcov_output_dod[tce_id]["Malaysia"]
                if malaysia == None: malaysia = "0.0000"
                popcov_info_str = world + "\t" + asian + "\t" + southeast_asia + "\t" + malaysia
            else:
                popcov_info_str = "0.0000\t0.0000\t0.0000\t0.0000"

            if report_stat == 1: print("popcov_info:", popcov_info_str)


            # ir_output_dolol[tce_id_str] = [bce_list]  # bce_list = [bce_id, tce_ovlp_bce, tce_prct_ovlp]
            # tce_ovlp_bce = sbjt_start + "\t" + tce_ovlp_bce[:int(sbjt_start)-1] + alligned_bce_seq + tce_ovlp_bce[int(sbjt_end):] + "\t" + sbjt_end
            fx.file_val(bce_ep)
            bce_ep_f = open(bce_ep, "r")
            bce_ep_lines = bce_ep_f.read().splitlines()
            bce_id_str = ""
            bce_seq_str = ""
            bce_all_dos = OrderedDict()  # to make a dictionary for bce_seq for population coverage input file
            for line in range(len(bce_ep_lines)):
                if re.match(r"^>(\d+)", bce_ep_lines[line]):
                    bce_id_str = re.findall(r"^>(\d+)", bce_ep_lines[line])[0]
                if not re.match(r"^>(\d+)", bce_ep_lines[line]):
                    bce_seq_str = bce_ep_lines[line].rstrip()
                    bce_all_dos[bce_id_str] = bce_seq_str
            if report_stat == 1: print("bce_all_dos:", bce_all_dos)

            bce_list_str = ""
            bce_lol = ir_output_dolol[tce_id]
            if len(bce_lol) > 0:
                ir_sum_table = [[], [], 0]  # ir_sum_table = [[tce_ovlp_bce], [ir_info_str], bce_counter]

                for bce_list in bce_lol:
                    bce_id, tce_ovlp_bce, tce_prct_ovlp = bce_list
                    lc_count = sum(map(str.islower, tce_ovlp_bce))
                    unmatch_bce_seq_len = str(len(bce_all_dos[bce_id]) - lc_count)
                    prcntge_bce = round(float(lc_count / len(bce_all_dos[bce_id]) * 100), 2)
                    unmatch_tce_seq_len = str(len(tce_ovlp_bce) - lc_count) #BCE_id:  match_lc_count (-unmatched_bce_in_tce) / total_tce (% bce_match_in_tce)
                    ir_info_str = "bceid_" + bce_id + ":" + str(lc_count) + "(-" + unmatch_bce_seq_len + ")/" + str(len(bce_all_dos[bce_id])) + "(" + str(prcntge_bce) + "%)" + "<=>tceid_" + tce_id + ":" + str(lc_count) + "(-" + unmatch_tce_seq_len + ")/" + str(len(tce_ovlp_bce)) + "(" + str(tce_prct_ovlp) + "%)"
                    if prcntge_bce >= 80.00:
                        ir_sum_table[0].append(tce_ovlp_bce)
                        ir_sum_table[1].append(ir_info_str)
                        ir_sum_table[2] += 1

                    #if bce_id == "429190" and tce_id == "58225":
                        #print("bce_id, tce_ovlp_bce, tce_prct_ovlp:" , bce_id, tce_ovlp_bce, tce_prct_ovlp)
                if ir_sum_table[2] > 0:
                    bce_list_str = ','.join(str(s) for s in ir_sum_table[0]) + "\t" + ','.join(str(s) for s in ir_sum_table[1]) + "\t" + str(ir_sum_table[2])
                else:
                    bce_list_str = "Nil\tNil\t0"
            else:
                bce_list_str = "Nil\tNil\t0"

            print("tce_id: " + tce_id + " ir_info:" + bce_list_str)

            #if bce_id == "429190" and tce_id == "58225":
                #print("tce_id: " + tce_id + " ir_info:" + bce_list_str)
                #exit()

            # hyperconservation
            # Structure: tce_strainprotein_dol[tce] = [#count of matched epitope] list follow strain_fn ordering

            hyp_op_str = ""
            if tce_id in tce_strainprotein_dol.keys():
                tce_id_hyp = tce_id
                tce_strain_list = tce_strainprotein_dol[tce_id_hyp]
                percent = round((len(tce_strain_list) - tce_strain_list.count(0)) / len(tce_strain_list), 4)  # (total list-list with value=0)/total list*100 = float in 2 decimal place
                if percent > 0.0000:
                    minimum = min(tce_strain_list)
                    maximum = max(tce_strain_list)
                    avg = round(sum(tce_strain_list) / len(tce_strain_list), 2)
                    sd = round(statistics.stdev(tce_strain_list), 2)
                    hyp_op_str = str(percent) + "\t" + str(minimum) + "," + str(maximum) + "," + str(avg) + "," + str(sd)
                    print("tce_id_hyp:", tce_id_hyp, "hyp_op_str:", hyp_op_str)

                else:
                    hyp_op_str = "0.0000\tNil"


            # Structure: mtb_smry_dolol = {Np_id:[seq_normal, seq_coded,[TCE_id, qu_seq,sub_seq, st,sp], [BCE_id, qu_seq, sub_seq, st,sp]]}
            # Structure: ort_smry_dolol = {TCE_id:[tce_qu_seq,[[prtn_id, tce_loc]], [[prtn_id, bce_base_ovlp]], tce_bce_ovlp_ratio }
            heg_op_str = ""
            if tce_id in tce_smry_dolol.keys():
                tce_id_heg = tce_id
                info_heg = tce_smry_dolol[tce_id_heg]
                heg_table_list = [[], [], set(), [], 0]  # [[heg_id], [location], heg_count, [species,condition] ,score for HEG expression]
                tce_qu_seq, tce_loc_lol, bce_loc_lol = info_heg
                if len(tce_loc_lol) > 0:
                    for tce_loc_list in tce_loc_lol:
                        rv_id, tce_loc, tce_bce_ovlp_ratio = tce_loc_list
                        heg_table_list[0].append(rv_id)
                        heg_table_list[1].append(tce_loc)
                        heg_table_list[2].add(rv_id)
                        species = []
                        HEG_score = 0
                        if rv_id in mtb_gene_lab_rvid:
                            species.append("in-vitro")
                            HEG_score += 1
                        if rv_id in mtb_gene_mice_rvid:
                            species.append("in-vivo(mouse)")
                            HEG_score += 2
                        if rv_id in mtb_gene_human_rvid:
                            species.append("in-vivo(human)")
                            HEG_score += 3
                        heg_table_list[3] = species
                        heg_table_list[4] = HEG_score
                    heg_op_str = ','.join(str(s) for s in heg_table_list[0]) + "\t" + ','.join(str(s) for s in heg_table_list[1]) \
                                 + "\t" + str(len(heg_table_list[2])) + "\t" + ','.join(str(s) for s in heg_table_list[3]) + "\t" \
                                 + str(heg_table_list[4])
                else:
                    heg_op_str = "Nil\tNil\t0\tNil\t0"

                print("tce_id_heg:", tce_id_heg, "heg_op_str:", heg_op_str)

            else:
                heg_op_str = "Nil\tNil\t0\tNil\t0"

            print(tce_id_heg + "\t" + sau_output_str + "\t" + mhc_op_str + "\t" + popcov_info_str + "\t" + bce_list_str)
            tce_output_file_handler.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (tce_id, tce_all_dos[tce_id], sau_output_str, mhc_op_str, popcov_info_str, bce_list_str, hyp_op_str, heg_op_str))