import os
import re
import shutil

#file validation
def file_val(file_name):
    try:
        with open(file_name) as file:
            pass
    except IOError as e:
        print("Unable to open file", file_name)


#make new directory
def makenewdir(foldername):
    try:
        os.makedirs(foldername)
    except OSError:
        pass

def countRVheader(mtb_NCBI_filename):
    NCBI_ctr = 0
    mtb_NCBI_list = open(mtb_NCBI_filename).readlines()
    for i in mtb_NCBI_list:
        if i.startswith("Rv"):
            NCBI_ctr += 1
#print("mtb_NCBI = ", mtb_NCBI_filename)
#print("NCBI_ctr = ", NCBI_ctr)
    return NCBI_ctr

def filterppe():
    #make a list of ppe gene_id
    file_val('list_ppe.txt')
    with open('list_ppe.txt') as ppefile:
        line=ppefile.read().splitlines()
        ppeid=[]  
        for i in range(len(line)):
            while not re.match(r'^(Rv\S+).*',line[i]):
                i+=1
            ppeid.append(re.findall(r'^(Rv\S+).*',line[i])[0])
        

    #convertion of gene_id of PPE/Pe gene
    file_val('GCF_000195955.2_ASM19595v2_feature_table.txt')
    with open('GCF_000195955.2_ASM19595v2_feature_table.txt')as featurefile:
        line=featurefile.read().splitlines()
        geneid={}
        for i in range(0,len(line)):
            x=line[i].split('\t')
            if x[0]=='CDS':
                if x[16] in ppeid:
                    geneid[x[10]]=x[16]  #protein_id =x[10] as keys and Rv_id = x[16] as value in geneid dict
                        
    return geneid
