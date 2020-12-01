#!/usr/bin/env python

"""
This script has function to create itol annotation label file for deletion

input
1. Need VCF extract file  same format as use for statistical test for finding marker
    VCF extract DATA will be in format like this
        1st row is header ==> POS|Sample1|Sample2|Sample3|...|SampleN
        column POS contain information of chr and position merge into single string EX "1:100-300" mean chr1, start=100, end=300
        other column contain representation number of genotype that got from GT in VCF. 0 => homoREF, 1=>hetero, 2=>homoALT

output
1. phylib binary file (Noted right now we have only convert to binary in phylip format)
2. rename index file (will be use to rename to the original name after don't create tre with script from P'cladia)

"""

import pandas as pd
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, to_tree
import stat_utility
import os
import numpy as np

__author__ = "Worawich Phornsiricharoenphant"
__copyright__ = ""
__credits__ = ["Worawich Phornsiricharonphant"]
__license__ = "GPL-3.0"
__version__ = "1.0"
__maintainer__ = "Worawich Phornsiricharoenphant"
__email__ = "worawich.ph@gmail.com"
__status__ = "Development"

variant_extract_file = "/Volumes/10TBSeagateBackupPlus/NBT/TB_project/1170_delprofiler/wasna_paper/snp_sv_dendogram/snp_del_extract_ready.txt"
outputPath = "/Volumes/10TBSeagateBackupPlus/NBT/TB_project/1170_delprofiler/wasna_paper/snp_sv_dendogram"
homo_only = True
binary = True

def encrypt(string, length):
    ## add space every end of number of length
    return ' '.join(string[i:i+length] for i in range(0,len(string),length))

def addSpaceToString(input_str,round):
    for i in range(round):
        input_str = input_str + " "
    return input_str

def write_data_to_phylip(save_path, input_dataframe):
    count_sample=1
    phylip_file = os.path.join(save_path, "variant.phy")
    sampleID_index = os.path.join(save_path, "sample_rename_idx.txt")
    sample_rename_dict = dict()

    num_sample = len(input_dataframe.axes[0])
    num_variant = len(input_dataframe.axes[1])

    with open(phylip_file, "w") as phylip:
        # write header
        phylip.write(str(num_sample) + " " + str(num_variant))
        phylip.write("\n")

        # loop through datafame and try
        variant_info = input_dataframe.to_csv(header=None, index=False, sep='\t').split('\n')
        count = 0
        for index, row in df_t.iterrows():
            sampleName = index
            new_sampleName = "PHY" + str(count_sample)
            sample_rename_dict[sampleName] = new_sampleName
            count_sample+=1
            num_add_space = 10 - len(new_sampleName)

            new_sampleName_final = addSpaceToString(new_sampleName,num_add_space)

            dummy = variant_info[count].split("\t")
            dummy2 = "".join(dummy)
            variant_info_final = encrypt(dummy2, 10)
            full_sequence_string = new_sampleName_final + variant_info_final

            phylip.write(full_sequence_string+"\n")
            count+1

    with open(sampleID_index, "w") as indexFile:
        for key,value in sample_rename_dict.items():
            old_name = key
            new_name = value

            indexFile.write(old_name + "\t" + new_name + "\n")

## main script
first_line = True
index = []
data = []
# Read variant extract file and save to dataframe
with open(variant_extract_file) as f:
    for line in f:
        data_line = line.splitlines()[0]
        if first_line == True :
            header = data_line.split("\t")
            column = header[1:]
            #print(header[-1])
            #print(column[-1])
            first_line = False
        else:
            dummy_data = data_line.split("\t")
            index.append(dummy_data[0])
            dummy_core_data = list(map(int,dummy_data[1:]))
            data.append(dummy_core_data)
            #print(line)
df = pd.DataFrame(data,columns=column,index=index)
df_t = df.transpose()

####################################################

# Check Homo flag, If True we will convert hetero value 1 to 0 (Will update better way to handle this later)
if homo_only == True:
    for (columnName,columnData) in df_t.iteritems():
        df_t.loc[df_t[columnName] < 2, columnName] = 0
########################################################

# Check binary flag, If True we will convert value that more than 0 to 1 (Will update better way to handle this later)
if binary == True:
    for (columnName,columnData) in df_t.iteritems():
        df_t.loc[df_t[columnName] > 0, columnName] = 1
########################################################

write_data_to_phylip(outputPath,df_t)