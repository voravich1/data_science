#!/usr/bin/env python

"""
This script has function to create itol annotation label file for deletion

input
1. Need VCF extract file  same format as use for statistical test for finding marker
2. Need variantion target file contain chr and start stop position
    The file has only one column and row contain any variation info
    Variation info must be in the same format as VCF extract file e.g. chr:start-stop

Noted: variation info can be other format But!! VCF extract file and variation target file must use the same info format

output
1. specify output folder and script will auto generate .itol.txt annotation file

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
target_variant_file = "/Volumes/10TBSeagateBackupPlus/NBT/TB_project/1170_delprofiler/wasna_paper/snp_sv_dendogram/deletion_target.txt"
outputfolder = "/Volumes/10TBSeagateBackupPlus/NBT/TB_project/1170_delprofiler/wasna_paper/snp_sv_dendogram/homo_only"
itol_ann = os.path.join(outputfolder,"variant_label.itol.txt")
homo_only = True


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

# Read target variant file
# Simultaniously map variant target and extract all sample that have variant
variant_hit_dict = dict()
variant_list = []
with open(target_variant_file) as f:
    for line in f:
        target_variant = line.splitlines()[0]
        variant_list.append(target_variant)
        #target_data = df_t[target_variant]
        #target_hit_data = target_data[target_data==2] # 2 represent homo alt that mean it has variant
        #hit_sample = target_hit_data.index
        #hit_sample_list = hit_sample.to_list()
        #variant_hit_dict[target_variant] = hit_sample_list

df_t_variant = df_t[variant_list]

# prepare data for ito annotation (header and info)
if homo_only == True:
    for (columnName,columnData) in df_t_variant.iteritems():
        df_t_variant.loc[df_t_variant[columnName] < 2, columnName] = 0
    for (columnName, columnData) in df_t_variant.iteritems():
        df_t_variant.loc[df_t_variant[columnName] == 2, columnName] = 1
else:
    for (columnName,columnData) in df_t_variant.iteritems():
        df_t_variant.loc[df_t_variant[columnName] > 0, columnName] = 1

variant_header = df_t_variant.columns.to_list()
varaint_info_for_itol = df_t_variant.to_csv(header=None, index=True, sep='\t').split('\n')

## Create itol label file
with open(itol_ann, "w") as text_file:
    text_file.write("DATASET_BINARY\n")
    text_file.write("SEPARATOR TAB\n")
    text_file.write("DATASET_LABEL\tVariant\n")
    text_file.write("COLOR\t#ffff00\n\n")

    text_file.write("SHOW_LABELS\t1\n")
    legend_shape = "FIELD_SHAPES"
    legend_color = "FIELD_COLORS"
    legend_label = "FIELD_LABELS"
    for variant in variant_header:
        legend_shape = legend_shape + "\t2"
        color = "black"
        legend_color = legend_color + "\t" + color
        legend_label = legend_label + "\t" + str(variant)

    text_file.write(legend_shape + "\n")
    text_file.write(legend_color + "\n")
    text_file.write(legend_label + "\n\n")

    text_file.write("DATA\n")

    for info in varaint_info_for_itol:
        text_file.write(info+"\n")
##############################################################################