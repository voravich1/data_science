#!/usr/bin/env python

"""
This script has function to create aonprofiler summary result file from two column meta data 

input
1. Need two column txt format (sampltname|lineage) with header
Noted lineage column must be in format L1.1.1 or lineage1.1.1 or l1.1.1 or Lineage1.1.1

output
1. summary result file in format aonprofiler


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

meta_data_file = "/Users/worawich/Downloads/TB_del_paper/indian_ocean/analysis/nextflow_sv/summary_x_meta_edit_for_compat_with_manta_res2.txt"
outputDir = "/Users/worawich/Downloads/TB_del_paper/indian_ocean/analysis/nextflow_sv"
header = True

"""
##### rewrite lineage
lineage="l1.2.1"
lin_split = lineage.split(".")
lin_tier = len(lin_split)
new_lineage = ""
for index in range(lin_tier):
    dummy_newlin = ""
    for i in range(index+1):
        if i == 0:
            dummy_newlin = lin_split[i]
        else:
            dummy_newlin = dummy_newlin+"."+lin_split[i]
    if index == 0:
        new_lineage = dummy_newlin
    else:
        new_lineage = new_lineage + "|" + dummy_newlin

print(new_lineage)
"""

output_summary_file = os.path.join(outputDir, "converted_meta_data_to_summary_res2.txt")
output = open(output_summary_file,"w")

with open(meta_data_file) as f:
    for line in f:
        data = line.splitlines()[0]

        if header == True:
            header = False
            output.write(data+"\n")
            continue
        
        sample_name = data.split("\t")[0]
        lineage = data.split("\t")[1]

        lin_split = lineage.split(".")
        lin_tier = len(lin_split)
        new_lineage = ""
        for index in range(lin_tier):
            dummy_newlin = ""
            for i in range(index+1):
                if i == 0:
                    dummy_newlin = lin_split[i]
                else:
                    dummy_newlin = dummy_newlin+"."+lin_split[i]
            if index == 0:
                new_lineage = dummy_newlin
            else:
                new_lineage = new_lineage + "|" + dummy_newlin
        
        output.write(sample_name + "\t" + new_lineage + "\n")

output.close()
        




