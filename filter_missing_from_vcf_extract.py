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

variant_extract_file = "/Users/worawich/Downloads/TB_del_paper/1170_x_manilla/analysis/sv/original_sample_plus_135_manilla/SNP_extract/snp_extract_inc_missing_gt.txt"
outputPath = "/Users/worawich/Downloads/TB_del_paper/1170_x_manilla/analysis/sv/original_sample_plus_135_manilla/SNP_extract"
log = "/Users/worawich/Downloads/TB_del_paper/1170_x_manilla/analysis/sv/original_sample_plus_135_manilla/SNP_extract/extract_log.txt"
missing_th = 10

## main script
output_filename = "vcf_extract_filter_missing_" + str(missing_th) + ".txt"
output_file = os.path.join(outputPath, output_filename)
output = open(output_file,"w")
output_log = open(log,"w")
first_line = True
# Read variant extract file and save to dataframe
with open(variant_extract_file) as f:
    for line in f:
        data_line = line.splitlines()[0]
        if first_line == True :
            output.write(data_line+"\n")
            first_line = False
        else:
            dummy_data = data_line.split("\t")
            num_sample = len(dummy_data) - 1 
            missing_count = 0
            for info in dummy_data:
                if info == "-1":
                    missing_count+=1
            if (missing_count * 100)/num_sample <= missing_th:
                output.write(data_line+"\n")
                output_log.write(str(missing_count)+"\n")
output.close()
output_log.close()