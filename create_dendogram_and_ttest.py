#!/usr/bin/env python

"""This program is use for 2 purpose
    1. create dendrogram from any VCF extract DATA
    What is "VCF extract DATA"?
    VCF extract DATA will be in format like this
        1st row is header ==> POS|Sample1|Sample2|Sample3|...|SampleN
        column POS contain information of chr and position merge into single string EX "1:100-300" mean chr1, start=100, end=300
        other column contain representation number of genotype that got from GT in VCF. 0 => homoREF, 1=>hetero, 2=>homoALT

    2. Do ttest statistic for finding significant varaint between two population. Do on all possible combination of population

    "color_th" is major factor that indicate how many poppulation will be for this analysis. High value of color_th will give low number
    population.

    An easy way to get VCF extract DATA for SV vcf such as manta SV is using
    "run_svtk_standard_cluster_extract_vcf_huge_sample_set_TB_only.sh" from sv_tb_pipeline
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

#input = "./test_data/master_extract.txt"
#input = r'G:\TB_BGI\all_sample_res\BGI_174_combine_extract.txt'
#input = r'C:\Users\vorav\Downloads\1188_mantaV1_4_combine_extract.txt'
#input = "/Users/worawich/Downloads/1170_delprofiler/del_analysis/lin1/svtk_batch500/sv_master_del_extract.txt"
input = "/Volumes/10TBSeagateBackupPlus/NBT/TB_project/1170_delprofiler/wasna_paper/snp_sv_dendogram/snp_del_extract_ready.txt"
outputfolder = "/Volumes/10TBSeagateBackupPlus/NBT/TB_project/1170_delprofiler/wasna_paper/snp_sv_dendogram"
outputDendogram = os.path.join(outputfolder,"dendogram.pdf")
treFile = os.path.join(outputfolder,"scipy_dendrogram.tre")
itol_ann = os.path.join(outputfolder,"scipy_dendrogram.itol.txt")
pvalue_csv_file = os.path.join(outputfolder,"pvalue_fisher_cluster.csv")
freq_csv_file = os.path.join(outputfolder,"freq_cluster.csv")

#color_th = 8.5  # this treshold can be adjust it will effect clustering and group coloring on dendrogram (old threshold when we do just del alone)
color_th = 25
homo_only = False # Homo flag ==> if True consider homo region by convert hetero value 1 to 0 (Will update better way to handle this later)
ttest = False
## function convert linkage resut to newick file
# credit https://stackoverflow.com/questions/28222179/save-dendrogram-to-newick-format
def getNewick(node, newick, parentdist, leaf_names):
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parentdist - node.dist, newick)
        else:
            newick = ");"
        newick = getNewick(node.get_left(), newick, node.dist, leaf_names)
        newick = getNewick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
        newick = "(%s" % (newick)
        return newick

first_line = True
column = list()
index = list()
data = list()

# read input file: This input file must be SV_std_combine_extract_vcf get from SV post process analysis
with open(input) as f:
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
############################################################

# Check Homo flag, If True we will convert hetero value 1 to 0 (Will update better way to handle this later)
if homo_only == True:
    for (columnName,columnData) in df_t.iteritems():
        df_t.loc[df_t[columnName] < 2, columnName] = 0
########################################################


plt.figure(num=None, figsize=(60, 120), dpi=80, facecolor='w', edgecolor='k')

# Do linkage function can be consider as calculate distance matrix (My understanding)
Z = linkage(df_t, 'ward')
#############################

# Create dedogram plot and save to pdf file
#color_th = 8.5 # this treshold can be ad just it will effect group coloring on dendrogram (old threshold when we do just del alone)
plt.title('Hierarchical Clustering Dendrogram')
plt.xlabel('sample index')
plt.ylabel('distance (Ward)')
dn = dendrogram(Z, labels=df_t.index, orientation='left',color_threshold=color_th)
plt.savefig(outputDendogram, dpi=1200)
######################################################################################


# generate newick string form linkage result and save to .tre file
tree = to_tree(Z,False)
newick_string = getNewick(tree, "", tree.dist, leaf_names=df_t.index)

with open(treFile, "w") as text_file:
    text_file.write(newick_string)
###############################################

# Extract cluster group. Make itol annotation for cluster.
max_d = color_th
cluster = fcluster(Z, max_d, criterion='distance')
index = df_t.index
idx_l = index.tolist()


index128colors = [
        "#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
        "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
        "#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
        "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
        "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
        "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
        "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
        "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
        "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
        "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
        "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
        "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
        "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C",
        "#83AB58", "#001C1E", "#D1F7CE", "#004B28", "#C8D0F6", "#A3A489", "#806C66", "#222800",
        "#BF5650", "#E83000", "#66796D", "#DA007C", "#FF1A59", "#8ADBB4", "#1E0200", "#5B4E51",
        "#C895C5", "#320033", "#FF6832", "#66E1D3", "#CFCDAC", "#D0AC94", "#7ED379", "#012C58"]


with open(itol_ann, "w") as text_file:
    text_file.write("DATASET_STYLE\n")
    text_file.write("SEPARATOR SPACE\n")
    text_file.write("DATASET_LABEL Cluster style\n")
    text_file.write("COLOR #ffff00\n")

    text_file.write("LEGEND_TITLE Cluster\n")
    max = cluster.max()
    legend_shape = "LEGEND_SHAPES"
    legend_color = "LEGEND_COLORS"
    legend_label = "LEGEND_LABELS"
    for group in range(1,max+1):
        legend_shape = legend_shape + " 1"
        color = index128colors[group]
        legend_color = legend_color + " " + color
        legend_label = legend_label + " group" + str(group)

    text_file.write(legend_shape + "\n")
    text_file.write(legend_color + "\n")
    text_file.write(legend_label + "\n\n")

    text_file.write("DATA\n")

    for idx,val in enumerate(cluster):
        sample = idx_l[idx]
        group = val
        color = index128colors[val]

        itol_str1 = sample + " label node #000000 1 bold-italic " + color + "\n"
        itol_str2 = sample + " branch node " + color + " 3 normal \n"

        text_file.write(itol_str1)
        text_file.write(itol_str2)
##############################################################################


## create cluster dict (mapping cluster group and sample name by index order of list)
cluster_dict = dict()
for idx, group in enumerate(cluster):
    sample_name = column[idx]

    if group in cluster_dict:
        dummy_list = cluster_dict[group]
        dummy_list.append(sample_name)

        cluster_dict[group] = dummy_list
    else:
        dummy_list = [sample_name]
        cluster_dict[group] = dummy_list

##################################################

# Tranform Data in to every possible combination (one group VS other group) and do Ttest
list_pvalue_df_res = []
list_score_df_res = []
list_freq_df_res = []

for group in cluster_dict:
    combination_name = "group " + str(group) + " vs other"
    sample_list_groupA = cluster_dict[group]
    exclude_key_dict = {group}
    other_group_dict = stat_utility.without_keys(cluster_dict,exclude_key_dict)
    list_other_group = other_group_dict.values()
    sample_list_groupB = []
    for l in list_other_group:
        for item in l:
            sample_list_groupB.append(item)

    dataframe_groupA = df_t.loc[sample_list_groupA]
    dataframe_groupB = df_t.loc[sample_list_groupB]

    # Do ttest or chi square for this combination
    if ttest == True:
        pvalue_res, score_res = stat_utility.multipleColumnTtest(dataframe_groupA, dataframe_groupB)
        pvalue_res.index = [combination_name]
        score_res.index = [combination_name]

        list_pvalue_df_res.append(pvalue_res)
        list_score_df_res.append(score_res) ## list of dataframe of tscore
    else:
        pvalue_res, score_res, frequency_res = stat_utility.multipleColumnFisherTest(dataframe_groupA, dataframe_groupB)

        pvalue_res.index = [combination_name]
        score_res.index = [combination_name]
        frequency_res.index = [combination_name]

        list_pvalue_df_res.append(pvalue_res)
        list_score_df_res.append(score_res) ## list of dataframe of odd ratio
        list_freq_df_res.append(frequency_res)

if ttest == True:
    pvalue_res_all = pd.concat(list_pvalue_df_res)
    score_res_all = pd.concat(list_score_df_res)

    pvalue_res_all.to_csv(pvalue_csv_file)
else:
    pvalue_res_all = pd.concat(list_pvalue_df_res)
    score_res_all = pd.concat(list_score_df_res)
    freq_res_all = pd.concat(list_freq_df_res)

    pvalue_res_all.to_csv(pvalue_csv_file)
    freq_res_all.to_csv(freq_csv_file)

########################################################




