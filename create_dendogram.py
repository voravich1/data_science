
import pandas as pd
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, to_tree
import numpy as np

#input = "./test_data/master_extract.txt"
#input = r'G:\TB_BGI\all_sample_res\BGI_174_combine_extract.txt'
#input = r'C:\Users\vorav\Downloads\1188_mantaV1_4_combine_extract.txt'
input = "/Users/worawich/Downloads/1170_delprofiler/del_analysis/lin1/sv_master_combine_extract.txt"
outputDendogram = "/Users/worawich/Downloads/1170_delprofiler/del_analysis/lin1/dendogram.pdf"
treFile = "/Users/worawich/Downloads/1170_delprofiler/del_analysis/lin1/scipy_dendrogram.tre"

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

plt.figure(num=None, figsize=(60, 120), dpi=80, facecolor='w', edgecolor='k')

# Do linkage function can be consider as calculate distance matrix (My understanding)
Z = linkage(df_t, 'ward')
#############################

# Create dedogram plot and save to pdf file
color_th = 8.5 # this treshold can be ad just it will effect group coloring on dendrogram
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






