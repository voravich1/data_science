
import pandas as pd
from scipy.spatial.distance import squareform, pdist
import seaborn as sns

#input = "./test_data/master_extract.txt"
#input = r'G:\TB_BGI\all_sample_res\BGI_174_combine_extract.txt'
#input = r'C:\Users\vorav\Downloads\1188_mantaV1_4_combine_extract.txt'
input = "/Users/worawich/Downloads/1170_delprofiler/del_analysis/lin1/sv_master_combine_extract.txt"
outputDendogram = "/Users/worawich/Downloads/1170_delprofiler/del_analysis/lin1/dendogram.pdf"
treFile = "/Users/worawich/Downloads/1170_delprofiler/del_analysis/lin1/test2.tre"

first_line = True
column = list()
index = list()
data = list()
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
#print(df)
#print(df_t)
#print(df_t.iloc[:, 1:])

#df_t.to_csv (r'C:\Users\vorav\Downloads\1188_mantaV1_4_combine_transform.csv', index=True, header=True)

#dist_matrix=pd.DataFrame(squareform(pdist(df_t)),columns=column,index=column)
#print(dist_matrix)
#sns.set(font_scale=1)
#sns.heatmap(dist_matrix.iloc[0:30,0:30], xticklabels=True, yticklabels=True, linewidths=.1, cmap="Greens")


aon = df_t.index

from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, to_tree
import numpy as np


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


#url = 'https://python-graph-gallery.com/wp-content/uploads/mtcars.csv'
#df = pd.read_csv(url)
#df = df.set_index('model')

plt.figure(num=None, figsize=(60, 120), dpi=80, facecolor='w', edgecolor='k')

Z = linkage(df_t, 'ward')



plt.title('Hierarchical Clustering Dendrogram')
plt.xlabel('sample index')
plt.ylabel('distance (Ward)')
dn = dendrogram(Z, labels=df_t.index, orientation='left', color_threshold=8.5)

ivl = dn['ivl']
ivl.reverse()
tree = to_tree(Z,False)
newick_string = getNewick(tree, "", tree.dist, leaf_names=df_t.index)

with open(treFile, "w") as text_file:
    text_file.write(newick_string)


color_new = dn['color_list']
dn['color_list'] = ['#0000ff','#0000ff','#0000ff','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','#000000']

#plt.show()

plt.savefig(outputDendogram, dpi=1200)

max_d = 8.5
cluster = fcluster(Z, max_d, criterion='distance')

dn_full = dendrogram(Z, labels=df_t.index, orientation='left')
label = dn_full['ivl']

new_color_label = list()
#for idx, value in enumerate(label):






