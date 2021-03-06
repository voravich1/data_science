
import pandas as pd
from scipy.spatial.distance import squareform, pdist
import seaborn as sns
import stat_utility

#input = "./test_data/master_extract.txt"
#input = r'G:\TB_BGI\all_sample_res\BGI_174_combine_extract.txt'
#input = r'C:\Users\vorav\Downloads\1188_mantaV1_4_combine_extract.txt'
input = "/Users/worawich/Downloads/1170_delprofiler/del_analysis/lin1/svtk_batch500/sv_master_del_extract.txt"
outputDendogram = "/Users/worawich/Downloads/1170_delprofiler/del_analysis/lin1/svtk_batch500/dendogram.pdf"
treFile = "/Users/worawich/Downloads/1170_delprofiler/del_analysis/lin1/svtk_batch500/test2.tre"
itolFile = "/Users/worawich/Downloads/1170_delprofiler/del_analysis/lin1/svtk_batch500/label_cluster.itol.txt"
pvalue_csv_file = "/Users/worawich/Downloads/1170_delprofiler/del_analysis/lin1/svtk_batch500/pvalue_cluster.csv"
freq_csv_file = "/Users/worawich/Downloads/1170_delprofiler/del_analysis/lin1/svtk_batch500/freq_cluster.csv"

homo_only = True
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



if homo_only == True:
    for (columnName,columnData) in df_t.iteritems():
        df_t.loc[df_t[columnName] < 2, columnName] = 0


    #df_t = (df_t.iloc[:,:] == 1).replace(int(0)).astype(int)
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

#with open(treFile, "w") as text_file:
 #   text_file.write(newick_string)


#plt.show()

#plt.savefig(outputDendogram, dpi=1200)

max_d = 8.5
cluster = fcluster(Z, max_d, criterion='distance')
index = df_t.index
idx_l = index.tolist()
dn_full = dendrogram(Z, labels=df_t.index, orientation='left')


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
a = dict()

def without_keys(input_dict,exclude_key_dict):
    return {k:v for k,v in input_dict.items() if k not in exclude_key_dict}


list_pvalue_df_res = []
list_tscore_df_res = []
list_freq_df_res = []

for group in cluster_dict:
    combination_name = "group " + str(group) + " vs other"
    sample_list_groupA = cluster_dict[group]
    exclude_key_dict = {group}
    other_group_dict = without_keys(cluster_dict,exclude_key_dict)
    list_other_group = other_group_dict.values()
    sample_list_groupB = []
    for l in list_other_group:
        for item in l:
            sample_list_groupB.append(item)

    dataframe_groupA = df_t.loc[sample_list_groupA]
    dataframe_groupB = df_t.loc[sample_list_groupB]

    # Do ttest for this combination
    pvalue_res, tscore_res, freq_res = stat_utility.multipleColumnFisherTest(dataframe_groupA, dataframe_groupB)

    pvalue_res.index = [combination_name]
    tscore_res.index = [combination_name]
    freq_res.index = [combination_name]

    list_pvalue_df_res.append(pvalue_res)
    list_tscore_df_res.append(tscore_res)
    list_freq_df_res.append(freq_res)

pvalue_res_all = pd.concat(list_pvalue_df_res)
tscore_res_all = pd.concat(list_tscore_df_res)
freq_res_all = pd.concat(list_freq_df_res)

########################################################


pvalue_res_all.to_csv(pvalue_csv_file)
freq_res_all.to_csv(freq_csv_file)

pvalue_res_all_t = pvalue_res_all.transpose()


for col_name, col_data in pvalue_res_all_t.iteritems():
    sample_n = col_name
    t = col_data[:] < 0.05

    f = col_data[t]
    np_data = f.to_numpy()

    print(f)
    print()

test = pvalue_res_all_t[:] < 0.05

filt = pvalue_res_all_t[test]

#group = df_t.loc[['ERR752267','ERR718415']]
#print(group)

#n_array = group.iloc[:,0].values
#print(n_array)

print()






