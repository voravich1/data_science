
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
itol_ann = "/Users/worawich/Downloads/1170_delprofiler/del_analysis/lin1/scipy_dendrogram.itol.txt"

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

# Extract cluster group. Make itol annotation foe cluster.
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



