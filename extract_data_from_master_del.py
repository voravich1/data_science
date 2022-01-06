import os

list_target_del_file="/Users/worawich/Downloads/TB_del_paper/1170_x_manilla/analysis/sv/original_sample_plus_135_manilla/target_del.txt"
master_del_extract_file="/Users/worawich/Downloads/TB_del_paper/1170_x_manilla/analysis/sv/original_sample_plus_135_manilla/sv_master_del_extract.txt"
outputDir="/Users/worawich/Downloads/TB_del_paper/1170_x_manilla/analysis/sv/original_sample_plus_135_manilla/"

output_summary_file = os.path.join(outputDir, "target_del_extract_data.txt")
output = open(output_summary_file,"w")

list_target_del = []

with open(list_target_del_file) as f:

    for line in f:
        info = line.splitlines()[0]
        list_target_del.append(info)


#output.write()

header = True

data_dict = dict()

with open(master_del_extract_file) as f:
    for line in f:
        #info = line.split("\t")
        #pos = info[0].split(":")[1]

        if header == True:
            header_info = line.splitlines()[0]
            output.write(header_info + "\n")
            header = False
        else:
            info = line.splitlines()[0]
            info_list = info.split("\t")
            pos = info_list[0].split(":")[1]
            if pos in list_target_del:
                data_dict[pos] = info
                #output.write(info + "\n")


for pos in list_target_del:
    if pos in data_dict:
        info = data_dict[pos]
        output.write(info + "\n")

output.close