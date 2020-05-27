import csv

#####
# script for compare master key result file with concatinate vcf file
# use master key file as main DB. Master file must be in tab separate format sample|lineage
# input must be concatinate vcf got from concat vcf python script (It in tbplatform pycharm project)
#####

master_file = r'D:\TB\Martin_new_list\data_for_analyze_1170_sample\mast_lineage_result_1170_sample.txt' #should be tab separate with 2 column sample|lineage
input_file = r'C:\Users\vorav\Downloads\summary_result.txt' #should be summary from tb_collate.py from tb_platform project
output_file = r'D:\TB\Martin_new_list\data_for_analyze_1170_sample\tbd1_found_compare_map.csv'
fail_output_file = r'D:\TB\Martin_new_list\data_for_analyze_1170_sample\tbd1_found_compare_fail.csv'
summary_output_file = r'D:\TB\Martin_new_list\data_for_analyze_1170_sample\tbd1_found_compare_summary.csv'

master_sample_lineage_dict = dict()
dup_list = list()
dup_dict = dict()
lineage_list = list()


def import_master_template(master_filename):
    # loop create dict for master check info
    with open(master_filename) as master:
        for line in master:
            info = line.split("\t")
            sampleID = info[0].split("_")[0]
            lineage = info[1].splitlines()[0]

            if sampleID not in master_sample_lineage_dict:
                master_sample_lineage_dict[sampleID] = lineage
            else:
                if sampleID in dup_dict:
                    lineage_list = dup_dict[sampleID]
                    lineage_list.append(lineage)
                    dup_dict[sampleID] = lineage_list
                    master_sample_lineage_dict[sampleID] = lineage_list
                else:
                    lineage_exist = master_sample_lineage_dict[sampleID]
                    lineage_list = [lineage_exist,lineage]
                    dup_dict[sampleID] = lineage_list
                    master_sample_lineage_dict[sampleID] = lineage_list

    return master_sample_lineage_dict


def import_input_summary_result(summary_input_file):
    # loop input for check with master
    header = True
    result_dict = dict()
    fail_result_dict = master_sample_lineage_dict.copy()
    result_count_dict = dict()
    fail_count_dict = dict()
    with open(input_file) as input:
        for line in input:
            if header == True: # skip header
                header = False
                continue

            info_list = line.split("\t")
            input_sampleID = info_list[0].split("_")[0]
            #input_lineage = info_list[2]

            if input_sampleID in master_sample_lineage_dict:
                master_lineage = master_sample_lineage_dict[input_sampleID]
                if input_sampleID in result_dict:
                    result_list = result_dict[input_sampleID]
                    result_list.append(master_lineage)
                    result_dict[input_sampleID] = result_list
                else:
                    #result_list = [input_lineage, master_lineage]
                    result_list = [master_lineage]
                    result_dict[input_sampleID] = result_list
                    del fail_result_dict[input_sampleID] # delete sampleID from fail dict [at first fail dict is copy of master dict]. the remaining member of this dict is fail sample that not map to input


                    if master_lineage in result_count_dict:
                        dummy_count = result_count_dict[master_lineage]
                        dummy_pass_count = dummy_count[0] + 1
                        result_count_dict[master_lineage] = [dummy_pass_count, 0]
                    else:
                        result_count_dict[master_lineage] = [1, 0]

    for key, val in fail_result_dict.items():

        if type(val) is str:
            lineage = val
            if lineage in result_count_dict:
                dummy_count = result_count_dict[lineage]
                pass_count = dummy_count[0]
                fail_count = dummy_count[1] + 1
                result_count_dict[lineage] = [pass_count, fail_count]
            else:
                result_count_dict[lineage] = [0, 1]
        else:
            for lineage in val:
                if lineage in result_count_dict:
                    dummy_count = result_count_dict[lineage]
                    pass_count = dummy_count[0]
                    fail_count = dummy_count[1] + 1
                    result_count_dict[lineage] = [pass_count, fail_count]
                else:
                    result_count_dict[lineage] = [0, 1]



# write result to csv file
w = csv.writer(open(output_file, "w", newline=''))
w.writerow(["sample","lineage"])
for key, val in result_dict.items():
    w.writerow([key, val])

w = csv.writer(open(fail_output_file, "w", newline=''))
w.writerow(["sample","lineage"])
for key, val in fail_result_dict.items():
    w.writerow([key, val])

summary = csv.writer(open(summary_output_file, "w", newline=''))
summary.writerow(["sample","map","fail"])
for key, val in result_count_dict.items():
    map_count = val[0]
    fail_count = val[1]
    summary.writerow([key, map_count, fail_count])



