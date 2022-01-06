#!/usr/bin/env python

# This script use for resolve range of number in position start and stop column
# It will find overlab range of position and then generate a new range value 
# the expect input file shoud have three column [chr,start,stop]
# 

import os

range_exclude_file = "/Users/worawich/Downloads/TB_del_paper/exclude_file/sum_range_Pbig_Pwas_exclude.txt" #file should be sort on start column. chr column must be the same
outputDir = "/Users/worawich/Downloads/TB_del_paper/exclude_file"

output_summary_file = os.path.join(outputDir, "exclude_resolve.txt")
output = open(output_summary_file,"w")

def group(L):
    first = last = L[0]
    for n in L[1:]:
        if n - 1 == last: # Part of the group, bump the end
            last = n
        else: # Not part of the group, yield current group and start a new
            yield first, last
            first = last = n
    yield first, last # Yield the last group


with open(range_exclude_file) as f:

    data_list = list()

    for line in f:
        info = line.split()
        new_chr = info[0]
        new_start = int(info[1])
        new_stop = int(info[2])

        #convert range to sequence data
        range_data = range(new_start,(new_stop+1))

        for i in range_data:
            data_list.append(i)

# sort and do uniq value
data_list.sort()
data_list_uniq = list(set(data_list))
data_list_uniq.sort()

# use function group to convert from sequnce data to range 
range_data = group(data_list_uniq)

for data in range_data:
    start = data[0]
    stop = data[1]
    output.write(new_chr + "\t" + str(start) + "\t" + str(stop) + "\n")


output.close()             
