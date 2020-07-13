import numpy as np
import pandas as pd
from scipy import stats


def multipleColumnTtest(dataframe_popA, dataframe_popB):
    p_value_res = list()
    t_score_res = list()

    num_A , feature_A = dataframe_popA.shape
    num_B, feature_B = dataframe_popB.shape

    column_name = dataframe_popA.columns

    if feature_A != feature_B:
        return None

    for i in range(feature_A):
        n_array_A = dataframe_popA.iloc[:,i].values
        n_array_B = dataframe_popB.iloc[:,i].values

        t_score, p_value = stats.ttest_ind(n_array_A,n_array_B)
        p_value_res.append(p_value)
        t_score_res.append(t_score)

    row_pvalue = [p_value_res]
    row_tscore = [t_score_res]
    dataframe_p_value_res = pd.DataFrame(row_pvalue, columns=column_name, index=['p_value'])
    dataframe_t_score_res = pd.DataFrame(row_tscore, columns=column_name, index=['t_score'])

    return dataframe_p_value_res, dataframe_t_score_res

def multipleColumnFisherTest(dataframe_popA, dataframe_popB):
    p_value_res = list()
    odd_ratio_res = list()
    frequency_res = list()

    num_A , feature_A = dataframe_popA.shape
    num_B, feature_B = dataframe_popB.shape

    column_name = dataframe_popA.columns

    if feature_A != feature_B:
        return None

    column_rename_dict = dict()
    for i in range(feature_A):
        n_array_A = dataframe_popA.iloc[:,i].values
        n_array_B = dataframe_popB.iloc[:,i].values

        ## create contingency table
        count_A_alt = 0
        count_A_ref = 0
        count_B_alt = 0
        count_B_ref = 0

        for val in n_array_A:
            if val == 0:
                count_A_ref = count_A_ref + 1
            elif val == 1:
                count_A_ref = count_A_ref + 1
            elif val == 2:
                count_A_alt = count_A_alt + 1

        for val in n_array_B:
            if val == 0:
                count_B_ref = count_B_ref + 1
            elif val == 1:
                count_B_ref = count_B_ref + 1
            elif val == 2:
                count_B_alt = count_B_alt + 1

        contingency_A = [count_A_alt, count_A_ref]
        contingency_B = [count_B_alt, count_B_ref]
        contingency_table = [contingency_A, contingency_B]
        #############################
        ## Do fisher exact test. If p value < 0.05 mean popA and popB has significant different.
        ## For more interpretation of the word significant different we may need to pbserve count of ALT and REF of both popA and B.
        ## In order to interprete that which population prefer ALT or REF
        odd_ratio, p_value = stats.fisher_exact(contingency_table)
        p_value_res.append(p_value)
        odd_ratio_res.append(odd_ratio)
        #############################
        count_data = " (altA:refA|altB:refB," + str(count_A_alt) + ":" + str(count_A_ref) + "|" + str(count_B_alt) + ":" + str(count_B_ref) + ")"
        frequency_res.append(count_data)
        ## Create Rename column dict
        #col_name = column_name[i]
        #new_col_name = col_name + " (altA:refA|altB:refB," + str(count_A_alt) + ":" + str(count_A_ref) + "|" + str(count_B_alt) + ":" + str(count_B_ref) + ")"
        #column_rename_dict[col_name] = new_col_name
        #############################

    row_pvalue = [p_value_res]
    row_tscore = [odd_ratio_res]
    row_frequency = [frequency_res]
    dataframe_p_value_res = pd.DataFrame(row_pvalue, columns=column_name, index=['p_value'])
    dataframe_odd_ratio_res = pd.DataFrame(row_tscore, columns=column_name, index=['odd_ratio'])
    dataframe_frequency_res = pd.DataFrame(row_frequency, columns=column_name, index=['REF_ALT_Frequency'])

    #dataframe_p_value_res.rename(columns=column_rename_dict, inplace = True)
    #dataframe_odd_ratio_res.rename(columns=column_rename_dict, inplace = True)


    return dataframe_p_value_res, dataframe_odd_ratio_res, dataframe_frequency_res

def extract_marker_from_pvalue(dataframe_p_value):
    dataframe_p_value_t = dataframe_p_value.transpose()




def without_keys(input_dict,exclude_key_dict):
    return {k:v for k,v in input_dict.items() if k not in exclude_key_dict}