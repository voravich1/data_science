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
    freq_ratio_res = list()
    eval_metrices_res = list()
    marker_res = list()

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
        count_A_alt = 0     # number of sample in A that have alt
        count_A_ref = 0     # number of sample in A that have ref
        count_B_alt = 0     # number of sample in B that have alt
        count_B_ref = 0     # number of sample in B that have ref

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
        ## For more interpretation of the word significant different we may need to observe count of ALT and REF of both popA and B.
        ## In order to interprete that which population prefer ALT or REF
        odd_ratio, p_value = stats.fisher_exact(contingency_table)
        #p_value_res.append(p_value) change to p-value full report format
        odd_ratio_res.append(odd_ratio)
        #############################
        count_data = " (altA:refA|altB:refB," + str(count_A_alt) + ":" + str(count_A_ref) + "|" + str(count_B_alt) + ":" + str(count_B_ref) + ")"
        frequency_res.append(count_data)

        ## ratio calculate
        if (count_A_alt + count_A_ref) != 0:
            altA_ratio = count_A_alt / (count_A_alt + count_A_ref)
        else:
            altA_ratio = 0
        
        if (count_B_alt + count_B_ref) != 0:
            altB_ratio = count_B_alt / (count_B_alt + count_B_ref)
        else:
            altB_ratio = 0
        
        if (count_A_alt + count_A_ref) != 0 :
            refA_ratio = count_A_ref / (count_A_alt + count_A_ref)
        else:
            refA_ratio = 0

        if (count_B_alt + count_B_ref) != 0:
            refB_ratio = count_B_ref / (count_B_alt + count_B_ref)
        else:
            refB_ratio = 0

        ratio_data = str(altA_ratio) + ":" + str(refA_ratio) + "|" + str(altB_ratio) + ":" + str(refB_ratio)
        freq_ratio_res.append(ratio_data)
        ## Evaluation metric calculation
        ## imagine that popA and popB is our dataset
        ## popA contain sample that should have deletion (should have all count as altA)
        ## popB contain sample taht shoul not have deletion (should have all count as refB)

        total_A = count_A_ref + count_A_alt     # consider A as class 1
        total_B = count_B_ref + count_B_alt     # consider B as class 0
        TP = count_A_alt    # number of sample in A that found alt (Can say in ML term as actual is alt and predict is also alt)
        FP = count_B_alt    # number of sample in B that found alt (Can say in ML term as actual is ref and predict is alt)
        FN = count_A_ref    # number of sample in A that found ref (Can say in ML term as actual is alt and predict is ref)
        TN = count_B_ref    # number of sample in B that found ref (Can say in ML term as actual is ref and predict is also ref)
        sensitivity = TP/(TP+FN)
        specificity = TN/(TN+FP+0.00000001)
        #if TP+FP == 0:
        #    precision = int(0.000001)
        #else:
        #    precision = TP/(TP+FP)

        #f_score = 2*((precision*sensitivity)/(precision+sensitivity))
        #evaluation_data = "(sen|spec|pre|f1)" + "(" + str(round(sensitivity),4) + "|" + str(round(specificity),4) + "|" + str(round(precision),4) + "|" + str(round(f_score),4) + ")"
        evaluation_data = "(sen|spec|pre|f1)" + "(" + str(round(sensitivity, 4)) + "|" + str(round(specificity,4)) + ")"
        eval_metrices_res.append(evaluation_data)

        ## marker judgement
        judgement = marker_judgement(p_value,count_A_alt,count_B_alt,count_A_ref,count_B_ref)
        judgement_data = 0
        if judgement == True:
            judgement_data = str(round(p_value, 6)) + "|" + str(round(altA_ratio, 6)) + "|" + str(round(altB_ratio, 6)) + "|" + str(count_A_alt) + ":" + str(count_A_ref) + "|" + str(count_B_alt) + ":" + str(count_B_ref)

        marker_res.append(judgement_data)

        ## new modified pvalue report file format
        ## change to full report format contain pvalue|frequency_ratio|frequency_count instead of only p-value
        pvalue_full_report = str(round(p_value, 6)) + "|" + str(round(altA_ratio, 6)) + "|" + str(round(altB_ratio, 6)) + "|" + str(count_A_alt) + ":" + str(count_A_ref) + "|" + str(count_B_alt) + ":" + str(count_B_ref)
        p_value_res.append(pvalue_full_report)
        ## Create Rename column dict
        #col_name = column_name[i]
        #new_col_name = col_name + " (altA:refA|altB:refB," + str(count_A_alt) + ":" + str(count_A_ref) + "|" + str(count_B_alt) + ":" + str(count_B_ref) + ")"
        #column_rename_dict[col_name] = new_col_name
        #############################

    row_pvalue = [p_value_res]
    row_tscore = [odd_ratio_res]
    row_frequency = [frequency_res]
    row_freq_ratio = [freq_ratio_res]
    row_marker = [marker_res]
    row_eval_metrices = [eval_metrices_res]
    dataframe_p_value_res = pd.DataFrame(row_pvalue, columns=column_name, index=['p_value'])
    dataframe_odd_ratio_res = pd.DataFrame(row_tscore, columns=column_name, index=['odd_ratio'])
    dataframe_frequency_res = pd.DataFrame(row_frequency, columns=column_name, index=['REF_ALT_Frequency'])
    dataframe_freq_ratio = pd.DataFrame(row_freq_ratio, columns=column_name, index=['REF_ALT_Frequency_Ratio'])
    dataframe_marker = pd.DataFrame(row_marker, columns=column_name, index=['marker_judgement'])
    dataframe_eval_metrices = pd.DataFrame(row_eval_metrices, columns=column_name, index=['eval_metrices'])

    #dataframe_p_value_res.rename(columns=column_rename_dict, inplace = True)
    #dataframe_odd_ratio_res.rename(columns=column_rename_dict, inplace = True)


    return dataframe_p_value_res, dataframe_odd_ratio_res, dataframe_frequency_res, dataframe_freq_ratio, dataframe_marker, dataframe_eval_metrices

def marker_judgement(p_value,count_A_alt,count_B_alt,count_A_ref,count_B_ref):
    altA_ratio = count_A_alt / (count_A_alt + count_A_ref)
    altB_ratio = count_B_alt / (count_B_alt + count_B_ref+0.00000000001)
    refA_ratio = count_A_ref / (count_A_alt + count_A_ref)
    refB_ratio = count_B_ref / (count_B_alt + count_B_ref+0.00000000001)

    if p_value <= 0.05 and altA_ratio > 0.9 and altB_ratio < 0.9:   # altA_ratio ==> we call target deletion ratio, altB_ratio ==> we call non-target deletion ratio
    #if p_value <= 0.05 and altA_ratio > 0.9 and refB_ratio > 0.9:
    #if p_value <= 0.05:
        true_marker = True
    else:
        true_marker = False

    return true_marker

def extract_marker_from_pvalue(dataframe_p_value):
    dataframe_p_value_t = dataframe_p_value.transpose()




def without_keys(input_dict,exclude_key_dict):
    return {k:v for k,v in input_dict.items() if k not in exclude_key_dict}