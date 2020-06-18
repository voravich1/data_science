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