import sys
import pandas as df

# Calculate Evaluate metrix
# input must be in format
# odd column is predict
# even column is actual value
# It will consider pair column for calculation
# value in actual column should have "R" for resist and "S" for sensitive and "U" for nodata
# Value in predict column must have "S" for sensitive but for resist can be any stirng. No missing data allow in column

input = "/Users/worawich/Downloads/drug_Pbum.csv"
output = "/Users/worawich/Downloads/evaluate_drug.csv"

data = df.read_csv(input)

sample_size, num_drug = data.shape


resist_predict = 0
sensitive_predict = 0
resist_actual = 0
sensitive_actual = 0
no_data_actual = 0

count = 0
actual_list = list()
predict_list = list()

report = open(output, "a")
report.write("Drug,Number of sample,Specificity, Sensitivity(Recall),Precission,F score,Accuracy\n")
for (columnName, columnData) in data.iteritems():

    #num_actual_data = 0
    count = count + 1
    #aon = columnData.value_counts().to_dict()
    if count % 2 == 0:
        #all_data = columnData.count()
        #actual_dict = columnData.value_counts().to_dict()
        #if "S" in actual_dict:
        #    sensitive_actual = actual_dict["S"]
        #if "U" in actual_dict:
        #   no_data_actual = actual_dict["U"]
        #if "R" in actual_dict:
        #    resist_actual = actual_dict["R"]

        #num_actual_data = all_data - no_data_actual
        #if num_actual_data != (resist_actual + sensitive_actual):
        #    print("Error")

        actual_list = columnData.tolist()
    else:
        #all_data = columnData.count()
        #predict_dict = columnData.value_counts().to_dict()
        #if "S" in predict_dict:
        #    sensitive_predict = predict_dict["S"]

        #resist_predict = all_data - sensitive_predict


        predict_list = columnData.tolist()

    # calculate evaluation value. cal every even count.
    if count % 2 == 0:
        TN=0
        FN=0
        TP=0
        FP=0
        P=0
        N=0

        for idx, value in enumerate(predict_list):

            predict_value = value
            actual_value = actual_list[idx]

            if actual_value == "S":
                N = N+1
            elif actual_value == "R":
                P = P+1


            if actual_value == "S" and predict_value == "S":
                TN = TN + 1
            elif actual_value == "R" and predict_value == "S":
                FN = FN + 1
            elif actual_value == "R" and predict_value != "S":
                TP = TP + 1
            elif actual_value == "S" and predict_value != "S":
                FP = FP + 1

        acc = ((TP+TN)/(P+N))*100
        sensitivity = TP/(TP+FN)
        specificity = TN/(TN+FP)
        if TP+FP == 0:
            precission = 0
        else:
            precision = TP/(TP+FP)

        f_score = 2*((precision * sensitivity)/(precision + sensitivity))
        num_all_sample = TN + FN + TP + FP
        drug = columnName.split(" ")[0]
        report_str = drug + "," + str(num_all_sample) + "," + str(round(specificity,2)) + "," + str(round(sensitivity,2)) + "," + str(round(precision,2)) + "," + str(round(f_score,2)) + "," + str(round(acc,2))

        report.write(report_str+"\n")

report.close()
