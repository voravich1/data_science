import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
import warnings
warnings.filterwarnings("ignore")
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn.metrics import confusion_matrix
from scipy.spatial.distance import squareform, pdist

np.random.seed(123)

#data = pd.read_csv('./test_data/data.csv')
data = pd.read_csv('./test_data/clean_cut_mcv_mch.csv')
print(data)

data = data.iloc[:,1:]

label_encoder = LabelEncoder()
data.iloc[:,0] = label_encoder.fit_transform(data.iloc[:,0]).astype('float64')

corr = data.corr()

sns.heatmap(corr)

##

from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

X = data.iloc[:,1:].values
y = data.iloc[:,0].values

clf = LinearDiscriminantAnalysis()
clf.fit(X, y)

coef = clf.coef_

abs_coef = coef.__abs__()

## get rid off feature that has coef lower than 1



##


## remove one of two features that have a correlation higher than 0.9

columns = np.full((corr.shape[0],), True, dtype=bool)
for i in range(corr.shape[0]):
    for j in range(i+1, corr.shape[0]):
        if corr.iloc[i,j] >= 0.9:
            if columns[j]:
                columns[j] = False
selected_columns = data.columns[columns]
data = data[selected_columns]



## Selecting columns based on p-value

selected_columns = selected_columns[1:].values
import statsmodels.api as sm

def backwardElimination(x, Y, sl, columns):
    numVars = len(x[0])
    for i in range(0, numVars):
        regressor_OLS = sm.OLS(Y, x).fit()
        maxVar = max(regressor_OLS.pvalues).astype(float)
        if maxVar > sl:
            for j in range(0, numVars - i):
                if (regressor_OLS.pvalues[j].astype(float) == maxVar):
                    x = np.delete(x, j, 1)
                    columns = np.delete(columns, j)

    regressor_OLS.summary()
    return x, columns

SL = 0.05
data_modeled, selected_columns = backwardElimination(data.iloc[:, 1:].values, data.iloc[:, 0].values, SL, selected_columns)


# move result to new dataframe
result = pd.DataFrame()
result['diagnosis'] = data.iloc[:,0]

data = pd.DataFrame(data = data_modeled, columns = selected_columns)



# visualize feature
fig = plt.figure(figsize = (20, 25))
j = 0
for i in data.columns:
    plt.subplot(2, 3, j+1)
    j += 1
    sns.distplot(data[i][result['diagnosis']==0], color='g', label = 'Alpha_thal_1_trait')
    sns.distplot(data[i][result['diagnosis']==1], color='r', label = 'Alpha_thal_2_trait')
    plt.legend(loc='best')
fig.suptitle('Alpha thal trait Analysis')
fig.tight_layout()
fig.subplots_adjust(top=0.95)
plt.show()




# split sample to train test
x_train, x_test, y_train, y_test = train_test_split(data.values, result.values, test_size = 0.2)

svc=SVC() # The default kernel used by SVC is the gaussian kernel
#svc=SVC(kernel='rbf')
#from sklearn.model_selection import cross_val_score
#scores = cross_val_score(clf, data.values, result.values, cv=5)
#print(scores)

svc.fit(x_train, y_train)

prediction = svc.predict(x_test)

#print("Accuracy: %0.2f (+/- %0.2f)" % (scores.mean(), scores.std() * 2))


#from sklearn import metrics
#scores = cross_val_score(clf, X, y, cv=5, scoring='f1_macro')
#print(scores)

cm = confusion_matrix(y_test, prediction)
sum = 0
for i in range(cm.shape[0]):
    sum += cm[i][i]

accuracy = sum / x_test.shape[0]
print(accuracy)

print("")
## Test create distance matrix ##
#data2 = pd.read_csv('./test_data/data.csv')
#data2_edited = data2.iloc[:,0:-1]
#label_encoder = LabelEncoder()
#data2_edited.iloc[:,1] = label_encoder.fit_transform(data2_edited.iloc[:,1]).astype('float64')
#dist_matrix=pd.DataFrame(squareform(pdist(data2_edited.iloc[:, 1:])), columns=data2.id.unique(), index=data2.id.unique())
#sns.set(font_scale=1)
#sns.heatmap(dist_matrix.iloc[0:30,0:30], xticklabels=True, yticklabels=True, linewidths=.1, cmap="Greens")
####################################


