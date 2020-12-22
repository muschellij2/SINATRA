import numpy as np
from scipy.stats import norm
import pandas

from models import ClassificationModel

# import data 
path = './Results/sinatra_ECT_dataROC.txt'
df = pandas.read_csv(path, sep = ",", header=None)
data = df.values
data = data.astype(int)

# interesting if preprocessing the data first would help. Standardize each feature to be standard normal.
Y = data[:,0]
X = data[:,1:]

# preprocess
X = X - np.mean(X,axis = 0)
sd = np.std(X,axis = 0)
sd[sd == 0] = 1
X = X/sd

# apply the sghmc procedure
model = ClassificationModel()
model.fit(X, Y)

#get posterior samples
m, v, ms = model.predict(X)

ms = ms.reshape((ms.shape[0],ms.shape[1]))
np.savetxt('./Results/posteriorsamples_sinatraROC.txt', ms)
np.savetxt('./Results/designmatrix_sinatraROC.txt', X)
