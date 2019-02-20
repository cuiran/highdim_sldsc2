import numpy as np
import sklearn
from sklearn import linear_model
import keras.backend as K
from keras.callbacks import Callback
from keras.models import Sequential
from keras.models import Model
from keras.layers.core import Dense
from keras import optimizers
from keras import regularizers

cov = np.array([[1,0.9,0,0,0.8],[0.9,1,0,0,0.7],[0,0,1,0.9,0.2],[0,0,0.9,1,0.1],[0.8,0.7,0.2,0.1,1]])
mu = np.zeros(5)
X = np.random.multivariate_normal(mu,cov,10000)
beta = [2,1.8,1,0,3]
epsilon = np.random.normal(0,1,10000)
y = X.dot(beta)+epsilon

def sklearn_fit(X,y):
    sklasso = linear_model.LassoCV(fit_intercept=False)
    sklasso.fit(X,y)
    return sklasso.coef_,sklasso.alpha_

def keras_lasso(X,y,rp,lr_hat):
    # rp is regularization parameter, here we use the alpha generated by sklearn LassoCV
    print('performing Lasso regression using keras')
    model = Sequential()
    p = X.shape[1]
    model.add(Dense(1,input_dim=p,use_bias=False,kernel_regularizer=regularizers.l1(lr_hat)))
    sgd = optimizers.SGD(lr=lr_hat)
    model.compile(loss='mse',optimizer=sgd)
    batch_size = 1000
    model.fit(X,y,batch_size=batch_size,epochs=800,verbose=1)
    coefs = model.get_weights()[0]
    return coefs.reshape((coefs.shape[0],))

def run(X,y):
    lassocoef,lassoalpha = sklearn_fit(X,y)
    kerascoef = keras_lasso(X,y,lassoalpha,0.001)
    print('LassoCV coef',lassocoef,'LassoCV alpha',lassoalpha)
    print('keras coef',kerascoef)
    print('true coef',beta)
    return

if __name__=='__main__':
    run(X,y)
