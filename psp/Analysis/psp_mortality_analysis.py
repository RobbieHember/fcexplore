
#%% Import modules

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy.io as spio
import xgboost as xgb
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import train_test_split
import statsmodels.api as sm
import statsmodels.formula.api as smf
import fcgadgets.macgyver.utilities_general as gu
from fcgadgets.psp import psp_utilities as utl

#%% Import data

d=gu.ipickle(meta['Paths']['DB'] + '\\Processed\\L2\\L2_BC.pkl')
tl=d['tobs'].copy()

#%% Filter

ikp=np.where( (tl['Mortality']>=0) & \
             (tl['ID Species']==meta['LUT']['Species']['BL']) & \
             (tl['DA t0']!=meta['LUT']['Damage Agents']['Insect']) & (tl['DA t1']!=meta['LUT']['Damage Agents']['Insect']) & \
             (tl['DBH t0']>=0) & (tl['H t0']>0) & (tl['Stand Age t0']>=0) & \
             (tl['Stand Cag L Larger t0']>=0) )[0]

#%% Set variables

vL=['Mortality','Delta t','H t0','Stand BA L Larger t0','Stand Age t0'] #,'DBH t0'
d={}
for v in vL:
    d[v]=tl[v][ikp]
#d['Mortality']=d['Mortality']/tl['Delta t'][ikp]

d['H t0^2']=d['H t0']**2
#d['DBH t0^2']=d['DBH t0']**2

df=pd.DataFrame(d)

#%%

X,y=df.iloc[:,1:],df.iloc[:,0]

ddmat=xgb.DMatrix(data=X,label=y)

X_train,X_test,y_train,y_test=train_test_split(X,y,test_size=0.2,random_state=123)

xg_reg=xgb.XGBRegressor(objective='reg:linear',colsample_bytree=0.3,learning_rate=0.1,max_depth=5,alpha=10,n_estimators=10)

xg_reg.fit(X_train,y_train)

preds=xg_reg.predict(X_test)

rmse=np.sqrt(mean_squared_error(y_test,preds))
print("RMSE: %f" % (rmse))

plt.close('all')
xgb.plot_importance(xg_reg)

#%%

x=X_test.copy()
for c in x:
    x[c]=np.nanmean(x[c])
x['Stand Age t0']=np.linspace(0,150,x.shape[0])
preds=xg_reg.predict(x)

plt.close('all')
plt.plot(x['Stand Age t0'],preds,'b.')

#%%

x=X_test.copy()
for c in x:
    x[c]=np.mean(x[c])
x['H t0^2']=np.linspace(0,35**2,x.shape[0])
preds=xg_reg.predict(x)

plt.close('all')
plt.plot(x['H t0^2'],preds,'b.')

#%%

x=X_test.copy()
for c in x:
    x[c]=np.mean(x[c])
x['H t0']=np.linspace(0,35,x.shape[0])
preds=xg_reg.predict(x)

plt.close('all')
plt.plot(x['H t0'],preds,'b.')

#%%

ikp=np.where( (tl['ID Species']==meta['LUT']['Species']['PL']) )[0]
x=tl['Stand Cag L Larger t0'][ikp]
y=tl['Mortality'][ikp]

bw=5; bin=np.arange(0,150,bw)
N,mu,sig,se=gu.discres(x,y,bw,bin)

plt.close('all')
plt.plot(bin,mu,'bo')

#%%

ikp=np.where( (tl['ID Species']==meta['LUT']['Species']['PL']) & (tl['DA t1']==meta['LUT']['Damage Agents']['Insect']) )[0]
x=tl['Stand Age t0'][ikp]
y=tl['Mortality'][ikp]
bw=1; bin=np.arange(0,50,bw)
N,mu,sig,se=gu.discres(x,y,bw,bin)
plt.plot(bin,mu,'rs')


#%%

da={}
u=np.unique(tl['DA t0'])
for i in range(u.size):
    da[i]=np.zeros(tl['DA t0'].size)
    ind=np.where( (tl['DA t0']==u[i]) )[0]
    da[i][ind]=1

ikp=np.where( (tl['DBH t0']>=0) & (tl['ID Species']==meta['LUT']['Species']['BL']) )[0]
x=tl['H t0'][ikp]
y=da[2][ikp]

bw=1; bin=np.arange(0,50,bw)
N,mu,sig,se=gu.discres(x,y,bw,bin)

plt.close('all')
plt.plot(bin,mu,'bo')

#%%

ikp=np.where(tl['ID Species']==meta['LUT']['Species']['FD'])[0]
x=tl['Csw t0'][ikp]
y=tl['Csw G'][ikp]

bw=1; bin=np.arange(0,50,bw)
N,mu,sig,se=gu.discres(x,y,bw,bin)

plt.close('all')
plt.plot(bin,mu,'bo')



