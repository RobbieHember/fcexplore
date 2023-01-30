#%% Import modules

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.io
import copy
import geopandas as gpd
import fiona
import rasterio
from rasterio import features
import scipy.io as spio
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats, linalg
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, confusion_matrix
from sklearn import metrics

import fcgadgets.macgyver.utilities_general as gu
import fcgadgets.macgyver.utilities_gis as gis

#%% Graphics parameters

gp=gu.SetGraphics('Manuscript')

#%% Import raster datasets

# Land mask
zMask=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\BC_Land_Mask.tif')

zSlope=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Terrain\slope.tif')
#zElev=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Terrain\elevation.tif')
zDTF=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Management\DistanceFromForestryFacility.tif')
zDTR=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Terrain\DistanceFromRoads.tif')

# Harvest mask
zH=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\VEG_CONSOLIDATED_CUT_BLOCKS_SP_All.tif')

#%% Filters and sampling

ssi=50

xMask=zMask['Data'][0::ssi,0::ssi].flatten()
xSlope=zSlope['Data'][0::ssi,0::ssi].flatten().astype(float)
#xElev=zElev['Data'][0::ssi,0::ssi].flatten().astype(float)
xDTR=zDTR['Data'][0::ssi,0::ssi].flatten().astype(float)
xDTF=zDTF['Data'][0::ssi,0::ssi].flatten().astype(float)
xH=zH['Data'][0::ssi,0::ssi].flatten().astype(float)

ikp=np.where(xMask>0)[0]

#%% Prepare variables

x=np.column_stack((xSlope[ikp],xDTR[ikp],xDTF[ikp]))
y=xH[ikp]

# z-score x variables
xz,mu,sig=gu.zscore(x)

#%% Regression analysis

# Is it balanced
print(np.sum(y)/y.size)

x_train,x_test,y_train,y_test=train_test_split(xz,y,test_size=0.5,random_state=0)

mod=LogisticRegression()#solver='liblinear',random_state=0)

mod.fit(x_train,y_train)

yhat=mod.predict(x_test)

beta=mod.coef_[0]
intercept=mod.intercept_[0]
print(intercept)
print(beta)

score=mod.score(x_test,y_test)
print(score)

cm=metrics.confusion_matrix(y_test,yhat)
print(cm)

#%% Predict map of annual probability of harvest (%/yr)

# Duration (Consolidated cutblocks goes from about 1955-2021)
Ivl=2021-1954

zSlope['Data']=np.maximum(0,zSlope['Data'])

zHhat=zMask.copy()
zHhat['Data']=zHhat['Data'].astype('float32')
zHhat['Data']=intercept + \
    (beta[0]*((zSlope['Data'].astype(float)-mu[0])/sig[0])) + \
    (beta[1]*((zDTR['Data'].astype(float)-mu[1])/sig[1])) + \
    (beta[2]*((zDTF['Data'].astype(float)-mu[2])/sig[2]))

zHhat['Data']=(np.exp(zHhat['Data'])/(1+np.exp(zHhat['Data'])))/Ivl*100

zHhat['Data'][zMask['Data']==0]=0

#%%

plt.close('all')
plt.hist(zHhat['Data'][0::ssi,0::ssi].flatten())

#%% Map annual probability of harvest (%/yr)

ssi=5

plt.close('all')
plt.matshow(zHhat['Data'][0::ssi,0::ssi],clim=[0,0.5])
plt.colorbar()

print(np.percentile(100*zHhat['Data'][0::ssi,0::ssi],99.9))

#%% Save Annal probability of harvest (%/yr)

z1=zMask.copy()
z1['Data']=zHhat['Data']*100
z1['Data']=z1['Data'].astype('int16')
gis.SaveGeoTiff(z1,r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\HarvestProbability.tif')