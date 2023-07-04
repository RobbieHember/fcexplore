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

import fcgadgets.bc1ha.bc1ha_utilities as u1ha
import fcgadgets.macgyver.utilities_general as gu
import fcgadgets.macgyver.utilities_gis as gis

#%% Import paths and look-up-tables

meta=u1ha.Init()
meta=u1ha.ImportLUTs(meta)
gp=gu.SetGraphics('Manuscript')

#%% Import raster datasets

# Land mask
zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])

z={}
z['Ref']=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])['Data']
z['Spc']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Species\\NTEMS_Spc1.tif')['Data']
z['RD']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\TA_REGIONAL_DISTRICTS_SVW\\REGIONAL_DISTRICT_NAME.tif')['Data']
z['Slope']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Terrain\\slope.tif')['Data']
z['Elev']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Terrain\elevation.tif')['Data']
z['T']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Climate\\Seasonal\\BC1ha_tmean_ann_norm_1971to2000_si_hist_v1.tif')['Data'].astype('float')/10
z['P']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Climate\\Seasonal\\BC1ha_prcp_ann_norm_1971to2000_si_hist_v1.tif')['Data']

#%% Filters and sampling

ssi=10

zf={}
for k in z.keys():
    zf[k]=z[k][0::ssi,0::ssi].flatten().astype(float)

#%% Prepare variables

z['Spc hat']=z['Spc'].copy()
z['Prob']=np.zeros(zRef['Data'].shape)

u=np.unique(zf['Spc'][zf['Spc']>0])

for iU in range(u.size):
    print(u[iU])

    flg=0
    if flg==1:
        y=np.zeros(z['Spc'].shape,dtype='int8')
        ind=np.where(z['Spc']==u[iU])
        y[ind]=1
        plt.matshow(y)

    y=np.zeros(zf['Spc'].size)
    ind=np.where(zf['Spc']==u[iU])[0]
    y[ind]=1

    ikp=np.where( (zf['Ref']>0) & (zf['T']>-20) & (zf['P']>=0) & (zf['RD']==meta['LUT']['TA_REGIONAL_DISTRICTS_SVW']['REGIONAL_DISTRICT_NAME']['CAPITAL']) )[0]

    y=y[ikp]
    print(np.sum(y))
    if np.sum(y)<100:
        continue

    x=np.column_stack((zf['Elev'][ikp],zf['T'][ikp],zf['T'][ikp]**2,zf['P'][ikp],zf['T'][ikp]*zf['P'][ikp]))

    # z-score x variables
    xz,mu,sig=gu.zscore(x)

    # Regression analysis

    # Is it balanced
    #print(np.sum(y)/y.size)
    x_train,x_test,y_train,y_test=train_test_split(x,y,test_size=0.5,random_state=0)
    mod=LogisticRegression()#solver='liblinear',random_state=0)
    mod.fit(x_train,y_train)
    yhat=mod.predict(x_test)
    beta=mod.coef_[0]
    intercept=mod.intercept_[0]
    #print(intercept)
    #print(beta)
    score=mod.score(x_test,y_test)
    #print(score)
    cm=metrics.confusion_matrix(y_test,yhat)
    #print(cm)
    Phat=intercept + \
        beta[0]*z['Elev'].astype(float) + \
        beta[1]*z['T'].astype(float) + \
        beta[2]*z['T'].astype(float)**2 + \
        beta[3]*z['P'].astype(float) + \
        beta[4]*z['P']*z['T']
    Phat=(np.exp(Phat)/(1+np.exp(Phat)))*100

    #plt.close('all')
    #plt.matshow(Phat,clim=[0,100])

    iFill=np.where( (z['Spc']==0) & (Phat>z['Prob']) )
    z['Spc hat'][iFill]=u[iU]
    z['Prob'][iFill]=Phat[iFill]

# Save
z1=zRef.copy()
z1['Data']=z['Spc hat']
gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\Species\\NTEMS_Spc1_Filled_RD_CAPITAL.tif')

#%% Site index from SPL

zSpc=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Species\\NTEMS_Spc1_Filled_RD_CAPITAL.tif')
u=np.unique(zSpc['Data'])
zSI=zRef.copy()
for iU in range(u.size):
    if u[iU]==0:
        continue
    ind0=np.where(zSpc['Data']==u[iU])
    ind1=np.where(meta['LUT']['Raw']['spc1_ntems']['Value']==u[iU])[0]
    cd1=meta['LUT']['Raw']['spc1_ntems']['Code'][ind1[0]]
    cd1=cd1[0] + cd1[1].lower()
    try:
        z1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\SPL\\Site_Prod_' + cd1 + '.tif')
    except:
        continue
    zSI['Data'][ind0]=z1['Data'][ind0]

plt.matshow(zSI['Data'],clim=[10,40])
gis.SaveGeoTiff(zSI,meta['Paths']['bc1ha'] + '\\SPL\\SI_From_NTEMS_Spc1_Filled_RD_CAPITAL.tif')

#%% Gap-fill with interpolation

from scipy.interpolate import griddata

zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
zSI=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\SPL\\SI_From_NTEMS_Spc1_Filled_RD_CAPITAL.tif')

ivl=5
iFill=np.where( (zRef['Data']==1) & (zSI['Data']<=3) )
iCal=np.where( (zRef['Data'][0::ivl,0::ivl]==1) & (zSI['Data'][0::ivl,0::ivl]>3) )
xy=np.column_stack([zRef['X'][0::ivl,0::ivl][iCal],zRef['Y'][0::ivl,0::ivl][iCal]])
vals=zSI['Data'][0::ivl,0::ivl][iCal]
zFill=griddata(xy,vals,(zRef['X'][iFill],zRef['Y'][iFill]),method='nearest')

zSI['Data'][iFill]=zFill

plt.close('all')
plt.matshow(zSI['Data'])
gis.SaveGeoTiff(zSI,meta['Paths']['bc1ha'] + '\\SPL\\SI_From_NTEMS_Spc1_Filled_RD_CAPITAL_GF.tif')

#%%



