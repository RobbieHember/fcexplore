#%% Import modules
import os
import numpy as np
import time
import matplotlib.pyplot as plt
import geopandas as gpd
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
import scipy.io as spio
import cv2
import copy
from sklearn.metrics import accuracy_score
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn import datasets
#from rasterio.transform import from_origin
#from shapely.geometry import Point,Polygon
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.macgyver.util_query_gdb as qgdb
import fcgadgets.bc1ha.bc1ha_util as u1ha
import fcgadgets.cbrunner.cbrun_util as cbu
gp=gu.SetGraphics('Manuscript')

#%% Import paths and look-up-tables
meta=u1ha.Init()

#%% Import data
vList=['refg','lc_comp1_2019','bsr_yr','bsr_sc','h_vri15','age_vri15','tdc_vri15','spc1_vri15','prcp_ann_n'] # ,'ibm_yr','harv_yr_comp1','geomorph','pdead_vri23'
z0=u1ha.Import_Raster(meta,[],vList,'Extract Grid')
#plt.matshow(z0['bsr_sc'],clim=[0,5])
#plt.matshow(z0['twi'],clim=[0,400])

#%% Regression
ikp=np.where( (z0['lc_comp1_2019']==meta['LUT']['Derived']['lc_comp1']['Forest']) & (z0['bsr_yr']>0) & (z0['bsr_sc']>=1) & (z0['bsr_sc']<=4) )
#plt.hist(z0['bsr_sc'][ikp][0::1000])
ivl=50
z1={}
for k in z0.keys():
	z1[k]=z0[k][ikp][0::ivl]

z1['tdc2']=(z1['tdc_vri15']==2).astype(int)
z1['tdc3']=(z1['tdc_vri15']==3).astype(int)

z1['SW']=(z1['spc1_vri15']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['SW']).astype(int)
z1['PL']=(z1['spc1_vri15']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['PL']).astype(int)
z1['BL']=(z1['spc1_vri15']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['BL']).astype(int)
z1['AT']=(z1['spc1_vri15']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['AT']).astype(int)
z1['EP']=(z1['spc1_vri15']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['EP']).astype(int)
#z1['AC']=(z1['spc1_vri15']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['AC']) | (z1['spc1_vri15']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['ACT']).astype(int)
z1['AC']=(z1['spc1_vri15']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['AC']).astype(int)

#%%
y=z1['bsr_sc']
#x=np.column_stack([z1['age_vri'],z1['height'],z1['tdc2'],z1['tdc3'],z1['gm2'],z1['gm3'],z1['gm5'],z1['gm6'],z1['gm7'],z1['gm9'],z1['gm10'] ])
x=np.column_stack([z1['age_vri15'],z1['h_vri15'],z1['prcp_ann_n'],z1['tdc2'],z1['tdc3'],z1['AT'],z1['BL'] ])
#x,mu,sig=gu.zscore(x)
x=sm.add_constant(x)

#%%
md=sm.MNLogit(y,x)
rs=md.fit()
print(rs.summary())

# Indices
iA=1
iH=2
iP=3
iD2=4
iD3=5
iAT=6
iBL=7

#%% Plot tree density class
plt.close('all')
fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(22,8));

x1=np.ones(x.shape[1])
x1[iA]=50
x1[iH]=15
x1[iP]=500
x1[iD2]=0
x1[iD3]=0
x1[iAT]=0
x1[iBL]=0
yhat0=rs.predict(x1)
ax[0].plot(yhat0[0,:],'-bo')

x1=np.ones(x.shape[1])
x1[iA]=50
x1[iH]=15
x1[iP]=2000
x1[iD2]=1
x1[iD3]=0
x1[iAT]=0
x1[iBL]=0
yhat1=rs.predict(x1)
ax[0].plot(yhat1[0,:],'-rs')

ax[1].plot(yhat1[0,:]-yhat0[0,:],'ob-')

#%% Plot tree density class
plt.close('all')
fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(22,8));

x1=np.ones(x.shape[1])
x1[iA]=50
x1[iH]=15
x1[iP]=800
#x1[iPD]=5
x1[iD2]=0
x1[iD3]=0
x1[iAT]=0
x1[iBL]=0
yhat0=rs.predict(x1)
ax[0].plot(yhat0[0,:],'-bo')

x1=np.ones(x.shape[1])
x1[iA]=50
x1[iH]=15
x1[iP]=800
#x1[iPD]=5
x1[iD2]=1
x1[iD3]=0
x1[iAT]=0
x1[iBL]=0
yhat1=rs.predict(x1)
ax[0].plot(yhat1[0,:],'-gs')

x1=np.ones(x.shape[1])
x1[iA]=50
x1[iH]=15
x1[iP]=800
#x1[iPD]=5
x1[iD2]=0
x1[iD3]=1
x1[iAT]=0
x1[iBL]=0
yhat2=rs.predict(x1)
ax[0].plot(yhat2[0,:],'-r^')

ax[1].plot(yhat2[0,:]-yhat0[0,:],'ob-')

#%% Plot species
iA=1;
iH=2;
iD2=3;
iD3=4
iAT=5
iBL=6

x1=np.ones(x.shape[1])
x1[iA]=50
x1[iH]=15
#x1[iPD]=5
x1[iD2]=0
x1[iD3]=0
x1[iAT]=0
x1[iBL]=1
yhat0=rs.predict(x1)
plt.close('all')
fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(22,8));
ax[0].plot(yhat0[0,:],'or-')

x1=np.ones(x.shape[1])
x1[iA]=50
x1[iH]=15
#x1[iPD]=5
x1[iD2]=0
x1[iD3]=0
x1[iAT]=1
x1[iBL]=0
yhat1=rs.predict(x1)
ax[0].plot(yhat1[0,:],'og--')

ax[1].plot(yhat1[0,:]-yhat0[0,:],'ob-')


#%%


#%% Compare distributions from national and provincial datasets
plt.close('all'); fig,ax=plt.subplots(2,1)
ikp=np.where( (z0['bsr_sc']>0) & (z0['bsr_sc']<5) & (z0['bsr_yr']<2016) )
ax[0].hist(z0['bsr_sc'][ikp][0::30])
ikp=np.where( (z0['bsr_sc']>0) & (z0['bsr_sc']<5) & (z0['bsr_yr']>=2017) )
ax[1].hist(z0['bsr_sc'][ikp][0::30])
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\BurnSeverityDistributionComparison','png',900)

#%% Analysis of variability
z1=copy.deepcopy(z0)
# meta['LUT']['Derived']['burnsev_comp1']
ikp=np.where( (z1['bsr_sc']>0) & (z1['bsr_sc']<5) & (z1['bsr_yr']>=2017) )
for k in z1.keys():
    z1[k]=z1[k][ikp]
z1['bm']=np.zeros(z1['bsr_yr'].size)
z1['bm'][np.where(z1['bsr_sc']>=3)[0]]=1
z1['ub']=np.zeros(z1['bsr_yr'].size)
z1['ub'][np.where(z1['bsr_sc']==1)[0]]=1

#ikp=np.where( z1['bgcz']==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['SBS'] )[0]

#%% Age
plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(7.8,6))
bw=10; bin=np.arange(5,265,bw)
N,mu1,med,sig,se=gu.discres(z1['age_vri']-(2022-z1['bsr_yr']),z1['bm'],bw,bin)
N,mu2,med,sig,se=gu.discres(z1['age_vri']-(2022-z1['bsr_yr']),z1['ub'],bw,bin)
ax.plot(bin,mu1*100,'-ro',ms=3,label='Moderate and high severity')
ax.plot(bin,(1-mu1-mu2)*100,'-ys',ms=3,label='Low severity')
ax.plot(bin,mu2*100,'-g^',ms=3,label='Unburned')
ax.set(xlabel='Stand age, years',ylabel='Frequency (%)',ylim=[0,100])
ax.legend(loc='upper left',frameon=False,fontsize=7)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=1.5)
#gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\BurnSeverity\BSR_Vs_AgeVRI','png',900)

#%% Height
plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(7.8,6))
bw=1; bin=np.arange(1,50,bw)
N,mu1,med,sig,se=gu.discres(z1['height']-(2022-z1['bsr_yr']),z1['bm'],bw,bin)
N,mu2,med,sig,se=gu.discres(z1['height']-(2022-z1['bsr_yr']),z1['ub'],bw,bin)
ax.plot(bin,mu1*100,'-ro',ms=3,label='Moderate and high severity')
ax.plot(bin,(1-mu1-mu2)*100,'-ys',ms=3,label='Low severity')
ax.plot(bin,mu2*100,'-g^',ms=3,label='Unburned')
ax.set(xlabel='Height (m)',ylabel='Frequency (%)',ylim=[0,100])
ax.legend(loc='upper left',frameon=False,fontsize=7)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=1.5)
#gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\BurnSeverity\BSR_Vs_AgeVRI','png',900)

#%% 
bw=5; bin=np.arange(2.5,50,bw)
plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(7.8,6))
N,mu,med,sig,se=gu.discres(z1['bsr_yr']-z1['harv_yr_comp1'],z1['bm'],bw,bin)
ax.plot(bin,mu*100,'-ro',ms=3,label='Moderate and high severity')
N,mu,med,sig,se=gu.discres(z1['bsr_yr']-z1['harv_yr_comp1'],z1['ub'],bw,bin)
ax.plot(bin,mu*100,'-bs',ms=3,label='Unburned')
ax.set(xlabel='Time since harvest, years',ylabel='Frequency (%)',ylim=[0,100],xlim=[0,50])
ax.legend(loc='upper left',frameon=False,fontsize=7)

ind=np.where(z1['harv_yr_comp1']==0)
np.mean(z1['bm'][ind])
ind=np.where( (z1['harv_yr_comp1']>0) & (z1['harv_yr_comp1']<z1['bsr_yr']) )
np.mean(z1['bm'][ind])

#%%
bw=2; bin=np.arange(0,50,bw)
plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(7.8,6))
N,mu,med,sig,se=gu.discres(z1['bsr_yr']-z1['ibm_yr'],z1['bm'],bw,bin)
ax.plot(bin,mu*100,'-ro',ms=3,label='Moderate and high severity')
N,mu,med,sig,se=gu.discres(z1['bsr_yr']-z1['ibm_yr'],z1['ub'],bw,bin)
ax.plot(bin,mu*100,'-bs',ms=3,label='Unburned')
ax.set(xlabel='Time since beetle outbreak, years',ylabel='Frequency (%)',ylim=[0,100],xlim=[0,50])
ax.legend(loc='upper left',frameon=False,fontsize=7)

ind=np.where(z1['ibm_yr']==0)
np.mean(z1['bm'][ind])
ind=np.where(z1['ibm_yr']>0)
np.mean(z1['bm'][ind])

#%%

bw=5; bin=np.arange(2.5,50,bw)
N,mu,med,sig,se=gu.discres(z1['bsr_yr']-z1['harv_yr_comp1'],z1['bm'],bw,bin)
plt.plot(bin,mu,'ko')

bw=5; bin=np.arange(2.5,50,bw)
N,mu,med,sig,se=gu.discres(z1['bsr_yr']-z1['ibm_yr'],z1['bm'],bw,bin)
plt.close('all'); plt.plot(bin,mu,'ko')

bw=1; bin=np.arange(0,5,bw)
N,mu,med,sig,se=gu.discres(z1['tdc'],z1['bm'],bw,bin)
plt.close('all'); plt.plot(bin,mu,'ko')

bw=5; bin=np.arange(1,55,bw)
N,mu,med,sig,se=gu.discres(z1['height'],z1['bm'],bw,bin)
plt.close('all'); plt.plot(bin,mu,'ko')