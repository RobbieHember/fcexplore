#%% Import modules
import os
import numpy as np
import time
import matplotlib.pyplot as plt
import geopandas as gpd
import pandas as pd
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

#%% Burn severity analysis
zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])    

vList=['height','ibm_yr','bsr_yr','bsr_sc','harv_yr_comp1','age_vri','tdc'] 
z0=u1ha.Import_Raster(meta,[],vList,'Extract Grid')

#plt.matshow(z0['bsr_sc'],clim=[0,5])
#plt.matshow(z0['twi'],clim=[0,400])

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