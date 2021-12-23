"""

"""

#%% 

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import geopandas as gpd
import openpyxl
import statistics
from scipy import stats
from numpy import matlib as mb
from shapely.geometry import Point
import gc
import fiona
import time

from fcgadgets.macgyver import utilities_general as gu
from fcgadgets.macgyver import utilities_gis as gis
from fcgadgets.macgyver import utilities_inventory as invu
from fcgadgets.cbrunner import cbrun_utilities as cbu
from fcgadgets.cbrunner import cbrun as cbr

#%% Import data

d=gu.ReadExcel(r'D:\Data\FCI_Projects\SummaryReforestation\Inputs\SummaryAttributes_BySXY.xlsx')

#%% Unique BGC zones/subzones

uBGC=np.unique(np.column_stack((d['BGCz'],d['BGCsz'])),axis=0)
uBGC.shape

aBGC=np.zeros(uBGC.shape[0])
for i in range(uBGC.shape[0]):
    ind=np.where( (d['BGCz']==uBGC[i,0]) & (d['BGCsz']==uBGC[i,1]) & (d['PL SPH']>0) )[0]
    aBGC[i]=aBGC[i]+np.sum(d['Area'][ind])

uBGC=np.column_stack((uBGC,aBGC))

idx=np.flip(np.argsort(aBGC))
uBGC=uBGC[idx,:]


#%%

tv=np.arange(1970,2021,1)

rs={}
for iBGC in range(uBGC.shape[0]):
    rs[iBGC]={}
    rs[iBGC]['A']=np.zeros(tv.size)
    rs[iBGC]['PD']=np.zeros(tv.size)

for iBGC in range(uBGC.shape[0]):    
    for iT in range(tv.size):
        ind=np.where( (d['BGCz']==uBGC[iBGC,0]) & (d['BGCsz']==uBGC[iBGC,1]) & (d['PL SPH']>0) & (np.floor(d['Year'])==tv[iT]) )[0]    
        rs[iBGC]['A'][iT]=rs[iBGC]['A'][iT]+np.sum(d['Area'][ind])
        rs[iBGC]['PD'][iT]=rs[iBGC]['PD'][iT]+np.sum(d['PL SPH'][ind]*d['Area'][ind])/np.sum(d['Area'][ind])

iMIRP=np.where( (tv>=2006) & (tv<=2015) )[0]
iPost=np.where( (tv>=2016) )[0]

for iBGC in range(uBGC.shape[0]): 
    rs[iBGC]['MIRP']=np.nanmean(rs[iBGC]['PD'][iMIRP])    
    rs[iBGC]['PD Post']=np.nanmean(rs[iBGC]['PD'][iPost])
    rs[iBGC]['PD Post Inc']=rs[iBGC]['PD Post']-rs[iBGC]['MIRP']
    rs[iBGC]['A Post']=np.nanmean(rs[iBGC]['A'][iPost])    
    rs[iBGC]['NT Post Inc']=(rs[iBGC]['PD Post']-rs[iBGC]['MIRP'])*rs[iBGC]['A Post']

plt.close('all')
for iBGC in range(30):
    plt.bar(iBGC,rs[iBGC]['A Post'],1,fc=[0.8,0.8,0.8])

plt.close('all')
for iBGC in range(30):
    plt.bar(iBGC,rs[iBGC]['MIRP'],1,fc=[0.8,0.8,0.8])    
    plt.plot(iBGC,rs[iBGC]['PD Post'],'ko',mfc='w')


y=np.zeros(uBGC.shape[0])
for iBGC in range(30):
    y[iBGC]=rs[iBGC]['NT Post Inc']

np.sum(y)

