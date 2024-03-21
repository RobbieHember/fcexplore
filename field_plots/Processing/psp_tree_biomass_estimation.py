
"""
TREE BIOMASS ESTIMATION
"""

#%% IMPORT MODULES

import sys
import numpy as np
import gc as garc
from osgeo import gdal
from osgeo import osr
import time
import matplotlib.pyplot as plt
import geopandas as gpd
import pandas as pd
import fiona
import rasterio
import pyproj
from rasterio import features
from shapely.geometry import Point, Polygon,box,shape
from shapely import geometry
from rasterio.transform import from_origin
from osgeo import osr
from scipy.interpolate import griddata
import statsmodels.formula.api as smf
from scipy.io import loadmat
import pyproj

import fcgadgets.macgyver.utilities_general as gu
import fcgadgets.macgyver.utilities_gis as gis
from fcexplore.psp.Processing import psp_utilities as utl

#%% Set figure properties
gp=gu.SetGraphics('Manuscript')

#%% Import data
meta={}
meta['Paths']={}
meta['Paths']['DB']=r'C:\Users\rhember\Documents\Data\GroundPlots\PSP-NADB2'
meta=utl.ImportParameters(meta)
d=gu.ipickle(meta['Paths']['DB'] + '\\Processed\\L2\\L2_BC.pkl')
sl=d['sobs']
tl=d['tobs']
del d

#%% Descriptive stats by species

u=np.unique(tl['ID Species'][np.isnan(tl['ID Species'])==False])

ds={}
ds['ID Species']=np.zeros(u.size)
ds['Name Species']=np.array(['' for _ in range(u.size)],dtype=object)
ds['N']=np.zeros(u.size)
ds['N DBH t0 Missing']=np.zeros(u.size)
ds['N H t0 Missing']=np.zeros(u.size)
ds['N Vws t0 Missing']=np.zeros(u.size)
for iS in range(u.size):
    ind=np.where( (tl['ID Species']==u[iS]) & (tl['Vital Status t1']==1) )[0]
    ds['ID Species'][iS]=u[iS]
    ind_allo=np.where(meta['Allo B']['ID']==u[iS])[0]
    if ind_allo.size>0:
        ds['Name Species'][iS]=meta['Allo B']['Common Name'][ind_allo[0]]
    ds['N'][iS]=ind.size
    ind1=np.where( (tl['ID Species']==u[iS]) & (tl['Vital Status t1']==1) & (tl['DBH t0']>0) & (tl['DBH t0']<1000) )[0]
    ds['N DBH t0 Missing'][iS]=ind.size-ind1.size
    ind1=np.where( (tl['ID Species']==u[iS]) & (tl['Vital Status t1']==1) & (tl['H t0']>0) & (tl['H t0']<1000) )[0]
    ds['N H t0 Missing'][iS]=ind.size-ind1.size
    ind1=np.where( (tl['ID Species']==u[iS]) & (tl['Vital Status t1']==1) & (tl['Vws t0']>0) & (tl['Vws t0']<1000) )[0]
    ds['N Vws t0 Missing'][iS]=ind.size-ind1.size

df=pd.DataFrame(ds)
df.to_excel(r'C:\Users\rhember\Documents\Data\GroundPlots\PSP-NADB2\Processed\Descriptive Stats TL.xlsx')

#%% QA by species

u=np.unique(tl['ID Species'][np.isnan(tl['ID Species'])==False])

for iS in range(u.size):

    #spc='SW'
    #meta['LUT']['Species'][spc]
    #cd=utl.lut_id2cd(meta,'Species',u[iS])
    ind_allo=np.where(meta['Allo B']['ID']==u[iS])[0]

    #indD=np.where( (tl['ID Species']==u[iS]) & (tl['Vital Status t1']!=meta['LUT']['Vital Status']['Live']) )[0]
    indL=np.where( (tl['ID Species']==u[iS]) & (tl['Vital Status t0']==meta['LUT']['Vital Status']['Live']) & (tl['DBH t0']>0) & (tl['H t0']>0) )[0]

    plt.close('all')
    fig,ax=plt.subplots(2,3,figsize=gu.cm2inch(15,6))
    ms=2.5

    #ax[0,0].plot(tl['DBH'][indD],tl['H'][indD],'r.')
    ax[0,0].plot(tl['DBH t0'][indL],tl['H t0'][indL],'b.',ms=ms)
    ax[0,0].set(xlabel='DBH',ylabel='H')
    ax[0,1].plot(tl['DBH t0'][indL],tl['Vws t0'][indL],'b.',ms=ms)
    ax[0,1].set(xlabel='DBH',ylabel='Volume whole stem')
    ax[0,2].plot(tl['Vws t0'][indL],tl['Csw t0'][indL],'b.',ms=ms)
    ax[0,2].set(xlabel='Volume',ylabel='C stemwood')
    ax[1,0].plot(tl['Vws t0'][indL],tl['Cbr t0'][indL],'b.',ms=ms)
    ax[1,0].set(xlabel='Volume',ylabel='C branch')
    ax[1,1].plot(tl['Vws t0'][indL],tl['Cbk t0'][indL],'b.',ms=ms)
    ax[1,1].set(xlabel='Volume',ylabel='C bark')
    ax[1,2].plot(tl['Vws t0'][indL],tl['Cf t0'][indL],'b.',ms=ms)
    ax[1,2].set(xlabel='Volume',ylabel='C foliage')

    gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Biomass\Ground Plot QA\2022-11\qa_' + meta['Allo B']['Common Name'][ind_allo[0]],'png',150)


#%% Branch biomass

ind=np.where( (tl['ID Species']==meta['LUT']['Species']['PL']) & (tl['Vital Status']==meta['LUT']['Vital Status']['Live']) )[0]
th=np.nanpercentile(tl['Cbr'][ind],99.75)
print(th)
#np.nanmean(tl['Cbr'][ind])
#3*np.nanstd(tl['Cbr'][ind])

plt.hist(tl['Cbr'][ind],np.arange(0,200,10))

bwD=2; binD=np.arange(0,100,bwD)
bwH=1; binH=np.arange(0,55,bwH)
z=np.zeros((binH.size,binD.size))
for i in range(binH.size):
    for j in range(binD.size):
        z[i,j]=0.5*(0.0285*binD[j]**3.3764*binH[i]**-1.4395)
plt.close('all')
plt.matshow(np.flip(z,1).T,extent=[binD[0],binD[-1],binH[0],binH[-1]],clim=[0,200])

#%%

y=tl['Cbr t0']
spc='SS'

ind=np.where( (tl['ID Species']==meta['LUT']['Species'][spc]) & (tl['Vital Status t0']==meta['LUT']['Vital Status']['Live']) )[0]
th=np.nanpercentile(y[ind],99.75)
print(th)

ms=6
plt.close('all')
ind=np.where( (tl['Vital Status t0']==meta['LUT']['Vital Status']['Live']) )[0]
plt.plot(tl['Vws t0'][ind],y[ind],'k.',mfc=[0.75,0.75,0.75],ms=ms,mec='none')

ind=np.where( (tl['ID Species']==meta['LUT']['Species'][spc]) & (tl['Vital Status t0']==meta['LUT']['Vital Status']['Live']) )[0]
plt.plot(tl['Vws t0'][ind],y[ind],'b.',ms=ms)

ind=np.where( (tl['ID Species']==meta['LUT']['Species'][spc]) & (tl['Vital Status t0']==meta['LUT']['Vital Status']['Live']) & (y<th) )[0]
plt.plot(tl['Vws t0'][ind],y[ind],'g.',ms=ms)

#%%

ind=np.where( (tl['ID Species']==meta['LUT']['Species']['PL']) )[0]
plt.close('all')
plt.plot(tl['H'][ind],tl['Cbr'][ind],'b.')

#%%

ind=np.where( (tl['ID Species']==meta['LUT']['Species']['PL']) & (tl['Cbr']>1000) )[0]
tl['H'][ind]


