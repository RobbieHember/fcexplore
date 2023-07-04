
#%% Import modules

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib.collections import PatchCollection
import geopandas as gpd
import pandas as pd
import copy
import fiona
import time
from scipy.interpolate import griddata
from shapely.geometry import Polygon,Point,box
import fcgadgets.macgyver.utilities_general as gu
import fcgadgets.macgyver.utilities_gis as gis
from fcgadgets.cbrunner import cbrun_utilities as cbu
from fcgadgets.bc1ha import bc1ha_utilities as u1ha

#%% Path management

meta={}
meta['Paths']={}
meta['Paths']['BC1ha']=r'C:\Users\rhember\Documents\Data\BC1ha'
meta['Paths']['Forest Inventory Disturbances']=r'C:\Users\rhember\Documents\Data\ForestInventory\Disturbances\20220422\Disturbances.gdb'
#meta['Paths']['Figures']=r'C:\Users\rhember\OneDrive - Government of BC\Figures\Reforestation'

#%% Plotting parameters

meta['Graphics']={}
meta['Graphics']['figwidth']=16

#meta['Graphics']['sidespace']=0.25
meta['Graphics']['sidespace']=0

meta['Graphics']['ax1 pos']=[0,0,1-meta['Graphics']['sidespace']-0.01,1]
meta['Graphics']['ax1 vis']='off'
meta['Graphics']['ax1 gridvis']=False
meta['Graphics']['ax2 pos']=[1-meta['Graphics']['sidespace']+meta['Graphics']['ax1 pos'][0]+0.01,0.6,0.03,0.35]
meta['Graphics']['ax2 pos long']=[1-meta['Graphics']['sidespace']+meta['Graphics']['ax1 pos'][0]+0.01,0.1,0.03,0.8]

gp=gu.SetGraphics('Manuscript')

meta['LUT']=u1ha.Import_BC1ha_LUTs()

#%% Import data

#gdf=u1ha.Import_GDBs_ProvinceWide()

zRef=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\BC_Land_Mask.tif')

zTSA=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\TSA.tif')
zLCC1=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\LandUseLandCover\LandCoverClass1.tif')
zProt=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\LandUseLandCover\PROTECTED_LANDS_DESIGNATION.tif')
#zFCR=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Results\RSLT_FOREST_COVER_RESERVE_SVW.tif')
zRET=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Results\SILV_RESERVE_CODE_Consolidated.tif')
zPL=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Results\Planting_FromRESULTS_MaskCount.tif')
#zST=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Results\stocktype.tif')
zHRT=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Results\Harvest_Regen_Type.tif')
zFYL=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\PROT_HISTORICAL_FIRE_POLYS_SP_YearLast.tif')
zHCC=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\VEG_CONSOLIDATED_CUT_BLOCKS_SP_MaskAll.tif')
zHYL=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\VEG_CONSOLIDATED_CUT_BLOCKS_SP_YearLast.tif')
zAgeVRI=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\age1.tif')

zVLiveVRI=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\v_live.tif')
zVLiveVRI['Data']=zVLiveVRI['Data'].astype('int16')

# Stemwood biomass from PFI
zB=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\ForestInventory\LiDAR\Boundary TSA\From GQ 20230302\BoundaryTSA_PFI_LnTotalWSV_bc1ha.tif')
zB['Data']=0.5*0.5*zB['Data']

dTIPSY=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\TIPSY_Boundary.xlsx')

#%% Filter non-study area

ind=np.where( (zTSA['Data']==meta['LUT']['tsa']['Boundary TSA']) & (zB['Data']>=0)  ) # & (zLCC1['Data']==1)

d={}
d['B']=zB['Data'][ind]
d['LCC1']=zLCC1['Data'][ind]
d['RET']=zRET['Data'][ind]
d['PL']=zPL['Data'][ind]
d['HCC']=zHCC['Data'][ind]
d['HYL']=zHYL['Data'][ind]
d['HRT']=zHRT['Data'][ind]
d['FYL']=zFYL['Data'][ind]
d['ASH']=2020-d['HYL']
d['ASF']=2020-d['FYL']
d['AgeVRI']=zAgeVRI['Data'][ind]
d['VLiveVRI']=zVLiveVRI['Data'][ind].astype('float')

#%% Age

plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(8,7))
ax.fill_between([0,6],-5,125,color=[1,0.875,0.775],lw=0,label='Harvest incomplete')
ax.fill_between([6,2020-1987],-5,125,color=[1,0.925,0.875],lw=0,label='FRPA')
ax.fill_between([2020-1987,100],-5,125,color=[1,0.975,0.95],lw=0,label='Pre-FRPA')

ind=np.where( (d['HCC']>=0) )[0]
bw=5; bin=np.arange(1,200,bw)
N,mu,med,sig,se=gu.discres(d['ASH'][ind],d['B'][ind],bw,bin)
ind=np.where(N>1000)[0]
ax.plot(bin[ind],mu[ind],'ko',ms=3,label='Predictive Forest Inventory\nchronosequence')
ax.errorbar(bin[ind],mu[ind],yerr=[3*se[ind],3*se[ind]],color='k',ls='',capsize=1.5,lw=0.25)

ax.plot(dTIPSY['Age'],dTIPSY['Tot C'],'k-',color=[0.27,0.49,0.79],lw=2,label='TIPSY (SI = 20m)')
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=2)
ax.set(xlabel='Age, years',ylabel='Stemwood carbon (MgC/ha)',xlim=[0,260],ylim=[0,110])

leg=ax.legend(loc='upper right',frameon=True,fontsize=7)
frame=leg.get_frame()
frame.set_facecolor([0.9,0.9,0.9])
frame.set_linewidth(0)

ind=np.where( (d['FYL']>0) )[0]
N,mu,med,sig,se=gu.discres(2020-d['FYL'][ind],d['B'][ind],bw,bin)
ind=np.where(N>100)[0]
ax.plot(bin[ind],mu[ind],'go',mec=[0.6,1,0],mfc=[0.6,1,0],ms=3,label='Predictive Forest Inventory\nchronosequence')

# ind=np.where( (d['HCC']==1) & (d['RET']!=6) )[0]
# bw=1; bin=np.arange(1,200,bw)
# N,mu,med,sig,se=gu.discres(d['ASH'][ind],d['B'][ind],bw,bin)
# ind=np.where(N>1000)[0]
# ax.plot(bin[ind],mu[ind],'gs',ms=3,label='Predictive Forest Inventory\nchronosequence')

# N,mu,med,sig,se=gu.discres(d['ASH'][ind],d['VLiveVRI'][ind],bw,bin)
# ind=np.where(N>1000)[0]
# ax.plot(bin[ind],mu[ind],'rs',ms=3,label='Predictive Forest Inventory\nchronosequence')

#%% Age VRI

plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(8,7))
ind=np.where( (d['HCC']>=0) )[0]
bw=10; bin=np.arange(1,200,bw)
N,mu,med,sig,se=gu.discres(d['AgeVRI'][ind],d['B'][ind],bw,bin)
ind=np.where(N>1000)[0]
ax.plot(bin[ind],mu[ind],'go',mec=[0.6,1,0],mfc=[0.6,1,0],ms=3,label='Predictive Forest Inventory\nchronosequence')
ax.errorbar(bin[ind],mu[ind],yerr=[3*se[ind],3*se[ind]],color='k',ls='',capsize=1.5,lw=0.25)

#%% Age

plt.close('all')
fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(15,8))
ind=np.where( (d['HCC']==1) & (d['ASH']>=7) & (d['ASH']<=30) )[0]
ax[0].hist(d['B'][ind],np.arange(0,220,10))
ind=np.where( (d['HCC']==0) )[0]
ax[1].hist(d['B'][ind],np.arange(0,220,10))

plt.close('all')
ind=np.where(d['HCC']==1)[0]
bw=5; bin=np.arange(5,200,bw)
N,mu,med,sig,se=gu.discres(d['ASH'][ind],d['B'][ind],bw,bin)
plt.plot(bin,mu,'bo',ms=5)

N,mu,med,sig,se=gu.discres(d['ASF'][ind],d['B'][ind],bw,bin)
plt.plot(bin,mu,'rs',ms=5)

#%% Effect of LCC1

u=np.unique(d['LCC1'])
u=u[u!=0]
mu=np.zeros(u.size)
err=np.zeros(u.size)
for i in range(u.size):
    ind=np.where( (d['LCC1']==u[i]) & (d['HCC']>=1) )[0]
    mu[i]=np.mean(d['B'][ind])
    err[i]=2*np.std(d['B'][ind])/np.sqrt(ind.size)
plt.close('all')
plt.bar(u,mu)

#%% Effect of LCC

u=np.unique(d['RET'])
u=u[u!=0]
mu=np.zeros(u.size)
err=np.zeros(u.size)
for i in range(u.size):
    ind=np.where( (d['ASH']>=7) & (d['ASH']<=30) & (d['HCC']==1) & (d['RET']==u[i]) )[0]
    mu[i]=np.mean(d['B'][ind])
    err[i]=2*np.std(d['B'][ind])/np.sqrt(ind.size)
plt.close('all')
plt.bar(u,(mu-mu[5])/mu[5]*100)

#%% Regen Type

u=np.unique(d['HRT'])
u=u[u!=0]
mu=np.zeros(u.size)
err=np.zeros(u.size)
for i in range(u.size):
    ind=np.where( (d['ASH']>=7) & (d['ASH']<=30) & (d['HCC']==1) & (d['HRT']==u[i]) )[0]
    mu[i]=np.mean(d['B'][ind])
    err[i]=2*np.std(d['B'][ind])/np.sqrt(ind.size)
plt.close('all')
plt.bar(u,(mu-mu[0])/mu[0]*100)

#%% Time series

print(np.mean(B))

tv=np.arange(1950,2022,1)
A_h=np.zeros(tv.size)
y_h=np.zeros(tv.size)
y_fcr=np.zeros(tv.size)
y_pl=np.zeros(tv.size)
y_npl=np.zeros(tv.size)
for iT in range(tv.size):
    zHy=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\VEG_CONSOLIDATED_CUT_BLOCKS_SP_' + str(tv[iT]) + '.tif')
    #zHy=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\Harvest_FromRESULTS_' + str(tv[iT]) + '.tif')
    hy=zHy['Data'][ind0]
    ind=np.where( (hy>0) & (fcr==0) )
    A_h[iT]=ind[0].size
    y_h[iT]=np.mean(B[ind])
    ind=np.where( (hy>0) & (fcr==0) & (pl==1) )
    y_pl[iT]=np.mean(B[ind])
    ind=np.where( (hy>0) & (fcr==0) & (pl==0) )
    y_npl[iT]=np.mean(B[ind])
    ind=np.where( (hy>0) & (fcr>0) )
    y_fcr[iT]=np.mean(B[ind])
    print(tv[iT])

#%%

plt.close('all')
plt.plot(tv,y_h,'-bo')
#plt.plot(tv,y_npl,'-bo')
#plt.plot(tv,y_pl,'-gs')

#%% Age



#%%

ind0=np.where( (zTSA['Data']==meta['LUT']['tsa']['Boundary TSA']) & (zB['Data']>=0) & (zProt['Data']>0)  ) # & (zLCC1['Data']==1)
B=zB['Data'][ind0]
print(np.mean(B))

ind0=np.where( (zTSA['Data']==meta['LUT']['tsa']['Boundary TSA']) & (zB['Data']>=0) & (zProt['Data']==0)  ) # & (zLCC1['Data']==1)
B=zB['Data'][ind0]
print(np.mean(B))

ind0=np.where( (zTSA['Data']==meta['LUT']['tsa']['Boundary TSA']) & (zB['Data']>=0) & (zFCR['Data']==1)  ) # & (zLCC1['Data']==1)
B=zB['Data'][ind0]
print(np.mean(B))