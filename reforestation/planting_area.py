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

meta['Graphics']={}
meta['Graphics']['figwidth']=16
#meta['Graphics']['sidespace']=0.25
meta['Graphics']['sidespace']=0
meta['Graphics']['ax1 pos']=[0,0,1-meta['Graphics']['sidespace']-0.01,1]
meta['Graphics']['ax1 vis']='off'
meta['Graphics']['ax1 gridvis']=False
meta['Graphics']['ax2 pos']=[1-meta['Graphics']['sidespace']+meta['Graphics']['ax1 pos'][0]+0.01,0.6,0.03,0.35]
meta['Graphics']['ax2 pos long']=[1-meta['Graphics']['sidespace']+meta['Graphics']['ax1 pos'][0]+0.01,0.1,0.03,0.8]

meta['LUT']=u1ha.Import_BC1ha_LUTs()

gp=gu.SetGraphics('Manuscript')

#%% Import base maps

gdf=u1ha.Import_GDBs_ProvinceWide()

zRef=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\BC_Land_Mask.tif')

zCrown=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\CrownForestLandMask.tif')
#zLCC1=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\LandUseLandCover\LandCoverClass1.tif')
zFCR=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Results\RSLT_FOREST_COVER_RESERVE_SVW.tif')
zPL=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Results\Planting_FromRESULTS_MaskCount.tif')
zHRE=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\Harvest_FromRESULTS_Mask.tif')
zHCC=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\VEG_CONSOLIDATED_CUT_BLOCKS_SP_All.tif')
zH=zRef.copy(); zH['Data']=np.maximum(zHRE['Data'],zHCC['Data'])

dH_ts=gu.ipickle(r'C:\Users\rhember\Documents\Data\Harvest\Harvest Area\HarvestAreaBC.pkl')
np.sum(dH_ts['Area Harvested CC'])/1e6
np.sum(dH_ts['Area Planted RESULTS'])/1e6


#%% Overall proportion of harvested area that was subsequently planted

ind0=np.where( (zCrown['Data']>0) & (zFCR['Data']==0) & (zH['Data']>0) ); print(ind0[0].size/1e6)

ind1=np.where( (zCrown['Data']>0) & (zFCR['Data']==0) & (zH['Data']>0) & (zPL['Data']>0) ); print(ind1[0].size/1e6)
ind1=np.where( (zPL['Data']>0) ); print(ind1[0].size/1e6)

#print(ind1[0].size/ind0[0].size)

#%% Area planted

d={}
d['tv']=np.arange(1950,2022,1)
d['A Planted First H']=np.zeros(d['tv'].size)
d['A Planted First']=np.zeros(d['tv'].size)
d['A Planted All']=np.zeros(d['tv'].size)
flg=np.zeros(zRef['Data'].shape,dtype='int8')
for iT in range(d['tv'].size):
    zPLy=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Results\Planting_FromRESULTS_' + str(d['tv'][iT]) + '.tif')

    ind=np.where( (zPLy['Data']>0) & (flg==0) & (zH['Data']>0) )
    d['A Planted First H'][iT]=ind[0].size

    ind=np.where( (zPLy['Data']>0) & (flg==0) )
    d['A Planted First'][iT]=ind[0].size

    ind=np.where(zPLy['Data']>0)
    d['A Planted All'][iT]=ind[0].size

    flg[ind]=1
    print(d['tv'][iT])

plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(12,12));
plt.plot(dH_ts['tv'],dH_ts['Area Harvested CC']/1e3,'-bo')
plt.plot(d['tv'],d['A Planted First']/1e3,'-gs')
plt.plot(d['tv'],d['A Planted First H']/1e3,'-k*')
plt.plot(dH_ts['tv'],dH_ts['Area Planted RESULTS']/1e3,'-cd')

fig,ax=plt.subplots(1,figsize=gu.cm2inch(12,12));
plt.plot(dH_ts['tv'],d['A Planted First']/dH_ts['Area Harvested CC'],'-bo')


fig,ax=plt.subplots(1,figsize=gu.cm2inch(12,12));
plt.plot(d['tv'],d['A Planted First']-d['A Planted First H'],'-gs')

#%%

yr=2005
zHREy=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\Harvest_FromRESULTS_' + str(yr) + '.tif')
zHCCy=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\VEG_CONSOLIDATED_CUT_BLOCKS_SP_' + str(yr) + '.tif')
zHy=zRef.copy(); zHy['Data']=np.maximum(zHREy['Data'],zHCCy['Data'])

z=np.zeros(zRef['Data'].shape,dtype='int8')
ind0=np.where( (zHy['Data']>0) & (zFCR['Data']==0) ); z[ind0]=1

N=np.zeros((15,2))
for iT in range(N.shape[0]):
    print(iT)
    try:
        zPLy=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Results\Planting_FromRESULTS_' + str(yr+iT) + '.tif')
    except:
        break
    ind=np.where( (zHy['Data']>0) & (zFCR['Data']==0) & (zPLy['Data']>0) )
    if ind[0].size>0:
        N[iT,0]=ind[0].size
    ind=np.where( (z==1) & (zHy['Data']>0) & (zFCR['Data']==0) & (zPLy['Data']>0) )
    if ind[0].size>0:
        z[ind]=0
        N[iT,1]=ind[0].size


plt.plot(np.cumsum(N[:,1]/ind0[0].size*100),'-rd')

plt.close('all')
plt.plot(np.cumsum(N[:,0]/ind0[0].size*100),'-bo')
plt.plot(np.cumsum(N[:,1]/ind0[0].size*100),'-gs')

#%% Map areas

z=np.zeros(zRef['Data'].shape,dtype='int8')
ind=np.where( (zH['Data']>0) & (zFCR['Data']==0) & (zCrown['Data']>0) ); z[ind]=1
ind=np.where( (zH['Data']>0) & (zFCR['Data']==0) & (zCrown['Data']>0) & (zPL['Data']==0) ); z[ind]=2
plt.matshow(z)

ind=np.where( (zH['Data']>0) & (zFCR['Data']==0) & (zCrown['Data']>0) ); print(ind[0].size/1e6)
ind=np.where( (zHRE['Data']>0) & (zFCR['Data']==0) & (zCrown['Data']>0) ); print(ind[0].size/1e6)
ind=np.where( (zHCC['Data']>0) & (zFCR['Data']==0) & (zCrown['Data']>0) ); print(ind[0].size/1e6)


ind1=np.where( (zH['Data']>0) & (zFCR['Data']==0) & (zCrown['Data']>0) );
ind2=np.where( (zHCC['Data']>0) & (zFCR['Data']==0) & (zCrown['Data']>0) );
ind1[0].size-ind2[0].size


#%% Compare harvest sources

ind=np.where( (zCrown['Data']==1) & (zFCR['Data']==0) & (zHRE['Data']>0) )
A_RE=ind[0].size

ind=np.where( (zCrown['Data']==1) & (zFCR['Data']==0) & (zHCC['Data']>0) )
A_CC=ind[0].size

ind=np.where( (zCrown['Data']==1) & (zFCR['Data']==0) & (zHCC['Data']>0) & (zHRE['Data']>0) )
A_wCC_wRE=ind[0].size

ind=np.where( (zCrown['Data']==1) & (zFCR['Data']==0) & (zHCC['Data']>0) & (zHRE['Data']==0) )
A_wCC_woRE=ind[0].size

ind=np.where( (zCrown['Data']==1) & (zFCR['Data']==0) & (zHCC['Data']==0) & (zHRE['Data']==0) )
A_woCC_woRE=ind[0].size

ind=np.where( (zCrown['Data']==1) & (zFCR['Data']==0) & (zHCC['Data']==0) & (zHRE['Data']>0) )
A_woCC_wRE=ind[0].size

#%%

print(A_CC/1e3)
print(A_RE/1e3)
print(A_wCC_wRE/1e3)
print(A_wCC_woRE/1e3)
#print(A_woCC_woRE/1e3)
print(A_woCC_wRE/1e3)



#%% Area of harvest that was planted

d={}
d['tv']=np.arange(1950,2022,1)
d['A Harvest']=np.zeros(d['tv'].size)
#d['A Harvest CC']=np.zeros(d['tv'].size)
d['A Harvest and Planted']=np.zeros(d['tv'].size)
for iT in range(d['tv'].size):
    #zHREy=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\Harvest_FromRESULTS_' + str(d['tv'][iT]) + '.tif')
    zHCCy=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\VEG_CONSOLIDATED_CUT_BLOCKS_SP_' + str(d['tv'][iT]) + '.tif')

    #zHy=zRef.copy();
    #zHy['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
    #ind=np.where( (zHREy['Data']>0) | (zHCCy['Data']>0) )
    #zHy['Data'][ind]=1

    ind=np.where( (zCrown['Data']>0) & (zFCR['Data']==0) & (zHCCy['Data']>0) )
    d['A Harvest'][iT]=ind[0].size

    #ind=np.where( (zCrown['Data']>0) & (zFCR['Data']==0) & (zHCCy['Data']>0) )
    #d['A Harvest CC'][iT]=ind[0].size

    ind=np.where( (zCrown['Data']>0) & (zFCR['Data']==0) & (zHCCy['Data']>0) & (zPL['Data']>0) )
    d['A Harvest and Planted'][iT]=ind[0].size
    print(d['tv'][iT])

gu.opickle(r'C:\Users\rhember\Documents\Data\Harvest\Harvest Area\HarvestAreaThatWasPlanted.pkl',d)


#%%



plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(12,12));
plt.plot(d2['tv'],d2['Area Harvested CC']/1e3,'-bo')
#plt.plot(d['tv'],d['A Harvest']/1e3,'-bo')
#plt.plot(d['tv'],d['A Harvest CC']/1e3,'-cd')
plt.plot(d['tv'],d['A Harvest and Planted']/1e3,'-gs')

fig,ax=plt.subplots(1,figsize=gu.cm2inch(12,12));
plt.plot(d['tv'],d['A Harvest and Planted']/d2['Area Harvested CC'],'-ko')

print(np.sum(d['A Harvest'])-np.sum(d['A Harvest CC']))

#%%

yr=2009
zHREy=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\Harvest_FromRESULTS_' + str(yr) + '.tif')
zHCCy=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\VEG_CONSOLIDATED_CUT_BLOCKS_SP_' + str(yr) + '.tif')
zHy=zRef.copy(); zHy['Data']=np.maximum(zHREy['Data'],zHCCy['Data'])

z=np.zeros(zRef['Data'].shape,dtype='int8')
ind0=np.where( (zHy['Data']>0) & (zFCR['Data']==0) ); z[ind0]=1

N=np.zeros((20,2))
for iT in range(N.shape[0]):
    print(iT)
    try:
        zPLy=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Results\Planting_FromRESULTS_' + str(yr+iT) + '.tif')
    except:
        break
    ind=np.where( (zHy['Data']>0) & (zFCR['Data']==0) & (zPLy['Data']>0) )
    if ind[0].size>0:
        N[iT,0]=ind[0].size
    ind=np.where( (z==1) & (zHy['Data']>0) & (zFCR['Data']==0) & (zPLy['Data']>0) )
    if ind[0].size>0:
        z[ind]=0
        N[iT,1]=ind[0].size


plt.plot(np.cumsum(N[:,1]/ind0[0].size*100),'-rd')

plt.close('all')
plt.plot(np.cumsum(N[:,0]/ind0[0].size*100),'-bo')
plt.plot(np.cumsum(N[:,1]/ind0[0].size*100),'-gs')
