
#%% Import Modules

import sys
import numpy as np
from osgeo import gdal
from osgeo import osr
import matplotlib.pyplot as plt
import geopandas as gpd
from rasterio.transform import from_origin
import pandas as pd
import scipy.io as spio
import fiona
import rasterio
from rasterio import features
from shapely.geometry import Point,Polygon
from shapely import geometry
import cv2

import fcgadgets.macgyver.utilities_general as gu
import fcgadgets.macgyver.utilities_gis as gis
from fcgadgets.macgyver import utilities_inventory as invu
import fcgadgets.macgyver.utilities_query_gdb as qgdb
import fcgadgets.bc1ha.bc1ha_utilities as bc1hau

gp=gu.SetGraphics('Manuscript')

#%% Import data

zRef=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\BC_Land_Mask.tif')
zCrown=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\CrownForestLandMask.tif')
zGFC=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\GlobalForestChange_LossYear_2021.tif')

zGFC2=zRef.copy()
zGFC2['Data']=np.zeros(zRef['Data'].shape,dtype='int8')

# Keep track of IBM so there is no duplication
ibm_ind=np.zeros(zRef['Data'].shape,dtype='int8')

tv=np.arange(2001,2022,1)
d={}
d['A IBM']=np.zeros(tv.size)
d['A WF']=np.zeros(tv.size)
d['A CC']=np.zeros(tv.size)
d['A GFC']=np.zeros(tv.size)
d['A GFC Priv']=np.zeros(tv.size)
for iT in range(tv.size):
    print(tv[iT])

    ind=np.where( (zCrown['Data']==0) & (zGFC['Data']+2000==tv[iT]) )
    d['A GFC Priv'][iT]=ind[0].size
    ind=np.where( (zCrown['Data']==1) & (zGFC['Data']+2000==tv[iT]) )
    d['A GFC'][iT]=ind[0].size

    z0=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\PEST_INFESTATION_POLY_IBM_SeverityClass_' + str(tv[iT]) + '.tif')
    ind=np.where( (zCrown['Data']==1) & (z0['Data']>2) & (ibm_ind==0) ) # Excluding trace
    ibm_ind[ind]=1
    d['A IBM'][iT]=ind[0].size
    ind=np.where( (zCrown['Data']==1) & (zGFC['Data']>0) & (z0['Data']>2) )
    zGFC2['Data'][ind]=0

    z0=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\VEG_CONSOLIDATED_CUT_BLOCKS_SP_' + str(tv[iT]) + '.tif')
    ind=np.where( (zCrown['Data']==1) & (z0['Data']==1) )
    d['A CC'][iT]=ind[0].size
    ind=np.where( (zCrown['Data']==1) & (zGFC['Data']>0) & (z0['Data']>0) )
    zGFC2['Data'][ind]=0

    z0=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\PROT_HISTORICAL_FIRE_POLYS_SP_' + str(tv[iT]) + '.tif')
    ind=np.where( (zCrown['Data']==1) & (z0['Data']>0) )
    d['A WF'][iT]=ind[0].size
    ind=np.where( (zCrown['Data']==1) & (zGFC['Data']>0) & (z0['Data']>0) )
    zGFC2['Data'][ind]=0

#%% Plot

wd=0.8
plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,7))
ax.bar(tv,d['A IBM']/1e6,wd,facecolor=[0.8,1,0.7],label='Mountian pine beetle')
ax.bar(tv,d['A WF']/1e6,wd,facecolor=[1,0.8,0.7],bottom=d['A IBM']/1e6,label='Wildfire')
ax.bar(tv,d['A CC']/1e6,wd,facecolor=[0.8,0.9,1],bottom=(d['A IBM']+d['A WF'])/1e6,label='Harvest')
ax.plot(tv,d['A GFC']/1e6,'-ko',mfc='k',ms=5,label='Tree cover loss (Global Forest Change 2021)')
#ax.plot(tv,d['A GFC Priv']/1e6,'-ks',mfc='w',ms=5,label='Tree cover loss (Global Forest Change 2021)')
ax.set(xlabel='Time, years',xticks=np.arange(2000,2026,1),ylabel='Area affected (million hectares)',xlim=[2000.25,2021.75])
ax.legend(loc='upper right',facecolor=[1,1,1],frameon=False);
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Global Forest Change\Comparison','png',900)

