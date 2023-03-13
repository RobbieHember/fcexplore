
#%% IMPORT MODULES

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
import fcgadgets.bc1ha.bc1ha_utilities as u1ha

#%% Import look-up-tables

gp=gu.SetGraphics('Manuscript')
lut=u1ha.Import_BC1ha_LUTs()
zRef=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\BC_Land_Mask.tif')

#%% Import PFI data

gdf=gpd.read_file(r'C:\Users\rhember\Documents\Data\Basemaps\Basemaps.gdb',layer='NRC_POLITICAL_BOUNDARIES_1M_SP')
pth1=r'C:\Users\rhember\Documents\Data\ForestInventory\LiDAR\Boundary TSA\From GQ 20230302\BoundaryTSA_PFI_LnTotalWSV.tif'
pth2=r'C:\Users\rhember\Documents\Data\ForestInventory\LiDAR\Boundary TSA\From GQ 20230302\BoundaryTSA_PFI_LnTotalWSV_bc1ha2.tif'
pth_ref=r'C:\Users\rhember\Documents\Data\BC1ha\Admin\BC_Land_Mask.tif'
zB=gis.ReprojectRasterAndClipToRaster(pth1,pth2,pth_ref,gdf.crs)

#%%

zB=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\ForestInventory\LiDAR\Boundary TSA\From GQ 20230302\BoundaryTSA_PFI_LnTotalWSV_bc1ha2.tif')
#z=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\ForestInventory\LiDAR\Boundary TSA\From GQ 20230302\BoundaryTSA_PFI_LnTotalWSV.tif')

#plt.hist(zB['Data'][(zB['Data']>0) & (zB['Data']<1000) ].flatten()[0::100])

zBpfi=zB.copy()
zBpfi['Data']=0.5*0.5*zBpfi['Data']

plt.close('all')
plt.matshow(zBpfi['Data'],clim=[0,160])
plt.colorbar()

plt.hist(zBpfi['Data'].flatten()[0::100])

#%% Import Global Biomass

zBglo=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\Biomass\GlobBiomass\N80W140_agb pc8.tif')
zBglo['Data']=0.5*zBglo['Data']

plt.matshow(zBglo['Data'])

#%%

dp=(zBglo['Data']-zBpfi['Data'])/zBpfi['Data']*100
ind=np.where( (zBpfi['Data']<=0) | (zBpfi['Data']>400) ); dp[ind]=0

plt.close('all')
plt.matshow(dp,clim=[-100,100])
plt.colorbar()

#%%

ivl=10
x=zBpfi['Data'][0::ivl,0::ivl].flatten()
y=zBglo['Data'][0::ivl,0::ivl].flatten()

ikp=np.where( (x>0) & (x<300) & (y>0) )[0]
x=x[ikp]
y=y[ikp]

rs,txt=gu.GetRegStats(x,y)

plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(11,11))
ax.plot([0,400],[0,400],'-k',lw=2,color=[0.75,0.75,0.75])
ax.plot(x,y,'b.')
ax.plot(rs['xhat'],rs['yhat'],'k-',lw=1,label='Best fit')
ax.text(390,15,txt,fontsize=10,color='k',ha='right')
ax.text(330,330,'1:1',fontsize=8,ha='center')


bw=10; bin=np.arange(10,400+bw,bw)
N,mu,med,sig,se=gu.discres(x,y,bw,bin)

fig,ax=plt.subplots(1,figsize=gu.cm2inch(11,11))
plt.plot(bin,mu,'bo')