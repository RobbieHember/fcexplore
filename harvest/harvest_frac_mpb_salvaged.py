
'''
AREA OF MPB OUTBREAK THAT HAS BEEN HARVESTED
'''
#%% Import modules
import os
import numpy as np
import time
import matplotlib.pyplot as plt
import geopandas as gpd
from rasterio.transform import from_origin
import pandas as pd
import pyproj
import scipy.io as spio
import fiona
import rasterio
from rasterio import features
from shapely.geometry import Point,Polygon
from shapely import geometry
from scipy.interpolate import griddata
import cv2
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.macgyver.util_query_gdb as qgdb
import fcgadgets.bc1ha.bc1ha_util as u1ha

#%% Import paths and look-up-tables
meta=u1ha.Init()
gp=gu.SetGraphics('Manuscript')

#%% Import data
zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
zLCC1_c=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\LandCoveruse\LandCoverClass1_Current.tif')

zY=np.zeros(zRef['Data'].shape,dtype='int16')
for iY in range(10):
    zY0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\PEST_SEVERITY_CODE_IBM_' + str(iY+1) + '_Year.tif')['Data']
    zS0=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\PEST_SEVERITY_CODE_IBM_' + str(iY+1) + '_Severity.tif')['Data']
    ind=np.where( (zS0>=meta['LUT']['PEST_INFESTATION_POLY']['PEST_SEVERITY_CODE']['L']) & (zY==0) )
    zY[ind]=zY0[ind]

zH=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\Harvest_Consol1_Year.tif')
zF=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\PROT_HISTORICAL_FIRE_POLYS_SP\FIRE_YEAR_YearLast.tif')

zFo=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\PROT_HISTORICAL_FIRE_POLYS_SP\FIRE_YEAR_2022.tif')
ind=np.where(zFo['Data']>0); zF['Data'][ind]=2022
zFo=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\PROT_HISTORICAL_FIRE_POLYS_SP\FIRE_YEAR_2023.tif')
ind=np.where(zFo['Data']>0); zF['Data'][ind]=2023

#%% Map

zMask=np.zeros(zRef['Data'].shape,dtype='int16')
ind1=np.where(zY>0); zMask[ind1]=1
ind2=np.where( (zY>0) & (zH['Data']>zY) ); zMask[ind2]=2
ind3=np.where( (zY>0) & (zF['Data']>zY) ); zMask[ind3]=3
ind4=np.where( (zY>0) & (zH['Data']>zY) & (zF['Data']>zY) ); zMask[ind4]=4
plt.close('all'); plt.matshow(zMask)
print(ind1[0].size)
print(ind2[0].size/ind1[0].size*100)
print(ind3[0].size/ind1[0].size*100)
print(ind4[0].size/ind1[0].size*100)

ind5=np.where( (zY>0) & (zH['Data']>zY) & (zF['Data']>zH['Data']) );
print(ind5[0].size/ind2[0].size)

#%%

n1=np.where( (zLCC1_c['Data']==1) & (zY>0) & (zF['Data']<zY) )[0].size
n2=np.where( (zLCC1_c['Data']==1) & (zY>0) & (zF['Data']>zY) )[0].size

n3=np.where( (zLCC1_c['Data']==1) & (zY==0) & (zF['Data']<2000) )[0].size
n4=np.where( (zLCC1_c['Data']==1) & (zY==0) & (zF['Data']>=2000) )[0].size

print(n2/n1)
print(n4/n3)
