#%% Import modules
import os
import numpy as np
import time
import matplotlib.pyplot as plt
import geopandas as gpd
import pandas as pd
import pyproj
import scipy.io as spio
import fiona
import rasterio
from rasterio import features
from shapely import geometry
from scipy.interpolate import griddata
import cv2
import copy
#from rasterio.transform import from_origin
#from shapely.geometry import Point,Polygon
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.bc1ha.bc1ha_util as u1ha
import fcgadgets.cbrunner.cbrun_util as cbu
gp=gu.SetGraphics('Manuscript')

#%% Import paths and look-up-tables
meta=u1ha.Init()
zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])

z=u1ha.Import_Raster(meta,[],['lc_comp1_2019','lu_comp1_2019','lu_comp1_2049s1','lu_comp1_2049s2','lu_comp1_2049s3','lu_comp1_2049s4'],'Extract Grid')

#%%

d0=copy.deepcopy(meta['LUT']['Derived']['lu_comp1'])
for k in meta['LUT']['Derived']['lu_comp1'].keys():
    ind=np.where( (z['lc_comp1_2019']==1) & (z['lu_comp1_2019']==meta['LUT']['Derived']['lu_comp1'][k]) )
    d0[k]=ind[0].size/1e6

d=copy.deepcopy(meta['LUT']['Derived']['lu_comp1'])
for k in meta['LUT']['Derived']['lu_comp1'].keys():
    ind=np.where( (z['lc_comp1_2019']==1) & (z['lu_comp1_2049s4']==meta['LUT']['Derived']['lu_comp1'][k]) )
    d[k]=ind[0].size/1e6

print(d['Conservation Natural']-d0['Conservation Natural'])
print(d['Conservation Consistent']-d0['Conservation Consistent'])
print(d['Timber Production Intense']-d0['Timber Production Intense'])
print(d['Wildfire Risk Management']-d0['Wildfire Risk Management'])
print(d['Bioenergy Production']-d0['Bioenergy Production'])


