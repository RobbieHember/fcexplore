
#%% Import modules

import os
import numpy as np
import time
import matplotlib.pyplot as plt
import geopandas as gpd
import pandas as pd
import whitebox
import scipy.io as spio
import rasterio
from rasterio import features
from shapely import geometry
from scipy.interpolate import griddata
import cv2
import copy
import whitebox
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.bc1ha.bc1ha_util as u1ha
gp=gu.SetGraphics('Manuscript')

#%%
import sys
sys.path.insert(1,r'C:\whitebox-python-master') #This is where whitebox tools is stored.
from whitebox_tools import WhiteboxTools
wbt = WhiteboxTools()

#%%

wbt=whitebox.WhiteboxTools()

pth=r'C:\Users\rhember\Documents\Data\BC1ha\Terrain'
wbt.set_working_dir(pth)

dem=pth + '\\elevation.tif'
wbt.convert_raster_format(dem,pth + '\\dem.dep')
dep=pth + '\\dem.dep'

#%% Ridges

fout=pth + '\\ridges.tif'
wbt.find_ridges(dem,fout,line_thin=True)

z=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Terrain\ridges.tif')
plt.matshow(z['Data'])

#%% Saddles

fout=pth + '\\shape_index.tif'
wbt.shape_index(dem,fout,zfactor=1.0) # ,callback=default_callback

#%% Geomorphons

fout=pth + '\\geomorphon.tif'
wbt.geomorphons(dem,fout,search=50,threshold=0.0,fdist=0,skip=0,forms=True,residuals=False)

z=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Terrain\geomorphons.tif')
plt.matshow(z['Data'])