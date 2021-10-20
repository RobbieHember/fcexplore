'''

BC Climate5K Grid

'''

#%% Import modules

import os
import numpy as np
import pandas as pd
import geopandas as gpd
import fiona
import matplotlib.pyplot as plt
import numpy.matlib as ml
import matplotlib.pyplot as plt
import time
from shapely.geometry import Polygon,Point
from fcgadgets.macgyver import utilities_general as gu
from fcgadgets.macgyver import utilities_gis as gis
from fcgadgets.macgyver import utilities_inventory as invu
from fcgadgets.cbrunner import cbrun_utilities

#%% Define paths

pth=r'C:\Users\rhember\Documents\Data\BC5k'

#%% Define sparse grid sample

# Define regular grid sampling frequency
#sfreq=100 # 10 km
sfreq=50 # 5 km

# Import TSA maps
zTSA=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\tsa.tif')
lut_tsa=pd.read_excel(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\lut_tsa.xlsx')
tsa_boundaries=gpd.read_file(r'C:\Users\rhember\Documents\Data\TSA\tsa_boundaries.shp')

# Extract subgrid
zTSA['X']=zTSA['X'][0::sfreq,0::sfreq]
zTSA['Y']=zTSA['Y'][0::sfreq,0::sfreq]
zTSA['m'],zTSA['n']=zTSA['X'].shape
zTSA['Data']=zTSA['Data'][0::sfreq,0::sfreq]

# Define additional inclusion criteria (inside political boundary of BC)
iIreg=np.where( (zTSA['Data']!=255) )

# Save grid
flg=0
if flg==1:
    z=zTSA.copy()
    z['Data']=np.zeros(zTSA['Data'].shape,dtype='int16')
    z['Data'][iIreg]=1 # Treed = 1
    gis.SaveGeoTiff(z,pth + '\\GridSXY.tiff')
    plt.matshow(z['Data'])

#%% Generate sparse grid

# Apply filters to BC1ha grid 
sxy={}
sxy['x']=zTSA['X'][iIreg]
sxy['y']=zTSA['Y'][iIreg]
sxy['ID_TSA']=zTSA['Data'][iIreg]

# Save to pickle file
gu.opickle(meta['Paths']['Geospatial'] + '\\sxy.pkl',sxy)

# Save as geojson
flg=1
if flg==1:
    points=[]
    for k in range(sxy['x'].size):
        points.append(Point(sxy['x'][k],sxy['y'][k]))
    gdf_sxy=gpd.GeoDataFrame({'geometry':points,'ID_TSA':sxy['ID_TSA']})
    gdf_sxy.crs=tsa_boundaries.crs  
    gdf_sxy.to_file(meta['Paths']['Geospatial'] + '\\sxy.geojson',driver='GeoJSON')



#%% CMIP6
    
pin=r'C:\Users\rhember\Documents\Data\CMIP6'

namM=['ACCESS-CM2']
namV=['pr']
namS=['historical','ssp245','ssp585']

iM=0
iV=0

zH=gis.OpenGeoTiff(pin + '\\' + namM[iM] + '\\' + '\\' + namV[iV] + '_Amon_' + namM[iM] + '_historical_r1i1p1f1_gn_185001-201412.nc')

zH=gis.OpenGeoTiff(pin + '\\' + namM[iM] + '\\' + '\\' + namV[iV] + '_Amon_' + namM[iM] + '_' + namS[iS] + '_r1i1p1f1_gn_185001-201412.nc')

pr_Amon_ACCESS-CM2_ssp245_r1i1p1f1_gn_201501-210012

tv=gu.tvec('m',1850,2014)


tv.shape








