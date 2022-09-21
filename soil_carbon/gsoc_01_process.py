
"""
Global Soil Carbon
"""

#%% IMPORT MODULES

import sys
import numpy as np
import gc as garc
from osgeo import gdal
from osgeo import osr
import matplotlib.pyplot as plt
import geopandas as gpd
import pandas as pd
import fiona
import rasterio
import pyproj
from rasterio import features
#import cv2
from shapely.geometry import Point, Polygon
from shapely import geometry
from rasterio.transform import from_origin
from osgeo import osr
from scipy.interpolate import griddata
import statsmodels.formula.api as smf


import fcgadgets.macgyver.utilities_general as gu
import fcgadgets.macgyver.utilities_gis as gis
from fcgadgets.macgyver import utilities_inventory as invu

#%% Import data

#  GLOSIS (https://www.fao.org/documents/card/en/c/I8891EN)
pth0=r'C:\Users\rhember\Documents\Data\Soils\GSOC\GSOCmap1.6.1.tif'
zS=gis.OpenGeoTiff(pth0)

# Clip to BC area
zS1=gis.ClipRasterByXYLimits(zS,[-142,-113],[47,61])

pth1=r'C:\Users\rhember\Documents\Data\Soils\GSOC\GSOCmap1.6.1 c.tif'
gis.SaveGeoTiff(zS1,pth1)

#%% Reproject

# Get CRS for BC
gdf=gpd.read_file(r'C:\Users\rhember\Documents\Data\Basemaps\Basemaps.gdb',layer='NRC_POLITICAL_BOUNDARIES_1M_SP')

pth2=r'C:\Users\rhember\Documents\Data\Soils\GSOC\GSOCmap1.6.1 cp.tif'
gis.ReprojectGeoTiff(pth1,pth2,gdf.crs)

#%%

pth_ref=r'C:\Users\rhember\Documents\Data\BC1ha\Admin\tsa.tif'
pth3=r'C:\Users\rhember\Documents\Data\Soils\GSOC\GSOCmap1.6.1 cpc.tif'
gis.ReviseRasterExtent(pth2,pth_ref,pth3)

#%% Clip to TSA grid

pth_ref=r'C:\Users\rhember\Documents\Data\BC1ha\Admin\tsa.tif'
zTSA=gis.OpenGeoTiff(pth_ref)

zIn=gis.OpenGeoTiff(pth3)

z=gis.ClipToRaster(zIn,zTSA)

ind=np.where(zTSA['Data']==255)
z['Data'][ind]=0
z['Data']=z['Data'].astype('int16')

zTSA['Data']=z['Data']

#plt.matshow(z['Data'],clim=[0,150])
pth4=r'C:\Users\rhember\Documents\Data\BC1ha\Soil\gsoc2010_bc1ha.tif'
gis.SaveGeoTiff(zTSA,pth4)


#%% Query by BGC zone

zSOC=gis.OpenGeoTiff(pth4)
zSOC['Data']=zSOC['Data'].astype('float')

pthBGC=r'C:\Users\rhember\Documents\Data\BC1ha\VRI\becz.tif'
zBGC=gis.OpenGeoTiff(pthBGC)
lutBGC=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\becz_lut.xlsx')

pthLC2=r'C:\Users\rhember\Documents\Data\BC1ha\VRI\lc2.tif'
zLC2=gis.OpenGeoTiff(pthLC2)

u=np.unique(zBGC['Data'])
u=u[u<255]
lab=np.array(['' for _ in range(u.size)],dtype=object)
soc_mu=np.zeros(u.size)
for i in range(u.size):
    ind=np.where( (zLC2['Data']==4) & (zBGC['Data']==u[i]) & (zSOC['Data']>0) & (zSOC['Data']<1000) )
    soc_mu[i]=np.mean(zSOC['Data'][ind])
    ind=np.where(lutBGC['VALUE']==u[i])[0]
    if ind.size>0:
        lab[i]=lutBGC['ZONE'][ind][0]

plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(14,8))
ax.bar(np.arange(u.size),soc_mu)
ax.set(xticks=np.arange(u.size),xticklabels=lab)

#%% Import Canadian Upland forest database

ufd=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\Soils\Shaw et al 2018 Database\SITES.xlsx')
ufd_p=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\Soils\Shaw et al 2018 Database\PROFILES.xlsx')

# Add depth from profiles to site
ufd['Depth']=np.zeros(ufd['LOCATION_ID'].size)
for i in range(ufd['LOCATION_ID'].size):
    ind=np.where(ufd_p['LOCATION_ID']==ufd['LOCATION_ID'][i])[0]
    if ind.size==0:
        continue
    ind2=np.where( ufd_p['UPPER_HZN_LIMIT'][ind]==np.max(ufd_p['UPPER_HZN_LIMIT'][ind]) )[0]
    ufd['Depth'][i]=ufd_p['UPPER_HZN_LIMIT'][ind[ind2]]+ufd_p['HZN_THICKNESS'][ind[ind2]]

# Keep BC
ind=np.where( (ufd['PROV_TERR']=='BC') )[0]
for k in ufd.keys():
    ufd[k]=ufd[k][ind]

# Dummy variables
uS=np.unique(ufd['LEAD_SPECIES_CODE'])
for iS in uS:
    ufd[iS]=np.zeros(ufd['LOCATION_ID'].size)
    ind=np.where(ufd['LEAD_SPECIES_CODE']==iS)[0]
    ufd[iS][ind]=1

uO=np.unique(ufd['ORDER'])
for iO in uO:
    ufd[iO]=np.zeros(ufd['ORDER'].size)
    ind=np.where(ufd['ORDER']==iO)[0]
    ufd[iO][ind]=1

# Get BC coords
srs=gis.ImportSRSs()
ufd['x']=np.zeros(ufd['LOCATION_ID'].size)
ufd['y']=np.zeros(ufd['LOCATION_ID'].size)
for i in range(ufd['x'].size):
    ufd['x'][i],ufd['y'][i]=srs['Proj']['BC1ha'](ufd['LONGITUDE'][i],ufd['LATITUDE'][i])

# Get BC1ha index
x=zTSA['X'][0,:]
y=zTSA['Y'][:,0]
ix=np.zeros(ufd['LOCATION_ID'].size,dtype=int)
iy=np.zeros(ufd['LOCATION_ID'].size,dtype=int)
for i in range(ufd['x'].size):
    ix[i]=np.where( np.abs(ufd['x'][i]-x)==np.min(np.abs(ufd['x'][i]-x)) )[0]
    iy[i]=np.where( np.abs(ufd['y'][i]-y)==np.min(np.abs(ufd['y'][i]-y)) )[0]
# Check that it worked
x=zTSA['X'][iy,ix]
y=zTSA['Y'][iy,ix]
plt.plot(x,y,'k.')

# Import BC1ha data

elev=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Terrain\elevation.tif')
ws=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\BC1ha_ws_gs_norm_1971to2000_comp_hist_v1.tif')
dwf=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\BC1ha_dwf_ann_norm_1971to2000_si_hist_v1.tif')
becz=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\becz.tif')
age1=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\age1.tif')

ufd['elev']=elev['Data'][iy,ix]
ufd['ws']=ws['Data'][iy,ix]
ufd['dwf']=dwf['Data'][iy,ix]
ufd['becz']=becz['Data'][iy,ix]
ufd['age1']=age1['Data'][iy,ix]
ufd['age1'][np.where(ufd['age1']<=0)[0]]=0

#%%

lutBGC=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\becz_lut.xlsx')

u=np.unique(ufd['becz'])
u=u[u<255]
lab=np.array(['' for _ in range(u.size)],dtype=object)
soc_mu=np.zeros(u.size)
for i in range(u.size):
    ind=np.where( (ufd['becz']==u[i]) & (ufd['TOT_C_THA']>0) )
    soc_mu[i]=np.mean(ufd['MIN_C_THA'][ind])
    ind=np.where(lutBGC['VALUE']==u[i])[0]
    if ind.size>0:
        lab[i]=lutBGC['ZONE'][ind][0]

plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(14,8))
ax.bar(np.arange(u.size),soc_mu)
ax.set(xticks=np.arange(u.size),xticklabels=lab)

#%%

u=np.unique(ufd['ORDER'])
lab=np.array(['' for _ in range(u.size)],dtype=object)
soc_mu=np.zeros(u.size)
for i in range(u.size):
    ind=np.where( (ufd['ORDER']==u[i]) & (ufd['TOT_C_THA']>0) )
    soc_mu[i]=np.mean(ufd['TOT_C_THA'][ind])


plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(14,8))
ax.bar(np.arange(u.size),soc_mu)
ax.set(xticks=np.arange(u.size),xticklabels=u)

#%%

df=pd.DataFrame.from_dict(ufd)

#x = df[['debt_ratio', 'industry']]
#y = df['cash_flow']
#%%

# NB. unlike sm.OLS, there is "intercept" term is included here
#form="TOT_C_THA ~ C(becz)"
#form="TOT_C_THA ~ C(ORDER)"
form="TOT_C_THA ~ PRECIP + MAT"
mr=smf.ols(formula=form, data=df).fit()

print(mr.summary())

#%%

x=ufd['Depth']; mn,mx=gu.minmax(x); r=mx-mn; bw=r/10; bin=np.arange(mn,mx,bw)
N,mu,sig,se=gu.discres(x,ufd['TOT_C_THA'],bw,bin)

plt.close('all')
plt.plot(bin,mu,'ko')

#%% Convert to dataframe

data = pd.DataFrame(boston.data)
data.columns = boston.feature_names










