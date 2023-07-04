
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

#%% Import data

gdf=u1ha.Import_GDBs_ProvinceWide_Simple()

zRef=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\BC_Land_Mask.tif')
zRD=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\RegionalDistrict.tif')
zH=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\NTEM_harv85to20_bc1ha.tif')
zF=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\PROT_HISTORICAL_FIRE_POLYS_SP_All.tif')
zIBM=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\IBM_Mask_ExcTrace.tif')

dCEC=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\BC1ha\LandUseLandCover\luc_cec_lc.xlsx')
zCEC10=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\LandUseLandCover\CEC_LandCover2010.tif')
zCEC20=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\LandUseLandCover\CEC_LandCover2020.tif')

#%% NALCMS CEC 2020
# Compress categories to remove tropics

meta,zCEC10c,zCEC20c=u1ha.NALCMS_Compress(meta,zRef,zCEC10,zCEC20)

#%% Land use change analysis -BC wide (based on compressed categories)

id=np.array(list(meta['LUT']['ceclc20 rev'].values()))

# Areas at start and end
v=np.zeros( (2,id.size) )
for i in range(id.size):
    ind=np.where( (zCEC10c['Data']==id[i]) & (zH['Data']==0) & (zF['Data']==0) & (zIBM['Data']==0) )
    v[0,i]=ind[0].size
    ind=np.where( (zCEC20c['Data']==id[i]) & (zH['Data']==0) & (zF['Data']==0) & (zIBM['Data']==0) )
    v[1,i]=ind[0].size
df=pd.DataFrame(v)
df.to_excel(r'C:\Users\rhember\OneDrive - Government of BC\Manuscripts\Land Use Change Assessment\CEC Change Analysis\BC Wide\cec_c0_nodist.xlsx')

# Transitions
v=np.zeros( (id.size,id.size) )
for i in range(id.size):
    for j in range(id.size):
        ind=np.where( (zCEC10c['Data']==id[i]) & (zCEC20c['Data']==id[j]) & (zH['Data']==0) & (zF['Data']==0) & (zIBM['Data']==0) )
        v[i,j]=ind[0].size
df=pd.DataFrame(v)
df.to_excel(r'C:\Users\rhember\OneDrive - Government of BC\Manuscripts\Land Use Change Assessment\CEC Change Analysis\BC Wide\cec_c1_nodist.xlsx')

#%% Land use change analysis - by regional district (based on compressed categories)

id=np.array(list(meta['LUT']['ceclc20 rev'].values()))

rs={}

for k in meta['LUT']['regdis'].keys():
    # Areas at start and end
    vI=np.zeros( (2,id.size) )
    for i in range(id.size):
        ind=np.where( (zCEC10c['Data']==id[i]) & (zH['Data']==0) & (zF['Data']==0) & (zIBM['Data']==0) & (zRD['Data']==meta['LUT']['regdis'][k]) )
        vI[0,i]=ind[0].size
        ind=np.where( (zCEC20c['Data']==id[i]) & (zH['Data']==0) & (zF['Data']==0) & (zIBM['Data']==0) & (zRD['Data']==meta['LUT']['regdis'][k]) )
        vI[1,i]=ind[0].size
    df=pd.DataFrame(vI)
    df.to_excel(r'C:\Users\rhember\OneDrive - Government of BC\Manuscripts\Land Use Change Assessment\CEC Change Analysis\cec_c0_nodist_' + k + '.xlsx')

    # Transitions
    vT=np.zeros( (id.size,id.size) )
    for i in range(id.size):
        for j in range(id.size):
            ind=np.where( (zCEC10c['Data']==id[i]) & (zCEC20c['Data']==id[j]) & (zH['Data']==0) & (zF['Data']==0) & (zIBM['Data']==0) & (zRD['Data']==meta['LUT']['regdis'][k]) )
            vT[i,j]=ind[0].size
    df=pd.DataFrame(vT)
    df.to_excel(r'C:\Users\rhember\OneDrive - Government of BC\Manuscripts\Land Use Change Assessment\CEC Change Analysis\cec_c1_nodist_' + k + '.xlsx')

    # Net deforestation summary
    ind=np.array([meta['LUT']['ceclc20 rev']['Cropland'],meta['LUT']['ceclc20 rev']['Barren Ground'],meta['LUT']['ceclc20 rev']['Urban']],dtype='int8')
    Ad=np.sum(vT[0,ind-1])
    Aa=np.sum(vT[ind-1,0])
    Pd=np.sum(vT[0,ind-1])/vI[0,0]*100
    Pa=np.sum(vT[ind-1,0])/vI[0,0]*100
    Pn=Pa-Pd

#%% Map deforestation
# Notes:
#   1) Forest to Grassland and Forest to Shrubland excluded due to suspicion of error
#   2) Areas with harvest excluded (harvest from NTEM)

z=zRef.copy()
z['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
ind=np.where( (zCEC10['Data']==meta['LUT']['Forest']) & (zCEC20['Data']==meta['LUT']['Urban']) & (zH['Data']==0) )
z['Data'][ind]=1
ind=np.where( (zCEC10['Data']==meta['LUT']['Forest']) & (zCEC20['Data']==meta['LUT']['Cropland']) & (zH['Data']==0) )
z['Data'][ind]=2
ind=np.where( (zCEC10['Data']==meta['LUT']['Forest']) & (zCEC20['Data']==meta['LUT']['Barren Ground']) & (zH['Data']==0) )
z['Data'][ind]=3
#plt.matshow(z['Data'])
gis.SaveGeoTiff(z,r'C:\Users\rhember\Documents\Data\BC1ha\LandUseLandCover\deforestation10to20_FromCEC.tif')

#%% Map Afforestation (from cropland, urban, or barren ground)

z=zRef.copy()
z['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
ind=np.where( (zCEC10['Data']==meta['LUT']['Urban']) & (zCEC20['Data']==meta['LUT']['Forest']) & (zH['Data']==0) )
z['Data'][ind]=1
ind=np.where( (zCEC10['Data']==meta['LUT']['Cropland']) & (zCEC20['Data']==meta['LUT']['Forest']) & (zH['Data']==0) )
z['Data'][ind]=2
ind=np.where( (zCEC10['Data']==meta['LUT']['Barren Ground']) & (zCEC20['Data']==meta['LUT']['Forest']) & (zH['Data']==0) )
z['Data'][ind]=3
gis.SaveGeoTiff(z,r'C:\Users\rhember\Documents\Data\BC1ha\LandUseLandCover\afforestation10to20_FromCEC.tif')

ind=np.where( (zH['Data']>0) )
z['Data'][ind]=2
plt.matshow(z['Data'])

ind=np.where(z['Data']==1)
ind[1].size/1000000


#%% Map Forest to Grassland and Forest to Shrubland
# *** Way too much to be real! ***

z=zRef.copy()
z['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
ind=np.where( (zCEC10['Data']==meta['LUT']['Forest']) & (zCEC20['Data']==meta['LUT']['Grassland']) )
z['Data'][ind]=1
ind=np.where( (zCEC10['Data']==meta['LUT']['Forest']) & (zCEC20['Data']==meta['LUT']['Shrubland']) )
z['Data'][ind]=1
ind=np.where( (zH['Data']>0) )
z['Data'][ind]=2
plt.matshow(z['Data'])







ind=np.where(zH['Data']>0)
ind[1].size/1000000