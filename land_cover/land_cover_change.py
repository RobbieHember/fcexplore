
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
import fcgadgets.bc1ha.bc1ha_utilities as u1ha

#%% Import parameters

meta=u1ha.Init()
meta=u1ha.ImportLUTs(meta)
gp=gu.SetGraphics('Manuscript')

#%% Import data

zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])

vList=['rd','lcc_cec_2010','lcc_cec_2020','harv_yr_con1','fire_yr','ibm_yr'] #'lcc1_c','gfcly','gfcly_filt',
z=u1ha.Import_Raster(meta,[],vList)

#%% Land use change analysis - BC wide

id=np.array(list(meta['LUT']['Derived']['lcc_cec_c'].values()))
# Areas at start and end
v=np.zeros( (2,id.size) )
for i in range(id.size):
    ind=np.where( (z['lcc_cec_2010']['Data']==id[i]) & (z['harv_yr_con1']['Data']==0) & (z['fire_yr']['Data']==0) & (z['ibm_yr']['Data']==0) )
    v[0,i]=ind[0].size
    ind=np.where( (z['lcc_cec_2020']['Data']==id[i]) & (z['harv_yr_con1']['Data']==0) & (z['fire_yr']['Data']==0) & (z['ibm_yr']['Data']==0) )
    v[1,i]=ind[0].size
df=pd.DataFrame(v)
df.to_excel(r'C:\Users\rhember\OneDrive - Government of BC\Manuscripts\Land Use Change Assessment\CEC Change Analysis\BC Wide\cec_c0_nodist.xlsx')

# Transitions
v=np.zeros( (id.size,id.size) )
for i in range(id.size):
    for j in range(id.size):
        ind=np.where( (z['lcc_cec_2010']['Data']==id[i]) & (z['lcc_cec_2020']['Data']==id[j]) & (z['harv_yr_con1']['Data']==0) & (z['fire_yr']['Data']==0) & (z['ibm_yr']['Data']==0) )
        v[i,j]=ind[0].size
df=pd.DataFrame(v)
df.to_excel(r'C:\Users\rhember\OneDrive - Government of BC\Manuscripts\Land Use Change Assessment\CEC Change Analysis\BC Wide\cec_c1_nodist.xlsx')

#%% Filter by regional district

ind=np.where(z['rd']['Data']==meta['LUT']['TA_REGIONAL_DISTRICTS_SVW']['REGIONAL_DISTRICT_NAME']['CAPITAL'])
for k in z.keys():
    z[k]['Data']=z[k]['Data'][ind]

#%% Land use change analysis - by regional district (based on compressed categories)

id=np.array(list(meta['LUT']['Derived']['lcc_cec_c'].values()))

rs={}
for k in meta['LUT']['TA_REGIONAL_DISTRICTS_SVW']['REGIONAL_DISTRICT_NAME'].keys():
    # Areas at start and end
    vI=np.zeros( (2,id.size) )
    for i in range(id.size):
        ind=np.where( (z['lcc_cec_2010']['Data']==id[i]) & (z['harv_yr_con1']['Data']==0) & (z['fire_yr']['Data']==0) & (z['ibm_yr']['Data']==0) & (z['rd']['Data']==meta['LUT']['TA_REGIONAL_DISTRICTS_SVW']['REGIONAL_DISTRICT_NAME'][k]) )
        vI[0,i]=ind[0].size
        ind=np.where( (z['lcc_cec_2020']['Data']==id[i]) & (z['harv_yr_con1']['Data']==0) & (z['fire_yr']['Data']==0) & (z['ibm_yr']['Data']==0) & (z['rd']['Data']==meta['LUT']['TA_REGIONAL_DISTRICTS_SVW']['REGIONAL_DISTRICT_NAME'][k]) )
        vI[1,i]=ind[0].size
    df=pd.DataFrame(vI)
    df.to_excel(r'C:\Users\rhember\OneDrive - Government of BC\Manuscripts\Land Use Change Assessment\CEC Change Analysis\cec_c0_nodist_' + k + '.xlsx')

    # Transitions
    vT=np.zeros( (id.size,id.size) )
    for i in range(id.size):
        for j in range(id.size):
            ind=np.where( (z['lcc_cec_2010']['Data']==id[i]) & (z['lcc_cec_2020']['Data']==id[j]) & (z['harv_yr_con1']['Data']==0) & (z['fire_yr']['Data']==0) & (z['ibm_yr']['Data']==0) & (z['rd']['Data']==meta['LUT']['TA_REGIONAL_DISTRICTS_SVW']['REGIONAL_DISTRICT_NAME'][k]) )
            vT[i,j]=ind[0].size
    df=pd.DataFrame(vT)
    df.to_excel(r'C:\Users\rhember\OneDrive - Government of BC\Manuscripts\Land Use Change Assessment\CEC Change Analysis\cec_c1_nodist_' + k + '.xlsx')

    # Net deforestation summary
    ind=np.array([meta['LUT']['Derived']['lcc_cec_c']['Cropland'],meta['LUT']['Derived']['lcc_cec_c']['Barren Ground'],meta['LUT']['Derived']['lcc_cec_c']['Urban']],dtype='int8')
    Ad=np.sum(vT[0,ind-1])
    Aa=np.sum(vT[ind-1,0])
    Pd=np.sum(vT[0,ind-1])/vI[0,0]*100
    Pa=np.sum(vT[ind-1,0])/vI[0,0]*100
    Pn=Pa-Pd

#%% Map deforestation
# Notes:
#   1) Forest to Grassland and Forest to Shrubland excluded due to suspicion of error
#   2) Areas with harvest excluded (harvest from NTEM)

z1=zRef.copy()
z1['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
ind=np.where( (z['lcc_cec_2010']['Data']==meta['LUT']['Derived']['lcc_cec_c']['Forest']) & (z['lcc_cec_2020']['Data']==meta['LUT']['Derived']['lcc_cec_c']['Urban']) & (z['harv_yr_con1']['Data']==0) )
z1['Data'][ind]=1
ind=np.where( (z['lcc_cec_2010']['Data']==meta['LUT']['Derived']['lcc_cec_c']['Forest']) & (z['lcc_cec_2020']['Data']==meta['LUT']['Derived']['lcc_cec_c']['Cropland']) & (z['harv_yr_con1']['Data']==0) )
z1['Data'][ind]=2
ind=np.where( (z['lcc_cec_2010']['Data']==meta['LUT']['Derived']['lcc_cec_c']['Forest']) & (z['lcc_cec_2020']['Data']==meta['LUT']['Derived']['lcc_cec_c']['Barren Ground']) & (z['harv_yr_con1']['Data']==0) )
z1['Data'][ind]=3
#plt.matshow(z1['Data'])
gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\Disturbances\\LandUseChange_Deforestation_10to20_CEC.tif')

#%% Map Afforestation (from cropland, urban, or barren ground)

z1=zRef.copy()
z1['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
ind=np.where( (z['lcc_cec_2010']['Data']==meta['LUT']['Derived']['lcc_cec_c']['Urban']) & (z['lcc_cec_2020']['Data']==meta['LUT']['Derived']['lcc_cec_c']['Forest']) & (z['harv_yr_con1']['Data']==0) )
z1['Data'][ind]=1
ind=np.where( (z['lcc_cec_2010']['Data']==meta['LUT']['Derived']['lcc_cec_c']['Cropland']) & (z['lcc_cec_2020']['Data']==meta['LUT']['Derived']['lcc_cec_c']['Forest']) & (z['harv_yr_con1']['Data']==0) )
z1['Data'][ind]=2
ind=np.where( (z['lcc_cec_2010']['Data']==meta['LUT']['Derived']['lcc_cec_c']['Barren Ground']) & (z['lcc_cec_2020']['Data']==meta['LUT']['Derived']['lcc_cec_c']['Forest']) & (z['harv_yr_con1']['Data']==0) )
z1['Data'][ind]=3
gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\Disturbances\\LandUseChange_Afforestation_10to20_CEC.tif')

# ind=np.where( (z['harv_yr_con1']['Data']>0) )
# z['Data'][ind]=2
# plt.matshow(z['Data'])

# ind=np.where(z['Data']==1)
# ind[1].size/1000000

#%% Map Forest to Grassland and Forest to Shrubland
# *** Way too much to be real! ***

z=zRef.copy()
z['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
ind=np.where( (z['lcc_cec_2010']['Data']==meta['LUT']['Forest']) & (z['lcc_cec_2020']['Data']==meta['LUT']['Grassland']) )
z['Data'][ind]=1
ind=np.where( (z['lcc_cec_2010']['Data']==meta['LUT']['Forest']) & (z['lcc_cec_2020']['Data']==meta['LUT']['Shrubland']) )
z['Data'][ind]=1
ind=np.where( (z['harv_yr_con1']['Data']>0) )
z['Data'][ind]=2
plt.matshow(z['Data'])







ind=np.where(z['harv_yr_con1']['Data']>0)
ind[1].size/1000000

#%%

tv=np.arange(2001,2022,1)
#List=[meta['LUT']['Derived']['lcc_cec_c']['Urban'],meta['LUT']['Derived']['lcc_cec_c']['Barren Ground'],meta['LUT']['Derived']['lcc_cec_c']['Cropland'],meta['LUT']['Derived']['lcc_cec_c']['Grassland']]
List=[meta['LUT']['Derived']['lcc_cec_c']['Urban'],meta['LUT']['Derived']['lcc_cec_c']['Cropland']]
A1=np.zeros(tv.size)
A2=np.zeros(tv.size)
A3=np.zeros(tv.size)
for iT in range(tv.size):
    print(tv[iT])
    ind=np.where( (z['gfcly']['Data']==tv[iT]) )[0]
    A1[iT]=ind.size
    ind=np.where( (z['gfcly_filt']['Data']==tv[iT]) )[0]
    A2[iT]=ind.size
    #ind=np.where( (z['gfcly']['Data']==tv[iT]) & (np.abs(z['gfcly']['Data']-z['harv_yr_con1']['Data'])>5) & (z['lcc1_c']['Data']==meta['LUT']['Derived']['lcc1']['Forest']) )[0]
    ind=np.where( (z['gfcly']['Data']==tv[iT]) & (np.isin(z['lcc_cec_2020']['Data'],List)==True) )[0]
    A3[iT]=ind.size

plt.close('all')
plt.plot(tv,A1,'-bo')
plt.plot(tv,A2,'-gs')
plt.plot(tv,A3,'-cd')

#%% Import data

# gdf=u1ha.Import_GDBs_ProvinceWide_Simple()

# zRef=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\BC_Land_Mask.tif')
# zRD=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\RegionalDistrict.tif')
# z['harv_yr_con1']=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\NTEM_harv85to20_bc1ha.tif')
# z['fire_yr']=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\PROT_HISTORICAL_FIRE_POLYS_SP_All.tif')
# z['ibm_yr']=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\IBM_Mask_ExcTrace.tif')

# dCEC=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\BC1ha\LandUseLandCover\luc_cec_lc.xlsx')
# z['lcc_cec_2010']=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\LandUseLandCover\CEC_LandCover2010.tif')
# z['lcc_cec_2020']=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\LandUseLandCover\CEC_LandCover2020.tif')