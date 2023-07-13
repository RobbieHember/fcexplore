
#%% Import modules

import numpy as np
import matplotlib.pyplot as plt
import fcgadgets.macgyver.utilities_general as gu
import fcgadgets.macgyver.utilities_gis as gis
import fcgadgets.bc1ha.bc1ha_utilities as u1ha

#%% Import parameters

meta=u1ha.Init()
meta=u1ha.ImportLUTs(meta)
gp=gu.SetGraphics('Manuscript')

#%% Import data

zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])

vList=['lc2','lc4','wbt','lcc_cec_2010','lcc_cec_2020','lcc_ntem_2019','harv_yr_cc','harv_yr_ntem']
z=u1ha.Import_Raster(meta,[],vList)

#%% Land cover class 1 (current)

# In the past, LC2 treed area has differed from that of LC4 (treatment of TFLs and private land), but that appears to be fixed in 2023

zLCC1=zRef.copy()
zLCC1['Data']=np.zeros(zRef['Data'].shape,dtype='int8')

ind=np.where( (np.isin(z['lcc_ntem_2019']['Data'],[meta['LUT']['Derived']['lcc_ntem']['Coniferous'],meta['LUT']['Derived']['lcc_ntem']['Broadleaf'],meta['LUT']['Derived']['lcc_ntem']['Mixedwood'],meta['LUT']['Derived']['lcc_ntem']['Wetland-treed']])==True) )
zLCC1['Data'][ind]=meta['LUT']['Derived']['lcc1']['Forest']

ind=np.where( (np.isin(z['lcc_ntem_2019']['Data'],[meta['LUT']['Derived']['lcc_ntem']['Shrubs']])==True) )
zLCC1['Data'][ind]=meta['LUT']['Derived']['lcc1']['Shrub']

ind=np.where( (np.isin(z['lcc_ntem_2019']['Data'],[meta['LUT']['Derived']['lcc_ntem']['Herbs']])==True) )
zLCC1['Data'][ind]=meta['LUT']['Derived']['lcc1']['Herb']

ind=np.where( (np.isin(z['lcc_ntem_2019']['Data'],[meta['LUT']['Derived']['lcc_ntem']['Wetland']])==True) )
zLCC1['Data'][ind]=meta['LUT']['Derived']['lcc1']['Wetland']

ind=np.where( (np.isin(z['lcc_ntem_2019']['Data'],[meta['LUT']['Derived']['lcc_ntem']['Bryoids']])==True) )
zLCC1['Data'][ind]=meta['LUT']['Derived']['lcc1']['Bryoid']

ind=np.where( (np.isin(z['lcc_ntem_2019']['Data'],[meta['LUT']['Derived']['lcc_ntem']['Rock and rubble'],meta['LUT']['Derived']['lcc_ntem']['Exposed barren land']])==True) )
zLCC1['Data'][ind]=meta['LUT']['Derived']['lcc1']['Earth and Rock']

ind=np.where( (np.isin(z['lcc_ntem_2019']['Data'],[meta['LUT']['Derived']['lcc_ntem']['Snow and ice']])==True) )
zLCC1['Data'][ind]=meta['LUT']['Derived']['lcc1']['Snow and Ice']

ind=np.where( (np.isin(z['lcc_ntem_2019']['Data'],[meta['LUT']['Derived']['lcc_ntem']['Water']])==True) )
zLCC1['Data'][ind]=meta['LUT']['Derived']['lcc1']['Water']

# If NTEM says harvest and CEC says not settlement, assume forest
ind=np.where( (z['harv_yr_ntem']['Data']>0) & (z['lcc_cec_2020']['Data']!=meta['LUT']['Derived']['lcc_cec_c']['Urban']) & (z['lcc_cec_2020']['Data']!=meta['LUT']['Derived']['lcc_cec_c']['Cropland']) )
zLCC1['Data'][ind]=meta['LUT']['Derived']['lcc1']['Forest']

# Override NTEM with FWA wetlands
ind=np.where( (z['wbt']['Data']==meta['LUT']['FWA_WETLANDS_POLY']['WATERBODY_TYPE']['W']) )
zLCC1['Data'][ind]=meta['LUT']['Derived']['lcc1']['Wetland']

# Override with forest from VRI LC4
ind=np.where( (z['lc4']['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['TC']) | (z['lc4']['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['TM']) | (z['lc4']['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['TB']) )
zLCC1['Data'][ind]=meta['LUT']['Derived']['lcc1']['Forest']

# Convert areas  where CEC map says urban
ind=np.where( (z['lcc_cec_2020']['Data']==meta['LUT']['Derived']['lcc_cec_c']['Urban']) )
zLCC1['Data'][ind]=meta['LUT']['Derived']['lcc1']['Settlement']

# CEC is more reliable indicator of cropland -> convert to shrubland
ind=np.where( (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Herb']) & (z['lcc_cec_2020']['Data']!=meta['LUT']['Derived']['lcc_cec_c']['Cropland']) )
zLCC1['Data'][ind]=meta['LUT']['Derived']['lcc1']['Shrub']

# Check for any cells left unclassified
# Very small number, likely water, reclassify as water
a=np.zeros(zLCC1['Data'].shape,dtype='int8');
ind=np.where( (zRef['Data']==1) & (zLCC1['Data']==0) ); a[ind]=1
print(ind[0].size)
zLCC1['Data'][ind]=meta['LUT']['Derived']['lcc1']['Water']

# Reclassify outside land mask as water
ind=np.where( (zRef['Data']==0) )
zLCC1['Data'][ind]=meta['LUT']['Derived']['lcc1']['Water']

plt.close('all')
plt.matshow(zLCC1['Data'])

# Save
gis.SaveGeoTiff(zLCC1,meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCoverClass1_Current.tif')

#%%

a=zRef['Data'].copy()
ind=np.where(zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Herb'])
a[ind]=2
plt.matshow(a)


a=zRef['Data'].copy()
ind=np.where(z['lcc_cec_2020']['Data']==meta['LUT']['Derived']['lcc_cec_c']['Cropland'])
a[ind]=2
plt.matshow(a)


#%% Land cover class 1 (pre-contact)

zLCC1_pc=zLCC1.copy()
ind=np.where( (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Settlement']) | (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Herb']) )
zLCC1_pc['Data'][ind]=meta['LUT']['Derived']['lcc1']['Forest']
gis.SaveGeoTiff(zLCC1_pc,meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCoverClass1_PreContact.tif')

#%% Summarize areas

labL=list(meta['LUT']['Derived']['lcc1'].keys())
A=np.zeros(len(labL))
cnt=0
for k in meta['LUT']['Derived']['lcc1'].keys():
    ind=np.where( (zRef['Data']==1) & (zLCC1['Data']==meta['LUT']['Derived']['lcc1'][k]) )
    A[cnt]=ind[0].size/1e6
    cnt=cnt+1

# Plot bar chart of area
x=np.arange(0,A.size,1)
plt.close('all');
fig,ax=plt.subplots(1,figsize=gu.cm2inch(10,6.5));
ax.bar(x,A)
ax.set(xticks=x,xticklabels=labL,ylabel='Area (Mha)')
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])

#%% Compare area treed from L2 and L4

ind=np.where( (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) )
A_LCC1=ind[0].size

ind=np.where( (z['lc2']['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_2']['T']) )
A_treed_lc2=ind[0].size

ind=np.where( (z['lc4']['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['TC']) | (z['lc4']['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['TM']) | (z['lc4']['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['TB']) )
A_treed_lc4=ind[0].size

ind=np.where( (np.isin(z['lcc_ntem_2019']['Data'],[meta['LUT']['Derived']['lcc_ntem']['Coniferous'],meta['LUT']['Derived']['lcc_ntem']['Broadleaf'],meta['LUT']['Derived']['lcc_ntem']['Mixedwood'],meta['LUT']['Derived']['lcc_ntem']['Wetland-treed']])==True) )
A_NTEM=ind[0].size

ind=np.where( (z['lcc_cec_2020']['Data']==meta['LUT']['Derived']['lcc_cec_c']['Forest']) )
A_CEC=ind[0].size

print(A_LCC1/1e6)
print(A_treed_lc2/1e6)
print(A_treed_lc4/1e6)
print(A_NTEM/1e6)
print(A_CEC/1e6)

#%% Disagreement in forest class

ind=np.where( (z['lc2']['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_2']['T']) & \
             (np.isin(z['lcc_ntem_2019']['Data'],[meta['LUT']['Derived']['lcc_ntem']['Coniferous'],meta['LUT']['Derived']['lcc_ntem']['Broadleaf'],meta['LUT']['Derived']['lcc_ntem']['Mixedwood'],meta['LUT']['Derived']['lcc_ntem']['Wetland-treed']])==True) )
print(ind[0].size/1e6)

ind=np.where( (z['lc2']['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_2']['T']) & \
             (np.isin(z['lcc_ntem_2019']['Data'],[meta['LUT']['Derived']['lcc_ntem']['Coniferous'],meta['LUT']['Derived']['lcc_ntem']['Broadleaf'],meta['LUT']['Derived']['lcc_ntem']['Mixedwood'],meta['LUT']['Derived']['lcc_ntem']['Wetland-treed']])==False) )
print(ind[0].size/1e6)

ind=np.where( (z['lc2']['Data']!=meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_2']['T']) & \
             (np.isin(z['lcc_ntem_2019']['Data'],[meta['LUT']['Derived']['lcc_ntem']['Coniferous'],meta['LUT']['Derived']['lcc_ntem']['Broadleaf'],meta['LUT']['Derived']['lcc_ntem']['Mixedwood'],meta['LUT']['Derived']['lcc_ntem']['Wetland-treed']])==True) )
print(ind[0].size/1e6)

#%%

a=np.zeros(zLCC1['Data'].shape,dtype='int8');
#ind=np.where( (zLCC1['Data']==1) & (z['lc2']['Data']!=meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_2']['T']) )
#a[ind]=1

ind=np.where( (z['lc2']['Data']!=meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_2']['T']) & \
             (z['lcc_cec_2020']['Data']==meta['LUT']['Derived']['lcc_cec_c']['Forest']) & \
             (np.isin(z['lcc_ntem_2019']['Data'],[meta['LUT']['Derived']['lcc_ntem']['Coniferous'],meta['LUT']['Derived']['lcc_ntem']['Broadleaf'],meta['LUT']['Derived']['lcc_ntem']['Mixedwood'],meta['LUT']['Derived']['lcc_ntem']['Wetland-treed']])==True) )
a[ind]=1

plt.close('all')
plt.matshow(a)

labL=list(meta['LUT']['Derived']['lcc_cec_c'])
d={}
for k,i in meta['LUT']['Derived']['lcc_cec_c'].items():
    ind=np.where( (z['lcc_cec_2020']['Data']==i) & (z['lc2']['Data']!=meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_2']['T']) & (np.isin(z['lcc_ntem_2019']['Data'],[meta['LUT']['Derived']['lcc_ntem']['Coniferous'],meta['LUT']['Derived']['lcc_ntem']['Broadleaf'],meta['LUT']['Derived']['lcc_ntem']['Mixedwood'],meta['LUT']['Derived']['lcc_ntem']['Wetland-treed']])==True) )
    d[k]=ind[0].size/1e6

#%% Generate a random sample of points from areas with disagreement in forest

import pyproj
import geopandas as gpd
from shapely.geometry import Polygon,Point

ind=np.where( (z['lc2']['Data']!=meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_2']['T']) & \
             (z['lcc_cec_2020']['Data']==meta['LUT']['Derived']['lcc_cec_c']['Forest']) & \
             (np.isin(z['lcc_ntem_2019']['Data'],[meta['LUT']['Derived']['lcc_ntem']['Coniferous'],meta['LUT']['Derived']['lcc_ntem']['Broadleaf'],meta['LUT']['Derived']['lcc_ntem']['Mixedwood'],meta['LUT']['Derived']['lcc_ntem']['Wetland-treed']])==True) )

srs=gis.ImportSRSs()

lon,lat=pyproj.transform(srs['Proj']['BC1ha'],srs['Proj']['Geographic'],zRef['X'][ind],zRef['Y'][ind])
lon_s=lon[0::30000]
lat_s=lat[0::30000]
print(lon_s.size)

flg=1
if flg==1:
    points=[]
    for k in range(lon_s.size):
        points.append(Point(lon_s[k],lat_s[k]))
    gdf_xy=gpd.GeoDataFrame({'geometry':points,'ID':np.arange(lon_s.size,dtype='int32')})
    #gdf_xy.crs=gdf_bc_boundary.crs
    gdf_xy.to_file(meta['Paths']['bc1ha'] + '\\forest.geojson',driver='GeoJSON')





