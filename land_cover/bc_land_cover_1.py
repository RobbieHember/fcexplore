
#%% Import modules

import numpy as np
import matplotlib.pyplot as plt
import fcgadgets.macgyver.utilities_general as gu
import fcgadgets.macgyver.utilities_gis as gis
from fcgadgets.bc1ha import bc1ha_utilities as u1ha

#%% Import parameters

meta=u1ha.Init()
meta=u1ha.ImportLUTs(meta)
gp=gu.SetGraphics('Manuscript')

#%% Import data

zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
zLC2=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\BCLCS_LEVEL_2.tif')
zLC4=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\BCLCS_LEVEL_4.tif')
zWet=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\FWA_WETLANDS_POLY\\WATERBODY_TYPE.tif')
zLCC_CEC10=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCoverClass_CEC_2010_Compressed.tif')
zLCC_CEC20=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCoverClass_CEC_2020_Compressed.tif')
zLCC_NTEM=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCoverClass_NTEMS_2019.tif')
zH=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\VEG_CONSOLIDATED_CUT_BLOCKS_SP_YearLast.tif')
zH_NTEM=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\Harvest_NTEM_Year.tif')
#zBTM=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\landuse.btm.tif')

#%% Classify land cover

# In the past, LC2 treed area has differed from that of LC4 (treatment of TFLs and private land), but that appears to be fixed in 2023

vNam='lcc1'

zLCC1=zRef.copy()
zLCC1['Data']=np.zeros(zLCC1['Data'].shape,dtype='int8')

ind=np.where( (np.isin(zLCC_NTEM['Data'],[meta['LUT']['Derived']['lcc_ntems']['Coniferous'],meta['LUT']['Derived']['lcc_ntems']['Broadleaf'],meta['LUT']['Derived']['lcc_ntems']['Mixedwood'],meta['LUT']['Derived']['lcc_ntems']['Wetland-treed']])==True) )
zLCC1['Data'][ind]=meta['LUT']['Derived'][vNam]['Forest']

ind=np.where( (np.isin(zLCC_NTEM['Data'],[meta['LUT']['Derived']['lcc_ntems']['Shrubs']])==True) )
zLCC1['Data'][ind]=meta['LUT']['Derived'][vNam]['Shrub']

ind=np.where( (np.isin(zLCC_NTEM['Data'],[meta['LUT']['Derived']['lcc_ntems']['Herbs']])==True) )
zLCC1['Data'][ind]=meta['LUT']['Derived'][vNam]['Herb']

ind=np.where( (np.isin(zLCC_NTEM['Data'],[meta['LUT']['Derived']['lcc_ntems']['Wetland']])==True) )
zLCC1['Data'][ind]=meta['LUT']['Derived'][vNam]['Wetland']

ind=np.where( (np.isin(zLCC_NTEM['Data'],[meta['LUT']['Derived']['lcc_ntems']['Bryoids']])==True) )
zLCC1['Data'][ind]=meta['LUT']['Derived'][vNam]['Bryoid']

ind=np.where( (np.isin(zLCC_NTEM['Data'],[meta['LUT']['Derived']['lcc_ntems']['Rock and rubble'],meta['LUT']['Derived']['lcc_ntems']['Exposed barren land']])==True) )
zLCC1['Data'][ind]=meta['LUT']['Derived'][vNam]['Earth and Rock']

ind=np.where( (np.isin(zLCC_NTEM['Data'],[meta['LUT']['Derived']['lcc_ntems']['Snow and ice']])==True) )
zLCC1['Data'][ind]=meta['LUT']['Derived'][vNam]['Snow and Ice']

ind=np.where( (np.isin(zLCC_NTEM['Data'],[meta['LUT']['Derived']['lcc_ntems']['Water']])==True) )
zLCC1['Data'][ind]=meta['LUT']['Derived'][vNam]['Water']

# If NTEM says harvest and CEC says not settlement, assume forest
ind=np.where( (zH_NTEM['Data']>0) & (zLCC_CEC20['Data']!=meta['LUT']['Derived']['lcc_cec_c']['Urban']) & (zLCC_CEC20['Data']!=meta['LUT']['Derived']['lcc_cec_c']['Cropland']) )
zLCC1['Data'][ind]=meta['LUT']['Derived'][vNam]['Forest']

# Override NTEM with FWA wetlands
ind=np.where( (zWet['Data']==meta['LUT']['FWA_WETLANDS_POLY']['WATERBODY_TYPE']['W']) )
zLCC1['Data'][ind]=meta['LUT']['Derived'][vNam]['Wetland']

# Override with forest from VRI LC4
ind=np.where( (zLC4['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['TC']) | (zLC4['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['TM']) | (zLC4['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['TB']) )
zLCC1['Data'][ind]=meta['LUT']['Derived'][vNam]['Forest']

# Convert areas  where CEC map says urban
ind=np.where( (zLCC_CEC20['Data']==meta['LUT']['Derived']['lcc_cec_c']['Urban']) )
zLCC1['Data'][ind]=meta['LUT']['Derived'][vNam]['Settlement']

# Check for any cells left unclassified
# Very small number, likely water, reclassify as water
a=np.zeros(zLCC1['Data'].shape,dtype='int8');
ind=np.where( (zRef['Data']==1) & (zLCC1['Data']==0) ); a[ind]=1
print(ind[0].size)
zLCC1['Data'][ind]=meta['LUT']['Derived'][vNam]['Water']

# Reclassify outside land mask as water
ind=np.where( (zRef['Data']==0) )
zLCC1['Data'][ind]=meta['LUT']['Derived'][vNam]['Water']

plt.close('all')
plt.matshow(zLCC1['Data'])

# Save
gis.SaveGeoTiff(zLCC1,meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCoverClass1.tif')

#%% Summarize areas

labL=list(meta['LUT']['Derived'][vNam].keys())
A=np.zeros(len(labL))
cnt=0
for k in meta['LUT']['Derived'][vNam].keys():
    ind=np.where( (zRef['Data']==1) & (zLCC1['Data']==meta['LUT']['Derived'][vNam][k]) )
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

ind=np.where( (zLCC1['Data']==meta['LUT']['Derived'][vNam]['Forest']) )
A_LCC1=ind[0].size

ind=np.where( (zLC2['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_2']['T']) )
A_treed_lc2=ind[0].size

ind=np.where( (zLC4['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['TC']) | (zLC4['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['TM']) | (zLC4['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_4']['TB']) )
A_treed_lc4=ind[0].size

ind=np.where( (np.isin(zLCC_NTEM['Data'],[meta['LUT']['Derived']['lcc_ntems']['Coniferous'],meta['LUT']['Derived']['lcc_ntems']['Broadleaf'],meta['LUT']['Derived']['lcc_ntems']['Mixedwood'],meta['LUT']['Derived']['lcc_ntems']['Wetland-treed']])==True) )
A_NTEM=ind[0].size

ind=np.where( (zLCC_CEC20['Data']==meta['LUT']['Derived']['lcc_cec_c']['Forest']) )
A_CEC=ind[0].size

print(A_LCC1/1e6)
print(A_treed_lc2/1e6)
print(A_treed_lc4/1e6)
print(A_NTEM/1e6)
print(A_CEC/1e6)

#%% Disagreement in forest class

ind=np.where( (zLC2['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_2']['T']) & \
             (np.isin(zLCC_NTEM['Data'],[meta['LUT']['Derived']['lcc_ntems']['Coniferous'],meta['LUT']['Derived']['lcc_ntems']['Broadleaf'],meta['LUT']['Derived']['lcc_ntems']['Mixedwood'],meta['LUT']['Derived']['lcc_ntems']['Wetland-treed']])==True) )
print(ind[0].size/1e6)

ind=np.where( (zLC2['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_2']['T']) & \
             (np.isin(zLCC_NTEM['Data'],[meta['LUT']['Derived']['lcc_ntems']['Coniferous'],meta['LUT']['Derived']['lcc_ntems']['Broadleaf'],meta['LUT']['Derived']['lcc_ntems']['Mixedwood'],meta['LUT']['Derived']['lcc_ntems']['Wetland-treed']])==False) )
print(ind[0].size/1e6)

ind=np.where( (zLC2['Data']!=meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_2']['T']) & \
             (np.isin(zLCC_NTEM['Data'],[meta['LUT']['Derived']['lcc_ntems']['Coniferous'],meta['LUT']['Derived']['lcc_ntems']['Broadleaf'],meta['LUT']['Derived']['lcc_ntems']['Mixedwood'],meta['LUT']['Derived']['lcc_ntems']['Wetland-treed']])==True) )
print(ind[0].size/1e6)

#%%

a=np.zeros(zLCC1['Data'].shape,dtype='int8');
#ind=np.where( (zLCC1['Data']==1) & (zLC2['Data']!=meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_2']['T']) )
#a[ind]=1

ind=np.where( (zLC2['Data']!=meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_2']['T']) & \
             (zLCC_CEC20['Data']==meta['LUT']['Derived']['lcc_cec_c']['Forest']) & \
             (np.isin(zLCC_NTEM['Data'],[meta['LUT']['Derived']['lcc_ntems']['Coniferous'],meta['LUT']['Derived']['lcc_ntems']['Broadleaf'],meta['LUT']['Derived']['lcc_ntems']['Mixedwood'],meta['LUT']['Derived']['lcc_ntems']['Wetland-treed']])==True) )
a[ind]=1

plt.close('all')
plt.matshow(a)

labL=list(meta['LUT']['Derived']['lcc_cec_c'])
d={}
for k,i in meta['LUT']['Derived']['lcc_cec_c'].items():
    ind=np.where( (zLCC_CEC20['Data']==i) & (zLC2['Data']!=meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_2']['T']) & (np.isin(zLCC_NTEM['Data'],[meta['LUT']['Derived']['lcc_ntems']['Coniferous'],meta['LUT']['Derived']['lcc_ntems']['Broadleaf'],meta['LUT']['Derived']['lcc_ntems']['Mixedwood'],meta['LUT']['Derived']['lcc_ntems']['Wetland-treed']])==True) )
    d[k]=ind[0].size/1e6

#%% Generate a random sample of points from areas with disagreement in forest

import pyproj
import geopandas as gpd
from shapely.geometry import Polygon,Point

ind=np.where( (zLC2['Data']!=meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_2']['T']) & \
             (zLCC_CEC20['Data']==meta['LUT']['Derived']['lcc_cec_c']['Forest']) & \
             (np.isin(zLCC_NTEM['Data'],[meta['LUT']['Derived']['lcc_ntems']['Coniferous'],meta['LUT']['Derived']['lcc_ntems']['Broadleaf'],meta['LUT']['Derived']['lcc_ntems']['Mixedwood'],meta['LUT']['Derived']['lcc_ntems']['Wetland-treed']])==True) )

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





