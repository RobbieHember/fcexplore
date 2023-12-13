'''
IMPACTS OF HARVEST RETENTION ON FOREST CARBON STOCKS IN THE BOUNDARY TIMBER SUPPLY AREA
Datasets:
    - whole stem volume from LiDAR-driven Predictive Forest Invnetory (from Geoff Quinn and Chris Butson)
    - Consolidated cutblocks database (year of harvest)
    - Harvest retention compilation 1 (see look-up table for definitions)
'''

#%% Import modules
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib.colors import LightSource
from shapely.geometry import Polygon,Point
import copy
import time
from shapely.geometry import Polygon,Point,box,shape
from rasterio import features
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.bc1ha.bc1ha_util as u1ha
import fcgadgets.bc1ha.bc1ha_plot as p1ha
import fcgadgets.macgyver.util_query_gdb as qgdb

#%% Import parameters
meta=u1ha.Init()
meta['Graphics']['Map']['RGSF']=1
meta['Graphics']['Map']['Fig Width']=15.5
meta['Graphics']['Map']['Side Space']=0.25
meta['Graphics']['Map']['Map Position']=[0,0,1-meta['Graphics']['Map']['Side Space']-0.01,1]
meta['Graphics']['Map']['Map Axis Vis']='off'
meta['Graphics']['Map']['Map Grid Vis']=False
meta['Graphics']['Map']['Legend X']=1-meta['Graphics']['Map']['Side Space']+meta['Graphics']['Map']['Map Position'][0]#+0.01,0.6,0.03,0.35]
meta['Graphics']['Map']['Legend Width']=0.0275
meta['Graphics']['Map']['Legend Font Size']=7
meta['Graphics']['Map']['Legend Text Space']=0.035
meta['Graphics']['Map']['Show Bound Land Mask']='On'
meta['Graphics']['Map']['Show Bound Within']='On'
meta['Graphics']['Map']['Show Lakes']='Off'
meta['Graphics']['Map']['Show Rivers']='Off'
meta['Graphics']['Map']['Show Roads']='On'
meta['Graphics']['Plot Style']='Manuscript'
meta['Graphics']['gp']=gu.SetGraphics(meta['Graphics']['Plot Style'])
meta['Graphics']['Print Figures']='On'
#meta['Graphics']['Print Figure Path']=r'C:\Users\rhember\OneDrive - Government of BC\Figures\\LICS'

meta['Paths']['Working']=r'C:\Data\BC20m\Boundary TSA'

#%% Define region of interest
roi={}
roi['Type']='ByTSA'
roi['Name']='Boundary TSA'
roi['List']=['Boundary TSA']
gdf=u1ha.Import_GDBs_ProvinceWide(meta)
roi=u1ha.DefineROI(meta,roi,gdf)

#%% Override the 1ha defaults with the Predictive Forest Inentory dataset from Geoff Quinn
zV=gis.OpenGeoTiff(r'C:\Data\LiDAR\Boundary TSA\From GQ 20230302\BoundaryTSA_PFI_LnTotalWSV.tif')
zV=gis.ClipRasterByXYLimits(zV,roi['grd']['xlim'],roi['grd']['ylim'])
roi['grd']=copy.deepcopy(zV)

# Convert to stemwood biomass (tC/ha)
zB=copy.deepcopy(zV)
zB['Data']=0.5*0.45*zB['Data']

#%% Rasterize TSA boundary at 20 m
df0=roi['gdf']['bound within']
df0=df0[df0.geometry!=None]
df0=df0.reset_index()
df0['ID']=1
shapes=((geom,value) for geom, value in zip(df0.geometry,df0['ID']))
z0=np.zeros(zV['Data'].shape,dtype=float)
burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=roi['grd']['Transform'])
roi['grd']['Data']=burned.astype('int8')

# Remove exterior from volume
ind=np.where(roi['grd']['Data']==0); zV['Data'][ind]=0
ind=np.where(zV['Data']<0); zV['Data'][ind]=-99
zV['Data']=zV['Data'].astype('int16')
gis.SaveGeoTiff(zV,meta['Paths']['Working'] + '\\VolumeWholeStem.tif')

#%% Import forest cover reserves
roi=u1ha.Import_GDB_Over_ROI(meta,roi,['fcres'])

#%% Plot ROI
cm=plt.cm.get_cmap('viridis',30)
plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
im=ax[0].matshow(zV['Data'],extent=roi['grd']['Extent'],clim=[0,700],cmap=cm)
roi=p1ha.PlotVectorBaseMaps(meta,roi,ax[0])
ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])
#roi['gdf']['H_CC'].plot(ax=ax[0],edgecolor=[0.85,0.35,1],facecolor='None',linewidth=1)
#gdf1['gdf']['fcres']['gdf'].plot(ax=ax[0],linestyle='--',edgecolor=[0.6,1,0],facecolor='None',linewidth=1)
cb=plt.colorbar(im,cax=ax[1],cmap=cm)
ax[1].set(position=[0.76,0.5,0.03,0.35])

#%% Rasterize harvest cutblocks
lNam='VEG_CONSOLIDATED_CUT_BLOCKS_SP'
vNam='HARVEST_YEAR'
ind=np.where( (meta['Geos']['Variable Info']['Layer Name']==lNam) & (meta['Geos']['Variable Info']['Variable Name']==vNam) )[0]
pthin=meta['Paths']['GDB'][ meta['Geos']['Variable Info']['Geodatabase Name'][ind[0]] ]
df=gpd.read_file(pthin,layer=lNam)
df=df[df.geometry!=None]
df=df.reset_index()

# Clip to ROI
roi['gdf']['H_CC']=df.cx[roi['grd']['xmin']:roi['grd']['xmax'],roi['grd']['ymin']:roi['grd']['ymax']]
roi['gdf']['H_CC']=roi['gdf']['H_CC'].reset_index(drop=True)
roi['gdf']['H_CC']=gpd.overlay(roi['gdf']['H_CC'],roi['gdf']['bound within'],how='intersection')

# Rasterize
zYearLast=copy.deepcopy(roi['grd'])
zYearLast['Data']=np.zeros(roi['grd']['Data'].shape,dtype='int16')
uYear=roi['gdf']['H_CC'][vNam].unique()
tv=np.arange(np.min(uYear),np.max(uYear),1)
for iT in range(tv.size):
    print(tv[iT])
    df0=roi['gdf']['H_CC'][roi['gdf']['H_CC'][vNam]==tv[iT]].copy()
    shapes=((geom,value) for geom, value in zip(df0.geometry,df0[vNam]))
    z0=np.zeros(roi['grd']['Data'].shape,dtype=float)
    if len(df0)>0:
        burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=roi['grd']['Transform'])
    z1=copy.deepcopy(roi['grd'])
    z1['Data']=z0.astype(meta['Geos']['Variable Info']['Precision'][ind[0]])
    #gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(tv[iT]) + '.tif')
    zYearLast['Data'][burned>0]=tv[iT]
gis.SaveGeoTiff(zYearLast,meta['Paths']['Working'] + '\\Harvest_CC_YearLast.tif')
#plt.close('all'); plt.matshow(zYearLast['Data'],clim=[1955,2022])

#%% Clip and resample harvest retention compilation 1
z1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Management\\HarvestRetentionComp1.tif') 
z2=gis.ClipRasterByXYLimits(z1,roi['grd']['xlim'],roi['grd']['ylim'])
pth=meta['Paths']['Working'] + '\\HarvestRetentionComp1.tif'
gis.SaveGeoTiff(z2,pth)
gis.ResampleRaster(pth,5)
z2=gis.OpenGeoTiff(pth) 
z2=gis.ClipToRaster(z2,roi['grd'])
gis.SaveGeoTiff(z2,pth)
plt.matshow(z2['Data'])

#%% Clip and resample land cover compilation 1 (2019)
z1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover_Comp1_2019.tif') 
z2=gis.ClipRasterByXYLimits(z1,roi['grd']['xlim'],roi['grd']['ylim'])
pth=meta['Paths']['Working'] + '\\LandCover_Comp1_2019.tif'
gis.SaveGeoTiff(z2,pth)
gis.ResampleRaster(pth,5)
z2=gis.OpenGeoTiff(pth) 
z2=gis.ClipToRaster(z2,roi['grd'])
gis.SaveGeoTiff(z2,pth)
#plt.matshow(z2['Data'])

#%% Clip and resample VRI volume
z1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\LIVE_STAND_VOLUME_125.tif') 
z2=gis.ClipRasterByXYLimits(z1,roi['grd']['xlim'],roi['grd']['ylim'])
pth=meta['Paths']['Working'] + '\\LIVE_STAND_VOLUME_125.tif'
gis.SaveGeoTiff(z2,pth)
gis.ResampleRaster(pth,5)
z2=gis.OpenGeoTiff(pth) 
z2=gis.ClipToRaster(z2,roi['grd'])
gis.SaveGeoTiff(z2,pth)
#plt.matshow(z2['Data'])

z1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\DEAD_STAND_VOLUME_125.tif') 
z2=gis.ClipRasterByXYLimits(z1,roi['grd']['xlim'],roi['grd']['ylim'])
pth=meta['Paths']['Working'] + '\\DEAD_STAND_VOLUME_125.tif'
gis.SaveGeoTiff(z2,pth)
gis.ResampleRaster(pth,5)
z2=gis.OpenGeoTiff(pth) 
z2=gis.ClipToRaster(z2,roi['grd'])
gis.SaveGeoTiff(z2,pth)
#plt.matshow(z2['Data'])

z1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\WHOLE_STEM_BIOMASS_PER_HA.tif') 
z2=gis.ClipRasterByXYLimits(z1,roi['grd']['xlim'],roi['grd']['ylim'])
pth=meta['Paths']['Working'] + '\\WHOLE_STEM_BIOMASS_PER_HA.tif'
gis.SaveGeoTiff(z2,pth)
gis.ResampleRaster(pth,5)
z2=gis.OpenGeoTiff(pth) 
z2=gis.ClipToRaster(z2,roi['grd'])
gis.SaveGeoTiff(z2,pth)
#plt.matshow(z2['Data'])

z1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\PROJ_AGE_1.tif') 
z2=gis.ClipRasterByXYLimits(z1,roi['grd']['xlim'],roi['grd']['ylim'])
pth=meta['Paths']['Working'] + '\\PROJ_AGE_1.tif'
gis.SaveGeoTiff(z2,pth)
gis.ResampleRaster(pth,5)
z2=gis.OpenGeoTiff(pth) 
z2=gis.ClipToRaster(z2,roi['grd'])
gis.SaveGeoTiff(z2,pth)
#plt.matshow(z2['Data'])

z1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\OGSR_TAP_PRIORITY_DEF_AREA_SP\\PRIORITY_DEFERRAL_ID.tif') 
z2=gis.ClipRasterByXYLimits(z1,roi['grd']['xlim'],roi['grd']['ylim'])
pth=meta['Paths']['Working'] + '\\PRIORITY_DEFERRAL_ID.tif'
gis.SaveGeoTiff(z2,pth)
gis.ResampleRaster(pth,5)
z2=gis.OpenGeoTiff(pth)
z2=gis.ClipToRaster(z2,roi['grd'])
gis.SaveGeoTiff(z2,pth)
#plt.matshow(z2['Data'])

#%% Import variables
zVL=gis.OpenGeoTiff(meta['Paths']['Working'] + '\\LIVE_STAND_VOLUME_125.tif') 
zVD=gis.OpenGeoTiff(meta['Paths']['Working'] + '\\DEAD_STAND_VOLUME_125.tif') 
zA_VRI=gis.OpenGeoTiff(meta['Paths']['Working'] + '\\PROJ_AGE_1.tif')
zB_VRI=gis.OpenGeoTiff(meta['Paths']['Working'] + '\\WHOLE_STEM_BIOMASS_PER_HA.tif')
zB_VRI['Data']=0.5*zB_VRI['Data']
zD=gis.OpenGeoTiff(meta['Paths']['Working'] + '\\PRIORITY_DEFERRAL_ID.tif')
zH=gis.OpenGeoTiff(meta['Paths']['Working'] + '\\Harvest_CC_YearLast.tif')
zR=gis.OpenGeoTiff(meta['Paths']['Working'] + '\\HarvestRetentionComp1.tif')

#%% Relationship between VRI volume and PFI volume

# ikp=np.where( (zV['Data']>1) & (zVL['Data']>1) )
# x=zVL['Data'][ikp][0::30]
# x2=zVL['Data'][ikp][0::30]+zVD['Data'][ikp][0::30]
# y=zV['Data'][ikp][0::30]

# plt.close('all')
# plt.plot([0,1200],[0,1200],'k-')
# plt.plot(x,y,'.')
# rs,txt=gu.GetRegStats(x,y)
# plt.plot(rs['xhat Line'],rs['yhat Line'],'r-')
# rs,txt=gu.GetRegStats(x2,y)
# plt.plot(rs['xhat Line'],rs['yhat Line'],'r--')

#%% Age responses 
ikp=np.where( (zB['Data']>1) & (zB_VRI['Data']>1) & (zA_VRI['Data']>0) & (zH['Data']>0) )
bw=5; bin=np.arange(5,260,bw)
x=zA_VRI['Data'][ikp][0::10]

plt.close('all'); 
fig,ax=plt.subplots(1,figsize=gu.cm2inch(11,8)); 
#plt.plot([0,300],[0,300],'k-')
y=zB['Data'][ikp][0::10]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
plt.plot(bin,mu,'o',label='Predictive Forest Inventory')
y=zB_VRI['Data'][ikp][0::10]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
plt.plot(bin,mu,'s',label='Vegetation Resource Inventory')
ax.set(ylabel='Stemwood biomass (tC/ha)',xlabel='Age, years',ylim=[0,150],xlim=[0,250])
ax.legend(loc='lower right',frameon=False)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=1.5)    
plt.tight_layout()

#%%

# ikp=np.where( (zD['Data']>0) & (zV['Data']>1) & (zVL['Data']>1) & (zA['Data']>0) )
# bw=10; bin=np.arange(0,260,10)
# x=zA['Data'][ikp][0::10]

# plt.close('all')
# y=0.45*zV['Data'][ikp][0::10]
# N,mu,med,sig,se=gu.discres(x,y,bw,bin)
# plt.plot(bin,mu,'o')
# bw=10; bin=np.arange(0,210,10)
# y=zB['Data'][ikp][0::10]
# N,mu,med,sig,se=gu.discres(x,y,bw,bin)
# plt.plot(bin,mu,'s')



#%% Plot chronosequence
a=zH['Data'][zH['Data']>0].flatten()[0::5]
plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(14,10))
ax.hist(a,np.arange(1950,2025,1),facecolor=[0.5,0.75,1])
ax.set(ylabel='Area harvested (hectares/year)',xlabel='Time, years')
ax.plot([2018,2018],[0,30000],'g--',color=[0.5,0,0],lw=3)
ax.legend(loc='lower right',frameon=False)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=1.5)    
plt.tight_layout()

#%% Plot chronosequence
fig,ax=plt.subplots(1,figsize=gu.cm2inch(14,10))
ax.hist(2018-a,np.arange(0,60,1),facecolor=[0.5,0.75,1])
ax.set(ylabel='Area harvested (hectares/year)',xlabel='Age at time of LiDAR (Years)')
ax.plot([0,0],[0,30000],'g--',color=[0.5,0,0],lw=3)

#%% Time since harvest
ikp=np.where( (zV['Data']>1) & (zVL['Data']>1) & (zH['Data']>0) )
bw=1; bin=np.arange(0,260,bw)
x=2018-zH['Data'][ikp][0::10]

dGY=gu.ReadExcel(r'C:\Users\rhember\OneDrive - Government of BC\Projects\MGEM 2023 Harvesting and Retention\TIPSY Age Response.xlsx')

plt.close('all'); 
fig,ax=plt.subplots(1,figsize=gu.cm2inch(11,8)); 
y=zB['Data'][ikp][0::10]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
plt.plot(bin,mu,'o',label='Predictive Forest Inventory')
y=zB_VRI['Data'][ikp][0::10]
#N,mu,med,sig,se=gu.discres(x,y,bw,bin)
#plt.plot(bin,mu,'s',label='Vegetation Resource Inventory')
#plt.plot(bin,1.*bin,'-',label='Model')
plt.plot(dGY['Age'],dGY['Stemwood biomass (tC/ha)'],'-',label='Model')
ax.set(ylabel='Stemwood biomass (tC/ha)',xlabel='Time since harvest, years',ylim=[0,90],xlim=[0,65])
ax.legend(loc='lower right',frameon=False)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=1.5)    
plt.tight_layout()

#%% OLD

# #%% Import modules

# import os
# import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.colors
# from matplotlib.collections import PatchCollection
# import geopandas as gpd
# import pandas as pd
# import copy
# import fiona
# import time
# from scipy.interpolate import griddata
# from shapely.geometry import Polygon,Point,box
# import fcgadgets.macgyver.util_general as gu
# import fcgadgets.macgyver.util_gis as gis
# from fcgadgets.cbrunner import cbrun_util as cbu
# from fcgadgets.bc1ha import bc1ha_util as u1ha

# #%%

# z=gis.OpenGeoTiff(r'C:\Data\LiDAR\Boundary TSA\From GQ 20230302\BoundaryTSA_PFI_LnTotalWSV.tif')

# z.keys()

# #%% Path management

# meta={}
# meta['Paths']={}
# meta['Paths']['BC1ha']=r'C:\Data\BC1ha'
# meta['Paths']['Forest Inventory Disturbances']=r'C:\Data\ForestInventory\Disturbances\20220422\Disturbances.gdb'
# #meta['Paths']['Figures']=r'C:\Users\rhember\OneDrive - Government of BC\Figures\Reforestation'

# #%% Plotting parameters

# meta['Graphics']={}
# meta['Graphics']['figwidth']=16

# #meta['Graphics']['sidespace']=0.25
# meta['Graphics']['sidespace']=0

# meta['Graphics']['ax1 pos']=[0,0,1-meta['Graphics']['sidespace']-0.01,1]
# meta['Graphics']['ax1 vis']='off'
# meta['Graphics']['ax1 gridvis']=False
# meta['Graphics']['ax2 pos']=[1-meta['Graphics']['sidespace']+meta['Graphics']['ax1 pos'][0]+0.01,0.6,0.03,0.35]
# meta['Graphics']['ax2 pos long']=[1-meta['Graphics']['sidespace']+meta['Graphics']['ax1 pos'][0]+0.01,0.1,0.03,0.8]

# gp=gu.SetGraphics('Manuscript')

# meta['LUT']=u1ha.Import_BC1ha_LUTs()

# #%% Import data

# #gdf=u1ha.Import_GDBs_ProvinceWide()

# zRef=gis.OpenGeoTiff(r'C:\Data\BC1ha\Admin\BC_Land_Mask.tif')

# zTSA=gis.OpenGeoTiff(r'C:\Data\BC1ha\Admin\TSA.tif')
# zLCC1=gis.OpenGeoTiff(r'C:\Data\BC1ha\LandUseLandCover\LandCoverClass1.tif')
# zProt=gis.OpenGeoTiff(r'C:\Data\BC1ha\LandUseLandCover\PROTECTED_LANDS_DESIGNATION.tif')
# #zFCR=gis.OpenGeoTiff(r'C:\Data\BC1ha\Results\RSLT_FOREST_COVER_RESERVE_SVW.tif')
# zRET=gis.OpenGeoTiff(r'C:\Data\BC1ha\Results\SILV_RESERVE_CODE_Consolidated.tif')
# zPL=gis.OpenGeoTiff(r'C:\Data\BC1ha\Results\Planting_FromRESULTS_MaskCount.tif')
# #zST=gis.OpenGeoTiff(r'C:\Data\BC1ha\Results\stocktype.tif')
# zHRT=gis.OpenGeoTiff(r'C:\Data\BC1ha\Results\Harvest_Regen_Type.tif')
# zFYL=gis.OpenGeoTiff(r'C:\Data\BC1ha\Disturbances\PROT_HISTORICAL_FIRE_POLYS_SP_YearLast.tif')
# zHCC=gis.OpenGeoTiff(r'C:\Data\BC1ha\Disturbances\VEG_CONSOLIDATED_CUT_BLOCKS_SP_MaskAll.tif')
# zHYL=gis.OpenGeoTiff(r'C:\Data\BC1ha\Disturbances\VEG_CONSOLIDATED_CUT_BLOCKS_SP_YearLast.tif')
# zAgeVRI=gis.OpenGeoTiff(r'C:\Data\BC1ha\VRI\age1.tif')

# zVLiveVRI=gis.OpenGeoTiff(r'C:\Data\BC1ha\VRI\v_live.tif')
# zVLiveVRI['Data']=zVLiveVRI['Data'].astype('int16')

# # Stemwood biomass from PFI
# zB=gis.OpenGeoTiff(r'C:\Data\ForestInventory\LiDAR\Boundary TSA\From GQ 20230302\BoundaryTSA_PFI_LnTotalWSV_bc1ha.tif')
# zB['Data']=0.5*0.5*zB['Data']

# dTIPSY=gu.ReadExcel(r'C:\Data\TIPSY_Boundary.xlsx')

# #%% Filter non-study area

# ind=np.where( (zTSA['Data']==meta['LUT']['tsa']['Boundary TSA']) & (zB['Data']>=0)  ) # & (zLCC1['Data']==1)

# d={}
# d['B']=zB['Data'][ind]
# d['LCC1']=zLCC1['Data'][ind]
# d['RET']=zRET['Data'][ind]
# d['PL']=zPL['Data'][ind]
# d['HCC']=zHCC['Data'][ind]
# d['HYL']=zHYL['Data'][ind]
# d['HRT']=zHRT['Data'][ind]
# d['FYL']=zFYL['Data'][ind]
# d['ASH']=2020-d['HYL']
# d['ASF']=2020-d['FYL']
# d['AgeVRI']=zAgeVRI['Data'][ind]
# d['VLiveVRI']=zVLiveVRI['Data'][ind].astype('float')

# #%% Age

# plt.close('all')
# fig,ax=plt.subplots(1,figsize=gu.cm2inch(8,7))
# ax.fill_between([0,6],-5,125,color=[1,0.875,0.775],lw=0,label='Harvest incomplete')
# ax.fill_between([6,2020-1987],-5,125,color=[1,0.925,0.875],lw=0,label='FRPA')
# ax.fill_between([2020-1987,100],-5,125,color=[1,0.975,0.95],lw=0,label='Pre-FRPA')

# ind=np.where( (d['HCC']>=0) )[0]
# bw=5; bin=np.arange(1,200,bw)
# N,mu,med,sig,se=gu.discres(d['ASH'][ind],d['B'][ind],bw,bin)
# ind=np.where(N>1000)[0]
# ax.plot(bin[ind],mu[ind],'ko',ms=3,label='Predictive Forest Inventory\nchronosequence')
# ax.errorbar(bin[ind],mu[ind],yerr=[3*se[ind],3*se[ind]],color='k',ls='',capsize=1.5,lw=0.25)

# ax.plot(dTIPSY['Age'],dTIPSY['Tot C'],'k-',color=[0.27,0.49,0.79],lw=2,label='TIPSY (SI = 20m)')
# ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=2)
# ax.set(xlabel='Age, years',ylabel='Stemwood carbon (MgC/ha)',xlim=[0,260],ylim=[0,110])

# leg=ax.legend(loc='upper right',frameon=True,fontsize=7)
# frame=leg.get_frame()
# frame.set_facecolor([0.9,0.9,0.9])
# frame.set_linewidth(0)

# ind=np.where( (d['FYL']>0) )[0]
# N,mu,med,sig,se=gu.discres(2020-d['FYL'][ind],d['B'][ind],bw,bin)
# ind=np.where(N>100)[0]
# ax.plot(bin[ind],mu[ind],'go',mec=[0.6,1,0],mfc=[0.6,1,0],ms=3,label='Predictive Forest Inventory\nchronosequence')

# # ind=np.where( (d['HCC']==1) & (d['RET']!=6) )[0]
# # bw=1; bin=np.arange(1,200,bw)
# # N,mu,med,sig,se=gu.discres(d['ASH'][ind],d['B'][ind],bw,bin)
# # ind=np.where(N>1000)[0]
# # ax.plot(bin[ind],mu[ind],'gs',ms=3,label='Predictive Forest Inventory\nchronosequence')

# # N,mu,med,sig,se=gu.discres(d['ASH'][ind],d['VLiveVRI'][ind],bw,bin)
# # ind=np.where(N>1000)[0]
# # ax.plot(bin[ind],mu[ind],'rs',ms=3,label='Predictive Forest Inventory\nchronosequence')

# #%% Age VRI

# plt.close('all')
# fig,ax=plt.subplots(1,figsize=gu.cm2inch(8,7))
# ind=np.where( (d['HCC']>=0) )[0]
# bw=10; bin=np.arange(1,200,bw)
# N,mu,med,sig,se=gu.discres(d['AgeVRI'][ind],d['B'][ind],bw,bin)
# ind=np.where(N>1000)[0]
# ax.plot(bin[ind],mu[ind],'go',mec=[0.6,1,0],mfc=[0.6,1,0],ms=3,label='Predictive Forest Inventory\nchronosequence')
# ax.errorbar(bin[ind],mu[ind],yerr=[3*se[ind],3*se[ind]],color='k',ls='',capsize=1.5,lw=0.25)

# #%% Age

# plt.close('all')
# fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(15,8))
# ind=np.where( (d['HCC']==1) & (d['ASH']>=7) & (d['ASH']<=30) )[0]
# ax[0].hist(d['B'][ind],np.arange(0,220,10))
# ind=np.where( (d['HCC']==0) )[0]
# ax[1].hist(d['B'][ind],np.arange(0,220,10))

# plt.close('all')
# ind=np.where(d['HCC']==1)[0]
# bw=5; bin=np.arange(5,200,bw)
# N,mu,med,sig,se=gu.discres(d['ASH'][ind],d['B'][ind],bw,bin)
# plt.plot(bin,mu,'bo',ms=5)

# N,mu,med,sig,se=gu.discres(d['ASF'][ind],d['B'][ind],bw,bin)
# plt.plot(bin,mu,'rs',ms=5)

# #%% Effect of LCC1

# u=np.unique(d['LCC1'])
# u=u[u!=0]
# mu=np.zeros(u.size)
# err=np.zeros(u.size)
# for i in range(u.size):
#     ind=np.where( (d['LCC1']==u[i]) & (d['HCC']>=1) )[0]
#     mu[i]=np.mean(d['B'][ind])
#     err[i]=2*np.std(d['B'][ind])/np.sqrt(ind.size)
# plt.close('all')
# plt.bar(u,mu)

# #%% Effect of LCC

# u=np.unique(d['RET'])
# u=u[u!=0]
# mu=np.zeros(u.size)
# err=np.zeros(u.size)
# for i in range(u.size):
#     ind=np.where( (d['ASH']>=7) & (d['ASH']<=30) & (d['HCC']==1) & (d['RET']==u[i]) )[0]
#     mu[i]=np.mean(d['B'][ind])
#     err[i]=2*np.std(d['B'][ind])/np.sqrt(ind.size)
# plt.close('all')
# plt.bar(u,(mu-mu[5])/mu[5]*100)

# #%% Regen Type

# u=np.unique(d['HRT'])
# u=u[u!=0]
# mu=np.zeros(u.size)
# err=np.zeros(u.size)
# for i in range(u.size):
#     ind=np.where( (d['ASH']>=7) & (d['ASH']<=30) & (d['HCC']==1) & (d['HRT']==u[i]) )[0]
#     mu[i]=np.mean(d['B'][ind])
#     err[i]=2*np.std(d['B'][ind])/np.sqrt(ind.size)
# plt.close('all')
# plt.bar(u,(mu-mu[0])/mu[0]*100)

# #%% Time series

# print(np.mean(B))

# tv=np.arange(1950,2022,1)
# A_h=np.zeros(tv.size)
# y_h=np.zeros(tv.size)
# y_fcr=np.zeros(tv.size)
# y_pl=np.zeros(tv.size)
# y_npl=np.zeros(tv.size)
# for iT in range(tv.size):
#     zHy=gis.OpenGeoTiff(r'C:\Data\BC1ha\Disturbances\VEG_CONSOLIDATED_CUT_BLOCKS_SP_' + str(tv[iT]) + '.tif')
#     #zHy=gis.OpenGeoTiff(r'C:\Data\BC1ha\Disturbances\Harvest_FromRESULTS_' + str(tv[iT]) + '.tif')
#     hy=zHy['Data'][ind0]
#     ind=np.where( (hy>0) & (fcr==0) )
#     A_h[iT]=ind[0].size
#     y_h[iT]=np.mean(B[ind])
#     ind=np.where( (hy>0) & (fcr==0) & (pl==1) )
#     y_pl[iT]=np.mean(B[ind])
#     ind=np.where( (hy>0) & (fcr==0) & (pl==0) )
#     y_npl[iT]=np.mean(B[ind])
#     ind=np.where( (hy>0) & (fcr>0) )
#     y_fcr[iT]=np.mean(B[ind])
#     print(tv[iT])

# #%%

# plt.close('all')
# plt.plot(tv,y_h,'-bo')
# #plt.plot(tv,y_npl,'-bo')
# #plt.plot(tv,y_pl,'-gs')

# #%% Age



# #%%

# ind0=np.where( (zTSA['Data']==meta['LUT']['tsa']['Boundary TSA']) & (zB['Data']>=0) & (zProt['Data']>0)  ) # & (zLCC1['Data']==1)
# B=zB['Data'][ind0]
# print(np.mean(B))

# ind0=np.where( (zTSA['Data']==meta['LUT']['tsa']['Boundary TSA']) & (zB['Data']>=0) & (zProt['Data']==0)  ) # & (zLCC1['Data']==1)
# B=zB['Data'][ind0]
# print(np.mean(B))

# ind0=np.where( (zTSA['Data']==meta['LUT']['tsa']['Boundary TSA']) & (zB['Data']>=0) & (zFCR['Data']==1)  ) # & (zLCC1['Data']==1)
# B=zB['Data'][ind0]
# print(np.mean(B))