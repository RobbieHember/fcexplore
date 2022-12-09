
"""
Observation-based estimates of total ecosystem carbon (TEC) across BC forests

Data sources:
- Biomass (Ground Plot Database - Forest Analysis and Inventory Branch)
- Soil (Canadian Upland Forest Soils in British Columbia, Shaw et al. 2018)

"""

#%% Import modules

import numpy as np
import gc as garc
from osgeo import osr
import geopandas as gpd
import matplotlib.pyplot as plt
import geopandas as gpd
import pandas as pd
import copy
import statsmodels.formula.api as smf
import cv2
from shapely.geometry import Point, Polygon
from shapely import geometry

import fcgadgets.macgyver.utilities_general as gu
import fcgadgets.macgyver.utilities_gis as gis
from fcexplore.psp.Processing import psp_utilities as utl_gp

#%% Set figure properties

#gp=gu.SetGraphics('Presentation Dark')
gp=gu.SetGraphics('Manuscript')

#%% Import Canadian Upland forest database

soc=gu.ipickle(r'C:\Users\rhember\Documents\Data\Soils\Shaw et al 2018 Database\SITES.pkl')

#%% Import ground plot data

metaGP={}
metaGP['Paths']={}
metaGP['Paths']['DB']=r'C:\Users\rhember\Documents\Data\GroundPlots\PSP-NADB2'
metaGP=utl_gp.ImportParameters(metaGP)
d=gu.ipickle(metaGP['Paths']['DB'] + '\\Processed\\L2\\L2_BC.pkl')
gplt=d['sobs']
del d

# Plot type indicator
gplt['pt_ind']=np.zeros(gplt['ID Plot'].size)
ind=np.where( (gplt['Plot Type']==metaGP['LUT']['Plot Type BC']['CMI']) |
             (gplt['Plot Type']==metaGP['LUT']['Plot Type BC']['NFI']) |
             (gplt['Plot Type']==metaGP['LUT']['Plot Type BC']['VRI']) )[0]
gplt['pt_ind'][ind]=1

#%% Import Raster grids

becz=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\becz.tif')
lutBGC=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\becz_lut.xlsx')

lc2=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\lc2.tif')

#%% Plot mean by BGC zone

u=np.unique(soc['becz'])
lab=np.array(['' for _ in range(u.size)],dtype=object)

ds={}
ds['N SOC tot']=np.zeros(u.size)
ds['mu SOC tot']=np.zeros(u.size)
ds['se SOC tot']=np.zeros(u.size)
ds['mu SOC min']=np.zeros(u.size)
ds['mu SOC org']=np.zeros(u.size)
ds['mu Ctot D t0']=np.zeros(u.size)
ds['mu Cbk L t0']=np.zeros(u.size)
ds['mu Cbr L t0']=np.zeros(u.size)
ds['mu Cf L t0']=np.zeros(u.size)
ds['mu Cr L t0']=np.zeros(u.size)
ds['mu Csw L t0']=np.zeros(u.size)
ds['Sum Area']=np.zeros(u.size)
ds['Sum TEC']=np.zeros(u.size)
for i in range(u.size):
    ind=np.where( (soc['becz']==u[i]) & (soc['TOT_C_THA']>0) )[0]
    ds['N SOC tot'][i]=ind.size
    ds['mu SOC tot'][i]=np.nanmean(soc['TOT_C_THA'][ind])
    ds['se SOC tot'][i]=np.nanstd(soc['TOT_C_THA'][ind])/np.sqrt(ind.size)
    ds['mu SOC min'][i]=np.nanmean(soc['MIN_C_THA'][ind])
    ds['mu SOC org'][i]=np.nanmean(soc['ORG_C_THA'][ind])
    ind=np.where(lutBGC['VALUE']==u[i])[0]
    if ind.size>0:
        lab[i]=lutBGC['ZONE'][ind][0]

    # Dead wood and biomass
    ind=np.where( (gplt['Ecozone BC L1']==metaGP['LUT']['Ecozone BC L1'][lab[i]]) & (gplt['pt_ind']==1) & (gplt['Ctot D t0']>=0) & (gplt['Ctot D t0']<2000) & (gplt['Ctot L t0']>=0) & (gplt['Ctot L t0']<2000) )[0]
    ds['mu Ctot D t0'][i]=np.nanmean(gplt['Ctot D t0'][ind])
    ds['mu Cbk L t0'][i]=np.nanmean(gplt['Cbk L t0'][ind])
    ds['mu Cbr L t0'][i]=np.nanmean(gplt['Cbr L t0'][ind])
    ds['mu Cf L t0'][i]=np.nanmean(gplt['Cf L t0'][ind])
    ds['mu Cr L t0'][i]=np.nanmean(gplt['Cr L t0'][ind])
    ds['mu Csw L t0'][i]=np.nanmean(gplt['Csw L t0'][ind])

# Total carbon

ds['mu C All']=ds['mu SOC tot']+ds['mu Ctot D t0']+ds['mu Cbk L t0']+ds['mu Cbr L t0']+ds['mu Cf L t0']+ds['mu Cr L t0']+ds['mu Csw L t0']

for i in range(lab.size):
    ind1=np.where(lutBGC['ZONE']==lab[i])[0]
    ind2=np.where( (becz['Data']==lutBGC['VALUE'][ind1[0]]) & (lc2['Data']==4) )
    ds['Sum TEC'][i]=ds['mu C All'][i]*ind2[0].size/1e9
    ds['Sum Area'][i]=ind2[0].size/1e6

tmp=ds.copy()
for k in tmp.keys():
    if k[0:2]=='mu':
        ds['Sum ' + k[3:]]=(ds['Sum Area']*1e6)*ds[k]/1e9

#%% Plot DOM

ord=np.argsort(ds['mu SOC tot']+ds['mu Ctot D t0'])
for k in ds:
    ds[k]=np.flip(ds[k][ord])
lab=np.flip(lab[ord])

plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,6))
ax.bar(np.arange(u.size),ds['mu SOC min'],facecolor=[0.45,0.3,0.3],label='Mineral soil horizon (to 100 cm depth)')
ax.bar(np.arange(u.size),ds['mu SOC org'],facecolor=[0.15,0.05,0.05],bottom=ds['mu SOC min'],label='Organic soil horizon')
ax.bar(np.arange(u.size),ds['mu Ctot D t0'],facecolor=[0.85,0.75,0.65],bottom=ds['mu SOC min']+ds['mu SOC org'],label='Standing + fallen dead wood')
#ax.errorbar(np.arange(u.size),ds['mu SOC tot'],yerr=ds['se'],color=gp['cla'],fmt='none',capsize=2)
for i in range(u.size):
    ax.text(i,10,str(ds['N SOC tot'][i].astype(int)),color=gp['cla'],ha='center',fontsize=8)
ax.set(position=[0.08,0.12,0.9,0.86],xticks=np.arange(u.size),
       xticklabels=lab,ylabel='Carbon stock (MgC ha$^{-1}$ yr$^{-1}$)',xlim=[-0.5,u.size-0.5],ylim=[0,375])
plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Soils\Mean Dead Carbon' ,'png',900)

#%% Plot Total Ecosystem Carbon

ord=np.argsort(ds['mu C All'])
for k in ds:
    ds[k]=np.flip(ds[k][ord])
lab=np.flip(lab[ord])

cl=np.array([ [0.45,0.3,0.3],[0.15,0.05,0.05],[0.85,0.75,0.65],[0.8,0,0],[0.6,1,0],[1,0.5,0],[0.45,1,1],[0.1,0.75,0] ])

plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,6))
ax.bar(np.arange(u.size),ds['mu SOC min'],facecolor=cl[0,:],label='Mineral soil horizon (100 cm depth)')
ax.bar(np.arange(u.size),ds['mu SOC org'],bottom=ds['mu SOC min'],facecolor=cl[1,:],label='Organic soil horizon')
ax.bar(np.arange(u.size),ds['mu Ctot D t0'],bottom=ds['mu SOC min']+ds['mu SOC org'],facecolor=cl[2,:],label='Standing + fallen dead wood')
y=ds['mu SOC min']+ds['mu SOC org']+ds['mu Ctot D t0']
ax.bar(np.arange(u.size),ds['mu Cr L t0'],bottom=y,facecolor=cl[3,:],label='Roots')
ax.bar(np.arange(u.size),ds['mu Cbk L t0'],bottom=y+ds['mu Cr L t0'],facecolor=cl[4,:],label='Bark')
ax.bar(np.arange(u.size),ds['mu Cbr L t0'],bottom=y+ds['mu Cr L t0']+ds['mu Cbk L t0'],facecolor=cl[5,:],label='Branches')
ax.bar(np.arange(u.size),ds['mu Cf L t0'],bottom=y+ds['mu Cr L t0']+ds['mu Cbk L t0']+ds['mu Cbr L t0'],facecolor=cl[6,:],label='Foliage')
ax.bar(np.arange(u.size),ds['mu Csw L t0'],bottom=y+ds['mu Cr L t0']+ds['mu Cbk L t0']+ds['mu Cbr L t0']+ds['mu Cf L t0'],facecolor=cl[7,:],label='Stemwood')

#ax.errorbar(np.arange(u.size),ds['mu SOC tot'],yerr=ds['se'],color=gp['cla'],fmt='none',capsize=2)
#for i in range(u.size):
#    ax.text(i,10,str(ds['N'][i].astype(int)),color=gp['cla'],ha='center',fontsize=8)
ax.set(position=[0.08,0.1,0.9,0.88],xticks=np.arange(u.size),
       xticklabels=lab,ylabel='Average stock (tC per hectare)',xlim=[-0.5,u.size-0.5],ylim=[0,550])
plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\TotalEcosystemCarbonBC\TEC Per-hectare Mean' ,'png',900)

#%% Plot Total Ecosystem Carbon

ord=np.argsort(ds['Sum C All'])
for k in ds:
    ds[k]=np.flip(ds[k][ord])
lab=np.flip(lab[ord])

cl=np.array([ [0.45,0.3,0.3],[0.15,0.05,0.05],[0.85,0.75,0.65],[0.8,0,0],[0.6,1,0],[1,0.5,0],[0.45,1,1],[0.1,0.75,0] ])

plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,6))
ax.bar(np.arange(u.size),ds['Sum SOC min'],facecolor=cl[0,:],label='Mineral soil horizon (100 cm depth)')
ax.bar(np.arange(u.size),ds['Sum SOC org'],bottom=ds['Sum SOC min'],facecolor=cl[1,:],label='Organic soil horizon')
ax.bar(np.arange(u.size),ds['Sum Ctot D t0'],bottom=ds['Sum SOC min']+ds['Sum SOC org'],facecolor=cl[2,:],label='Standing + fallen dead wood')
y=ds['Sum SOC min']+ds['Sum SOC org']+ds['Sum Ctot D t0']
ax.bar(np.arange(u.size),ds['Sum Cr L t0'],bottom=y,facecolor=cl[3,:],label='Roots')
ax.bar(np.arange(u.size),ds['Sum Cbk L t0'],bottom=y+ds['Sum Cr L t0'],facecolor=cl[4,:],label='Bark')
ax.bar(np.arange(u.size),ds['Sum Cbr L t0'],bottom=y+ds['Sum Cr L t0']+ds['Sum Cbk L t0'],facecolor=cl[5,:],label='Branches')
ax.bar(np.arange(u.size),ds['Sum Cf L t0'],bottom=y+ds['Sum Cr L t0']+ds['Sum Cbk L t0']+ds['Sum Cbr L t0'],facecolor=cl[6,:],label='Foliage')
ax.bar(np.arange(u.size),ds['Sum Csw L t0'],bottom=y+ds['Sum Cr L t0']+ds['Sum Cbk L t0']+ds['Sum Cbr L t0']+ds['Sum Cf L t0'],facecolor=cl[7,:],label='Stemwood')

ax.set(position=[0.08,0.1,0.9,0.88],xticks=np.arange(u.size),
       xticklabels=lab,ylabel='Total stock (Billion tonnes C)',xlim=[-0.5,u.size-0.5],ylim=[0,4])
plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\TotalEcosystemCarbonBC\TEC Sum' ,'png',900)

#%% Map

plt.close('all')
plt.plot(gplt['X'],gplt['Y'],'ro')
plt.plot(soc['x'],soc['y'],'bs')


# Get X and Y coords
#srs=gis.ImportSRSs()
#sl['x'],sl['y']=srs['Proj']['BC1ha'](sl['Lon'],sl['Lat'])


#%% Import vector files

gdf_ogsr=gpd.read_file(r'C:\Users\rhember\Documents\Data\ForestInventory\LandUse\20220422\LandUse.gdb',layer='OGSR_TAP_PRIORITY_DEF_AREA_SP')
gdf_bm=gpd.read_file(r'C:\Users\rhember\Documents\Data\Basemaps\Basemaps.gdb',layer='NRC_POLITICAL_BOUNDARIES_1M_SP')

#%% Create geodatabase

points=[]
for i in range(gplt['X'].size):
    points.append( Point(gplt['X'][i],gplt['Y'][i]) )
d={'geometry':points}
for k in gplt.keys():
    d[k]=gplt[k]
gdf_gplt=gpd.GeoDataFrame(d)
gdf_gplt.crs=gdf_bm.crs
#gdf_gp.to_file(r'C:\Users\rhember\Documents\Data\GroundPlots\PSP-NADB\ground_plots.geojson',driver='GeoJSON')
gdf_gplt_og=gpd.overlay(gdf_gplt,gdf_ogsr,how='intersection')

points=[]
for i in range(soc['x'].size):
    points.append( Point(soc['x'][i],soc['y'][i]) )
d={'geometry':points}
for k in soc.keys():
    d[k]=soc[k]
gdf_soc=gpd.GeoDataFrame(d)
gdf_soc.crs=gdf_bm.crs
#gdf_gp.to_file(r'C:\Users\rhember\Documents\Data\GroundPlots\PSP-NADB\ground_plots.geojson',driver='GeoJSON')
gdf_soc_og=gpd.overlay(gdf_soc,gdf_ogsr,how='intersection')

#%% Stats inside deferral areas

ds={}
ds['N SOC tot']=np.zeros(u.size)
ds['mu SOC tot']=np.zeros(u.size)
ds['se SOC tot']=np.zeros(u.size)
ds['mu SOC min']=np.zeros(u.size)
ds['mu SOC org']=np.zeros(u.size)
ds['mu Ctot D t0']=np.zeros(u.size)
ds['mu Cbk L t0']=np.zeros(u.size)
ds['mu Cbr L t0']=np.zeros(u.size)
ds['mu Cf L t0']=np.zeros(u.size)
ds['mu Cr L t0']=np.zeros(u.size)
ds['mu Csw L t0']=np.zeros(u.size)
ds['Sum Area']=np.zeros(u.size)
ds['Sum TEC']=np.zeros(u.size)
for i in range(u.size):
    ind=np.where( (gdf_soc_og['becz']==u[i]) & (gdf_soc_og['TOT_C_THA']>0) )[0]
    ds['N SOC tot'][i]=ind.size
    ds['mu SOC tot'][i]=np.nanmean(gdf_soc_og['TOT_C_THA'][ind])
    ds['se SOC tot'][i]=np.nanstd(gdf_soc_og['TOT_C_THA'][ind])/np.sqrt(ind.size)
    ds['mu SOC min'][i]=np.nanmean(gdf_soc_og['MIN_C_THA'][ind])
    ds['mu SOC org'][i]=np.nanmean(gdf_soc_og['ORG_C_THA'][ind])

    # Dead wood and biomass
    ind=np.where( (gdf_gplt_og['Ecozone BC L1']==metaGP['LUT']['Ecozone BC L1'][lab[i]]) & (gdf_gplt_og['pt_ind']==1) & (gdf_gplt_og['Ctot D t0']>=0) & (gdf_gplt_og['Ctot D t0']<2000) & (gdf_gplt_og['Ctot L t0']>=0) & (gdf_gplt_og['Ctot L t0']<2000) )[0]
    ds['mu Ctot D t0'][i]=np.nanmean(gdf_gplt_og['Ctot D t0'][ind])
    ds['mu Cbk L t0'][i]=np.nanmean(gdf_gplt_og['Cbk L t0'][ind])
    ds['mu Cbr L t0'][i]=np.nanmean(gdf_gplt_og['Cbr L t0'][ind])
    ds['mu Cf L t0'][i]=np.nanmean(gdf_gplt_og['Cf L t0'][ind])
    ds['mu Cr L t0'][i]=np.nanmean(gdf_gplt_og['Cr L t0'][ind])
    ds['mu Csw L t0'][i]=np.nanmean(gdf_gplt_og['Csw L t0'][ind])

# Total carbon
ds['mu C All']=ds['mu SOC tot']+ds['mu Ctot D t0']+ds['mu Cbk L t0']+ds['mu Cbr L t0']+ds['mu Cf L t0']+ds['mu Cr L t0']+ds['mu Csw L t0']

for k in ds.keys():
    ds[k]=np.nan_to_num(ds[k])

#%% Plot Total Ecosystem Carbon

ord=np.argsort(ds['mu C All'])
for k in ds:
    ds[k]=np.flip(ds[k][ord])
lab=np.flip(lab[ord])

cl=np.array([ [0.45,0.3,0.3],[0.15,0.05,0.05],[0.85,0.75,0.65],[0.8,0,0],[0.6,1,0],[1,0.5,0],[0.45,1,1],[0.1,0.75,0] ])

plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,6))
ax.bar(np.arange(u.size),ds['mu SOC min'],facecolor=cl[0,:],label='Mineral soil horizon (100 cm depth)')
ax.bar(np.arange(u.size),ds['mu SOC org'],bottom=ds['mu SOC min'],facecolor=cl[1,:],label='Organic soil horizon')
ax.bar(np.arange(u.size),ds['mu Ctot D t0'],bottom=ds['mu SOC min']+ds['mu SOC org'],facecolor=cl[2,:],label='Standing + fallen dead wood')
y=ds['mu SOC min']+ds['mu SOC org']+ds['mu Ctot D t0']
ax.bar(np.arange(u.size),ds['mu Cr L t0'],bottom=y,facecolor=cl[3,:],label='Roots')
ax.bar(np.arange(u.size),ds['mu Cbk L t0'],bottom=y+ds['mu Cr L t0'],facecolor=cl[4,:],label='Bark')
ax.bar(np.arange(u.size),ds['mu Cbr L t0'],bottom=y+ds['mu Cr L t0']+ds['mu Cbk L t0'],facecolor=cl[5,:],label='Branches')
ax.bar(np.arange(u.size),ds['mu Cf L t0'],bottom=y+ds['mu Cr L t0']+ds['mu Cbk L t0']+ds['mu Cbr L t0'],facecolor=cl[6,:],label='Foliage')
ax.bar(np.arange(u.size),ds['mu Csw L t0'],bottom=y+ds['mu Cr L t0']+ds['mu Cbk L t0']+ds['mu Cbr L t0']+ds['mu Cf L t0'],facecolor=cl[7,:],label='Stemwood')

#ax.errorbar(np.arange(u.size),ds['mu SOC tot'],yerr=ds['se'],color=gp['cla'],fmt='none',capsize=2)
#for i in range(u.size):
#    ax.text(i,10,str(ds['N'][i].astype(int)),color=gp['cla'],ha='center',fontsize=8)
ax.set(position=[0.08,0.1,0.9,0.88],xticks=np.arange(u.size),
       xticklabels=lab,ylabel='Average stock (tC per hectare)',xlim=[-0.5,u.size-0.5],ylim=[0,600])
plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\TotalEcosystemCarbonBC\TEC Per-hectare Mean (OGSR Priority Areas)' ,'png',900)