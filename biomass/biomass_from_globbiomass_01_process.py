
"""
Global Biomass (https://doi.pangaea.de/10.1594/PANGAEA.894711?format=html#download)

Santoro, Maurizio; Cartus, Oliver; Mermoz, Stephane; Bouvet, Alexandre; Le Toan, Thuy; Carvalhais, Nuno; Rozendaal, Danae; Herold, Martin; Avitabile, Valerio; Quegan, Shaun; Carreiras, Joao; Rauste, Yrjö; Balzter, Heiko; Schmullius, Christiane C; Seifert, Frank Martin (2018): A detailed portrait of the forest aboveground biomass pool for the year 2010 obtained from multiple remote sensing observations. Geophysical Research Abstracts, 20, EGU2018-18932, https://meetingorganizer.copernicus.org/EGU2018/EGU2018-18932.pdf

Metadata:
https://hs.pangaea.de/Maps/GlobBiomass/README_GLOBBIOMASS_global_20180531_web.pdf

- above ground biomass (AGB, unit: tons/ha i.e., Mg/ha) for the year 2010
Definition: the mass, expressed as oven-dry weight of the woody parts (stem, bark, branches and
twigs) of all living trees excluding stump and roots.

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
from shapely.geometry import Point, Polygon
from shapely import geometry
from rasterio.transform import from_origin
from osgeo import osr
from scipy.interpolate import griddata
import statsmodels.formula.api as smf

import fcgadgets.macgyver.utilities_general as gu
import fcgadgets.macgyver.utilities_gis as gis
from fcgadgets.macgyver import utilities_inventory as invu

#%% Plotting parameters

gp=gu.SetGraphics('Presentation Dark')

# fs=7
# params={'font.sans-serif':'Arial','font.size':fs,'axes.edgecolor':'black','axes.labelsize':fs,'axes.labelcolor':'black','axes.titlesize':fs,'axes.linewidth':0.5,'lines.linewidth':0.5,
#         'text.color':'black','xtick.color':'black','xtick.labelsize':fs,'xtick.major.width':0.5,'xtick.major.size':3,'xtick.direction':'in','ytick.color':'black','ytick.labelsize':fs,
#         'ytick.major.width':0.5,'ytick.major.size':3,'ytick.direction':'in','legend.fontsize':fs,'savefig.dpi':300,'savefig.transparent':True}
# plt.rcParams.update(params)

#%% Import data

flg=0
if flg==1:
    # Do all the processing
    pth0=r'C:\Users\rhember\Documents\Data\Biomass\GlobBiomass\N80W140_agb.tif'
    z0=gis.OpenGeoTiff(pth0)

    # Clip to BC area
    z0=gis.ClipRasterByXYLimits(z0,[-142,-113],[47,61])

    pth1=r'C:\Users\rhember\Documents\Data\Biomass\GlobBiomass\N80W140_agb c.tif'
    z0['Data']=z0['Data'].astype('int16')
    gis.SaveGeoTiff(z0,pth1)

    # Reproject

    # Get CRS for BC
    gdf=gpd.read_file(r'C:\Users\rhember\Documents\Data\Basemaps\Basemaps.gdb',layer='NRC_POLITICAL_BOUNDARIES_1M_SP')

    pth2=r'C:\Users\rhember\Documents\Data\Biomass\GlobBiomass\N80W140_agb pc8.tif'
    pth_ref=r'C:\Users\rhember\Documents\Data\BC1ha\Admin\tsa.tif'
    zB=gis.ReprojectRasterAndClipToRaster(pth1,pth2,pth_ref,gdf.crs)

else:
    # Skip processing
    zB=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\Biomass\GlobBiomass\N80W140_agb pc8.tif')

#%% Convert to carbon

zB['Data']=0.5*zB['Data']

#%% Plot map

plt.close('all')
plt.matshow(zB['Data'],clim=[0,150])

# Forest mask
btm=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\LandUseLandCover\landuse.btm2.tif')
iMask=np.where( (btm['Data']==8) | (btm['Data']==11) | (btm['Data']==12) | (btm['Data']==16) | (btm['Data']==21) )

print(np.sum(zB['Data'][iMask])/1e9)

#%% Import BEC

becz=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\becz.tif')
lutBGC=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\becz_lut.xlsx')

#%% Plot by BGC zone

u=np.unique(becz['Data'][becz['Data']!=255])

lab=np.array(['' for _ in range(u.size)],dtype=object)
B={}
B['mu']=np.zeros(u.size)
B['sd']=np.zeros(u.size)
B['se']=np.zeros(u.size)
for i in range(u.size):
    ind=np.where( (becz['Data']==u[i]) & (zB['Data']>0) & (zB['Data']<2000) )

    B['mu'][i]=np.nanmean(zB['Data'][ind])
    B['sd'][i]=np.nanstd(zB['Data'][ind])
    B['se'][i]=np.nanstd(zB['Data'][ind])/np.sqrt(ind[0].size)

    ind=np.where(lutBGC['VALUE']==u[i])[0]
    if ind.size>0:
        lab[i]=lutBGC['ZONE'][ind][0]

# Put in order
ord=np.argsort(B['mu'])
for k in B:
    B[k]=np.flip(B[k][ord])
lab=np.flip(lab[ord])
uo=u[ord]


#%% Stratified by BGC zone

plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,6))
ax.bar(np.arange(u.size),B['mu'],facecolor=[0.75,0.9,0])
ax.errorbar(np.arange(u.size),B['mu'],yerr=B['sd'],color='k',fmt='none',capsize=2)
ax.set(position=[0.08,0.12,0.9,0.86],xlim=[-0.5,u.size-0.5],ylim=[0,160],xticks=np.arange(u.size),
       xticklabels=lab,ylabel='Aboveground biomass (MgC ha$^{-1}$ yr$^{-1}$)')
#plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Biomass\GlobBiomass_ByGBCZone','png',900)


#%%

# plt.close('all')
# fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,6))
# ax.bar(np.arange(u.size)-0.225,B['mu'],0.45,facecolor=[0.29,0.47,0.79],label='Global Biomass')
# ax.errorbar(np.arange(u.size)-0.225,B['mu'],yerr=2*B['se'],color='k',fmt='none',capsize=2)
# ax.bar(np.arange(u.size)+0.225,gps['mu'],0.45,facecolor=[0.75,0.9,0],label='BC Plot Sample')
# ax.errorbar(np.arange(u.size)+0.225,gps['mu'],yerr=2*gps['se'],color='k',fmt='none',capsize=2)
# for i in range(u.size):
#     ax.text(i+0.225,5,str(N[i].astype(int)),color='k',ha='center',fontsize=7)
# ax.set(position=[0.08,0.12,0.9,0.86],xlim=[-0.5,u.size-0.5],ylim=[0,160],xticks=np.arange(u.size),
#        xticklabels=lab,ylabel='Aboveground biomass (MgC ha$^{-1}$ yr$^{-1}$)')
# plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
# ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')
# gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Biomass\GlobBiomass_CompWithGPs','png',900)


