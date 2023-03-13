
"""
PSP - DESCRIPTIVE STATISTICS
"""

#%% Import modules

import numpy as np
import gc as garc
import time
import matplotlib.pyplot as plt
import geopandas as gpd
import pandas as pd
from shapely.geometry import Polygon,Point
import statsmodels.formula.api as smf
import fcgadgets.macgyver.utilities_general as gu
import fcgadgets.macgyver.utilities_gis as gis
import fcexplore.psp.Processing.psp_utilities as utl_gp

#%% Import data

# Graphics
gp=gu.SetGraphics('Manuscript')

# Import BC boundary
gdf_bc_boundary=gpd.read_file(r'C:\Users\rhember\Documents\Data\Basemaps\Basemaps.gdb',layer='NRC_POLITICAL_BOUNDARIES_1M_SP')

# Ground plots
meta={}
meta['Paths']={}
meta['Paths']['DB']=r'C:\Users\rhember\Documents\Data\GroundPlots\PSP-NADB2'
meta['Paths']['Figs']=r'C:\Users\rhember\OneDrive - Government of BC\Figures\Biomass'
meta=utl_gp.ImportParameters(meta)
d=gu.ipickle(meta['Paths']['DB'] + '\\Processed\\L2\\L2_BC.pkl')
sl=d['sobs'].copy()
del d

#%% Plot

points=[]
for k in range(sl['X'].size):
    points.append(Point(sl['X'][k],sl['Y'][k]))

gdf=gpd.GeoDataFrame({'geometry':points,'ID_TSA':np.ones(sl['Y'].size)})
gdf.crs=gdf_bc_boundary.crs

plt.close('all')
fig,ax=plt.subplots(figsize=gu.cm2inch(14,14)); ms=2
gdf_bc_boundary.plot(ax=ax,facecolor=[0.8,0.8,0.8],edgecolor=[0,0,0],label='Political Boundary',linewidth=0.25,alpha=1)
gdf.plot(ax=ax,markersize=ms,facecolor=[0.75,0,0],edgecolor=None,linewidth=0.75,alpha=1)
ax.grid(color='k',linestyle='-',linewidth=0.25)
ax.set(position=[0.01,0.01,0.98,0.98],xticks=[],yticks=[])
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Ground Plots\Map_GroundPlots_All','png',900)

#%% VRI

a=sl.copy()

ind=np.where( (a['Plot Type']==meta['LUT']['Plot Type BC']['VRI']) )[0]
for k in a.keys():
    a[k]=a[k][ind]

points=[]
for k in range(a['X'].size):
    points.append(Point(a['X'][k],a['Y'][k]))

gdf=gpd.GeoDataFrame({'geometry':points,'ID_TSA':np.ones(a['Y'].size)})
gdf.crs=gdf_bc_boundary.crs

plt.close('all')
fig,ax=plt.subplots(figsize=gu.cm2inch(14,14)); ms=2
gdf_bc_boundary.plot(ax=ax,facecolor=[0.8,0.8,0.8],edgecolor=[0,0,0],label='Political Boundary',linewidth=0.25,alpha=1)
gdf.plot(ax=ax,markersize=ms,facecolor=[0.75,0,0],edgecolor=None,linewidth=0.75,alpha=1)
ax.grid(color='k',linestyle='-',linewidth=0.25)
ax.set(position=[0.01,0.01,0.98,0.98],xticks=[],yticks=[])
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Ground Plots\Map_GroundPlots_VRI','png',900)

#%% CMI + NFI

a=sl.copy()

ind=np.where( (a['Plot Type']==meta['LUT']['Plot Type BC']['CMI']) | (a['Plot Type']==meta['LUT']['Plot Type BC']['NFI']) )[0]
for k in a.keys():
    a[k]=a[k][ind]

points=[]
for k in range(a['X'].size):
    points.append(Point(a['X'][k],a['Y'][k]))

gdf=gpd.GeoDataFrame({'geometry':points,'ID_TSA':np.ones(a['Y'].size)})
gdf.crs=gdf_bc_boundary.crs

plt.close('all')
fig,ax=plt.subplots(figsize=gu.cm2inch(14,14)); ms=2
gdf_bc_boundary.plot(ax=ax,facecolor=[0.8,0.8,0.8],edgecolor=[0,0,0],label='Political Boundary',linewidth=0.25,alpha=1)
gdf.plot(ax=ax,markersize=ms,facecolor=[0.75,0,0],edgecolor=None,linewidth=0.75,alpha=1)
ax.grid(color='k',linestyle='-',linewidth=0.25)
ax.set(position=[0.01,0.01,0.98,0.98],xticks=[],yticks=[])
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Ground Plots\Map_GroundPlots_CMI_plus_NFI','png',900)

#%% CMI + NFI (remeasurements only)

a=sl.copy()

ind=np.where( (a['Plot Type']==meta['LUT']['Plot Type BC']['CMI']) & (a['N L t1']>0) | (a['Plot Type']==meta['LUT']['Plot Type BC']['NFI']) & (a['N L t1']>0) )[0]
for k in a.keys():
    a[k]=a[k][ind]

points=[]
for k in range(a['X'].size):
    points.append(Point(a['X'][k],a['Y'][k]))

gdf=gpd.GeoDataFrame({'geometry':points,'ID_TSA':np.ones(a['Y'].size)})
gdf.crs=gdf_bc_boundary.crs

plt.close('all')
fig,ax=plt.subplots(figsize=gu.cm2inch(14,14)); ms=2
gdf_bc_boundary.plot(ax=ax,facecolor=[0.8,0.8,0.8],edgecolor=[0,0,0],label='Political Boundary',linewidth=0.25,alpha=1)
gdf.plot(ax=ax,markersize=ms,facecolor=[0.75,0,0],edgecolor=None,linewidth=0.75,alpha=1)
ax.grid(color='k',linestyle='-',linewidth=0.25)
ax.set(position=[0.01,0.01,0.98,0.98],xticks=[],yticks=[])
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Ground Plots\Map_GroundPlots_CMI_plus_NFI_remes_only','png',900)

#%% YSM

a=sl.copy()

ind=np.where( (a['Plot Type']==meta['LUT']['Plot Type BC']['YSM']) | (a['Plot Type']==meta['LUT']['Plot Type BC']['YSM']) )[0]
for k in a.keys():
    a[k]=a[k][ind]

points=[]
for k in range(a['X'].size):
    points.append(Point(a['X'][k],a['Y'][k]))

gdf=gpd.GeoDataFrame({'geometry':points,'ID_TSA':np.ones(a['Y'].size)})
gdf.crs=gdf_bc_boundary.crs

plt.close('all')
fig,ax=plt.subplots(figsize=gu.cm2inch(14,14)); ms=2
gdf_bc_boundary.plot(ax=ax,facecolor=[0.8,0.8,0.8],edgecolor=[0,0,0],label='Political Boundary',linewidth=0.25,alpha=1)
gdf.plot(ax=ax,markersize=ms,facecolor=[0.75,0,0],edgecolor=None,linewidth=0.75,alpha=1)
ax.grid(color='k',linestyle='-',linewidth=0.25)
ax.set(position=[0.01,0.01,0.98,0.98],xticks=[],yticks=[])
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Ground Plots\Map_GroundPlots_YSM','png',900)

#%% YSM (remes only)

a=sl.copy()

ind=np.where( (a['Plot Type']==meta['LUT']['Plot Type BC']['YSM']) & (a['N L t1']>0) | (a['Plot Type']==meta['LUT']['Plot Type BC']['YSM']) & (a['N L t1']>0) )[0]
for k in a.keys():
    a[k]=a[k][ind]

points=[]
for k in range(a['X'].size):
    points.append(Point(a['X'][k],a['Y'][k]))

gdf=gpd.GeoDataFrame({'geometry':points,'ID_TSA':np.ones(a['Y'].size)})
gdf.crs=gdf_bc_boundary.crs

plt.close('all')
fig,ax=plt.subplots(figsize=gu.cm2inch(14,14)); ms=2
gdf_bc_boundary.plot(ax=ax,facecolor=[0.8,0.8,0.8],edgecolor=[0,0,0],label='Political Boundary',linewidth=0.25,alpha=1)
gdf.plot(ax=ax,markersize=ms,facecolor=[0.75,0,0],edgecolor=None,linewidth=0.75,alpha=1)
ax.grid(color='k',linestyle='-',linewidth=0.25)
ax.set(position=[0.01,0.01,0.98,0.98],xticks=[],yticks=[])
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Ground Plots\Map_GroundPlots_YSM_remes_only','png',900)
