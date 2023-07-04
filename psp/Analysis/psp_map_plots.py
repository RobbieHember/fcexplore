
"""
PSP - DESCRIPTIVE STATISTICS
"""

#%% Import modules

import numpy as np
import gc as garc
import time
import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib.collections import PatchCollection
import geopandas as gpd
import pandas as pd
from shapely.geometry import Polygon,Point
import statsmodels.formula.api as smf
import fcgadgets.macgyver.utilities_general as gu
import fcgadgets.macgyver.utilities_gis as gis
import fcexplore.psp.Processing.psp_utilities as ugp
from fcgadgets.bc1ha import bc1ha_utilities as u1ha

#%% Import data

# Ground plots
meta={}
meta['Paths']={}
meta['Paths']['DB']=r'C:\Users\rhember\Documents\Data\GroundPlots\PSP-NADB2'
meta['Paths']['Figs']=r'C:\Users\rhember\OneDrive - Government of BC\Figures\Biomass'

meta['Graphics']={}
meta['Graphics']['figwidth']=16

#meta['Graphics']['sidespace']=0.25
meta['Graphics']['sidespace']=0

meta['Graphics']['ax1 pos']=[0,0,1-meta['Graphics']['sidespace']-0.01,1]
meta['Graphics']['ax1 vis']='off'
meta['Graphics']['ax1 gridvis']=False
meta['Graphics']['ax2 pos']=[1-meta['Graphics']['sidespace']+meta['Graphics']['ax1 pos'][0]+0.01,0.6,0.03,0.35]
meta['Graphics']['ax2 pos long']=[1-meta['Graphics']['sidespace']+meta['Graphics']['ax1 pos'][0]+0.01,0.1,0.03,0.8]

gp=gu.SetGraphics('Manuscript')

meta=ugp.ImportParameters(meta)
d=gu.ipickle(meta['Paths']['DB'] + '\\Processed\\L2\\L2_BC.pkl')
sl=d['sobs'].copy()
del d

meta['LUT BC1ha']=u1ha.Import_BC1ha_LUTs()

gdf_bc_boundary=gpd.read_file(r'C:\Users\rhember\Documents\Data\Basemaps\Basemaps.gdb',layer='NRC_POLITICAL_BOUNDARIES_1M_SP')

zRef=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\BC_Land_Mask.tif')
zLCC1=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\LandUseLandCover\LandCoverClass1.tif')

ind=np.where(zLCC1['Data']==meta['LUT BC1ha']['lcc1']['Forest'])
zRef['Data'][ind]=2

cm=np.vstack( ((1,1,1,1),(0.95,0.95,0.95,1),(0.8,0.8,0.8,1)) )
cm=matplotlib.colors.ListedColormap(cm)

#%% Plot all

points=[]
for k in range(sl['X'].size):
    points.append(Point(sl['X'][k],sl['Y'][k]))

gdf=gpd.GeoDataFrame({'geometry':points,'ID_TSA':np.ones(sl['Y'].size)})
gdf.crs=gdf_bc_boundary.crs

plt.close('all')
fig,ax=plt.subplots(figsize=gu.cm2inch(14,14)); ms=2
gdf_bc_boundary.plot(ax=ax,facecolor=[0.8,0.8,0.8],edgecolor=[0,0,0],label='Political Boundary',linewidth=0.25,alpha=1)
gdf.plot(ax=ax,markersize=ms,facecolor=[0.75,0,0],edgecolor=None,linewidth=0.75,alpha=1)
#ax.grid(color='k',linestyle='-',linewidth=0.25)
ax.set(position=[0,0,1,1],xticks=[],yticks=[])
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Ground Plots\Map_GroundPlots_All','png',900)

#%% VRI

ms=2
plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(meta['Graphics']['figwidth'],(1-meta['Graphics']['sidespace'])*meta['Graphics']['figwidth']*zRef['yxrat']))
mp=ax.matshow(zRef['Data'],extent=zRef['Extent'],cmap=cm,label='Forest mask')
gdf_bc_boundary.plot(ax=ax,facecolor='None',edgecolor=[0,0,0],label='Political Boundary',linewidth=0.25,alpha=1)

a=sl.copy()
ind=np.where( (a['Plot Type']==meta['LUT']['Plot Type BC']['VRI']) )[0]
for k in a.keys():
    a[k]=a[k][ind]
points=[]
for k in range(a['X'].size):
    points.append(Point(a['X'][k],a['Y'][k]))
gdf=gpd.GeoDataFrame({'geometry':points,'ID_TSA':np.ones(a['Y'].size)})
gdf.crs=gdf_bc_boundary.crs
gdf.plot(ax=ax,markersize=ms,facecolor=[0.75,0,0],edgecolor=None,linewidth=0.75,alpha=1,label='CMI + NFI (All)')

ax.set(position=meta['Graphics']['ax1 pos'],xlim=zRef['xlim'],ylim=zRef['ylim'],aspect='auto')
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')
ax.grid(meta['Graphics']['ax1 gridvis'])
ax.axis(meta['Graphics']['ax1 vis'])
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Ground Plots\Map_GroundPlots_VRI','png',900)

#%% CMI + NFI

ms=2
plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(meta['Graphics']['figwidth'],(1-meta['Graphics']['sidespace'])*meta['Graphics']['figwidth']*zRef['yxrat']))
mp=ax.matshow(zRef['Data'],extent=zRef['Extent'],cmap=cm,label='Forest mask')
gdf_bc_boundary.plot(ax=ax,facecolor='None',edgecolor=[0,0,0],label='Political Boundary',linewidth=0.25,alpha=1)

a=sl.copy()
ind=np.where( (a['Plot Type']==meta['LUT']['Plot Type BC']['CMI']) | (a['Plot Type']==meta['LUT']['Plot Type BC']['NFI']) )[0]
for k in a.keys():
    a[k]=a[k][ind]
points=[]
for k in range(a['X'].size):
    points.append(Point(a['X'][k],a['Y'][k]))
gdf=gpd.GeoDataFrame({'geometry':points,'ID_TSA':np.ones(a['Y'].size)})
gdf.crs=gdf_bc_boundary.crs
gdf.plot(ax=ax,markersize=ms,facecolor=[0.75,0,0],edgecolor=None,linewidth=0.75,alpha=1,label='CMI + NFI (All)')

a=sl.copy()
ind=np.where( (a['Plot Type']==meta['LUT']['Plot Type BC']['CMI']) & (a['N L t1']>0) | (a['Plot Type']==meta['LUT']['Plot Type BC']['NFI']) & (a['N L t1']>0) )[0]
for k in a.keys():
    a[k]=a[k][ind]
points=[]
for k in range(a['X'].size):
    points.append(Point(a['X'][k],a['Y'][k]))
gdf=gpd.GeoDataFrame({'geometry':points,'ID_TSA':np.ones(a['Y'].size)})
gdf.crs=gdf_bc_boundary.crs
gdf.plot(ax=ax,marker='s',markersize=ms+1,facecolor=[0.25,0.5,1],edgecolor=None,linewidth=0.75,alpha=1,label='CMI + NFI (With remeasurements)')

ax.set(position=meta['Graphics']['ax1 pos'],xlim=zRef['xlim'],ylim=zRef['ylim'],aspect='auto')
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')
ax.grid(meta['Graphics']['ax1 gridvis'])
ax.axis(meta['Graphics']['ax1 vis'])
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Ground Plots\Map_GroundPlots_CMI_plus_NFI','png',900)

#%% YSM

ms=2
plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(meta['Graphics']['figwidth'],(1-meta['Graphics']['sidespace'])*meta['Graphics']['figwidth']*zRef['yxrat']))
mp=ax.matshow(zRef['Data'],extent=zRef['Extent'],cmap=cm,label='Forest mask')
gdf_bc_boundary.plot(ax=ax,facecolor='None',edgecolor=[0,0,0],label='Political Boundary',linewidth=0.25,alpha=1)

a=sl.copy()
ind=np.where( (a['Plot Type']==meta['LUT']['Plot Type BC']['YSM']) )[0]
for k in a.keys():
    a[k]=a[k][ind]
points=[]
for k in range(a['X'].size):
    points.append(Point(a['X'][k],a['Y'][k]))
gdf=gpd.GeoDataFrame({'geometry':points,'ID_TSA':np.ones(a['Y'].size)})
gdf.crs=gdf_bc_boundary.crs
gdf.plot(ax=ax,markersize=ms,facecolor=[0.75,0,0],edgecolor=None,linewidth=0.75,alpha=1,label='CMI + NFI (All)')

a=sl.copy()
ind=np.where( (a['Plot Type']==meta['LUT']['Plot Type BC']['YSM']) & (a['N L t1']>0) )[0]
for k in a.keys():
    a[k]=a[k][ind]
points=[]
for k in range(a['X'].size):
    points.append(Point(a['X'][k],a['Y'][k]))
gdf=gpd.GeoDataFrame({'geometry':points,'ID_TSA':np.ones(a['Y'].size)})
gdf.crs=gdf_bc_boundary.crs
gdf.plot(ax=ax,marker='s',markersize=ms+1,facecolor=[0.25,0.5,1],edgecolor=None,linewidth=0.75,alpha=1,label='CMI + NFI (With remeasurements)')

ax.set(position=meta['Graphics']['ax1 pos'],xlim=zRef['xlim'],ylim=zRef['ylim'],aspect='auto')
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')
ax.grid(meta['Graphics']['ax1 gridvis'])
ax.axis(meta['Graphics']['ax1 vis'])
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Ground Plots\Map_GroundPlots_YSM','png',900)

#%% Soils

dSOC=gu.ipickle(r'C:\Users\rhember\Documents\Data\Soils\Shaw et al 2018 Database\SITES.pkl')

ms=2
plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(meta['Graphics']['figwidth'],(1-meta['Graphics']['sidespace'])*meta['Graphics']['figwidth']*zRef['yxrat']))
mp=ax.matshow(zRef['Data'],extent=zRef['Extent'],cmap=cm,label='Forest mask')
gdf_bc_boundary.plot(ax=ax,facecolor='None',edgecolor=[0,0,0],label='Political Boundary',linewidth=0.25,alpha=1)

a=sl.copy()
ind=np.where( (a['Plot Type']==meta['LUT']['Plot Type BC']['NFI']) )[0]
for k in a.keys():
    a[k]=a[k][ind]
points=[]
for k in range(a['X'].size):
    points.append(Point(a['X'][k],a['Y'][k]))
gdf=gpd.GeoDataFrame({'geometry':points,'ID_TSA':np.ones(a['Y'].size)})
gdf.crs=gdf_bc_boundary.crs
gdf.plot(ax=ax,marker='s',markersize=ms+4,facecolor=[0.25,0.5,1],edgecolor=None,linewidth=0.75,alpha=1,label='CMI + NFI (With remeasurements)')

points=[]
for k in range(dSOC['x'].size):
    points.append(Point(dSOC['x'][k],dSOC['y'][k]))

gdf=gpd.GeoDataFrame({'geometry':points,'ID_TSA':np.ones(dSOC['y'].size)})
gdf.crs=gdf_bc_boundary.crs
gdf.plot(ax=ax,markersize=ms,facecolor=[0.75,0,0],edgecolor=None,linewidth=0.75,alpha=1,label='CMI + NFI (All)')

ax.set(position=meta['Graphics']['ax1 pos'],xlim=zRef['xlim'],ylim=zRef['ylim'],aspect='auto')
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')
ax.grid(meta['Graphics']['ax1 gridvis'])
ax.axis(meta['Graphics']['ax1 vis'])
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Soils\Map_Shawetal2018','png',900)

#%% Map for sampling power

ms=2
plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(meta['Graphics']['figwidth'],(1-meta['Graphics']['sidespace'])*meta['Graphics']['figwidth']*zRef['yxrat']))
mp=ax.matshow(zRef['Data'],extent=zRef['Extent'],cmap=cm,label='Forest mask')
gdf_bc_boundary.plot(ax=ax,facecolor='None',edgecolor=[0,0,0],label='Political Boundary',linewidth=0.25,alpha=1)

ind=np.where(zRef['Data']==2)
A=ind[0].size/1e3

ivl=10
ind=np.where(zRef['Data'][0::ivl,0::ivl]==2)
x=zRef['X'][0::ivl,0::ivl][ind]
y=zRef['Y'][0::ivl,0::ivl][ind]
rho1=x.size/A
print(rho1)

points=[]
for k in range(x.size):
    points.append(Point(x[k],y[k]))
gdf=gpd.GeoDataFrame({'geometry':points,'ID':np.ones(x.size)})
gdf.crs=gdf_bc_boundary.crs
gdf.plot(ax=ax,markersize=ms,facecolor=[0.75,0,0],edgecolor=None,linewidth=0.75,alpha=1,label='Proposed sample')

ax.set(position=meta['Graphics']['ax1 pos'],xlim=zRef['xlim'],ylim=zRef['ylim'],aspect='auto')
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')
ax.grid(meta['Graphics']['ax1 gridvis'])
ax.axis(meta['Graphics']['ax1 vis'])
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Biomass\SamplingPowerMap','png',900)