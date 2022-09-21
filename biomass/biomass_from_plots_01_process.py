
"""


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
from shapely.geometry import Point, Polygon,box,shape
from shapely import geometry
from rasterio.transform import from_origin
from osgeo import osr
from scipy.interpolate import griddata
import statsmodels.formula.api as smf
from scipy.io import loadmat
import pyproj

import fcgadgets.macgyver.utilities_general as gu
import fcgadgets.macgyver.utilities_gis as gis
from fcgadgets.macgyver import utilities_inventory as invu

#%% Set figure properties

fs=7
params={'font.sans-serif':'Arial','font.size':fs,'axes.edgecolor':'black','axes.labelsize':fs,'axes.labelcolor':'black','axes.titlesize':fs,'axes.linewidth':0.5,'lines.linewidth':0.5,
        'text.color':'black','xtick.color':'black','xtick.labelsize':fs,'xtick.major.width':0.5,'xtick.major.size':3,'xtick.direction':'in','ytick.color':'black','ytick.labelsize':fs,
        'ytick.major.width':0.5,'ytick.major.size':3,'ytick.direction':'in','legend.fontsize':fs,'savefig.dpi':300,'savefig.transparent':True}
plt.rcParams.update(params)

#%% Import data

mat=loadmat(r'C:\Users\rhember\Documents\Data\GroundPlots\PSP-NADB\FI_PSP_L3_SL_R19b.mat')
sl={n: mat['q'][n][0, 0] for n in  mat['q'].dtype.names}

for k in sl.keys():
    sl[k]=sl[k].flatten()

# Get X and Y coords
srs=gis.ImportSRSs()
sl['x'],sl['y']=srs['Proj']['BC1ha'](sl['Lon'],sl['Lat'])

#%% Filter

ind1=np.where( (sl['ID_DB']==1) & (sl['Lat']>0) & (sl['Lon']!=0) )[0]

for k in sl.keys():
    sl[k]=sl[k][ind1]

#%% Save geodatabase

points=[]
for i in range(sl['x'].size):
    points.append( Point(sl['x'][i],sl['y'][i]) )
d={'geometry':points}
for k in sl.keys():
    d[k]=sl[k]
gdf=gpd.GeoDataFrame(d)
gdf_bm=gpd.read_file(r'C:\Users\rhember\Documents\Data\Basemaps\Basemaps.gdb',layer='NRC_POLITICAL_BOUNDARIES_1M_SP')
gdf.crs=gdf_bm.crs
#gdf.to_file(r'C:\Users\rhember\Documents\Data\GroundPlots\PSP-NADB\ground_plots.geojson',driver='GeoJSON')

ogsr=gpd.read_file(r'C:\Users\rhember\Documents\Data\ForestInventory\LandUse\20220422\LandUse.gdb',layer='OGSR_TAP_PRIORITY_DEF_AREA_SP')

slo=gpd.overlay(gdf,ogsr,how='intersection')

#%% Stats inside deferral areas

ind=np.where( (slo['EcoZone_BC_L1']==9) & (slo['PlotType']==1) & (slo['Age_t0']>200) & (slo['Cag_L_t0']>0) & (slo['Cag_L_t0']<2000) |
             (slo['EcoZone_BC_L1']==9) & (slo['PlotType']==9) & (slo['Age_t0']>200) & (slo['Cag_L_t0']>0) & (slo['Cag_L_t0']<2000))[0]
print(ind.size)
print(np.mean(np.real(slo['Cag_L_t0'][ind])))
print(np.std(slo['Cag_L_t0'][ind])/np.sqrt(ind.size))

#%% Stats all

ind=np.where( (sl['EcoZone_BC_L1']==9) & (sl['PlotType']==1) & (sl['Age_t0']>200) & (sl['Cag_L_t0']>0) & (sl['Cag_L_t0']<2000) |
             (sl['EcoZone_BC_L1']==9) & (sl['PlotType']==9) & (sl['Age_t0']>200) & (sl['Cag_L_t0']>0) & (sl['Cag_L_t0']<2000))[0]
print(ind.size)
print(np.mean(np.real(sl['Cag_L_t0'][ind])))
print(np.std(sl['Cag_L_t0'][ind])/np.sqrt(ind.size))

#%%

fig,ax=plt.subplots(1,figsize=gu.cm2inch(10,10))
plt.plot(sl['Age_t0'][ind],sl['Ctot_L_t0'][ind],'.')

ind=np.where( (sl['EcoZone_BC_L1']==9) & (sl['Age_t0']>125) & (sl['Csw_L_t0']>0) & (sl['Csw_L_t0']<1500) )[0]
np.mean(sl['Cag_L_t0'][ind])


#%% Get VRI polygons that intersect plots

pth=r'C:\Users\rhember\Documents\Data\ForestInventory\VRI\20220404\VRI.gdb'
List=[None]*int(1e5)
cnt=0
with fiona.open(pth,layer='VEG_COMP_LYR_R1_POLY') as source:
    for feat in source:
        if feat['geometry']==None:
            continue
        shp=shape(feat['geometry'])
        #break
        ind=np.where(gdf.within(shp)==True)[0]
        if ind.size>0:
            List[cnt]=feat
            cnt=cnt+1
            #print('working')
List=List[0:cnt-1]
gu.opickle(r'C:\Users\rhember\Documents\Data\GroundPlots\PSP-NADB\ground_plots_vri.pkl',List)

gdf_vri=gpd.GeoDataFrame.from_features(List,crs=gdf_bm.crs)

gdf_vri.to_file(r'C:\Users\rhember\Documents\Data\GroundPlots\PSP-NADB\ground_plots_vri.geojson',driver='GeoJSON')


#%% Spatial join of VRI polygons and PSPs

#gdf2=gpd.sjoin(gdf,gdf_vri,how='left')
#gdf2=gdf.sjoin(gdf_vri,how='inner')

gdf2=gdf_vri.sjoin(gdf,how='inner')

# Plot to make sure it worked
plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(12,12))
gdf2.plot(ax=ax,edgecolor='r',facecolor='y',linewidth=0.5,label='Road',alpha=1)
gdf.plot(ax=ax,edgecolor='g',linewidth=0.5,markersize=14,label='Road',alpha=1)

#a=gpd.overlay(gdf_vri,gdf,how='intersection')

#%% Comparison of plot and VRI polygon

df2=pd.DataFrame(gdf2.drop(columns='geometry'))
d2={}
for k in df2.columns:
    d2[k]=df2[k].values

d2['AGB_VRI']=d2['WHOLE_STEM_BIOMASS_PER_HA']+d2['BRANCH_BIOMASS_PER_HA']+d2['FOLIAGE_BIOMASS_PER_HA']+d2['BARK_BIOMASS_PER_HA']

#%%

ind=np.where( (d2['EcoZone_BC_L1']>0) & (d2['PlotType']==1) & (d2['Age_t0']>140) & (d2['Cag_L_t0']>0) & (d2['Cag_L_t0']<2000) & (d2['AGB_VRI']>0) |
             (d2['EcoZone_BC_L1']>0) & (d2['PlotType']==9) & (d2['Age_t0']>140) & (d2['Cag_L_t0']>0) & (d2['Cag_L_t0']<2000) & (d2['AGB_VRI']>0) )[0]

x=np.real(d2['Cag_L_t0'][ind])
y=0.5*d2['AGB_VRI'][ind]
rs,txt=gu.GetRegStats(x,y);

plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(10,10))
ax.plot([0,1000],[0,1000],'-k',lw=3,color=[0.8,0.8,0.8])
ax.plot(x,y,'ko',markersize=3,markerfacecolor='w',markeredgecolor='b',lw=0.5,label='Provincial Ground Plots (VRI + CMI)')
ax.plot(rs['xhat'],rs['yhat'],'c-',lw=1,label='Best fit (Provincial Ground Plots)')
ax.text(125,525,txt,fontsize=10,color='c')
ax.text(770,770,'1:1',fontsize=8,ha='center')
ax.set(position=[0.12,0.12,0.84,0.84],xlim=[0,1000],ylim=[0,1000],xlabel='Ground Sample (MgC/ha)',ylabel='VRI polygon (MgC/ha)')
z=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\Dellasala et al 2022 Fig7 dig.xlsx')
ax.plot(z['Sample Plot'],z['VRI polygon'],'gs',markersize=3,markerfacecolor='w',markeredgecolor='g',lw=0.5,label='Dellasala et al. (2022)')
ax.legend()
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\SummaryBC20k_HAR\InventoryComparison','png',300)

