
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
import fiona
import statsmodels.formula.api as smf
import fcgadgets.macgyver.utilities_general as gu
import fcgadgets.macgyver.utilities_gis as gis
import fcexplore.psp.Processing.psp_utilities as ugp
from fcgadgets.bc1ha import bc1ha_utilities as u1ha

#%% Set figure properties

gp=gu.SetGraphics('Manuscript')

#%% Import data

metaGP,gplt=ugp.ImportPSPs(type='Stand')

#%% Import Raster grids

lut_1ha=u1ha.Import_BC1ha_LUTs()

zBGCZ=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\becz.tif')

zLCC1=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\LandUseLandCover\LandCoverClass1.tif')

#%% Compare plot AGB with VRI polygon estimates

# Create geodatabase for plot locations
points=[]
for i in range(gplt['X'].size):
    points.append( Point(gplt['X'][i],gplt['Y'][i]) )
d={'geometry':points}
for k in sl.keys():
    d[k]=gplt[k]
gdf=gpd.GeoDataFrame(d)
gdf_bm=gpd.read_file(r'C:\Users\rhember\Documents\Data\Basemaps\Basemaps.gdb',layer='NRC_POLITICAL_BOUNDARIES_1M_SP')
gdf.crs=gdf_bm.crs
#gdf.to_file(r'C:\Users\rhember\Documents\Data\GroundPlots\PSP-NADB\ground_plots.geojson',driver='GeoJSON')

#ogsr=gpd.read_file(r'C:\Users\rhember\Documents\Data\ForestInventory\LandUse\20220422\LandUse.gdb',layer='OGSR_TAP_PRIORITY_DEF_AREA_SP')
#slo=gpd.overlay(gdf,ogsr,how='intersection')

# Import VRI
flg=0
if flg==1:
    t0=time.time()
    pth=r'C:\Users\rhember\Documents\Data\ForestInventory\VRI\20220404\VRI.gdb'
    List=[None]*int(1e5)
    cnt=0
    with fiona.open(pth,layer='VEG_COMP_LYR_R1_POLY') as source:
        for feat in source:
            if feat['geometry']==None:
                continue
            shp=shape(feat['geometry'])
            ind=np.where(gdf.within(shp)==True)[0]
            if ind.size>0:
                List[cnt]=feat
                cnt=cnt+1
                #print('working')
    List=List[0:cnt-1]
    #gu.opickle(r'C:\Users\rhember\Documents\Data\GroundPlots\PSP-NADB\ground_plots_vri.pkl',List)
    gdf_vri=gpd.GeoDataFrame.from_features(List,crs=gdf_bm.crs)
    print((time.time()-t0)/60)
    #gdf_vri.to_file(meta['Paths']['DB'] + '\\ground_plots_vri.geojson',driver='GeoJSON')
else:
    gdf_vri=gpd.read_file(metaGP['Paths']['DB'] + '\\ground_plots_vri.geojson')

# Spatial join of VRI polygons and PSPs
gdf_vri_j=gdf_vri.sjoin(gdf,how='inner')
#gdf_vri_j=gpd.sjoin(gdf,gdf_vri,how='left'); #gdf_vri_j=gdf.sjoin(gdf_vri,how='inner')

# Plot to make sure it worked
flg=0
if flg==1:
    plt.close('all')
    fig,ax=plt.subplots(1,figsize=gu.cm2inch(12,12))
    gdf_vri_j.plot(ax=ax,edgecolor='r',facecolor='y',linewidth=0.5,label='Road',alpha=1)
    gdf.plot(ax=ax,edgecolor='g',linewidth=0.5,markersize=14,label='Road',alpha=1)

# Comparison of plot and VRI polygon
df_vri2=pd.DataFrame(gdf_vri_j.drop(columns='geometry'))
d_vri2={}
for k in df_vri2.columns:
    d_vri2[k]=df_vri2[k].values
d_vri2['AGB_VRI']=d_vri2['WHOLE_STEM_BIOMASS_PER_HA']+d_vri2['BRANCH_BIOMASS_PER_HA']+d_vri2['FOLIAGE_BIOMASS_PER_HA']+d_vri2['BARK_BIOMASS_PER_HA']

#%% Plot comparison between ground plot AGB and VRI polygon AGB

# Isolate a sample
ind=np.where( (d_vri2['Ecozone BC L1']>0) & (d_vri2['Flag Plot Types To Keep']==1) & (d_vri2['Age t0']>0) & (d_vri2['Cag L t0']>=0) & (d_vri2['Cag L t0']<2500) & (d_vri2['AGB_VRI']>0) )[0]

x=np.real(d_vri2['Cag L t0'][ind])
y=0.5*d_vri2['AGB_VRI'][ind]
rs,txt=gu.GetRegStats(x,y);

# Import Dellasala et al. 2022
zDS22=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\Dellasala et al 2022 Fig7 dig.xlsx')

plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(10,10))
ax.plot([0,1000],[0,1000],'-k',lw=3,color=[0.8,0.8,0.8])
ax.plot(x,y,'ko',markersize=2.5,markerfacecolor=[0.25,0.25,0.25],markeredgecolor='w',markeredgewidth=0.25,label='Provincial Ground Plots (VRI + CMI)')
ax.plot(rs['xhat'],rs['yhat'],'k-',lw=1,label='Best fit (Provincial Ground Plots)')
ax.text(125,525,txt,fontsize=10,color='k')
ax.text(770,770,'1:1',fontsize=8,ha='center')
ax.set(position=[0.12,0.12,0.84,0.84],xlim=[0,1000],ylim=[0,1000],xlabel='AGB ground sample (MgC/ha)',ylabel='AGB VRI polygon (MgC/ha)')
ax.plot(zDS22['Sample Plot'],zDS22['VRI polygon'],'gs',markersize=3,markerfacecolor='w',markeredgecolor=[0.4,0.8,0],markeredgewidth=0.65,label='Dellasala et al. (2022)')
ax.legend(frameon=False)
#gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Biomass\Biomass_GroundPlotsVsVRIPolygons','png',300)

#%% Plot by BGC zone

vL=['Cag L t0','Csw L t0','Cbr L t0','Cbk L t0','Cf L t0','AGB_VRI','WHOLE_STEM_BIOMASS_PER_HA','BRANCH_BIOMASS_PER_HA','BARK_BIOMASS_PER_HA','FOLIAGE_BIOMASS_PER_HA']
u=np.unique(d_vri2['Ecozone BC L1'])
lab=np.array(['' for _ in range(u.size)],dtype=object)
d={}
for v in vL:
    d[v]={}
    d[v]['N']=np.zeros(u.size)
    d[v]['mu']=np.zeros(u.size)
    d[v]['sd']=np.zeros(u.size)
    d[v]['se']=np.zeros(u.size)
for i in range(u.size):
    lab[i]=ugp.lut_id2cd(metaGP,'Ecozone BC L1',u[i])
    for v in vL:
        ind=np.where( (d_vri2['Ecozone BC L1']==u[i]) & (d_vri2['Flag Plot Types To Keep']==1) & (d_vri2[v]>=0) & (d_vri2[v]<2000) )[0]
        d[v]['N'][i]=ind.size
        d[v]['mu'][i]=np.nanmean(d_vri2[v][ind])
        d[v]['sd'][i]=np.nanstd(d_vri2[v][ind])
        #d[v]['se'][i]=np.nanstd(d_vri2[v][ind])/np.sqrt(ind[0].size)

# Put in order
v='Cag L t0'
ord=np.argsort(d[v]['mu'])
uo=u[ord]
lab=np.flip(lab[ord])
for v in d:
    for k in d[v].keys():
        d[v][k]=np.flip(d[v][k][ord])

# Plot
cl=np.array([[0.22,0.49,0.77],[0.6,0.9,0.25]])
plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,6))
ax.bar(np.arange(u.size)-0.17,d['Cag L t0']['mu'],0.3,facecolor=cl[0,:],label='Ground plots')
ax.bar(np.arange(u.size)+0.17,0.5*d['AGB_VRI']['mu'],0.3,facecolor=cl[1,:],label='VRI polygons')
#ax.errorbar(np.arange(u.size),d['Ctot L t0']['mu'],yerr=d['Ctot L t0']['se'],color='k',fmt='none',capsize=2,lw=0.5)
ax.set(position=[0.08,0.12,0.9,0.86],xlim=[-0.5,u.size-0.5],ylim=[0,300],xticks=np.arange(u.size),
       xticklabels=lab,ylabel='AGB (MgC ha$^{-1}$ yr$^{-1}$)')
plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Biomass\Biomass_ComparingPlotsAndVRIPolygons_ByGBCZone','png',900)


# # Plot
# plt.close('all')
# fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,6))
# ax.bar(np.arange(u.size)-0.17,d['Csw L t0']['mu'],0.3,facecolor=cl[0,:],label='Ground plots')
# ax.bar(np.arange(u.size)+0.17,0.5*d['WHOLE_STEM_BIOMASS_PER_HA']['mu'],0.3,facecolor=cl[1,:],label='VRI polygons')
# #ax.errorbar(np.arange(u.size),d['Ctot L t0']['mu'],yerr=d['Ctot L t0']['se'],color='k',fmt='none',capsize=2,lw=0.5)
# ax.set(position=[0.08,0.12,0.9,0.86],xlim=[-0.5,u.size-0.5],ylim=[0,300],xticks=np.arange(u.size),
#        xticklabels=lab,ylabel='AGB (MgC ha$^{-1}$ yr$^{-1}$)')
# plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
# ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')
# #gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Biomass\Biomass_ComparingPlotsAndVRIPolygons_ByGBCZone','png',900)

#%%

plt.close('all')
fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(16,16))
ax[0,0].plot([0,200],[0,200],'k--')
ax[0,0].plot(d['Csw L t0']['mu'],0.5*d['WHOLE_STEM_BIOMASS_PER_HA']['mu'],'k.')
ax[0,1].plot([0,200],[0,200],'k--')
ax[0,1].plot(d['Cbk L t0']['mu'],0.5*d['BARK_BIOMASS_PER_HA']['mu'],'k.')
ax[1,0].plot([0,200],[0,200],'k--')
ax[1,0].plot(d['Cbr L t0']['mu'],0.5*d['BRANCH_BIOMASS_PER_HA']['mu'],'k.')
ax[1,1].plot([0,200],[0,200],'k--')
ax[1,1].plot(d['Cf L t0']['mu'],0.5*d['FOLIAGE_BIOMASS_PER_HA']['mu'],'k.')
#ax[1,1].plot(d['Cbk L t0']['mu'],0.5*d['BARK_BIOMASS_PER_HA']['mu'],'k.')
#ax[1,2].plot(d['Cbr L t0']['mu'],0.5*d['BRANCH_BIOMASS_PER_HA']['mu'],'k.')
