
#%% Import modules

import os
import numpy as np
import gc
import gdal
import matplotlib.pyplot as plt
import matplotlib.colors
import geopandas as gpd
import pandas as pd
import fiona
import time
from shapely.geometry import Polygon,Point,box
import fcgadgets.macgyver.utilities_general as gu
import fcgadgets.macgyver.utilities_gis as gis
#import fcgadgets.macgyver.query_vector_db as qv
from fcgadgets.cbrunner import cbrun_utilities as cbu
from fcgadgets.bc1ha import bc1ha_utilities as bc1ha

#%% Set figure properties

pthout=r'C:\Users\rhember\OneDrive - Government of BC\Figures\Reforestation\Reforestation Monitoring'

params_graphic=cbu.Import_GraphicsParameters('bc1ha_1')
plt.rcParams.update(params_graphic)

figsize1=[900,700]
pos1=[0.04,0.02,0.74,0.95]
pos2=[0.79,0.6,0.03,0.35]

#%% Site Selection Maps

def ClimateSpaceProvWide():

    zT=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\BC1ha_tmin_ann_norm_1971to2000_si_hist_v1.tif')
    zT['Data']=zT['Data'].astype('float')/10
    
    zW=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\BC1ha_ws_gs_norm_1971to2000_comp_hist_v1.tif')
    zW['Data']=zW['Data'].astype('float')
    gc.collect()
    #plt.matshow(zW['Data']);plt.colorbar()

    bsrF={}
    bsrF['grd']=gis.OpenGeoTiff(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Disturbances\VEG_BURN_SEVERITY_SP_2017.tif')
    bsrF['key']=gu.ReadExcel(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Disturbances\VEG_BURN_SEVERITY_SP.xlsx') 
   
    iBurn=np.where( (bsrF['grd']['Data']==bsrF['key']['ID'][ np.where(bsrF['key']['Code']=='High')[0] ]) | (bsrF['grd']['Data']==bsrF['key']['ID'][ np.where(bsrF['key']['Code']=='Medium')[0] ]) )

    # Climate 
    clm={}
    clm['Tmin']=[-10,-10,-10,-8.25,-8.25,-8.25,-5,-5,-5]
    clm['W']=[150,65+(150-65)/2,65, 120,25+(120-27)/2,27, 120,25+(120-27)/2,27]
    clm['Tmin buffer']=0.25
    clm['W buffer']=10
    
    ivl=100
    
    plt.close('all')
    fig,ax=plt.subplots(1,figsize=gu.cm2inch(7.8,7))
    ax.plot(zT['Data'][0::ivl,0::ivl].flatten(),zW['Data'][0::ivl,0::ivl].flatten(),'.',markerfacecolor=[0.75,0.75,0.75],markeredgecolor='None')
    ax.plot(zT['Data'][iBurn].flatten()[0::10],zW['Data'][iBurn].flatten()[0::10],'.',markerfacecolor=[0.8,0.7,0.4],markeredgecolor='None')
    for i in range(len(clm['Tmin'])):
        ax.plot(clm['Tmin'][i],clm['W'][i],'s',markersize=15,markeredgecolor='k',mfc='None')
    ax.set(position=[0.14,0.14,0.8,0.8],xlim=[-14,2],ylim=[0,200],xlabel='Minimum monthly temperature (\circC)',ylabel='Soil water content (mm)')
    
    

#%% Import base maps

bm,tsa,road,district=bc1ha.Import_BaseMaps()

#%% Define region of interest
#
## By TSA
#t0=time.time()
#roi={}
#roi['Type']='ByTSA'
## Pick the TSAs to include
##roi['TSA List']=['Soo TSA']
##roi['TSA List']=['Kamloops TSA','100 Mile House TSA','Williams Lake TSA']
#roi['TSA List']=['Williams Lake TSA']
##roi['TSA List']=list(tsa['key']['Name'])
#roi=bc1ha.DefineROI(roi,tsa,bm,road)
#t1=time.time()
#print((t1-t0)/60)

# By Lat and Long
roi={}
roi['Type']='ByLatLon'
#roi['Centre']=[-123.1308, 51.8157]
#roi['Radius']=2500
# Hanceville Fire
roi['Centre']=[-122.85,51.964553]
roi['Radius']=45*1000
roi=bc1ha.DefineROI(roi,tsa,bm,road)

#%% Reforestation Monitoring Site Selection Maps

def ClimateSelection():

    gdf_wf=qv.GetWildfirePerimiter(2017,2018)
    
    zT_tmp=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\BC1ha_tmin_ann_norm_1971to2000_si_hist_v1.tif')
    zT=tsa['grd'].copy()
    zT['Data']=zT_tmp['Data']
    del zT_tmp
    zT=gis.ClipRaster(zT,roi['xlim'],roi['ylim'])  
    zT['Data']=zT['Data'].astype('float')/10
    gc.collect()    
    
    zW_tmp=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\BC1ha_ws_gs_norm_1971to2000_comp_hist_v1.tif')
    zW=tsa['grd'].copy()
    zW['Data']=zW_tmp['Data']
    del zW_tmp
    zW=gis.ClipRaster(zW,roi['xlim'],roi['ylim'])  
    zW['Data']=zW['Data'].astype('float')
    gc.collect()
    #plt.matshow(zW['Data']);plt.colorbar()

    bsr=bc1ha.Import_Raster_Over_ROI('bsr',roi)
    iBurn=np.where( (bsr['grd']['Data']==bsr['key']['ID'][ np.where(bsr['key']['Code']=='High')[0] ]) | (bsr['grd']['Data']==bsr['key']['ID'][ np.where(bsr['key']['Code']=='Medium')[0] ]) )
    iNoBurn=np.where( (bsr['grd']['Data']!=bsr['key']['ID'][ np.where(bsr['key']['Code']=='High')[0] ]) & (bsr['grd']['Data']!=bsr['key']['ID'][ np.where(bsr['key']['Code']=='Medium')[0] ]) )

    # Climate 
    clm={}
    clm['Tmin']=[-10,-10, -8.8,-7.8, -5,-5]
    clm['W']=   [150,65,  105,40,    120,27]
    clm['Tmin buffer']=0.25
    clm['W buffer']=15    
        
    plt.close('all'); ivl=10
    fig,ax=plt.subplots(1,figsize=gu.cm2inch(7.5,7))
    ax.plot(zT['Data'][0::ivl,0::ivl].flatten(),zW['Data'][0::ivl,0::ivl].flatten(),'.',markerfacecolor=[0.85,0.85,0.85],markeredgecolor='None')
    ax.plot(zT['Data'][iBurn].flatten()[0::10],zW['Data'][iBurn].flatten()[0::10],'.',markerfacecolor=[0.8,0.7,0.4],markeredgecolor='None')
    for i in range(2,4):
        ax.plot(clm['Tmin'][i],clm['W'][i],'s',markersize=8,markeredgecolor='k',mfc='None')
    ax.set(position=[0.14,0.14,0.8,0.8],xlim=[-14,2],ylim=[0,200],xlabel='Minimum monthly temperature (\circC)',ylabel='Soil water content (mm)')
    gu.PrintFig(pthout + '\\ClimateSpace_Hanceville','png',900)
    

    # Grid
    N_CB=6; N_Hidden=1; N_C=N_CB+N_Hidden
    z1=L*np.ones(zT['Data'].shape)
    for i in range(len(clm['Tmin'])):
        ind=np.where( (np.abs(zT['Data']-clm['Tmin'][i])<clm['Tmin buffer']) & (np.abs(zW['Data']-clm['W'][i])<clm['W buffer']) )
        z1[ind]=i
    z1[iNoBurn]=i+1    
    lab=['Cold/Wet','Cold/Dry','Mod/Wet','Mod/Dry','Warm/Wet','Warm/Dry','NA']
    #cm=plt.cm.get_cmap('viridis',6)    #cm=plt.cm.get_cmap('plasma',i)
    cm=np.vstack( ((0,0,0.5,1),(0.5,0.75,1,1), (0.25,0.5,0.25,1),(0.5,1,0.5,1), (1,0.75,0.25,1),(0.5,0,0,1), (0.95,0.95,0.95,1)) )
    cm=matplotlib.colors.ListedColormap(cm)
    
    plt.close('all'); fig,ax=plt.subplots(1,2);ivl=1
    mngr=plt.get_current_fig_manager(); mngr.window.setGeometry(100,100,950,750)
    im=ax[0].matshow(z1[0::ivl,0::ivl],clim=(0,N_C),extent=zT['Extent'],cmap=cm)
    gdf_wf.plot(ax=ax[0],fc='None',ec=[0,0,0],linewidth=0.25,label='Wildfire')
    ax[0].set(position=[0.04,0.02,0.82,0.96],xlim=roi['xlim'],ylim=roi['ylim'],aspect='auto'); ax[0].grid(False)    
    cb=plt.colorbar(im,cax=ax[1],boundaries=np.arange(0,N_C,1),ticks=np.arange(0.5,L+1.5,1))    
    cb.ax.set(yticklabels=lab)    
    cb.ax.tick_params(labelsize=N_CB,length=0)
    for i in range(0,N_C):
        ax[1].plot([0,100],[i/(N_C-1),i/(N_C-1)],'k-',linewidth=0.5)
    ax[1].set(position=[0.87,0.04,0.025,0.25])
    gu.PrintFig(pthout + '\\MapCandidateAreas_Hanceville','png',500)
    
    return

