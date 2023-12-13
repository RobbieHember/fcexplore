#%% Import modules
import os
import numpy as np
import time
import matplotlib.pyplot as plt
import geopandas as gpd
import pandas as pd
import pyproj
import scipy.io as spio
import fiona
import rasterio
from rasterio import features
from shapely import geometry
from scipy.interpolate import griddata
import cv2
import copy
#from rasterio.transform import from_origin
#from shapely.geometry import Point,Polygon
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.bc1ha.bc1ha_util as u1ha
import fcgadgets.cbrunner.cbrun_util as cbu
gp=gu.SetGraphics('Manuscript')

#%% Import paths and look-up-tables
meta=u1ha.Init()
zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])

z=u1ha.Import_Raster(meta,[],['lc_comp1_2019','sp_me','bgcz'],'Extract Grid')

#%%

N_pl_NoSP=0
N_pl_WithSP=0
N_rp_NoSP=0
N_rp_WithSP=0
N_fp_NoSP=0
N_fp_WithSP=0
N_bbp_NoSP=0
N_bbp_WithSP=0
bz='SBPS'
for i in range(6):
    rt=gis.OpenGeoTiff(r'C:\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\PL_All_' + str(i+1) + '_RegenType.tif')['Data']
    yr=gis.OpenGeoTiff(r'C:\Data\BC1ha\RSLT_ACTIVITY_TREATMENT_SVW\PL_All_' + str(i+1) + '_Year.tif')['Data']
    ind1=np.where( (z['lc_comp1_2019']==1) & (z['bgcz']==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'][bz]) & (z['sp_me']==0) & (yr>0) )
    ind2=np.where( (z['lc_comp1_2019']==1) & (z['bgcz']==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'][bz]) & (z['sp_me']==0) & (yr>0) & (rt==meta['LUT']['Derived']['RegenType']['Replanting']) )
    ind3=np.where( (z['lc_comp1_2019']==1) & (z['bgcz']==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'][bz]) & (z['sp_me']>0) & (yr>=z['sp_me']) )
    ind4=np.where( (z['lc_comp1_2019']==1) & (z['bgcz']==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'][bz]) & (z['sp_me']>0) & (yr>=z['sp_me']) & (rt==meta['LUT']['Derived']['RegenType']['Replanting']) )
    ind5=np.where( (z['lc_comp1_2019']==1) & (z['bgcz']==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'][bz]) & (z['sp_me']==0) & (yr>0) & (rt==meta['LUT']['Derived']['RegenType']['Fill Planting']) )
    ind6=np.where( (z['lc_comp1_2019']==1) & (z['bgcz']==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'][bz]) & (z['sp_me']>0) & (yr>=z['sp_me']) & (rt==meta['LUT']['Derived']['RegenType']['Fill Planting']) )
    ind5=np.where( (z['lc_comp1_2019']==1) & (z['bgcz']==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'][bz]) & (z['sp_me']==0) & (yr>0) & (rt==meta['LUT']['Derived']['RegenType']['Back-to-back Planting']) )
    ind6=np.where( (z['lc_comp1_2019']==1) & (z['bgcz']==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'][bz]) & (z['sp_me']>0) & (yr>=z['sp_me']) & (rt==meta['LUT']['Derived']['RegenType']['Back-to-back Planting']) )
    
    N_pl_NoSP=N_pl_NoSP+ind1[0].size
    N_pl_WithSP=N_pl_WithSP+ind3[0].size
    N_rp_NoSP=N_rp_NoSP+ind2[0].size    
    N_rp_WithSP=N_rp_WithSP+ind4[0].size
    N_fp_NoSP=N_fp_NoSP+ind5[0].size    
    N_fp_WithSP=N_fp_WithSP+ind6[0].size
    N_bbp_NoSP=N_bbp_NoSP+ind5[0].size    
    N_bbp_WithSP=N_bbp_WithSP+ind6[0].size

#%%
plt.close('all'); fig,ax=plt.subplots(3,1,figsize=gu.cm2inch(8.5,12)); barw=0.35;
ax[0].bar([0,1],[N_rp_NoSP/N_pl_NoSP*100,N_rp_WithSP/N_pl_WithSP*100])
ax[1].bar([0,1],[N_fp_NoSP/N_pl_NoSP*100,N_fp_WithSP/N_pl_WithSP*100])
ax[2].bar([0,1],[N_bbp_NoSP/N_pl_NoSP*100,N_bbp_WithSP/N_pl_WithSP*100])
ax[0].set(xticks=[0,1],xticklabels=['W/O MSP','With MSP'],ylabel='Frequency (%)',ylim=[0,10])
ax[1].set(xticks=[0,1],xticklabels=['W/O MSP','With MSP'],ylabel='Frequency (%)',ylim=[0,10])
ax[2].set(xticks=[0,1],xticklabels=['W/O MSP','With MSP'],ylabel='Frequency (%)',ylim=[0,10])
ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].tick_params(length=gp['tickl'])
ax[1].yaxis.set_ticks_position('both'); ax[1].xaxis.set_ticks_position('both'); ax[1].tick_params(length=gp['tickl'])
ax[2].yaxis.set_ticks_position('both'); ax[2].xaxis.set_ticks_position('both'); ax[2].tick_params(length=gp['tickl'])
plt.tight_layout()
gu.axletters(ax,plt,0.02,0.89,FontColor=meta['Graphics']['gp']['cla'],LetterStyle='NoPar',FontWeight='Bold',Labels=['Replanting','Fill-planting','Back-to-back planting'],LabelSpacer=0.025)
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Site Prep\SitePrepAndPlantingFrequency_SBPS','png',900)

#%%



