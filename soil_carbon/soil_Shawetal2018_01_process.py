
"""
Analyiss of Canadian Upland Forest Soils in British Columbia
Shaw et al. (2018)
"""

#%% Import modules

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
import copy
from shapely.geometry import Point, Polygon
from shapely import geometry
from rasterio.transform import from_origin
from osgeo import osr
from scipy.interpolate import griddata
import statsmodels.formula.api as smf
import cv2

import fcgadgets.macgyver.utilities_general as gu
import fcgadgets.macgyver.utilities_gis as gis
import fcgadgets.bc1ha.bc1ha_utilities as u1ha

#%% Import data
meta=u1ha.Init()
gp=gu.SetGraphics('Manuscript')

#%% Import Canadian Upland forest database

flg=0
if flg==1:

    ufd=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\Soils\Shaw et al 2018 Database\SITES.xlsx')
    ufd_p=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\Soils\Shaw et al 2018 Database\PROFILES.xlsx')

    # Add depth from profiles to site
    ufd['Depth']=np.zeros(ufd['LOCATION_ID'].size)
    for i in range(ufd['LOCATION_ID'].size):
        ind=np.where(ufd_p['LOCATION_ID']==ufd['LOCATION_ID'][i])[0]
        if ind.size==0:
            continue
        ind2=np.where( ufd_p['UPPER_HZN_LIMIT'][ind]==np.max(ufd_p['UPPER_HZN_LIMIT'][ind]) )[0]
        ufd['Depth'][i]=ufd_p['UPPER_HZN_LIMIT'][ind[ind2]]+ufd_p['HZN_THICKNESS'][ind[ind2]]

    # Just keep BC
    ind=np.where( (ufd['PROV_TERR']=='BC') )[0]
    for k in ufd.keys():
        ufd[k]=ufd[k][ind]

    # Get BC1ha grid coords
    srs=gis.ImportSRSs()
    ufd['x']=np.zeros(ufd['LOCATION_ID'].size)
    ufd['y']=np.zeros(ufd['LOCATION_ID'].size)
    for i in range(ufd['x'].size):
        ufd['x'][i],ufd['y'][i]=srs['Proj']['BC1ha'](ufd['LONGITUDE'][i],ufd['LATITUDE'][i])

    # Import BC1ha data
    z=u1ha.Import_Raster(meta,[],['bgcz'])
    ind=gis.GetGridIndexToPoints(z['bgcz'],ufd['x'],ufd['y'])
    ufd['bgcz']=z['bgcz']['Data'][ind]

    # Delete no BGC zone
    # Only two
    ind=np.where( (ufd['bgcz']==0) )[0]
    print(ind.size)

    ind=np.where( (ufd['bgcz']!=0) )[0]
    for k in ufd.keys():
        ufd[k]=ufd[k][ind]

    # Save to pickle
    gu.opickle(r'C:\Users\rhember\Documents\Data\Soils\Shaw et al 2018 Database\SITES.pkl',ufd)

else:

    ufd=gu.ipickle(r'C:\Users\rhember\Documents\Data\Soils\Shaw et al 2018 Database\SITES.pkl')

#%% Model 1 - BGC zone only

df=pd.DataFrame.from_dict(ufd)
form="TOT_C_THA ~ C(bgcz)"
mr1=smf.ols(formula=form, data=df).fit()
print(mr1.summary())
dP1=mr1.params.to_dict()

#%% Model 2 - soil type (Subgroup.Greatgroup)

df=pd.DataFrame.from_dict(ufd)
form="TOT_C_THA ~ C(CSSC_CODE)"
mr2=smf.ols(formula=form, data=df).fit()
print(mr2.summary())
dP2=mr2.params.to_dict()

#%% Model 3 - BGC zone + soil type (Subgroup.Greatgroup)

df=pd.DataFrame.from_dict(ufd)
form="TOT_C_THA ~ C(bgcz) + C(CSSC_CODE)"
mr3=smf.ols(formula=form,data=df).fit()
print(mr3.summary())
dP3=mr3.params.to_dict()

#%% Import BC Soil survey data

bcss=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Soil\BC Soil Surveys\soil_dev1a.tif')
bcss=gis.ClipToRaster(bcss,zRef)
lut_bcss=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\BC1ha\Soil\BC Soil Surveys\soil_dev1.tif.vat.dbf.xlsx')

plt.matshow(bcss['Data'])

#%% Import SLC
# Dead end - they don't include any way to reconstruct great group/subgroup outside the soil surveys

# slc=gpd.read_file(r'C:\Users\rhember\Documents\Data\Soils\Soil Landscapes of Canada\Version32\ca_all_slc_v3r2.shp')

# # Clip
# slc=slc.cx[-142:-112,47:61]
# slc=slc.reset_index(drop=True)

# # Reproject
# rst=rasterio.open(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\tsa.tif')
# slc=slc.to_crs(rst.crs.to_dict())

# out_arr=np.zeros(rst.shape,dtype=float)
# # This is where we create a generator of geom, value pairs to use in rasterizing
# shapes=((geom,value) for geom, value in zip(slc.geometry,slc.POLY_ID))
# z_slc=features.rasterize(shapes=shapes,fill=0,out=out_arr,transform=zRef['Transform'])

# # Delete info overlapping BC soil surveys
# # z_slc2=z_slc.copy()
# # ind=np.where( (bcss['Data']>=1) & (bcss['Data']<128) )
# # z_slc2[ind]=0

# lut_slc1=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\Soils\Soil Landscapes of Canada\Version32\ca_all_slc_v3r2_cmp.dbf.xlsx')
# # ind=np.where(lut_slc1['CMP']==1)[0]
# # for k in lut_slc1.keys():
# #     lut_slc1[k]=lut_slc1[k][ind]

# lut_slc2=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\Soils\Soil Landscapes of Canada\soil_name_bc_v2r20150610.xlsx')

# u_pid=np.unique(z_slc[z_slc>0])
# nam_slc=np.array(['' for _ in range(u_pid.size)],dtype=object)

# ind=np.where(lut_slc1['PROVINCE']=='BC')[0]
# np.unique(lut_slc1['POLY_ID'][ind]).size

# for i in range(u_pid.size):
#     ind2=np.where( lut_slc1['POLY_ID']==u_pid[i] )
#     if ind2[0].size==0:
#         ind=np.where(z_slc==u_pid[i])
#         a[ind]=1
#         continue
#     soil_name=lut_slc1['SOIL_ID'][ind2]

#     ind3=np.where( lut_slc2['SOIL_ID']==soil_name )
#     if ind3[0].size==0:
#         continue
#     gg=lut_slc2['G_GROUP2'][ind3]
#     sg=lut_slc2['S_GROUP2'][ind3]
#     nam_slc[i]=sg[0] + '.' + gg[0]
# u_nam_slc=np.unique(nam_slc)
# u_nam_slc=u_nam_slc[u_nam_slc!='']
# u_nam_slc=u_nam_slc[u_nam_slc!='-.-']

#%% Map predictions

flg=0
if flg==1:

    # Model 1
    zSOC1=copy.deepcopy(zRef)
    zSOC1['Data']=dP1['Intercept']+0*zSOC1['Data']

    # Add BEC
    u=np.unique(bgcz['Data'][bgcz['Data']<255])
    for i in range(u.size):
        ind1=np.where( bgcz['Data']==u[i] )
        for k in dP1:
            if k=='C(bgcz)[T.' + str(u[i]) + ']':
                break
        zSOC1['Data'][ind1]=zSOC1['Data'][ind1]+dP1[k]

    # Model 2
    zSOC2=copy.deepcopy(zRef)
    zSOC2['Data']=0*zSOC2['Data']
    for i in range(u_nam_slc.size):
        ind1=np.where( nam_slc==u_nam_slc[i] )[0]
        list=[]
        for j in range(ind1.size):
            list.append(u_pid[ind1[j]])
        ind2=np.where(np.isin(z_slc,list)==True)
        #zSOC2['Data'][ind2]=i+1
        #plt.matshow(zSOC2['Data'],clim=[0,40])
        nam=u_nam_slc[i]
        for k in dP2:
            if k=='C(CSSC_CODE)[T.' + nam + ']':
                zSOC2['Data'][ind2]=dP2['Intercept']+dP2[k]

    # # Add soil type from BC Soil Surveys
    # u=np.unique(bcss['Data'])
    # for i in range(u.size):
    #     ind1=np.where( bcss['Data']==u[i])
    #     ind2=np.where(lut_bcss['VALUE']==u[i])[0]
    #     if ind2.size==0:
    #         continue
    #     nam=lut_bcss['DEV_1'][ind2][0]
    #     for k in dP2:
    #         if k=='C(CSSC_CODE)[T.' + nam + ']':
    #             #zSOC2['Data'][ind1]=zSOC2['Data'][ind1]+dP2[k]
    #             zSOC2['Data'][ind1]=dP2[k]

    # Model 3
    zSOC3=copy.deepcopy(zRef)
    zSOC3['Data']=dP3['Intercept']+0*zSOC3['Data']

    # Add BEC
    u=np.unique(bgcz['Data'][bgcz['Data']<255])
    for i in range(u.size):
        ind1=np.where( bgcz['Data']==u[i] )
        for k in dP3:
            if k=='C(bgcz)[T.' + str(u[i]) + ']':
                break
        zSOC3['Data'][ind1]=zSOC3['Data'][ind1]+dP3[k]

    # Add soil type from BC Soil Surveys
    u=np.unique(bcss['Data'])
    for i in range(u.size):
        ind1=np.where( bcss['Data']==u[i])
        ind2=np.where(lut_bcss['VALUE']==u[i])[0]
        if ind2.size==0:
            continue
        nam=lut_bcss['DEV_1'][ind2][0]
        for k in dP3:
            if k=='C(CSSC_CODE)[T.' + nam + ']':
                zSOC3['Data'][ind1]=zSOC3['Data'][ind1]+dP3[k]

    # Combine model 1 and model 3
    zSOC=copy.deepcopy(zSOC1)
    ind=np.where( (bcss['Data']>0) & (bcss['Data']<=95) & (bcss['Data']!=5) & (bcss['Data']!=32) )
    zSOC['Data'][ind]=zSOC3['Data'][ind]

    # Save map
    zSOC['Data']=zSOC['Data'].astype('int16')
    gis.SaveGeoTiff(zSOC,r'C:\Users\rhember\Documents\Data\BC1ha\Soil\soc_tot_forest_Shawetal2018.tif')

else:
    # Import predicted SOC map
    zSOC=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Soil\soc_tot_forest_Shawetal2018.tif')

#plt.close('all')
#plt.matshow(zSOC['Data'],clim=[0,450])

#%% Calculate mean SOC

# Forest mask
btm=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\LandUseLandCover\landuse.btm.tif')
iMask=np.where( (btm['Data']==8) | (btm['Data']==11) | (btm['Data']==12) | (btm['Data']==16) | (btm['Data']==21) )

np.sum(zSOC['Data'][iMask])/1e9

#%% Plot by BGC zone

u=np.unique(ufd['bgcz'])

lab=np.array(['' for _ in range(u.size)],dtype=object)
soc={}
soc['N']=np.zeros(u.size)
soc['mu']=np.zeros(u.size)
soc['se']=np.zeros(u.size)
soc['min_mu']=np.zeros(u.size)
soc['org_mu']=np.zeros(u.size)
soc['Model A']=np.zeros(u.size)
soc['Model Soil C']=np.zeros(u.size)
for i in range(u.size):
    ind=np.where( (ufd['bgcz']==u[i]) & (ufd['TOT_C_THA']>0) )[0]
    soc['N'][i]=ind.size
    soc['mu'][i]=np.nanmean(ufd['TOT_C_THA'][ind])
    soc['se'][i]=np.nanstd(ufd['TOT_C_THA'][ind])/np.sqrt(ind.size)
    soc['min_mu'][i]=np.nanmean(ufd['MIN_C_THA'][ind])
    soc['org_mu'][i]=np.nanmean(ufd['ORG_C_THA'][ind])
    #soc['Model A'][i]=np.nanmean(mos['Scenarios'][0]['Mean']['A']['Ensemble Mean'][-1,0,ind])
    #soc['Model Soil C'][i]=np.nanmean(mos['Scenarios'][0]['Mean']['C_Soil_Tot']['Ensemble Mean'][-1,0,ind])
    ind=np.where(lutBGC['VALUE']==u[i])[0]
    if ind.size>0:
        lab[i]=lutBGC['ZONE'][ind][0]

ord=np.argsort(soc['mu'])
for k in soc:
    soc[k]=np.flip(soc[k][ord])
lab=np.flip(lab[ord])

#%%

plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,6))
ax.bar(np.arange(u.size),soc['min_mu'],facecolor=[0.45,0.3,0.3],label='Mineral horizons')
ax.bar(np.arange(u.size),soc['org_mu'],facecolor=[0.15,0.05,0.05],bottom=soc['min_mu'],label='Organic horizon')
ax.errorbar(np.arange(u.size),soc['mu'],yerr=soc['se'],color='k',fmt='none',capsize=2)
#ax.plot(np.arange(u.size),soc['Model Soil C'],'bo',lw=1,ms=8,mfc='w')
for i in range(u.size):
    ax.text(i,10,str(soc['N'][i].astype(int)),color='w',ha='center',fontsize=8)
ax.set(position=[0.08,0.12,0.9,0.86],xlim=[-0.5,u.size-0.5],ylim=[0,375],xticks=np.arange(u.size),
       xticklabels=lab,ylabel='Soil organic carbon (MgC ha$^{-1}$ yr$^{-1}$)')
plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')
#gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Soils\Bars_BGCSpinupRI' ,'png',900)

#%% Scatterplot

x=soc['min_mu']
y=soc['Model Soil C']
rs,txt=gu.GetRegStats(x,y)

plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(11,11))
ax.plot([0,1000],[0,1000],'-k',lw=3,color=[0.8,0.8,0.8])
ax.plot(x,y,'ko',mfc='k',mec='w',lw=0.5,ms=6)
ax.plot(rs['xhat'],rs['yhat'],'r-',lw=1,label='Best fit')
ax.text(350,30,txt,fontsize=10,color='r',ha='right')
ax.text(300,300,'1:1',fontsize=8,ha='center')
ax.set(position=[0.1,0.1,0.86,0.86],xlabel='Observed SOC (MgC ha$^{-1}$)',ylabel='Predicted SOC (MgC ha$^{-1}$)',xlim=[0,375],ylim=[0,375])
#plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Soils\Scatter_BGCSpinupRI' ,'png',900)

#%%

plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,6))
ax.plot(np.arange(u.size),soc['Model A'],'bo',lw=1,ms=8,mfc='w')

#%% Confirm that map predictions agree with point sample

soc['map_mu']=np.zeros(u.size)
for i in range(u.size):
    ind=np.where( bgcz['Data']==u[i] )
    soc['map_mu'][i]=np.nanmean(zSOC['Data'][ind])

# Reorder
soc['map_mu']=np.flip(soc['map_mu'][ord])

#%% Plot obs vs predicted comparison

plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,6))
ax.bar(np.arange(u.size)-0.2,soc['mu'],0.4,facecolor=[0.27,0.49,0.77],label='Observed sample')
ax.bar(np.arange(u.size)+0.2,soc['map_mu'],0.4,facecolor=[0.65,0.95,0.2],label='Predictions from map')
ax.errorbar(np.arange(u.size)-0.2,soc['mu'],yerr=soc['se'],color='k',fmt='none',capsize=2)
ax.set(position=[0.08,0.12,0.9,0.86],xlim=[-0.5,u.size-0.5],ylim=[0,375],xticks=np.arange(u.size),
       xticklabels=lab,ylabel='Soil organic carbon (MgC ha$^{-1}$ yr$^{-1}$)')
plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Soils\SOC_By_BGCZone_ObsVsPred' ,'png',900)

#%%

Mask=0*zRef['Data'].copy()
ind=np.where( (bcss['Data']>0) & (bcss['Data']<=95) & (bcss['Data']!=5) & (bcss['Data']!=32) )
Mask[ind]=1
plt.matshow(Mask)

# Create binary image
z=np.zeros((m,n,3),dtype=np.uint8)
z[tsa==id,:]=255
z=cv2.cvtColor(z,cv2.COLOR_BGR2GRAY) # Convert to grey scale

# Calculate contour of object
cont=cv2.findContours(image=Mask,mode=cv2.RETR_LIST,method=cv2.CHAIN_APPROX_SIMPLE)







