
"""
Analyiss of Canadian Upland Forest Soils in British Columbia
Shaw et al. (2018)
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
import copy
from shapely.geometry import Point, Polygon
from shapely import geometry
from rasterio.transform import from_origin
from osgeo import osr
from scipy.interpolate import griddata
import statsmodels.formula.api as smf

import fcgadgets.macgyver.utilities_general as gu
import fcgadgets.macgyver.utilities_gis as gis
from fcgadgets.macgyver import utilities_inventory as invu

#%% Set figure properties

fs=7
params={'font.sans-serif':'Arial','font.size':fs,'axes.edgecolor':'black','axes.labelsize':fs,'axes.labelcolor':'black','axes.titlesize':fs,'axes.linewidth':0.5,'lines.linewidth':0.5,
        'text.color':'black','xtick.color':'black','xtick.labelsize':fs,'xtick.major.width':0.5,'xtick.major.size':3,'xtick.direction':'in','ytick.color':'black','ytick.labelsize':fs,
        'ytick.major.width':0.5,'ytick.major.size':3,'ytick.direction':'in','legend.fontsize':fs,'savefig.dpi':300,'savefig.transparent':True}
plt.rcParams.update(params)

#%% Import Canadian Upland forest database

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

#%% Get BC1ha grid coords

zTSA=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\tsa.tif')

srs=gis.ImportSRSs()
ufd['x']=np.zeros(ufd['LOCATION_ID'].size)
ufd['y']=np.zeros(ufd['LOCATION_ID'].size)
for i in range(ufd['x'].size):
    ufd['x'][i],ufd['y'][i]=srs['Proj']['BC1ha'](ufd['LONGITUDE'][i],ufd['LATITUDE'][i])

x=zTSA['X'][0,:]
y=zTSA['Y'][:,0]
ix=np.zeros(ufd['LOCATION_ID'].size,dtype=int)
iy=np.zeros(ufd['LOCATION_ID'].size,dtype=int)
for i in range(ufd['x'].size):
    ix[i]=np.where( np.abs(ufd['x'][i]-x)==np.min(np.abs(ufd['x'][i]-x)) )[0]
    iy[i]=np.where( np.abs(ufd['y'][i]-y)==np.min(np.abs(ufd['y'][i]-y)) )[0]

# Check that it worked
x=zTSA['X'][iy,ix]
y=zTSA['Y'][iy,ix]
plt.plot(x,y,'k.')

# Import BC1ha data
becz=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\becz.tif')
lutBGC=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\becz_lut.xlsx')
ufd['becz']=becz['Data'][iy,ix]

#elev=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Terrain\elevation.tif')
#ws=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\BC1ha_ws_gs_norm_1971to2000_comp_hist_v1.tif')
#dwf=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\BC1ha_dwf_ann_norm_1971to2000_si_hist_v1.tif')
#age1=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\age1.tif')
#gsoc=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Soil\gsoc2010_bc1ha.tif')

#ufd['elev']=elev['Data'][iy,ix]
#ufd['gsoc']=gsoc['Data'][iy,ix]
#ufd['ws']=ws['Data'][iy,ix]
#ufd['dwf']=dwf['Data'][iy,ix]
#ufd['age1']=age1['Data'][iy,ix]
#ufd['age1'][np.where(ufd['age1']<=0)[0]]=0

#del elev,ws,age1,gsoc,becz,dwf
#del gsoc
garc.collect()

#%% Plot

# ind=np.where(ufd['gsoc']>0)[0]
# plt.close('all')
# plt.plot(ufd['MIN_C_THA'][ind],ufd['gsoc'][ind],'.')

#%% Delete no BGC zone

# Only two
ind=np.where( (ufd['becz']==255) )[0]

ind=np.where( (ufd['becz']!=255) )[0]
for k in ufd.keys():
    ufd[k]=ufd[k][ind]

#%% Model 1 - BGC zone only

df=pd.DataFrame.from_dict(ufd)
form="TOT_C_THA ~ C(becz)"
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
form="TOT_C_THA ~ C(becz) + C(CSSC_CODE)"
mr3=smf.ols(formula=form,data=df).fit()
print(mr3.summary())
dP3=mr3.params.to_dict()

#%% Import BC Soil survey data

bcss=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Soil\soil_dev1a.tif')
bcss=gis.ClipToRaster(bcss,zTSA)
lut_bcss=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\BC1ha\Soil\soil_dev1.tif.vat.dbf.xlsx')

#%% Import SLC

slc=gpd.read_file(r'C:\Users\rhember\Documents\Data\Soils\Soil Landscapes of Canada\Version32\ca_all_slc_v3r2.shp')

# Clip
slc=slc.cx[-142:-112,47:61]
slc=slc.reset_index(drop=True)

# Reproject
rst=rasterio.open(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\tsa.tif')
slc=slc.to_crs(rst.crs.to_dict())

out_arr=np.zeros(rst.shape,dtype=float)
# This is where we create a generator of geom, value pairs to use in rasterizing
shapes=((geom,value) for geom, value in zip(slc.geometry,slc.POLY_ID))
z_slc=features.rasterize(shapes=shapes,fill=0,out=out_arr,transform=out.transform)

# Delete info overlapping BC soil surveys
# z_slc2=z_slc.copy()
# ind=np.where( (bcss['Data']>=1) & (bcss['Data']<128) )
# z_slc2[ind]=0

lut_slc1=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\Soils\Soil Landscapes of Canada\Version32\ca_all_slc_v3r2_cmp.dbf.xlsx')
# ind=np.where(lut_slc1['CMP']==1)[0]
# for k in lut_slc1.keys():
#     lut_slc1[k]=lut_slc1[k][ind]

lut_slc2=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\Soils\Soil Landscapes of Canada\soil_name_bc_v2r20150610.xlsx')

u_pid=np.unique(z_slc[z_slc>0])
nam_slc=np.array(['' for _ in range(u_pid.size)],dtype=object)

ind=np.where(lut_slc1['PROVINCE']=='BC')[0]
np.unique(lut_slc1['POLY_ID'][ind]).size

for i in range(u_pid.size):
    ind2=np.where( lut_slc1['POLY_ID']==u_pid[i] )
    if ind2[0].size==0:
        ind=np.where(z_slc==u_pid[i])
        a[ind]=1
        continue
    soil_name=lut_slc1['SOIL_ID'][ind2]

    ind3=np.where( lut_slc2['SOIL_ID']==soil_name )
    if ind3[0].size==0:
        continue
    gg=lut_slc2['G_GROUP2'][ind3]
    sg=lut_slc2['S_GROUP2'][ind3]
    nam_slc[i]=sg[0] + '.' + gg[0]
u_nam_slc=np.unique(nam_slc)
u_nam_slc=u_nam_slc[u_nam_slc!='']
u_nam_slc=u_nam_slc[u_nam_slc!='-.-']

#%% Import bec zone

becz=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\becz.tif')

#%% Map predictions

flg=0
if flg==1:

    # Model 1
    zSOC1=copy.deepcopy(zTSA)
    zSOC1['Data']=dP1['Intercept']+0*zSOC1['Data']

    # Add BEC
    u=np.unique(becz['Data'][becz['Data']<255])
    for i in range(u.size):
        ind1=np.where( becz['Data']==u[i] )
        for k in dP1:
            if k=='C(becz)[T.' + str(u[i]) + ']':
                break
        zSOC1['Data'][ind1]=zSOC1['Data'][ind1]+dP1[k]

    # Model 2
    zSOC2=copy.deepcopy(zTSA)
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

    # Add soil type from BC Soil Surveys
    u=np.unique(bcss['Data'])
    for i in range(u.size):
        ind1=np.where( bcss['Data']==u[i])
        ind2=np.where(lut_bcss['VALUE']==u[i])[0]
        if ind2.size==0:
            continue
        nam=lut_bcss['DEV_1'][ind2][0]
        for k in dP2:
            if k=='C(CSSC_CODE)[T.' + nam + ']':
                #zSOC2['Data'][ind1]=zSOC2['Data'][ind1]+dP2[k]
                zSOC2['Data'][ind1]=dP2[k]

    # Model 3
    zSOC2=copy.deepcopy(zTSA)
    zSOC2['Data']=dP2['Intercept']+0*zSOC2['Data']

    # Add BEC
    u=np.unique(becz['Data'][becz['Data']<255])
    for i in range(u.size):
        ind1=np.where( becz['Data']==u[i] )
        for k in dP2:
            if k=='C(becz)[T.' + str(u[i]) + ']':
                break
        zSOC2['Data'][ind1]=zSOC2['Data'][ind1]+dP2[k]

    # Add soil type from Soil Landscapes of Canada
    for i in range(u_nam_slc.size):
        ind1=np.where( nam_slc==u_nam_slc[i] )[0]
        #print(ind1.size)
        list=[]
        for j in range(ind1.size):
            list.append(u_pid[ind1[j]])
        ind2=np.where(np.isin(z_slc,list)==True)
        nam=u_nam_slc[i]
        for k in dP2:
            if k=='C(CSSC_CODE)[T.' + nam + ']':
                zSOC2['Data'][ind1]=zSOC2['Data'][ind1]+dP2[k]

    # Add soil type from BC Soil Surveys
    u=np.unique(bcss['Data'])
    for i in range(u.size):
        ind1=np.where( bcss['Data']==u[i])
        ind2=np.where(lut_bcss['VALUE']==u[i])[0]
        if ind2.size==0:
            continue
        nam=lut_bcss['DEV_1'][ind2][0]
        for k in dP2:
            if k=='C(CSSC_CODE)[T.' + nam + ']':
                zSOC2['Data'][ind1]=zSOC2['Data'][ind1]+dP2[k]

    # Combined
    zSOC=copy.deepcopy(zSOC1)
    ind=np.where( (bcss['Data']>0) & (bcss['Data']<=95) & (bcss['Data']!=5) & (bcss['Data']!=32) )
    #ind=np.where( (zSOC2['Data']>0) & (zSOC2['Data']<250) )
    zSOC['Data'][ind]=zSOC2['Data'][ind]

    # Save map
    zSOC['Data']=zSOC['Data'].astype('int16')
    gis.SaveGeoTiff(zSOC,r'C:\Users\rhember\Documents\Data\BC1ha\Soil\soc_tot_forest_Shawetal2018.tif')

#%% Import predicted SOC map

zSOC=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Soil\soc_tot_forest_Shawetal2018.tif')

#plt.close('all')
#plt.matshow(zSOC['Data'],clim=[0,250])

#%% Calculate mean SOC

# Forest mask
btm=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\LandUseLandCover\landuse.btm2.tif')
iMask=np.where( (btm['Data']==8) | (btm['Data']==11) | (btm['Data']==12) | (btm['Data']==16) | (btm['Data']==21) )

np.sum(zSOC['Data'][iMask])/1e9

#%% Plot by BGC zone

u=np.unique(ufd['becz'])

lab=np.array(['' for _ in range(u.size)],dtype=object)
soc={}
soc['N']=np.zeros(u.size)
soc['mu']=np.zeros(u.size)
soc['se']=np.zeros(u.size)
soc['min_mu']=np.zeros(u.size)
soc['org_mu']=np.zeros(u.size)
for i in range(u.size):
    ind=np.where( (ufd['becz']==u[i]) & (ufd['TOT_C_THA']>0) )[0]
    soc['N'][i]=ind.size
    soc['mu'][i]=np.nanmean(ufd['TOT_C_THA'][ind])
    soc['se'][i]=np.nanstd(ufd['TOT_C_THA'][ind])/np.sqrt(ind.size)
    soc['min_mu'][i]=np.nanmean(ufd['MIN_C_THA'][ind])
    soc['org_mu'][i]=np.nanmean(ufd['ORG_C_THA'][ind])
    ind=np.where(lutBGC['VALUE']==u[i])[0]
    if ind.size>0:
        lab[i]=lutBGC['ZONE'][ind][0]

ord=np.argsort(soc['mu'])
for k in soc:
    soc[k]=np.flip(soc[k][ord])
lab=np.flip(lab[ord])

plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,6))
ax.bar(np.arange(u.size),soc['min_mu'],facecolor=[0.45,0.3,0.3],label='Mineral horizons')
ax.bar(np.arange(u.size),soc['org_mu'],facecolor=[0.15,0.05,0.05],bottom=soc['min_mu'],label='Organic horizon')
ax.errorbar(np.arange(u.size),soc['mu'],yerr=soc['se'],color='k',fmt='none',capsize=2)
for i in range(u.size):
    ax.text(i,10,str(soc['N'][i].astype(int)),color='w',ha='center',fontsize=8)
ax.set(position=[0.08,0.12,0.9,0.86],xlim=[-0.5,u.size-0.5],ylim=[0,375],xticks=np.arange(u.size),
       xticklabels=lab,ylabel='Soil organic carbon (MgC ha$^{-1}$ yr$^{-1}$)')
plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Soils\SOC_By_BGCZone_PointSample' ,'png',900)

#%% Confirm that map predictions agree with point sample

soc['map_mu']=np.zeros(u.size)
for i in range(u.size):
    ind=np.where( becz['Data']==u[i] )
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

u=np.unique(ufd['ORDER'])
lab=np.array(['' for _ in range(u.size)],dtype=object)
soc_mu=np.zeros(u.size)
for i in range(u.size):
    ind=np.where( (ufd['ORDER']==u[i]) & (ufd['TOT_C_THA']>0) )
    soc_mu[i]=np.mean(ufd['TOT_C_THA'][ind])


plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(14,8))
ax.bar(np.arange(u.size),soc_mu)
ax.set(xticks=np.arange(u.size),xticklabels=u)







