
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
import fcgadgets.macgyver.utilities_inventory as invu
import fcexplore.psp.Processing.psp_utilities as utl_gp

#%% Set figure properties

gp=gu.SetGraphics('Manuscript')

#%% Import data

# Ground plots
meta={}
meta['Paths']={}
meta['Paths']['DB']=r'C:\Users\rhember\Documents\Data\GroundPlots\PSP-NADB2'
meta['Paths']['Figs']=r'C:\Users\rhember\OneDrive - Government of BC\Figures\Biomass'
meta=utl_gp.ImportParameters(meta)
d=gu.ipickle(meta['Paths']['DB'] + '\\Processed\\L2\\L2_BC.pkl')
sl=d['sobs'].copy()
del d

#%% Import Raster grids

becz=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\becz.tif')
lutBGC=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\becz_lut.xlsx')
lc2=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\lc2.tif')

#%% Remove plots with bad coordinates

ind=np.where( (sl['Lat']>30) & (sl['Lon']!=0) )[0]
for k in sl.keys():
    sl[k]=sl[k][ind]

#%% Plot type filter (CMI, NFI, and VRI)

sl['Flag Plot Types To Keep']=np.zeros(sl['ID Plot'].size)
ind=np.where( (sl['Plot Type']==meta['LUT']['Plot Type BC']['CMI']) |
             (sl['Plot Type']==meta['LUT']['Plot Type BC']['NFI']) |
             (sl['Plot Type']==meta['LUT']['Plot Type BC']['VRI']) )[0]
#ind=np.where( (sl['Plot Type']==meta['LUT']['Plot Type BC']['CMI']) | (sl['Plot Type']==meta['LUT']['Plot Type BC']['NFI']) )[0]
sl['Flag Plot Types To Keep'][ind]=1
np.sum(sl['Flag Plot Types To Keep'])

sl['Flag Keep CMI+NFI']=np.zeros(sl['ID Plot'].size)
ind=np.where( (sl['Plot Type']==meta['LUT']['Plot Type BC']['CMI']) | (sl['Plot Type']==meta['LUT']['Plot Type BC']['NFI']) )[0]
sl['Flag Keep CMI+NFI'][ind]=1

#%% Export tabular summary

# Summary (Just PSPs)
ind=np.where( (sl['Flag Keep CMI+NFI']==1) & (sl['N L t1']>=0) )[0]
d={}
for k in sl.keys():
    d[k]=np.round(np.nanmean(sl[k][ind]),decimals=2)
df=pd.DataFrame(d,index=[0])
df.to_excel(r'C:\Users\rhember\Documents\Data\GroundPlots\PSP-NADB2\Processed\SummaryPSPs.xlsx')

# Summary by Plot Type
d={}
ind=np.where( (sl['Plot Type']==meta['LUT']['Plot Type BC']['VRI']) )[0]
for k in sl.keys():
    d[k]=np.round(np.nanmean(sl[k][ind]),decimals=2)
ind=np.where( (sl['Plot Type']==meta['LUT']['Plot Type BC']['YSM']) )[0]
for k in sl.keys():
    d[k]=np.append(d[k],np.round(np.nanmean(sl[k][ind]),decimals=2))
ind=np.where( (sl['Plot Type']==meta['LUT']['Plot Type BC']['CMI']) | (sl['Plot Type']==meta['LUT']['Plot Type BC']['NFI']) )[0]
for k in sl.keys():
    d[k]=np.append(d[k],np.round(np.nanmean(sl[k][ind]),decimals=2))
df=pd.DataFrame(d,index=[0,1,2])
df.to_excel(r'C:\Users\rhember\Documents\Data\GroundPlots\PSP-NADB2\Processed\SummarySL_ByPlotType.xlsx')

#%% Assess whether samples are representative of population crown closure

plt.close('all')

ikp=np.where(sl['Flag Plot Types To Keep']==1)[0]
kd=gu.ksdensity(sl['Crown Closure (VRI)'][ikp])
plt.plot(kd,'b-')
print(np.mean(sl['Crown Closure (VRI)'][ikp]))

ikp=np.where(sl['Flag Keep CMI+NFI']==1)[0]
kd=gu.ksdensity(sl['Crown Closure (VRI)'][ikp])
plt.plot(kd,'g-')
print(np.mean(sl['Crown Closure (VRI)'][ikp]))

zLC2=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\lc2.tif')
zCC=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\crownc.tif')
ikp=np.where( (zLC2['Data'][0::10,0::10].flatten()==4) )[0]
kd=gu.ksdensity(zCC['Data'][0::10,0::10].flatten()[ikp])
plt.plot(kd,'c--')
print(np.mean(zCC['Data'][0::10,0::10].flatten()[ikp]))

#%% Plot Stocks by BGC zone

vL=['Age t0','Cbk L t0','Cbr L t0','Cf L t0','Csw L t0','Cr L t0','Cag L t0','Ctot L t0','Ctot G Surv','Ctot G Recr','Ctot Mort','Ctot Mort Harv','Ctot Net']

u=np.unique(sl['Ecozone BC L1'])
u=u[u>0]
lab=np.array(['' for _ in range(u.size)],dtype=object)

d={}
for v in vL:
    d[v]={}
    d[v]['N']=np.zeros(u.size)
    d[v]['mu']=np.zeros(u.size)
    d[v]['sd']=np.zeros(u.size)
    d[v]['se']=np.zeros(u.size)
data=[None]*u.size
for i in range(u.size):
    lab[i]=utl_gp.lut_id2cd(meta,'Ecozone BC L1',u[i])
    for v in vL:
        ind=np.where( (sl['Ecozone BC L1']==u[i]) &
                     (sl['Flag Plot Types To Keep']==1) &
                     (sl['Cbk L t0']>=0) & (sl['Cbk L t0']<2000) &
                     (sl['Cbr L t0']>=0) & (sl['Cbr L t0']<2000) &
                     (sl['Cf L t0']>=0) & (sl['Cf L t0']<2000) &
                     (sl['Cr L t0']>=0) & (sl['Cr L t0']<2000) &
                     (sl['Csw L t0']>=0) & (sl['Csw L t0']<2000) &
                     (sl['Ctot L t0']>=0) & (sl['Ctot L t0']<10000))[0]
        d[v]['N'][i]=ind.size
        d[v]['mu'][i]=np.nanmean(sl[v][ind])
        d[v]['sd'][i]=np.nanstd(sl[v][ind])
        #d[v]['se'][i]=np.nanstd(sl[v][ind])/np.sqrt(ind[0].size)
    ind=np.where( (sl['Ecozone BC L1']==u[i]) & (sl['Flag Plot Types To Keep']==1) & (sl['Ctot L t0']>=0) & (sl['Ctot L t0']<10000))[0]
    data[i]=sl['Ctot L t0'][ind]

# Put in order
d['Ctot L t0']['mu']=d['Cbk L t0']['mu']+d['Cbr L t0']['mu']+d['Cf L t0']['mu']+d['Cr L t0']['mu']+d['Csw L t0']['mu']
ord=np.argsort(d['Ctot L t0']['mu'])
uo=u[ord]
lab=np.flip(lab[ord])
for v in d:
    for k in d[v].keys():
        d[v][k]=np.flip(d[v][k][ord])
data2=[None]*len(data)
for i in range(len(data)):
    data2[i]=data[np.flip(ord)[i]]

# Plot
cl=np.array([[1,0.75,0.55],[0.25,0.14,0.05],[1,0.5,0],[0.45,0.75,1],[0.6,1,0]])
plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,9))
ax.bar(np.arange(u.size),d['Csw L t0']['mu'],facecolor=cl[0,:],label='Stemwood')
ax.bar(np.arange(u.size),d['Cbk L t0']['mu'],bottom=d['Csw L t0']['mu'],facecolor=cl[2,:],label='Bark')
ax.bar(np.arange(u.size),d['Cbr L t0']['mu'],bottom=d['Csw L t0']['mu']+d['Cbk L t0']['mu'],facecolor=cl[3,:],label='Branches')
ax.bar(np.arange(u.size),d['Cf L t0']['mu'],bottom=d['Csw L t0']['mu']+d['Cbk L t0']['mu']+d['Cbr L t0']['mu'],facecolor=cl[4,:],label='Foliage')
ax.bar(np.arange(u.size),d['Cr L t0']['mu'],bottom=d['Csw L t0']['mu']++d['Cbk L t0']['mu']+d['Cbr L t0']['mu']+d['Cf L t0']['mu'],facecolor=cl[1,:],label='Roots')
ax.set(position=[0.08,0.065,0.9,0.92],yticks=np.arange(0,550,20),xticks=np.arange(u.size),
       xticklabels=lab,ylabel='Biomass (MgC ha$^{-1}$ yr$^{-1}$)',xlim=[-0.5,u.size-0.5],ylim=[0,200])
plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])

flg=0
if flg==1:
    vio=ax.violinplot(data2,np.arange(u.size),widths=0.7,showmeans=False,showextrema=True,showmedians=True)
    for pc in vio['bodies']:
        pc.set_facecolor('none')
        pc.set_edgecolor([0.25,0.25,0.25])
        pc.set_linewidth(0.5)
        #pc.set_alpha(1)
    for partname in ('cbars','cmins','cmaxes','cmedians'):
        vp=vio[partname]
        vp.set_edgecolor([0.75,0.75,0.75])
        vp.set_linewidth(0.5)
        vp.set_alpha(0.75)

    vp=vio['cbars']
    vp.set_alpha(0)

    vp=vio['cmedians']
    vp.set_edgecolor([0.75,0.75,0.75])
    vp.set_linewidth(2)
    #vp.set_alpha(0.75)

for i in range(u.size):
    ax.text(i,6,str(d['Csw L t0']['N'][i].astype(int)),color='k',ha='center',fontsize=7,fontweight='normal')

#gu.PrintFig(meta['Paths']['Figs'] + '\\BiomassFromPlots_ByGBCZone','png',900)

#%% Plot fluxes by BGC zone

vL=['Ctot G Surv','Ctot G Recr','Ctot Mort','Ctot Mort Harv','Ctot Net']

u=np.unique(sl['Ecozone BC L1'])
u=u[u>0]
lab=np.array(['' for _ in range(u.size)],dtype=object)

d={}
for v in vL:
    d[v]={}
    d[v]['N']=np.zeros(u.size)
    d[v]['mu']=np.zeros(u.size)
    d[v]['sd']=np.zeros(u.size)
    d[v]['se']=np.zeros(u.size)

for i in range(u.size):
    lab[i]=utl_gp.lut_id2cd(meta,'Ecozone BC L1',u[i])
    for v in vL:
        ind=np.where( (sl['Ecozone BC L1']==u[i]) &
                     (sl['Flag Keep CMI+NFI']==1) &
                     (sl['Ctot Mort Harv']==0) &
                     (sl['Ctot Net']>=-1000) & (sl['Ctot Net']<1000) )[0]
        d[v]['N'][i]=ind.size
        d[v]['mu'][i]=np.nanmean(sl[v][ind])
        d[v]['sd'][i]=np.nanstd(sl[v][ind])
        d[v]['se'][i]=np.nanstd(sl[v][ind])/np.sqrt(ind.size)

# Remove classes with inadequate data
ind=np.where(d['Ctot Net']['N']>=3)[0]
for v in vL:
    for k in d[v].keys():
        d[v][k]=d[v][k][ind]
u=u[ind]
lab=lab[ind]

# Put in order
ord=np.argsort(d['Ctot Net']['mu'])
uo=u[ord]
lab=np.flip(lab[ord])
for v in d:
    for k in d[v].keys():
        d[v][k]=np.flip(d[v][k][ord])

# Area weighting
d['Area']=np.zeros(lab.size)
for i in range(lab.size):
    ind1=np.where(lutBGC['ZONE']==lab[i])[0]
    ind2=np.where( (becz['Data']==lutBGC['VALUE'][ind1[0]]) & (lc2['Data']==4) )
    d['Area'][i]=ind2[0].size/1e6

wa=np.sum(d['Ctot Net']['mu']*d['Area'])/np.sum(d['Area'])

lab2=np.append(lab,'Area\nweighted')

# Plot

plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,9)); cl=np.array([[0.75,0.75,0.75],[0.5,0.5,0.5],[1,0.5,0],[0.45,0.75,1],[0.6,1,0]])
ax.bar(np.arange(u.size),d['Ctot Net']['mu'],facecolor=cl[0,:],label='')
for i in range(u.size):
    ax.errorbar(i,d['Ctot Net']['mu'][i],yerr=d['Ctot Net']['se'][i],color=gp['cla'],fmt='none',capsize=2,lw=0.5)
ax.bar(u.size,wa,facecolor=cl[1,:],label='Weighted average')
ax.errorbar(u.size,wa,yerr=np.mean(d['Ctot Net']['se']),color=gp['cla'],fmt='none',capsize=2,lw=0.5)
ax.plot([-1,20],[0,0],'-k',color=gp['cla'],lw=0.5)
ax.set(position=[0.08,0.075,0.9,0.91],yticks=np.arange(-5,5,0.2),xticks=np.arange(u.size+1),xticklabels=lab2,
       ylabel='Net biomass production (MgC ha$^{-1}$ yr$^{-1}$)',
       xlim=[-0.5,u.size+1-0.5],ylim=[-1,1.8])
#plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])

for i in range(u.size):
    if d['Ctot Net']['mu'][i]>0:
        adj=0.08
    else:
        adj=-0.08
    ax.text(i+0.2,d['Ctot Net']['mu'][i]+adj,str(d['Ctot Net']['N'][i].astype(int)),color='k',ha='center',va='center',fontsize=7,fontweight='normal')

gu.PrintFig(meta['Paths']['Figs'] + '\\NetBiomassProductionFromPlots_ByGBCZone','png',900)

#%% Net Biomass Production

#gp=gu.SetGraphics('Manuscript')
gp=gu.SetGraphics('Presentation Light')

# Import ground plot data
metaGP={}
metaGP['Paths']={}
metaGP['Paths']['DB']=r'C:\Users\rhember\Documents\Data\GroundPlots\PSP-NADB2'
metaGP=utl_gp.ImportParameters(metaGP)
d=gu.ipickle(metaGP['Paths']['DB'] + '\\Processed\\L2\\L2_BC.pkl')
sl=d['sobs']
del d

ikp=np.where( (sl['Flag Keep CMI+NFI']==1) & (sl['Year t1']>0) & (sl['Lat']>0) & (sl['Lon']!=0) )[0]

A=57000000

vL=['Year t0','Year t1','Ctot G Surv','Ctot G Recr','Ctot Mort','Ctot Mort Harv','Ctot Net']
sts={}
for v in vL:
    sts['mu ' + v]=A*np.nanmean(sl[v][ikp])/1e6*3.667
    sts['se ' + v]=A*np.nanstd(sl[v][ikp])/np.sqrt(ikp.size)/1e6*3.667

print( str(np.nanpercentile(sl['Year t0'],25)) + ' ' + str(np.nanpercentile(sl['Year t1'],75)) )
print(sts['mu Ctot Net'])
print(sts['mu Ctot Mort Harv'])

cl=np.array([[0.24,0.49,0.77],[0.6,1,0]])
cle=[0.05,0.2,0.45]
barw=0.5
lab=['Survivor\ngrowth','Recruitment\ngrowth','Natural\nmortality','Harvest\nmortality','Net\ngrowth']

plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(14.5,8))
ax.plot([0,6],[0,0],'-k',color=gp['cla'],lw=0.5)
ax.bar(1,sts['mu Ctot G Surv'],barw,facecolor=cl[0,:],label='Growth survivors')
ax.bar(2,sts['mu Ctot G Recr'],barw,facecolor=cl[0,:],label='Growth recruitment')
ax.bar(3,-sts['mu Ctot Mort']+sts['mu Ctot Mort Harv'],barw,facecolor=cl[0,:],label='Mortality')
ax.bar(4,-sts['mu Ctot Mort Harv'],barw,facecolor=cl[0,:],label='Mortality')
ax.bar(5,sts['mu Ctot Net'],barw,facecolor=cl[0,:],label='Mortality')
ax.errorbar(1,sts['mu Ctot G Surv'],yerr=sts['se Ctot G Surv'],color=cle,fmt='none',capsize=2,lw=0.5)
ax.errorbar(2,sts['mu Ctot G Recr'],yerr=sts['se Ctot G Recr'],color=cle,fmt='none',capsize=2,lw=0.5)
ax.errorbar(3,-sts['mu Ctot Mort']+sts['mu Ctot Mort Harv'],yerr=sts['se Ctot Mort'],color=cle,fmt='none',capsize=2,lw=0.5)
ax.errorbar(4,-sts['mu Ctot Mort Harv'],yerr=sts['se Ctot Mort Harv'],color=cle,fmt='none',capsize=2,lw=0.5)
ax.errorbar(5,sts['mu Ctot Net'],yerr=sts['se Ctot Net'],color=cle,fmt='none',capsize=2,lw=0.5)
ax.set(position=[0.14,0.12,0.84,0.86],xticks=np.arange(1,6),xticklabels=lab,ylabel='Carbon balance of trees (MtCO$_{2}$e yr$^{-1}$)',xlim=[0.5,5.5],ylim=[-300,300])
#plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Biomass\BC Tree Carbon Dynamics from Ground Plots','png',900)

#%% Net Volume Production

#gp=gu.SetGraphics('Manuscript')
gp=gu.SetGraphics('Presentation Light')

# Import ground plot data
metaGP={}
metaGP['Paths']={}
metaGP['Paths']['DB']=r'C:\Users\rhember\Documents\Data\GroundPlots\PSP-NADB2'
metaGP=utl_gp.ImportParameters(metaGP)
d=gu.ipickle(metaGP['Paths']['DB'] + '\\Processed\\L2\\L2_BC.pkl')
sl=d['sobs']
del d

ikp=np.where( (sl['Flag Keep CMI+NFI']==1) & (sl['Year t1']>0) & (sl['Lat']>0) & (sl['Lon']!=0) )[0]

vL=['Year t0','Year t1','Vws G Surv','Vws G Recr','Vws Mort','Vws Mort Harv','Vws Net']
sts={}
for v in vL:
    sts['mu ' + v]=np.nanmean(1000*sl[v][ikp])
    sts['se ' + v]=np.nanstd(1000*sl[v][ikp])/np.sqrt(ikp.size)

print( str(np.nanpercentile(sl['Year t0'],25)) + ' ' + str(np.nanpercentile(sl['Year t1'],75)) )
print(sts['mu Vws Net'])
print(sts['mu Vws Mort Harv'])

cl=np.array([[0.24,0.49,0.77],[0.6,1,0]])
cle=[0.05,0.2,0.45]
barw=0.5
lab=['Survivor\ngrowth','Recruitment\ngrowth','Natural\nmortality','Net\ngrowth'] # ,'Harvest\nmortality'

plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(14.5,8))
ax.plot([0,6],[0,0],'-k',color=gp['cla'],lw=0.5)
ax.bar(1,sts['mu Vws G Surv'],barw,facecolor=cl[0,:],label='Growth survivors')
ax.bar(2,sts['mu Vws G Recr'],barw,facecolor=cl[0,:],label='Growth recruitment')
ax.bar(3,-sts['mu Vws Mort'],barw,facecolor=cl[0,:],label='Mortality')
#ax.bar(4,-sts['mu Vws Mort Harv'],barw,facecolor=cl[0,:],label='Mortality')
ax.bar(4,sts['mu Vws Net'],barw,facecolor=cl[0,:],label='Mortality')
ax.errorbar(1,sts['mu Vws G Surv'],yerr=sts['se Vws G Surv'],color=cle,fmt='none',capsize=2,lw=0.5)
ax.errorbar(2,sts['mu Vws G Recr'],yerr=sts['se Vws G Recr'],color=cle,fmt='none',capsize=2,lw=0.5)
ax.errorbar(3,-sts['mu Vws Mort'],yerr=sts['se Vws Mort'],color=cle,fmt='none',capsize=2,lw=0.5)
#ax.errorbar(4,-sts['mu Vws Mort Harv'],yerr=sts['se Vws Mort Harv'],color=cle,fmt='none',capsize=2,lw=0.5)
ax.errorbar(4,sts['mu Vws Net'],yerr=sts['se Vws Net'],color=cle,fmt='none',capsize=2,lw=0.5)
ax.set(position=[0.14,0.12,0.84,0.86],xticks=np.arange(1,5),xticklabels=lab,ylabel='Whole stem volume chagne (m$^{3}$ ha$^{-1}$ yr$^{-1}$)',xlim=[0.5,4.5],ylim=[-3,3])
#plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Biomass\BC Tree Volume Dynamics from Ground Plots','png',900)


#%% Compare plot AGB with VRI polygon estimates

# Create geodatabase for plot locations
points=[]
for i in range(sl['X'].size):
    points.append( Point(sl['X'][i],sl['Y'][i]) )
d={'geometry':points}
for k in sl.keys():
    d[k]=sl[k]
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
    gdf_vri=gpd.read_file(meta['Paths']['DB'] + '\\ground_plots_vri.geojson')

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
    lab[i]=utl_gp.lut_id2cd(meta,'Ecozone BC L1',u[i])
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
