
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
from scipy import stats
import statsmodels.formula.api as smf
import fcgadgets.macgyver.utilities_general as gu
import fcgadgets.macgyver.utilities_gis as gis
import fcgadgets.bc1ha.bc1ha_utilities as u1ha
import fcexplore.psp.Processing.psp_utilities as ugp

gp=gu.SetGraphics('Manuscript')

#%% Import data

meta=u1ha.Init()
meta=u1ha.ImportLUTs(meta)
#lut_1ha=u1ha.Import_BC1ha_LUTs()
#metaGP,gpt=ugp.ImportPSPs(type='Stand')
meta,gpt=ugp.ImportPlotData(meta,type='Stand')

#%% Import Raster grids

zLCC1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCoverClass1.tif')
zBGCZ=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\BEC_BIOGEOCLIMATIC_POLY\\ZONE.tif')
zAge=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\PROJ_AGE_1.tif')

#%% Export summary by Plot Type

d={}
ind=np.where( (gpt['Plot Type']==meta['LUT']['GP']['Plot Type BC']['VRI']) )[0]
for k in gpt.keys():
    d[k]=np.round(np.nanmean(gpt[k][ind]),decimals=2)
ind=np.where( (gpt['Plot Type']==meta['LUT']['GP']['Plot Type BC']['YSM']) )[0]
for k in gpt.keys():
    d[k]=np.append(d[k],np.round(np.nanmean(gpt[k][ind]),decimals=2))
ind=np.where( (gpt['PTF CN']==1) )[0]
for k in gpt.keys():
    d[k]=np.append(d[k],np.round(np.nanmean(gpt[k][ind]),decimals=2))
df=pd.DataFrame(d,index=[0,1,2])
df.to_excel(meta['Paths']['GP']['DB'] + '\\Processed\\SummarySL_ByPlotType.xlsx')

#%% QA - compare standing biomass with standing volume whole stem

plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(7.8,7.2))
ind=np.where( (gpt['PTF CN']==1) & (gpt['Vws L t0']>0) & (gpt['Ctot L t0']>0) )[0]
ax.plot(gpt['Vws L t0'][ind],gpt['Ctot L t0'][ind],'b.',ms=4,mec='w',mfc=[0.27,0.44,0.79],markeredgewidth=0.25)
rs,txt=gu.GetRegStats(gpt['Vws L t0'][ind],gpt['Ctot L t0'][ind])
ax.plot(rs['xhat'],rs['yhat'],'-k')
ax.set(ylabel='Tree biomass (MgC ha$^{-1}$)',xlabel='Whole-stem volume (m$^{3}$ ha$^{-1}$)',xlim=[0,1800],ylim=[0,600])
ax.text(900,100,txt,fontsize=8)
plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])
gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_QA_BiomassVsVolume','png',900)

#%% Plot Stocks by BGC zone

vL=['Age VRI t0','Cbk L t0','Cbr L t0','Cf L t0','Csw L t0','Cr L t0','Cag L t0','Ctot L t0','Ctot G Surv','Ctot G Recr','Ctot Mort+Lost','Ctot Net']

u=np.unique(gpt['Ecozone BC L1'])
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
    lab[i]=ugp.lut_id2cd(meta['LUT']['GP']['Ecozone BC L1'],u[i])
    for v in vL:
        ind=np.where( (gpt['Ecozone BC L1']==u[i]) &
                     (gpt['PTF CNV']==1) &
                     (gpt['Cbk L t0']>=0) & (gpt['Cbk L t0']<2000) &
                     (gpt['Cbr L t0']>=0) & (gpt['Cbr L t0']<2000) &
                     (gpt['Cf L t0']>=0) & (gpt['Cf L t0']<2000) &
                     (gpt['Cr L t0']>=0) & (gpt['Cr L t0']<2000) &
                     (gpt['Csw L t0']>=0) & (gpt['Csw L t0']<2000) &
                     (gpt['Ctot L t0']>=0) & (gpt['Ctot L t0']<10000))[0]
        d[v]['N'][i]=ind.size
        d[v]['mu'][i]=np.nanmean(gpt[v][ind])
        d[v]['sd'][i]=np.nanstd(gpt[v][ind])
        #d[v]['se'][i]=np.nanstd(gpt[v][ind])/np.sqrt(ind[0].size)
    ind=np.where( (gpt['Ecozone BC L1']==u[i]) & (gpt['PTF CNV']==1) & (gpt['Ctot L t0']>=0) & (gpt['Ctot L t0']<10000))[0]
    data[i]=gpt['Ctot L t0'][ind]

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
cl=np.array([[0.65,0.85,0.05],[0.75,0.5,0.95],[0.6,0.85,1],[0,0.6,0],[0.75,0.5,0]])
plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,8))
ax.bar(np.arange(u.size),d['Csw L t0']['mu'],facecolor=cl[0,:],label='Stemwood')
ax.bar(np.arange(u.size),d['Cbk L t0']['mu'],bottom=d['Csw L t0']['mu'],facecolor=cl[1,:],label='Bark')
ax.bar(np.arange(u.size),d['Cbr L t0']['mu'],bottom=d['Csw L t0']['mu']+d['Cbk L t0']['mu'],facecolor=cl[2,:],label='Branches')
ax.bar(np.arange(u.size),d['Cf L t0']['mu'],bottom=d['Csw L t0']['mu']+d['Cbk L t0']['mu']+d['Cbr L t0']['mu'],facecolor=cl[3,:],label='Foliage')
ax.bar(np.arange(u.size),d['Cr L t0']['mu'],bottom=d['Csw L t0']['mu']++d['Cbk L t0']['mu']+d['Cbr L t0']['mu']+d['Cf L t0']['mu'],facecolor=cl[4,:],label='Roots')
ax.set(position=[0.08,0.065,0.9,0.92],yticks=np.arange(0,550,20),xticks=np.arange(u.size),
       xticklabels=lab,ylabel='Biomass (MgC ha$^{-1}$)',xlim=[-0.5,u.size-0.5],ylim=[0,180])
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

gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_Biomass_ByGBCZone_CNV','png',900)

#%% Plot net biomass production by BGC zone (CN)

vL=['Ctot G Surv','Ctot G Recr','Ctot Mort+Lost','Ctot Mort+Lost Harv','Ctot Net']

u=np.unique(gpt['Ecozone BC L1'])
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
    lab[i]=ugp.lut_id2cd(meta['LUT']['GP']['Ecozone BC L1'],u[i])
    for v in vL:
        ind=np.where( (gpt['Ecozone BC L1']==u[i]) &
                     (gpt['PTF CN']==1) &
                     (gpt['Ctot Net']>=-1000) & (gpt['Ctot Net']<1000) )[0]

                     #(gpt['Ctot Mort Harv']==0) &
        d[v]['N'][i]=ind.size
        d[v]['mu'][i]=np.nanmean(gpt[v][ind])
        d[v]['sd'][i]=np.nanstd(gpt[v][ind])
        d[v]['se'][i]=np.nanstd(gpt[v][ind])/np.sqrt(ind.size)

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
    ind2=np.where( (zBGCZ['Data']==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'][lab[i]]) & (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest']) )
    d['Area'][i]=ind2[0].size/1e6

wa=np.sum(d['Ctot Net']['mu']*d['Area'])/np.sum(d['Area'])

lab2=np.append(lab,'Area\nweighted')

plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,9)); cl=np.array([[0.75,0.75,0.75],[0.5,0.5,0.5],[1,0.5,0],[0.45,0.75,1],[0.6,1,0]])
ax.bar(np.arange(u.size),d['Ctot Net']['mu'],facecolor=cl[0,:],label='')
for i in range(u.size):
    ax.errorbar(i,d['Ctot Net']['mu'][i],yerr=d['Ctot Net']['se'][i],color=gp['cla'],fmt='none',capsize=2,lw=0.5)
ax.bar(u.size,wa,facecolor=cl[1,:],label='Weighted average')
ax.errorbar(u.size,wa,yerr=np.mean(d['Ctot Net']['se']),color=gp['cla'],fmt='none',capsize=2,lw=0.5)
ax.plot([-1,20],[0,0],'-k',color=gp['cla'],lw=0.5)
ax.set(position=[0.08,0.075,0.9,0.91],yticks=np.arange(-5,3.5,0.5),xticks=np.arange(u.size+1),xticklabels=lab2,
       ylabel='Net biomass production (MgC ha$^{-1}$ yr$^{-1}$)',
       xlim=[-0.5,u.size+1-0.5],ylim=[-1.5,3])
#plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])
for i in range(u.size):
    if d['Ctot Net']['mu'][i]>0:
        adj=0.08
    else:
        adj=-0.08
    ax.text(i+0.2,d['Ctot Net']['mu'][i]+adj,str(d['Ctot Net']['N'][i].astype(int)),color='k',ha='center',va='center',fontsize=7,fontweight='normal')
gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_NetGrowth_ByGBCZone_CN','png',900)

#%% Plot net biomass production by BGC zone (YSM)

vL=['Ctot G Surv','Ctot G Recr','Ctot Mort+Lost','Ctot Mort+Lost Harv','Ctot Net']

u=np.unique(gpt['Ecozone BC L1'])
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
    lab[i]=ugp.lut_id2cd(meta['LUT']['GP']['Ecozone BC L1'],u[i])
    for v in vL:
        ind=np.where( (gpt['Ecozone BC L1']==u[i]) &
                     (gpt['PTF YSM']==1) &
                     (gpt['Ctot Net']>=-1000) & (gpt['Ctot Net']<1000) )[0]

                     #(gpt['Ctot Mort Harv']==0) &
        d[v]['N'][i]=ind.size
        d[v]['mu'][i]=np.nanmean(gpt[v][ind])
        d[v]['sd'][i]=np.nanstd(gpt[v][ind])
        d[v]['se'][i]=np.nanstd(gpt[v][ind])/np.sqrt(ind.size)

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
    ind2=np.where( (zBGCZ['Data']==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'][lab[i]]) & (zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest Land']) )
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
ax.set(position=[0.08,0.075,0.9,0.91],yticks=np.arange(-5,10,0.5),xticks=np.arange(u.size+1),xticklabels=lab2,
       ylabel='Net biomass production (MgC ha$^{-1}$ yr$^{-1}$)',
       xlim=[-0.5,u.size+1-0.5],ylim=[-0.5,5.25])
#plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])
for i in range(u.size):
    if d['Ctot Net']['mu'][i]>0:
        adj=0.08
    else:
        adj=-0.08
    ax.text(i+0.2,d['Ctot Net']['mu'][i]+adj,str(d['Ctot Net']['N'][i].astype(int)),color='k',ha='center',va='center',fontsize=7,fontweight='normal')
gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_NetGrowth_ByGBCZone_YSM','png',900)

#%% Average biomass dynammics

ikp=np.where( (gpt['PTF CN']==1) & (gpt['Year t1']>0) )[0]
vL=['Year t0','Year t1','Ctot G Surv','Ctot G Recr','Ctot Mort+Lost','Ctot Net']
sts={}
for v in vL:
    sts['mu ' + v]=np.nanmean(gpt[v][ikp])
    sts['se ' + v]=np.nanstd(gpt[v][ikp])/np.sqrt(ikp.size)

print( str(np.nanpercentile(gpt['Year t0'],25)) + ' ' + str(np.nanpercentile(gpt['Year t1'],75)) )
print(sts['mu Ctot Net'])

cl=np.array([[0.75,0.75,0.75],[0.24,0.49,0.77],[0.6,1,0]])
cle=[0,0,0]#[0.05,0.2,0.45]
barw=0.6
lab=['Survivor\ngrowth','Recruitment\ngrowth','Natural\nmortality','Net\ngrowth'] #,'Harvest\nmortality'

plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(14.5,8))
ax.plot([0,6],[0,0],'-k',color=gp['cla'],lw=0.5)
ax.bar(1,sts['mu Ctot G Surv'],barw,facecolor=cl[0,:],label='Growth survivors')
ax.bar(2,sts['mu Ctot G Recr'],barw,facecolor=cl[0,:],label='Growth recruitment')
ax.bar(3,-sts['mu Ctot Mort+Lost'],barw,facecolor=cl[0,:],label='Mortality')
#ax.bar(4,-sts['mu Ctot Mort Harv'],barw,facecolor=cl[0,:],label='Mortality')
ax.bar(4,sts['mu Ctot Net'],barw,facecolor=cl[0,:],label='Mortality')
ax.errorbar(1,sts['mu Ctot G Surv'],yerr=sts['se Ctot G Surv'],color=cle,fmt='none',capsize=2,lw=0.5)
ax.errorbar(2,sts['mu Ctot G Recr'],yerr=sts['se Ctot G Recr'],color=cle,fmt='none',capsize=2,lw=0.5)
ax.errorbar(3,-sts['mu Ctot Mort+Lost'],yerr=sts['se Ctot Mort+Lost'],color=cle,fmt='none',capsize=2,lw=0.5)
#ax.errorbar(4,-sts['mu Ctot Mort Harv'],yerr=sts['se Ctot Mort Harv'],color=cle,fmt='none',capsize=2,lw=0.5)
ax.errorbar(4,sts['mu Ctot Net'],yerr=sts['se Ctot Net'],color=cle,fmt='none',capsize=2,lw=0.5)
ax.set(position=[0.14,0.12,0.84,0.86],xticks=np.arange(1,len(lab)+1),xticklabels=lab,ylabel='Carbon balance of trees (MgC ha${^-1}$ yr$^{-1}$)',xlim=[0.5,4.5],ylim=[-1.5,2])
#plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])
gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_BiomassDynamics_Mean_CN','png',900)

#%% Average biomass dynammics (totals in CO2e)

# Area of forest
ind=np.where(zLCC1['Data']==meta['LUT']['Derived']['lcc1']['Forest Land'])
A=ind[0].size

ikp=np.where( (gpt['PTF CN']==1) & (gpt['Year t1']>0) )[0]
vL=['Year t0','Year t1','Ctot G Surv','Ctot G Recr','Ctot Mort+Lost','Ctot Net']
sts={}
for v in vL:
    sts['mu ' + v]=A*np.nanmean(gpt[v][ikp])/1e6*3.667
    sts['se ' + v]=A*np.nanstd(gpt[v][ikp])/np.sqrt(ikp.size)/1e6*3.667

print( str(np.nanpercentile(gpt['Year t0'],25)) + ' ' + str(np.nanpercentile(gpt['Year t1'],75)) )
print(sts['mu Ctot Net'])
#print(sts['mu Ctot Mort Harv'])

cl=np.array([[0.75,0.75,0.75],[0.24,0.49,0.77],[0.6,1,0]])
cle=[0,0,0]#[0.05,0.2,0.45]
barw=0.6
lab=['Survivor\ngrowth','Recruitment\ngrowth','Natural\nmortality','Net\ngrowth'] #,'Harvest\nmortality'

plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(14.5,8))
ax.plot([0,6],[0,0],'-k',color=gp['cla'],lw=0.5)
ax.bar(1,sts['mu Ctot G Surv'],barw,facecolor=cl[0,:],label='Growth survivors')
ax.bar(2,sts['mu Ctot G Recr'],barw,facecolor=cl[0,:],label='Growth recruitment')
ax.bar(3,-sts['mu Ctot Mort+Lost'],barw,facecolor=cl[0,:],label='Mortality')
#ax.bar(4,-sts['mu Ctot Mort Harv'],barw,facecolor=cl[0,:],label='Mortality')
ax.bar(4,sts['mu Ctot Net'],barw,facecolor=cl[0,:],label='Mortality')
ax.errorbar(1,sts['mu Ctot G Surv'],yerr=sts['se Ctot G Surv'],color=cle,fmt='none',capsize=2,lw=0.5)
ax.errorbar(2,sts['mu Ctot G Recr'],yerr=sts['se Ctot G Recr'],color=cle,fmt='none',capsize=2,lw=0.5)
ax.errorbar(3,-sts['mu Ctot Mort+Lost'],yerr=sts['se Ctot Mort+Lost'],color=cle,fmt='none',capsize=2,lw=0.5)
#ax.errorbar(4,-sts['mu Ctot Mort Harv'],yerr=sts['se Ctot Mort Harv'],color=cle,fmt='none',capsize=2,lw=0.5)
ax.errorbar(4,sts['mu Ctot Net'],yerr=sts['se Ctot Net'],color=cle,fmt='none',capsize=2,lw=0.5)
ax.set(position=[0.14,0.12,0.84,0.86],xticks=np.arange(1,len(lab)+1),xticklabels=lab,ylabel='Carbon balance of trees (MtCO$_{2}$e yr$^{-1}$)',xlim=[0.5,4.5],ylim=[-350,350])
#plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])
gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_BiomassDynamics_Sum_CN','png',900)

#%% Net Volume Production

ikp=np.where( (gpt['PTF CN']==1) & (gpt['Year t1']>0) )[0]

vL=['Year t0','Year t1','Vws G Surv','Vws G Recr','Vws Mort','Vws Net'] #,'Vws Mort Harv'
sts={}
for v in vL:
    sts['mu ' + v]=np.nanmean(gpt[v][ikp])
    sts['se ' + v]=np.nanstd(gpt[v][ikp])/np.sqrt(ikp.size)

print( str(np.nanpercentile(gpt['Year t0'],25)) + ' ' + str(np.nanpercentile(gpt['Year t1'],75)) )
print(sts['mu Vws Net'])
#print(sts['mu Vws Mort Harv'])

cl=np.array([[0.8,0.8,0.8],[0.6,1,0]])
cle=[0,0,0]
barw=0.5
lab=['Survivor\ngrowth','Recruitment\ngrowth','Mortality','Net\ngrowth'] # ,'Harvest\nmortality'

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
ax.set(position=[0.14,0.12,0.84,0.86],xticks=np.arange(1,5),xticklabels=lab,ylabel='Whole stem volume chagne (m$^{3}$ ha$^{-1}$ yr$^{-1}$)',xlim=[0.5,4.5],ylim=[-3,5])
#plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])
gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_VolumeDynamics_Mean_CN','png',900)


#%% Age distributions

x=np.arange(0,401,1)
y=zAge['Data'][0::5,0::5].flatten()
y=y[y>=0]
kde=stats.gaussian_kde(y)
p1=kde(x)

plt.close('all'); lw=1.25
fig,ax=plt.subplots(1,figsize=gu.cm2inch(14.5,7))
#plt.plot(x,p1/np.sum(p1)*100,'b-',lw=lw,color=[0.27,0.49,0.77],label='VRI R1 Comp Polygons')

ind=np.where( (gpt['Plot Type']==meta['LUT']['GP']['Plot Type BC']['VRI']) & (gpt['Age Mean t0']>=0) )[0]
kde=stats.gaussian_kde(gpt['Age Mean t0'][ind])
p=kde(x)
plt.plot(x,p/np.sum(p)*100,'r-',lw=lw,color=[0.27,0.49,0.77],label='VRI ground plots')

ind=np.where( (gpt['PTF CN']==1) & (gpt['Age Mean t0']>=0) )[0]
kde=stats.gaussian_kde(gpt['Age Mean t0'][ind])
p=kde(x)
plt.plot(x,p/np.sum(p)*100,'r--',lw=lw,color=[1,0.5,0],label='CMI + NFI cores (mean)')

kde=stats.gaussian_kde(gpt['Age Min t0'][ind])
p=kde(x)
plt.plot(x,p/np.sum(p)*100,'r-.',lw=lw,color=[0.5,1,0],label='CMI + NFI cores (min)')

kde=stats.gaussian_kde(gpt['Age Max t0'][ind])
p=kde(x)
plt.plot(x,p/np.sum(p)*100,'r:',lw=lw,color=[0,0.4,0],label='CMI + NFI cores (max)')

# ind=np.where( (gpt['PTF YSM']==1) & (gpt['Age Mean t0']>=0) )[0]
# kde=stats.gaussian_kde(gpt['Age Mean t0'][ind])
# p=kde(x)
# plt.plot(x,p,'g--')

ax.set(ylabel='Frequency (%)',xlabel='Stand age, years',xlim=[0,400],ylim=[0,0.85])
plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])
plt.tight_layout()
gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_AgeDistribution','png',900)
