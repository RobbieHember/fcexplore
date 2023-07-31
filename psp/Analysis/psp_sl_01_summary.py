
"""
PSP - SUMMARY
"""

#%% Import modules
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import copy
from scipy import stats
import statsmodels.formula.api as smf
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.bc1ha.bc1ha_util as u1ha
import fcexplore.psp.Processing.psp_util as ugp
gp=gu.SetGraphics('Manuscript')

#%% Import data
meta=u1ha.Init()
meta,gpt=ugp.ImportPlotData(meta,type='Stand')

#%% Import Raster grids
vList=['tdc','refg','lcc1_c','bgcz','harv_yr_con1','fire_yr','ibm_yr','PROJ_AGE_1'] #'lcc1_c','gfcly','gfcly_filt',
roi={'points':{'x':gpt['X'],'y':gpt['Y']}}
z=u1ha.Import_Raster(meta,roi,vList)

zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
iGrd=gis.GetGridIndexToPoints(zRef,gpt['X'],gpt['Y'])

#%% Save age for input to model
flg=0
if flg==1:
    z1=zRef.copy()
    z1['Data']=90*np.ones(zRef['Data'].shape,dtype='int16')
    z1['Data'][iGrd]=np.maximum(gpt['Age VRI t0'],gpt['Age Med t0'])+(2022-gpt['Year t0']).astype('int16')
    #plt.hist(z1['Data'][iGrd])
    gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\Age\\CustomAge_GroundPlots.tif')

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

#%% Sample frequency per BGC zone

vList2=['lcc1_c','bgcz']
z2=u1ha.Import_Raster(meta,[],vList2)

u=np.unique(gpt['Ecozone BC L1']); u=u[u>0]
lab=np.array(['' for _ in range(u.size)],dtype=object)
d={'Area':np.zeros(lab.size),'N First':np.zeros(lab.size),'N Rem':np.zeros(lab.size)}
for iU in range(u.size):
    lab[iU]=ugp.lut_id2cd(meta['LUT']['GP']['Ecozone BC L1'],u[iU])
    ind0=np.where( (z2['bgcz']['Data']==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'][lab[iU]]) & (z2['lcc1_c']['Data']==meta['LUT']['Derived']['lcc1']['Forest']) )
    ind1=np.where( (gpt['PTF CN']==1) & (z['bgcz']==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'][lab[iU]]) & (z['lcc1_c']==meta['LUT']['Derived']['lcc1']['Forest']) & (gpt['Year t0']>0) )
    ind2=np.where( (gpt['PTF CN']==1) & (z['bgcz']==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'][lab[iU]]) & (z['lcc1_c']==meta['LUT']['Derived']['lcc1']['Forest']) & (gpt['Year t0']>0) & (gpt['Year t1']>=0) )
    d['Area'][iU]=ind0[0].size
    d['N First'][iU]=ind1[0].size
    d['N Rem'][iU]=ind2[0].size

plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(14,5)); barw=0.35;
ax.bar(np.arange(u.size)-barw/2,d['N First']/(d['Area']/1e6),barw,fc=[0.7,0.8,1],ec=None,label='Single measurement')
ax.bar(np.arange(u.size)+barw/2,d['N Rem']/(d['Area']/1e6),barw,fc=[0.8,1,0.7],ec=None,label='With remeasurement')
ax.set(xticks=np.arange(0,len(lab)),xticklabels=lab,ylabel='Sampling frequency (plots Mha$^{-1}$)',xlim=[-0.5,u.size-1+0.5],ylim=[0,50])
for iU in range(u.size):
    ax.text(iU-barw/2,d['N First'][iU]/(d['Area'][iU]/1e6)+1.5,str(int(d['N First'][iU])),color=0.5*np.array([0.7,0.8,1]),fontsize=6,ha='center')
    ax.text(iU+barw/2,d['N Rem'][iU]/(d['Area'][iU]/1e6)+1.5,str(int(d['N Rem'][iU])),color=0.5*np.array([0.8,1,0.7]),fontsize=6,ha='center')
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])
plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
plt.tight_layout()
gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_SamplingFrequencyByBGC_CN','png',900)

#%% QA - Biomass vs. volume whole stem

plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(7.8,7))
ax.plot([0,2000],[0,2000],'-',lw=2,color=[0.8,0.8,0.8])
ind=np.where( (gpt['PTF CN']==1) & (gpt['Vws L t0']>0) & (gpt['Ctot L t0']>0) )[0]
ax.plot(gpt['Vws L t0'][ind],gpt['Ctot L t0'][ind],'b.',ms=4,mec='w',mfc=[0.27,0.44,0.79],markeredgewidth=0.25)
rs,txt=gu.GetRegStats(gpt['Vws L t0'][ind],gpt['Ctot L t0'][ind])
ax.plot(rs['xhat'],rs['yhat'],'-k')
ax.set(ylabel='Tree biomass (tC ha$^{-1}$)',xlabel='Whole-stem volume (m$^{3}$ ha$^{-1}$)',xlim=[0,1700],ylim=[0,1700])
ax.text(925,480,txt,fontsize=8)
ax.text(1200,1200,'1:1',fontsize=8,ha='center')
plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])
gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_QA_BiomassVsVolume','png',900)

#%% Age class distribution

x=np.arange(0,401,1)
y=z['PROJ_AGE_1']
y=y[y>=0]
kde=stats.gaussian_kde(y)
p1=kde(x)

plt.close('all'); lw=1.25
fig,ax=plt.subplots(1,figsize=gu.cm2inch(14,6))
#plt.plot(x,p1/np.sum(p1)*100,'b-',lw=lw,color=[0.27,0.49,0.77],label='VRI R1 Comp Polygons')
ind=np.where( (gpt['Plot Type']==meta['LUT']['GP']['Plot Type BC']['VRI']) & (gpt['Age Mean t0']>=0) )[0]
kde=stats.gaussian_kde(gpt['Age Mean t0'][ind])
p=kde(x)
plt.plot(x,p/np.sum(p)*100,'k-',lw=2.25,color=[0.75,0.75,0.75],label='VRI ground plots')
ind=np.where( (gpt['PTF CN']==1) & (gpt['Age Mean t0']>=0) )[0]
kde=stats.gaussian_kde(gpt['Age Mean t0'][ind])
p=kde(x)
plt.plot(x,p/np.sum(p)*100,'k-',lw=0.75,color=[0,0,0],label='CMI + NFI cores (mean)')
kde=stats.gaussian_kde(gpt['Age Min t0'][ind])
p=kde(x)
plt.plot(x,p/np.sum(p)*100,'k-.',lw=0.75,color=[0,0,0],label='CMI + NFI cores (min)')
kde=stats.gaussian_kde(gpt['Age Max t0'][ind])
p=kde(x)
plt.plot(x,p/np.sum(p)*100,'k:',lw=0.75,color=[0,0,0],label='CMI + NFI cores (max)')
#ind=np.where( (gpt['PTF YSM']==1) & (gpt['Age Mean t0']>=0) )[0]
#kde=stats.gaussian_kde(gpt['Age Mean t0'][ind])
#p=kde(x)
#plt.plot(x,p/np.sum(p)*100,'g--')
ax.set(ylabel='Frequency (%)',xlabel='Stand age, years',xlim=[0,400],ylim=[0,0.85])
plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])
plt.tight_layout()
gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_AgeDistribution','png',900)

#%% Stats by BGC zone

vL=['Ctot L t0','Cbk L t0','Cbr L t0','Cf L t0','Cr L t0','Csw L t0','Ctot G Surv','Ctot G Recr','Ctot Mort Nat','Ctot Mort Harv','Ctot Net',
    'Vws G Surv','Vws G Recr','Vws Mort Nat','Vws Mort Harv','Vws Net','Area']

uZ=np.unique(gpt['Ecozone BC L1'])
uZ=uZ[uZ>0]
lab=np.array(['' for _ in range(uZ.size)],dtype=object)

d={'CN':{},'CNV':{},'YSM':{}}
for k in d.keys():
    for v in vL:
        d[k][v]={'N':np.zeros(uZ.size),'mu':np.zeros(uZ.size),'sd':np.zeros(uZ.size),'se':np.zeros(uZ.size),'sum':np.zeros(uZ.size)}
    for iU in range(u.size):
        lab[iU]=ugp.lut_id2cd(meta['LUT']['GP']['Ecozone BC L1'],uZ[iU])
        for v in vL:
            if np.isin(v,['Ctot L t0','Cbk L t0','Cbr L t0','Cf L t0','Cr L t0','Csw L t0'])==True:
                ind=np.where( (gpt['Ecozone BC L1']==uZ[iU]) & (gpt['PTF ' + k]==1) &
                         (gpt['Cbk L t0']>=0) & (gpt['Cbk L t0']<2000) &
                         (gpt['Cbr L t0']>=0) & (gpt['Cbr L t0']<2000) &
                         (gpt['Cf L t0']>=0) & (gpt['Cf L t0']<2000) &
                         (gpt['Cr L t0']>=0) & (gpt['Cr L t0']<2000) &
                         (gpt['Csw L t0']>=0) & (gpt['Csw L t0']<2000) &
                         (gpt['Ctot L t0']>=0) & (gpt['Ctot L t0']<2000))[0]
            else:
                ind=np.where( (gpt['Ecozone BC L1']==uZ[iU]) & (gpt['PTF ' + k]==1) & (gpt['Ctot Net']>=-1000) & (gpt['Ctot Net']<1000) )[0]
            d[k][v]['N'][iU]=ind.size
            d[k][v]['mu'][iU]=np.nanmean(gpt[v][ind])
            d[k][v]['sd'][iU]=np.nanstd(gpt[v][ind])
            d[k][v]['se'][iU]=np.nanstd(gpt[v][ind])/np.sqrt(ind.size)
        ind2=np.where( (z['bgcz']==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'][lab[iU]]) & (z['lcc1_c']==meta['LUT']['Derived']['lcc1']['Forest']) )
        d[k]['Area']['sum'][iU]=ind2[0].size/1e6

#%% Biomass by BGC zone

typ='CNV'

# Put in order
d[typ]['Ctot L t0 2']={}
d[typ]['Ctot L t0 2']['mu']=d[typ]['Cbk L t0']['mu']+d[typ]['Cbr L t0']['mu']+d[typ]['Cf L t0']['mu']+d[typ]['Cr L t0']['mu']+d[typ]['Csw L t0']['mu']
ord=np.argsort(d[typ]['Ctot L t0 2']['mu'])
uo=u[ord]
lab=np.flip(lab[ord])
for v in d[typ]:
    for k in d[typ][v].keys():
        d[typ][v][k]=np.flip(d[typ][v][k][ord])

lab2=np.append(lab,'Area\nweighted')

cl=np.array([[0.65,0.85,0.05],[0.75,0.5,0.95],[0.6,0.85,1],[0,0.6,0],[0.85,0.65,0.15]])
plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,7))
ax.bar(np.arange(uZ.size),d[typ]['Csw L t0']['mu'],facecolor=cl[0,:])
ax.bar(np.arange(uZ.size),d[typ]['Cbk L t0']['mu'],bottom=d[typ]['Csw L t0']['mu'],facecolor=cl[1,:])
ax.bar(np.arange(uZ.size),d[typ]['Cbr L t0']['mu'],bottom=d[typ]['Csw L t0']['mu']+d[typ]['Cbk L t0']['mu'],facecolor=cl[2,:])
ax.bar(np.arange(uZ.size),d[typ]['Cf L t0']['mu'],bottom=d[typ]['Csw L t0']['mu']+d[typ]['Cbk L t0']['mu']+d[typ]['Cbr L t0']['mu'],facecolor=cl[3,:])
ax.bar(np.arange(uZ.size),d[typ]['Cr L t0']['mu'],bottom=d[typ]['Csw L t0']['mu']+d[typ]['Cbk L t0']['mu']+d[typ]['Cbr L t0']['mu']+d[typ]['Cf L t0']['mu'],facecolor=cl[4,:])

dW={}
for k in d[typ].keys():
    if k=='Area':
        continue
    ikp=np.where(np.isnan(d[typ][k]['mu'])==False)[0]
    dW[k]=np.sum(d[typ][k]['mu'][ikp]*d[typ]['Area']['sum'][ikp])/np.sum(d[typ]['Area']['sum'][ikp])
totW=dW['Csw L t0']+dW['Cbk L t0']+dW['Cbr L t0']+dW['Cf L t0']+dW['Cr L t0']
ax.bar(u.size,dW['Csw L t0'],facecolor=cl[0,:],label='Stemwood (' + str(int(dW['Csw L t0']/totW*100)) + '%)')
ax.bar(u.size,dW['Cbk L t0'],bottom=dW['Csw L t0'],facecolor=cl[1,:],label='Bark (' + str(int(dW['Cbk L t0']/totW*100)) + '%)')
ax.bar(u.size,dW['Cbr L t0'],bottom=dW['Csw L t0']+dW['Cbk L t0'],facecolor=cl[2,:],label='Branches (' + str(int(dW['Cbr L t0']/totW*100)) + '%)')
ax.bar(u.size,dW['Cf L t0'],bottom=dW['Csw L t0']+dW['Cbk L t0']+dW['Cbr L t0'],facecolor=cl[3,:],label='Foliage (' + str(int(dW['Cf L t0']/totW*100)) + '%)')
ax.bar(u.size,dW['Cr L t0'],bottom=dW['Csw L t0']+dW['Cbk L t0']+dW['Cbr L t0']+dW['Cf L t0'],facecolor=cl[4,:],label='Roots (' + str(int(dW['Cr L t0']/totW*100)) + '%)')

ax.set(yticks=np.arange(0,550,20),xticks=np.arange(uZ.size+1),xticklabels=lab2,ylabel='Biomass (tC ha$^{-1}$)',xlim=[-0.5,uZ.size+1-0.5],ylim=[0,180])
plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])

#for iU in range(u.size):
#    ax.text(i,6,str(d[typ]['Csw L t0']['N'][iU].astype(int)),color='k',ha='center',fontsize=7,fontweight='normal')
gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_Biomass_ByGBCZone_CNV','png',900)

#%% Net growth by BGC zone (CN)

typ='CN'
d2=copy.deepcopy(d[typ])

# Remove classes with inadequate data
ind=np.where(d2['Ctot Net']['N']>=3)[0]
for v in vL:
    for k in d2[v].keys():
        d2[v][k]=d2[v][k][ind]
uZ2=uZ[ind]
lab2=lab[ind]

# Put in order
ord=np.argsort(d2['Ctot Net']['mu'])
uZo=uZ2[ord]
lab2=np.flip(lab2[ord])
for v in d2:
    for k in d2[v].keys():
        d2[v][k]=np.flip(d2[v][k][ord])

wa=np.sum(d2['Ctot Net']['mu']*d2['Area']['sum'])/np.sum(d2['Area']['sum'])
lab2=np.append(lab2,'Area\nweighted')

plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(10,6)); cl=np.array([[0.75,0.75,0.75],[0.5,0.5,0.5],[1,0.5,0],[0.45,0.75,1],[0.6,1,0]])
ax.bar(np.arange(uZ2.size),d2['Ctot Net']['mu'],facecolor=cl[0,:],label='')
for i in range(uZ2.size):
    ax.errorbar(i,d2['Ctot Net']['mu'][i],yerr=d2['Ctot Net']['se'][i],color=gp['cla'],fmt='none',capsize=2,lw=0.5)
ax.bar(uZ2.size,wa,facecolor=cl[1,:],label='Weighted average')
ax.errorbar(uZ2.size,wa,yerr=np.mean(d2['Ctot Net']['se']),color=gp['cla'],fmt='none',capsize=2,lw=0.5)
ax.plot([-1,20],[0,0],'-k',color=gp['cla'],lw=0.5)
ax.set(yticks=np.arange(-5,3.5,0.5),xticks=np.arange(uZ2.size+1),xticklabels=lab2,ylabel='Net growth (tC ha$^{-1}$ yr$^{-1}$)',xlim=[-0.5,uZ2.size+1-0.5],ylim=[-1.5,3])
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])

for i in range(uZ2.size):
    if d2['Ctot Net']['mu'][i]>0:
        adj=0.1
    else:
        adj=-0.1
    ax.text(i+0.24,d2['Ctot Net']['mu'][i]+adj,str(d2['Ctot Net']['N'][i].astype(int)),color='k',ha='center',va='center',fontsize=5,fontweight='normal')
plt.tight_layout()
gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_NetGrowth_ByGBCZone_CN','png',900)

#%% Net growth by BGC zone (YSM)

typ='YSM'
d2=copy.deepcopy(d[typ])

# Remove classes with inadequate data
ind=np.where(d2['Ctot Net']['N']>=3)[0]
for v in vL:
    for k in d2[v].keys():
        d2[v][k]=d2[v][k][ind]
uZ2=uZ[ind]
lab2=lab[ind]

# Put in order
ord=np.argsort(d2['Ctot Net']['mu'])
uZo=uZ2[ord]
lab2=np.flip(lab2[ord])
for v in d2:
    for k in d2[v].keys():
        d2[v][k]=np.flip(d2[v][k][ord])

wa=np.sum(d2['Ctot Net']['mu']*d2['Area']['sum'])/np.sum(d2['Area']['sum'])
lab2=np.append(lab2,'Area\nweighted')

plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(10,6)); cl=np.array([[0.75,0.75,0.75],[0.5,0.5,0.5],[1,0.5,0],[0.45,0.75,1],[0.6,1,0]])
ax.bar(np.arange(uZ2.size),d2['Ctot Net']['mu'],facecolor=cl[0,:],label='')
for i in range(uZ2.size):
    ax.errorbar(i,d2['Ctot Net']['mu'][i],yerr=d2['Ctot Net']['se'][i],color=gp['cla'],fmt='none',capsize=2,lw=0.5)
ax.bar(uZ2.size,wa,facecolor=cl[1,:],label='Weighted average')
ax.errorbar(uZ2.size,wa,yerr=np.mean(d2['Ctot Net']['se']),color=gp['cla'],fmt='none',capsize=2,lw=0.5)
ax.plot([-1,20],[0,0],'-k',color=gp['cla'],lw=0.5)
ax.set(yticks=np.arange(-5,5.5,0.5),xticks=np.arange(uZ2.size+1),xticklabels=lab2,ylabel='Net growth (tC ha$^{-1}$ yr$^{-1}$)',xlim=[-0.5,uZ2.size+1-0.5],ylim=[0,5.25])
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])

for i in range(uZ2.size):
    if d2['Ctot Net']['mu'][i]>0:
        adj=0.1
    else:
        adj=-0.1
    ax.text(i+0.24,d2['Ctot Net']['mu'][i]+adj,str(d2['Ctot Net']['N'][i].astype(int)),color='k',ha='center',va='center',fontsize=5,fontweight='normal')
plt.tight_layout()
gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_NetGrowth_ByGBCZone_YSM','png',900)

#%% Average biomass dynammics

typ='CN'
d2=copy.deepcopy(d[typ])

mu={}
se={}
for v in d2.keys():
    mu[v]=np.nansum(d2[v]['mu']*d2['Area']['sum'])/np.sum(d2['Area']['sum'])
    se[v]=np.nansum(d2[v]['se']*d2['Area']['sum'])/np.sum(d2['Area']['sum'])

cl=np.array([[0.75,0.75,0.75],[0.24,0.49,0.77],[0.6,1,0]])
cle=[0,0,0]#[0.05,0.2,0.45]
barw=0.6
lab=['Survivor\ngrowth','Recruitment\ngrowth','Natural\nmortality','Harvest\nmortality','Net\ngrowth'] #

plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(7.8,6)); yt=np.arange(-1.5,2.5,0.5)
ax.plot([0,6],[0,0],'-k',color=gp['cla'],lw=0.5)
ax.bar(1,mu['Ctot G Surv'],barw,facecolor=cl[0,:],label='Growth survivors')
ax.bar(2,mu['Ctot G Recr'],barw,facecolor=cl[0,:],label='Growth recruitment')
ax.bar(3,-mu['Ctot Mort Nat'],barw,facecolor=cl[0,:],label='Natural\nmortality')
ax.bar(4,-mu['Ctot Mort Harv'],barw,facecolor=cl[0,:],label='Harvest')
ax.bar(5,mu['Ctot Net'],barw,facecolor=cl[0,:],label='Net')
ax.errorbar(1,mu['Ctot G Surv'],yerr=se['Ctot G Surv'],color=cle,fmt='none',capsize=2,lw=0.5)
ax.errorbar(2,mu['Ctot G Recr'],yerr=se['Ctot G Recr'],color=cle,fmt='none',capsize=2,lw=0.5)
ax.errorbar(3,-mu['Ctot Mort Nat'],yerr=se['Ctot Mort Nat'],color=cle,fmt='none',capsize=2,lw=0.5)
ax.errorbar(4,-mu['Ctot Mort Harv'],yerr=se['Ctot Mort Harv'],color=cle,fmt='none',capsize=2,lw=0.5)
ax.errorbar(5,mu['Ctot Net'],yerr=se['Ctot Net'],color=cle,fmt='none',capsize=2,lw=0.5)
ax.set(xticks=np.arange(1,len(lab)+1),xticklabels=lab,ylabel='Carbon balance of trees (tC ha$^{-1}$ yr$^{-1}$)',xlim=[0.5,5.5],yticks=yt,ylim=[-1.5,2])
#plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])

ax2=ax.twinx()
ax2.set_ylabel('Carbon balance of trees (MtCO$_{2}$e yr$^{-1}$)')
A=62000000 # Consistent with LC2 layer from BC LSCS 2023
yt2=A*yt/1e6*3.667; #yt2.astype('int16')
ax2.set(yticks=np.arange(-300,500,100),ylim=[yt2[0],yt2[-1]])
plt.tight_layout()
gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_BiomassDynamics_Mean_CN','png',900)

#%% Average volume dynammics

typ='CN'
d2=copy.deepcopy(d[typ])

mu={}
se={}
for v in d2.keys():
    mu[v]=np.nansum(d2[v]['mu']*d2['Area']['sum'])/np.sum(d2['Area']['sum'])
    se[v]=np.nansum(d2[v]['se']*d2['Area']['sum'])/np.sum(d2['Area']['sum'])

cl=np.array([[0.75,0.75,0.75],[0.24,0.49,0.77],[0.6,1,0]])
cle=[0,0,0]#[0.05,0.2,0.45]
barw=0.6
lab=['Survivor\ngrowth','Recruitment\ngrowth','Natural\nmortality','Harvest\nmortality','Net\ngrowth'] #

plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(7.8,6)); yt=np.arange(-1.5,2.5,0.5)
ax.plot([0,6],[0,0],'-k',color=gp['cla'],lw=0.5)
ax.bar(1,mu['Vws G Surv'],barw,facecolor=cl[0,:],label='Growth survivors')
ax.bar(2,mu['Vws G Recr'],barw,facecolor=cl[0,:],label='Growth recruitment')
ax.bar(3,-mu['Vws Mort Nat'],barw,facecolor=cl[0,:],label='Natural\nmortality')
ax.bar(4,-mu['Vws Mort Harv'],barw,facecolor=cl[0,:],label='Harvest')
ax.bar(5,mu['Vws Net'],barw,facecolor=cl[0,:],label='Net')
ax.errorbar(1,mu['Vws G Surv'],yerr=se['Vws G Surv'],color=cle,fmt='none',capsize=2,lw=0.5)
ax.errorbar(2,mu['Vws G Recr'],yerr=se['Vws G Recr'],color=cle,fmt='none',capsize=2,lw=0.5)
ax.errorbar(3,-mu['Vws Mort Nat'],yerr=se['Vws Mort Nat'],color=cle,fmt='none',capsize=2,lw=0.5)
ax.errorbar(4,-mu['Vws Mort Harv'],yerr=se['Vws Mort Harv'],color=cle,fmt='none',capsize=2,lw=0.5)
ax.errorbar(5,mu['Vws Net'],yerr=se['Vws Net'],color=cle,fmt='none',capsize=2,lw=0.5)
ax.set(xticks=np.arange(1,len(lab)+1),xticklabels=lab,ylabel='Whole stem volume periodic increment (m$^{3}$ ha$^{-1}$ yr$^{-1}$)',xlim=[0.5,5.5],ylim=[-4,6])
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])
plt.tight_layout()
gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_VolumDynamics_Mean_CN','png',900)

#%% Age responses

bw=25; bin=np.arange(bw,250+bw,bw)
lw=0.5; ms=3; mec=[0.27,0.49,0.77]; mfc=[1,1,1]; cl=[0.27,0.49,0.77];
plt.close('all'); fig,ax1=plt.subplots(1,2,figsize=gu.cm2inch(15,6))
# Coast
ind=np.where( (gpt['PTF CN']==1) & (gpt['Ecozone BC L1']==meta['LUT']['GP']['Ecozone BC L1']['CWH']) | (gpt['PTF CN']==1) & (gpt['Ecozone BC L1']==meta['LUT']['GP']['Ecozone BC L1']['MH']) )[0]
x=gpt['Age Med t0'][ind]
y=gpt['Ctot L t0'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
ax1[0].plot(bin,mu,'-ko',ms=ms,lw=lw,mew=lw,color='k',mfc='w',mec='k',label='Biomass',zorder=1)

ind=np.where( (gpt['PTF CNV']==1) & (gpt['Ecozone BC L1']==meta['LUT']['GP']['Ecozone BC L1']['CWH']) | (gpt['PTF CN']==1) & (gpt['Ecozone BC L1']==meta['LUT']['GP']['Ecozone BC L1']['MH']) )[0]
x=gpt['Age Med t0'][ind]
y=gpt['Ctot L t0'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
ax1[0].plot(bin,mu,'--ks',ms=ms,lw=lw,mew=lw,color=[0.5,0.5,0.5],mfc='w',mec=[0.5,0.5,0.5],label='Biomass',zorder=1)
ax1[0].set(ylabel='Biomass (tC ha$^-$$^1$)',xlabel='Age, years',xticks=np.arange(0,400,bw),yticks=np.arange(0,300,20),xlim=[0,250+bw],ylim=[0,240])
ax1[0].yaxis.set_ticks_position('both'); ax1[0].xaxis.set_ticks_position('both'); ax1[0].tick_params(length=gp['tickl'])
ax2=ax1[0].twinx()
y=gpt['Ctot Net'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
#yhat=np.interp(xhat,bin[np.isnan(mu)==False],mu[np.isnan(mu)==False])
ax2.bar(bin,mu,0.9*bw,fc=[0.75,0.75,0.75],ec='none',label='Net growth',zorder=-1)
ax2.plot([0,500],[0,0],'-k',lw=lw)
ax2.set(ylabel='Net growth (tC ha$^-$$^1$ yr$^-$$^1$)',xlabel='',ylim=[-1.75,5.5])
ax1[0].set_zorder(ax2.get_zorder()+1)
ax1[0].patch.set_visible(False)
ax2.tick_params(length=gp['tickl'])

# Interior
#ind=np.where( (gpt['PTF CN']==1) & (gpt['Ecozone BC L1']!=meta['LUT']['GP']['Ecozone BC L1']['CWH']) & (gpt['Ecozone BC L1']!=meta['LUT']['GP']['Ecozone BC L1']['MH']) & (gpt['Ecozone BC L1']!=meta['LUT']['GP']['Ecozone BC L1']['ICH']) )[0]
ind=np.where( (gpt['PTF CN']==1) & (gpt['Ecozone BC L1']!=meta['LUT']['GP']['Ecozone BC L1']['CWH']) & (gpt['Ecozone BC L1']!=meta['LUT']['GP']['Ecozone BC L1']['MH']) )[0]
x=gpt['Age Med t0'][ind]
y=gpt['Ctot L t0'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
ax1[1].plot(bin,mu,'-ko',ms=ms,lw=lw,mew=lw,color='k',mfc='w',mec='k',label='CMI+NFI',zorder=1)

ind=np.where( (gpt['PTF CNV']==1) & (gpt['Ecozone BC L1']!=meta['LUT']['GP']['Ecozone BC L1']['CWH']) & (gpt['Ecozone BC L1']!=meta['LUT']['GP']['Ecozone BC L1']['MH']) )[0]
x=gpt['Age Med t0'][ind]
y=gpt['Ctot L t0'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
ax1[1].plot(bin,mu,'--ks',ms=ms,lw=lw,mew=lw,color=[0.5,0.5,0.5],mfc='w',mec=[0.5,0.5,0.5],label='CMI+NFI+VRI',zorder=1)
ax1[1].set(ylabel='Biomass (tC ha$^-$$^1$)',xlabel='Age, years',xticks=np.arange(0,400,bw),yticks=np.arange(0,300,20),xlim=[0,250+bw],ylim=[0,240])
ax1[1].yaxis.set_ticks_position('both'); ax1[1].xaxis.set_ticks_position('both'); ax1[1].tick_params(length=gp['tickl'])
ax3=ax1[1].twinx()
y=gpt['Ctot Net'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)

ax3.bar(bin,mu,0.9*bw,fc=[0.75,0.75,0.75],ec='none',label='Net growth CMI+NFI',zorder=-1)
ax3.plot([0,500],[0,0],'-k',lw=lw)
ax3.set(ylabel='Net growth (tC ha$^-$$^1$ yr$^-$$^1$)',xlabel='',ylim=[-1.75,5.5])
ax1[1].set_zorder(ax3.get_zorder()+1)
ax1[1].patch.set_visible(False)
ax3.tick_params(length=gp['tickl'])
leg=ax1[1].legend(loc='upper right',frameon=False,facecolor=None,edgecolor='w')
gu.axletters(ax1,plt,0.04,0.93,FontColor=meta['Graphics']['gp']['cla'],LetterStyle='Caps',FontWeight='Bold')
plt.tight_layout()
gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_AgeResponsesByRegion','png',900)

#%% Profiles - import data

#bw=5; bin=np.arange(-12.5,12.5,bw)
bw=5; bin=np.arange(-42.5,142.5,bw)

# Fire
fy=z['fire_yr']
tsf=gpt['Year t0']-fy
ikp=np.where(gpt['Year IBM']==0)[0]
prfl_wf={}
for k in gpt.keys():
    prfl_wf[k]={}
    prfl_wf[k]['N'],prfl_wf[k]['mu'],prfl_wf[k]['med'],prfl_wf[k]['sig'],prfl_wf[k]['se']=gu.discres(tsf[ikp],gpt[k][ikp],bw,bin)

# IBM
tso=-100*np.ones(gpt['Year t0'].size)
for iY in range(10):
    zY=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\PEST_SEVERITY_CODE_IBM_' + str(iY+1) + '_Year.tif')['Data'][iGrd]
    zS=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\PEST_SEVERITY_CODE_IBM_' + str(iY+1) + '_Severity.tif')['Data'][iGrd]
    ind=np.where( (zS==meta['LUT']['PEST_INFESTATION_POLY']['PEST_SEVERITY_CODE']['S']) & (gpt['Year t0']-zY>tso) )[0]
    tso[ind]=gpt['Year t0'][ind]-zY[ind]
ikp=np.where(gpt['Year Wildfire']<1995)[0]
prfl_ibm={}
for k in gpt.keys():
    prfl_ibm[k]={}
    prfl_ibm[k]['N'],prfl_ibm[k]['mu'],prfl_ibm[k]['med'],prfl_ibm[k]['sig'],prfl_ibm[k]['se']=gu.discres(tso[ikp],gpt[k][ikp],bw,bin)

# Harvest
hy=z['harv_yr_con1']
tsh=gpt['Year t0']-hy
ikp=np.where(gpt['Year Wildfire']<2200)[0]
prfl_h={}
for k in gpt.keys():
    prfl_h[k]={}
    prfl_h[k]['N'],prfl_h[k]['mu'],prfl_h[k]['med'],prfl_h[k]['sig'],prfl_h[k]['se']=gu.discres(tsh[ikp],gpt[k][ikp],bw,bin)

#%% Plot profiles

it0=np.where(bin<0)[0]; it1=np.where( (bin>0) )[0]

plt.close('all');
fig,ax=plt.subplots(3,3,figsize=gu.cm2inch(15.5,12)); ms=2
v='Ctot L t0'
ax[0,0].plot([0,0],[0,300],'k-',lw=2,color=[0.8,0.8,0.8])
ax[0,0].plot(bin,prfl_wf[v]['mu'],'ko',ms=ms,mfc='k',mec='k',label='Live')
ax[0,0].errorbar(bin,prfl_wf[v]['mu'],yerr=prfl_wf[v]['se'],color='k',fmt='none',capsize=1.5,lw=0.5)
mu0=np.nanmean(prfl_wf[v]['mu'][it0])
mu1=np.nanmean(prfl_wf[v]['mu'][it1])
ax[0,0].plot(bin[it0],mu0*np.ones(it0.size),'k--',color='k',lw=0.5)
ax[0,0].plot(bin[it1],mu1*np.ones(it1.size),'k-.',color='k',lw=0.5)
#ax[0,0].text(np.mean(bin[it1]),1.82*mu1,'Change = ' + str(np.round(mu1-mu0,decimals=1)) + '\n(' + str(np.round( (mu1-mu0)/mu0*100 ,decimals=1)) + '%)',ha='center',fontsize=7,style='normal',color='k')
ax[0,0].set(ylabel='Biomass (tC ha$^1$)',xlabel='Years since wildfire',ylim=[0,240],xlim=[-40,80])
ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=gp['tickl'])

ax[0,1].plot([0,0],[0,300],'k-',lw=2,color=[0.8,0.8,0.8])
ax[0,1].plot(bin,prfl_ibm[v]['mu'],'ko',ms=ms,mfc='k',mec='k',label='Live')
ax[0,1].errorbar(bin,prfl_ibm[v]['mu'],yerr=prfl_ibm[v]['se'],color='k',fmt='none',capsize=1.5,lw=0.5)
mu0=np.nanmean(prfl_ibm[v]['mu'][it0])
mu1=np.nanmean(prfl_ibm[v]['mu'][it1])
ax[0,1].plot(bin[it0],mu0*np.ones(it0.size),'k--',color='k',lw=0.5)
ax[0,1].plot(bin[it1],mu1*np.ones(it1.size),'k-.',color='k',lw=0.5)
#ax[0,1].text(np.mean(bin[it1]),1.3*mu1,'Change = ' + str(np.round(mu1-mu0,decimals=1)) + '\n(' + str(np.round( (mu1-mu0)/mu0*100 ,decimals=1)) + '%)',ha='center',fontsize=7,style='normal',color='k')
ax[0,1].set(ylabel='Biomass (tC ha$^1$)',xlabel='Years since severe outbreak',ylim=[0,240],xlim=[-40,80])
#ax[0].legend(loc='upper right',frameon=False,facecolor='w',edgecolor='w')
ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=gp['tickl'])

ax[0,2].plot([0,0],[0,300],'k-',lw=2,color=[0.8,0.8,0.8])
ax[0,2].plot(bin,prfl_h[v]['mu'],'ko',ms=ms,mfc='k',mec='k',label='Live')
ax[0,2].errorbar(bin,prfl_h[v]['mu'],yerr=prfl_h[v]['se'],color='k',fmt='none',capsize=1.5,lw=0.5)
mu0=np.nanmean(prfl_h[v]['mu'][it0])
mu1=np.nanmean(prfl_h[v]['mu'][it1])
ax[0,2].plot(bin[it0],mu0*np.ones(it0.size),'k--',color='k',lw=0.5)
ax[0,2].plot(bin[it1],mu1*np.ones(it1.size),'k-.',color='k',lw=0.5)
#ax[0,2].text(np.mean(bin[it1]),1.3*mu1,'Change = ' + str(np.round(mu1-mu0,decimals=1)) + '\n(' + str(np.round( (mu1-mu0)/mu0*100 ,decimals=1)) + '%)',ha='center',fontsize=7,style='normal',color='k')
ax[0,2].set(ylabel='Biomass (tC ha$^1$)',xlabel='Years since harvest',ylim=[0,240],xlim=[-40,80])
#ax[0].legend(loc='upper right',frameon=False,facecolor='w',edgecolor='w')
ax[0,2].yaxis.set_ticks_position('both'); ax[0,2].xaxis.set_ticks_position('both'); ax[0,2].tick_params(length=gp['tickl'])

v='Ctot D t0'
ax[1,0].plot([0,0],[0,200],'k-',lw=2,color=[0.8,0.8,0.8])
ax[1,0].plot(bin,prfl_wf[v]['mu'],'ko',ms=ms,mfc='k',mec='k',label='Live')
ax[1,0].errorbar(bin,prfl_wf[v]['mu'],yerr=prfl_wf[v]['se'],color='k',fmt='none',capsize=1.5,lw=0.5)
mu0=np.nanmean(prfl_wf[v]['mu'][it0])
mu1=np.nanmean(prfl_wf[v]['mu'][it1])
ax[1,0].plot(bin[it0],mu0*np.ones(it0.size),'k--',color='k',lw=0.5)
ax[1,0].plot(bin[it1],mu1*np.ones(it1.size),'k-.',color='k',lw=0.5)
#ax[1,0].text(np.mean(bin[it1]),1.55*mu1,'Change = ' + str(np.round(mu1-mu0,decimals=1)) + '\n(' + str(np.round( (mu1-mu0)/mu0*100 ,decimals=1)) + '%)',ha='center',fontsize=7,style='normal',color='k')
ax[1,0].set(ylabel='Dead wood (tC ha$^1$)',xlabel='Years since wildfire',ylim=[0,40],xlim=[-40,80])
ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=gp['tickl'])

ax[1,1].plot([0,0],[0,200],'k-',lw=2,color=[0.8,0.8,0.8])
ax[1,1].plot(bin,prfl_ibm[v]['mu'],'ko',ms=ms,mfc='k',mec='k',label='Live')
ax[1,1].errorbar(bin,prfl_ibm[v]['mu'],yerr=prfl_ibm[v]['se'],color='k',fmt='none',capsize=1.5,lw=0.5)
mu0=np.nanmean(prfl_ibm[v]['mu'][it0])
mu1=np.nanmean(prfl_ibm[v]['mu'][it1])
ax[1,1].plot(bin[it0],mu0*np.ones(it0.size),'k--',color='k',lw=0.5)
ax[1,1].plot(bin[it1],mu1*np.ones(it1.size),'k-.',color='k',lw=0.5)
#ax[1,1].text(np.mean(bin[it1]),0.7*mu1,'Change = ' + str(np.round(mu1-mu0,decimals=1)) + '\n(' + str(np.round( (mu1-mu0)/mu0*100 ,decimals=1)) + '%)',ha='center',fontsize=7,style='normal',color='k')
ax[1,1].set(ylabel='Dead wood (tC ha$^1$)',xlabel='Years since severe outbreak',ylim=[0,40],xlim=[-40,80])
#ax[0].legend(loc='upper right',frameon=False,facecolor='w',edgecolor='w')
ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=gp['tickl'])

ax[1,2].plot([0,0],[0,200],'k-',lw=2,color=[0.8,0.8,0.8])
ax[1,2].plot(bin,prfl_h[v]['mu'],'ko',ms=ms,mfc='k',mec='k',label='Live')
ax[1,2].errorbar(bin,prfl_h[v]['mu'],yerr=prfl_h[v]['se'],color='k',fmt='none',capsize=1.5,lw=0.5)
mu0=np.nanmean(prfl_h[v]['mu'][it0])
mu1=np.nanmean(prfl_h[v]['mu'][it1])
ax[1,2].plot(bin[it0],mu0*np.ones(it0.size),'k--',color='k',lw=0.5)
ax[1,2].plot(bin[it1],mu1*np.ones(it1.size),'k-.',color='k',lw=0.5)
#ax[1,2].text(np.mean(bin[it1]),0.7*mu1,'Change = ' + str(np.round(mu1-mu0,decimals=1)) + '\n(' + str(np.round( (mu1-mu0)/mu0*100 ,decimals=1)) + '%)',ha='center',fontsize=7,style='normal',color='k')
ax[1,2].set(ylabel='Dead wood (tC ha$^1$)',xlabel='Years since harvest',ylim=[0,40],xlim=[-40,80])
ax[1,2].yaxis.set_ticks_position('both'); ax[1,2].xaxis.set_ticks_position('both'); ax[1,2].tick_params(length=gp['tickl'])

v='Ctot Net'
ax[2,0].plot([0,0],[-200,200],'k-',lw=2,color=[0.8,0.8,0.8])
ax[2,0].plot(bin,prfl_wf[v]['mu'],'ko',ms=ms,mfc='k',mec='k',label='Live')
ax[2,0].errorbar(bin,prfl_wf[v]['mu'],yerr=prfl_wf[v]['se'],color='k',fmt='none',capsize=1.5,lw=0.5)
mu0=np.nanmean(prfl_wf[v]['mu'][it0])
mu1=np.nanmean(prfl_wf[v]['mu'][it1])
ax[2,0].plot(bin[it0],mu0*np.ones(it0.size),'k--',color='k',lw=0.5)
ax[2,0].plot(bin[it1],mu1*np.ones(it1.size),'k-.',color='k',lw=0.5)
#ax[2,0].text(np.mean(bin[it1]),1.55*mu1,'Change = ' + str(np.round(mu1-mu0,decimals=1)) + '\n(' + str(np.round( (mu1-mu0)/mu0*100 ,decimals=1)) + '%)',ha='center',fontsize=7,style='normal',color='k')
ax[2,0].set(ylabel='Net growth (tC ha$^{-1}$ yr$^{-1}$)',xlabel='Years since wildfire',ylim=[-2,4],xlim=[-40,80])
ax[2,0].yaxis.set_ticks_position('both'); ax[2,0].xaxis.set_ticks_position('both'); ax[2,0].tick_params(length=gp['tickl'])

ax[2,1].plot([0,0],[-200,200],'k-',lw=2,color=[0.8,0.8,0.8])
ax[2,1].plot(bin,prfl_ibm[v]['mu'],'ko',ms=ms,mfc='k',mec='k',label='Live')
ax[2,1].errorbar(bin,prfl_ibm[v]['mu'],yerr=prfl_ibm[v]['se'],color='k',fmt='none',capsize=1.5,lw=0.5)
mu0=np.nanmean(prfl_ibm[v]['mu'][it0])
mu1=np.nanmean(prfl_ibm[v]['mu'][it1])
ax[2,1].plot(bin[it0],mu0*np.ones(it0.size),'k--',color='k',lw=0.5)
ax[2,1].plot(bin[it1],mu1*np.ones(it1.size),'k-.',color='k',lw=0.5)
#ax[2,1].text(np.mean(bin[it1]),1.55*mu1,'Change = ' + str(np.round(mu1-mu0,decimals=1)) + '\n(' + str(np.round( (mu1-mu0)/mu0*100 ,decimals=1)) + '%)',ha='center',fontsize=7,style='normal',color='k')
ax[2,1].set(ylabel='Net growth (tC ha$^{-1}$ yr$^{-1}$)',xlabel='Years since severe outbreak',ylim=[-2,5],xlim=[-40,80])
ax[2,1].yaxis.set_ticks_position('both'); ax[2,1].xaxis.set_ticks_position('both'); ax[2,1].tick_params(length=gp['tickl'])

ax[2,2].plot([0,0],[-200,200],'k-',lw=2,color=[0.8,0.8,0.8])
ax[2,2].plot(bin,prfl_h[v]['mu'],'ko',ms=ms,mfc='k',mec='k',label='Live')
ax[2,2].errorbar(bin,prfl_h[v]['mu'],yerr=prfl_h[v]['se'],color='k',fmt='none',capsize=1.5,lw=0.5)
mu0=np.nanmean(prfl_h[v]['mu'][it0])
mu1=np.nanmean(prfl_h[v]['mu'][it1])
ax[2,2].plot(bin[it0],mu0*np.ones(it0.size),'k--',color='k',lw=0.5)
ax[2,2].plot(bin[it1],mu1*np.ones(it1.size),'k-.',color='k',lw=0.5)
#ax[2,2].text(np.mean(bin[it1]),1.55*mu1,'Change = ' + str(np.round(mu1-mu0,decimals=1)) + '\n(' + str(np.round( (mu1-mu0)/mu0*100 ,decimals=1)) + '%)',ha='center',fontsize=7,style='normal',color='k')
ax[2,2].set(ylabel='Net growth (tC ha$^{-1}$ yr$^{-1}$)',xlabel='Years since harvest',ylim=[-2,5],xlim=[-40,80])
ax[2,2].yaxis.set_ticks_position('both'); ax[2,2].xaxis.set_ticks_position('both'); ax[2,2].tick_params(length=gp['tickl'])

gu.axletters(ax,plt,0.04,0.88,FontColor=meta['Graphics']['gp']['cla'],LetterStyle='Caps',FontWeight='Bold')#,Labels=['Wildfire','MPB','Wildfire','MPB'])
plt.tight_layout()
gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\QA_Disturbance_Profiles','png',900)

#%% Mortality associations

# Estimate harvest mortality
gpt['Ctot Mort+Lost Harv']=np.zeros(gpt['PTF CN'].size)
ind=np.where( (gpt['Occ_Harv']==1) )[0]
gpt['Ctot Mort+Lost Harv'][ind]=gpt['Ctot Mort+Lost'][ind]

gpt['N Mort+Lost Harv']=np.zeros(gpt['PTF CN'].size)
ind=np.where( (gpt['Occ_Harv']==1) )[0]
gpt['N Mort+Lost Harv'][ind]=gpt['N Mort+Lost'][ind]

daL=['No damage','Insect','Disease','Plant competition','Animal browsing','Fire','Frost, snow, ice, hail','Water stress','Wind','Silviculture','Defects','Flooding, lightning, slides','Unknown','Harv']
#daL=['Insect', 'Disease', 'Plant Competition', 'Animal Browsing', 'Fire', 'Frost', 'Drought', 'Snow and Ice', 'Wind', 'Silviculture', 'Other', 'Unknown','Harv']
daL_lab=['No\ndamage\n(competition)','Insects', 'Diseases', 'Plant\ncompetition', 'Animal\nbrowsing', 'Fire', 'Frost, snow\nice and\nhail', 'Water\nstress','Wind', 'Silvi-\nculture','Defects\nand\nbreakage','Flooding\nlightning\nslides','Species\ndecline','Harvest']

# Make sure that all mortality is accounted for by damage agents, assign unaccounted mortality to no damage
ind=np.where( (gpt['PTF CN']==1) & (gpt['Ctot Mort+Lost']>0) & (gpt['Ctot Mort+Lost No damage']==0) & (gpt['Ctot Mort+Lost Insect']==0) & \
              (gpt['Ctot Mort+Lost Disease']==0) & (gpt['Ctot Mort+Lost Plant competition']==0) & (gpt['Ctot Mort+Lost Animal browsing']==0) & \
              (gpt['Ctot Mort+Lost Fire']==0) & (gpt['Ctot Mort+Lost Frost, snow, ice, hail']==0) & (gpt['Ctot Mort+Lost Water stress']==0) & \
              (gpt['Ctot Mort+Lost Silviculture']==0) & (gpt['Ctot Mort+Lost Defects']==0) & (gpt['Ctot Mort+Lost Flooding, lightning, slides']==0) & \
              (gpt['Ctot Mort+Lost Unknown']==0) & (gpt['Ctot Mort+Lost Harv']==0) )[0]
gpt['Ctot Mort+Lost No damage'][ind]=gpt['Ctot Mort+Lost'][ind]

# Make sure that all mortality is accounted for by damage agents, assign unaccounted mortality to no damage
ind=np.where( (gpt['PTF CN']==1) & (gpt['N Mort+Lost']>0) & (gpt['N Mort+Lost No damage']==0) & (gpt['N Mort+Lost Insect']==0) & \
              (gpt['N Mort+Lost Disease']==0) & (gpt['N Mort+Lost Plant competition']==0) & (gpt['N Mort+Lost Animal browsing']==0) & \
              (gpt['N Mort+Lost Fire']==0) & (gpt['N Mort+Lost Frost, snow, ice, hail']==0) & (gpt['N Mort+Lost Water stress']==0) & \
              (gpt['N Mort+Lost Silviculture']==0) & (gpt['N Mort+Lost Defects']==0) & (gpt['N Mort+Lost Flooding, lightning, slides']==0) & \
              (gpt['N Mort+Lost Unknown']==0) & (gpt['N Mort+Lost Harv']==0) )[0]
gpt['N Mort+Lost No damage'][ind]=gpt['N Mort+Lost'][ind]

d={}
ikp=np.where( (gpt['PTF CN']==1) )[0]
d['Mean All']=np.nanmean(gpt['Ctot Mort+Lost'][ikp])
d['Sum All']=np.nansum(gpt['Ctot Mort+Lost'][ikp])
d['Mean']=np.zeros(len(daL))
d['Sum']=np.zeros(len(daL))
d['N Mean All']=np.nanmean(gpt['N Mort+Lost'][ikp])
d['N Sum All']=np.nansum(gpt['N Mort+Lost'][ikp])
d['N Mean']=np.zeros(len(daL))
d['N Sum']=np.zeros(len(daL))
for i in range(len(daL)):
    d['Mean'][i]=np.nanmean(np.nan_to_num(gpt['Ctot Mort+Lost ' + daL[i]][ikp]))
    d['Sum'][i]=np.nansum(np.nan_to_num(gpt['Ctot Mort+Lost ' + daL[i]][ikp]))
    d['N Mean'][i]=np.nanmean(np.nan_to_num(gpt['N Mort+Lost ' + daL[i]][ikp]))
    d['N Sum'][i]=np.nansum(np.nan_to_num(gpt['N Mort+Lost ' + daL[i]][ikp]))
print(np.sum(d['Sum'])/d['Sum All'])

mup=d['Mean']/np.sum(d['Mean'])*100
mup2=d['N Mean']/np.sum(d['N Mean'])*100
#mup=d['Sum']#/np.sum(d['Mean'])*100
ord=np.flip(np.argsort(mup))
mup=mup[ord]
mup2=mup2[ord]
lab=np.array(daL_lab)[ord]

plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(15.5,5.5)); barw=0.35
ax.bar(np.arange(len(daL))-barw/1.85,mup,barw,fc=[0.7,0.9,0.6],ec=None,label='Gravimetric (biomass loss due to mortality)')
ax.bar(np.arange(len(daL))+barw/1.85,mup2,barw,fc=[0.8,0.85,1],ec=None,label='Demographic (number of trees)')
ax.set(xticks=np.arange(0,len(lab)),xticklabels=lab,ylabel='Frequency (%)',xlim=[-0.5,len(lab)-1+0.5])# ,ylim=[0,35]
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])
plt.legend(loc='upper right',frameon=False,facecolor=[1,1,1],labelspacing=0.25)
plt.tight_layout()
gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_MortalityAssociations_CN','png',900)

#%% Age response of mortality

# Age=gpt['Age Med t0']
# ind=np.where( (np.isnan(gpt['Age Med t0'])==True) | (gpt['Age Med t0']==0) )[0]
# Age[ind]=gpt['Age VRI t0'][ind]

# bwA=25; binA=np.arange(bwA,250+bwA,bwA)
# mu=np.zeros((binA.size,len(daL)))
# for iD in range(len(daL)):
#     for iA in range(binA.size):
#         ind=np.where( (gpt['PTF CN']==1) & (np.abs(Age-binA[iA])<=bwA/2) )[0]
#         mu[iA,iD]=np.nanmean(np.nan_to_num(gpt['Ctot Mort+Lost ' + daL[iD]][ind]))

# ord=np.flip(np.argsort(np.nanmean(mu,axis=0)))
# mu2=mu[:,ord]
# lab=np.array(daL)[ord]

# ind=np.where(np.nanmean(mu2,axis=0)>0)[0]
# mu2=mu2[:,ind]
# lab=lab[ind]

# plt.close('all'); fig,ax=plt.subplots(5,2,figsize=gu.cm2inch(15,12)); cnt=0
# for r in range(5):
#     for c in range(2):
#         ax[r,c].bar(binA,mu[:,cnt],0.75*bwA,fc=[0.75,0.75,0.75],ec='none')
#         #ax[r,c].plot([0,500],[0,0],'-k',lw=0.25)
#         ax[r,c].set(ylabel='Mortality\n(tC ha$^-$$^1$ yr$^-$$^1$)',xlabel='Age, years',xticks=np.arange(0,400,bwA),xlim=[0,250+bwA],ylim=[0,1.05*np.max(mu[:,cnt])])
#         ax[r,c].yaxis.set_ticks_position('both'); ax[r,c].xaxis.set_ticks_position('both'); ax[r,c].tick_params(length=meta['Graphics']['gp']['tickl'])
#         cnt=cnt+1
# plt.tight_layout()
# gu.axletters(ax,plt,0.03,0.73,FontColor=meta['Graphics']['gp']['cla'],LetterStyle='Caps',FontWeight='Bold',Labels=lab)
# # #gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_AgeResponsesByRegion','png',900)


#%% Initial tree density

u=np.unique(gpt['Ecozone BC L1']); u=u[u>0]
lab=np.array(['' for _ in range(u.size)],dtype=object)
d={'N L t0':np.zeros(lab.size)}
for iU in range(u.size):
    lab[iU]=ugp.lut_id2cd(meta['LUT']['GP']['Ecozone BC L1'],u[iU])
    ikp=np.where( (gpt['PTF CNV']==1) & (z['bgcz']==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'][lab[iU]]) & (gpt['Year t0']>0) & (gpt['Age Mean t0']>20) & (gpt['Age Mean t0']<60) )
    d['N L t0'][iU]=np.nanmean(gpt['N L t0'][ikp])
# Remove nans
ind=np.where(np.isnan(d['N L t0'])==False)[0]
lab=lab[ind]
u=u[ind]
d['N L t0']=d['N L t0'][ind]
# Put in order
ord=np.flip(np.argsort(d['N L t0']))
lab=lab[ord]
d['N L t0']=d['N L t0'][ord]

cl=np.array([[0.75,0.75,0.75]])
plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,5.5))
ax.bar(np.arange(u.size),d['N L t0'],facecolor=cl[0,:],label='Stemwood')
ax.set(position=[0.09,0.15,0.9,0.8],yticks=np.arange(0,3200,500),xticks=np.arange(u.size),
       xticklabels=lab,ylabel='Tree density (stems ha$^{-1}$)',xlabel='Biogeoclimatic zone',xlim=[-0.5,u.size-0.5],ylim=[0,3600])
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])
gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_Initial Tree Density By GBC Zone','png',900)

#%% Biomass by tree density class

d={}
for k in gpt.keys():
    d[k]=np.zeros(3)
    cnt=0
    for c in meta['LUT']['Derived']['tdc'].keys():
        ind=np.where( (gpt['PTF CN']==1) & (z['tdc']==meta['LUT']['Derived']['tdc'][c]) )[0]
        d[k][cnt]=np.nanmean(gpt[k][ind])
        cnt=cnt+1
d['Ctot L t0']
d['Ctot Mort Harv']

# #%% Modelling

# #fAg=inline(['max(0,b(1).*(1+((b(2).*(x(:,1)./b(3)).^b(4)-1)./(exp(x(:,1)./b(3))))))'],'b','x'); % Weibull
# #fAm=inline(['max(0,b(1).*(1+((b(2).*((x-b(5))./b(3)).^b(4)-1)./(exp((x-b(5))./b(3))))))'],'b','x');

# def fun(x,a,b,c,d):
#     yhat=d+a*((1-np.exp(-b*x))**(1/(1-c))) # From Dzierzon and Mason 2006
#     return yhat
# xhat=np.arange(0,100,1)

# indM=np.where( (tl['DBH']>0) & (np.isnan(tl['H'])==True) | (tl['DBH']>0) & (tl['H']<=0) )[0]
# indM.size/tl['DBH'].size

# indG=np.where( (tl['DBH']>0) & (tl['H']>0) )[0]
# indG.size/tl['DBH'].size

# # Global model
# ikp=np.where( (tl['DBH']>0) & (tl['H']>0) & (tl['Vital Status']==meta['LUT']['Vital Status']['Live']) & (tl['ID DA1']==meta['LUT']['Damage Agents']['None']) & (tl['Stature']==meta['LUT']['Stature']['Standing']) )[0]
# x=tl['DBH'][ikp]
# y=tl['H'][ikp]
# popt_glob0=[26,0.1,0.66,2]
# popt_glob,pcov=curve_fit(fun,x,y,popt_glob0)
# #yhat=fun(xhat,popt_glob[0],popt_glob[1],popt_glob[2],popt_glob[3])
# rs_glob,txt=gu.GetRegStats(x,y)



# #%% Aridity

# #ind=np.where( (gpt['PTF CN']==1) & (gpt['Ecozone BC L1']==meta['LUT']['GP']['Ecozone CA L1']['Pacific_Maritime']) )[0]
# ind=np.where( (gpt['PTF CN']==1) )[0]
# #ind=np.where( (gpt['PTF CN']==1) & (gpt['Ecozone BC L1']==meta['LUT']['GP']['Ecozone BC L1']['SBPS']) )[0]

# x=gpt['Climate Class'][ind]
# y=gpt['Ctot G Surv'][ind]

# bw=1; bin=np.arange(1,6,bw)
# N,mu,med,sig,se=gu.discres(x,y,bw,bin)

# plt.close('all')
# plt.plot(bin,mu,'bo')

