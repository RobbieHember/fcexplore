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
import fcexplore.field_plots.Processing.psp_util as ugp
import fcexplore.field_plots.Processing.psp_plot as pgp
gp=gu.SetGraphics('Manuscript')

#%% Import data
meta=u1ha.Init()
meta,gpt,soc=ugp.ImportGroundPlotData(meta,type='Stand',include_soil='True')

#%% Import Raster grids
vList=['tdc_vri15','refg','lc_comp1_2019','bgcz','harv_yr_comp1','fire_yr','ibm_yr','age_vri15','bsr_yr','bsr_sc'] #'lc_comp1_2019','gfcly','gfcly_filt',
roi={'points':{'x':gpt['X'],'y':gpt['Y']}}
z=u1ha.Import_Raster(meta,roi,vList)

zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
iGrd=gis.GetGridIndexToPoints(zRef,gpt['X'],gpt['Y'])

vList=['tdc_vri15','lc_comp1_2019','bgcz']
roi={'points':{'x':soc['x'],'y':soc['y']}}
zS=u1ha.Import_Raster(meta,roi,vList)

#%% Save age for input to model
flg=0
if flg==1:
    z1=zRef.copy()
    z1['Data']=90*np.ones(zRef['Data'].shape,dtype='int16')
    z1['Data'][iGrd]=np.maximum(gpt['Age VRI t0'],gpt['Age Med t0'])+(2022-gpt['Year t0']).astype('int16')
    #plt.hist(z1['Data'][iGrd])
    gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\Age\\CustomAge_GroundPlots.tif')

#%% Export summary by Plot Type
df=pgp.GP_ExportTableStatsByPlotType(meta,gpt,export=False)

#%% Sample frequency per BGC zone
pgp.GP_SamplingFrequencyByBGC(meta,gpt,soc)

#%% QA - Biomass vs. volume whole stem
pgp.GP_Biomass_vs_VolumeWholeStem(meta,gpt)

#%% Relationship between age and age
pgp.GP_QA_Age_vs_Age(meta,gpt)

#%% Age class distribution
pgp.GP_AgeClassDistribution(meta,gpt)

#%% Calculate statisics for BGC zone
dBGC=ugp.CalcStatsByBGCZone(meta,gpt,soc)
#dZ,dBGC['id'],dBGC['code']=ugp.StatsByBGCZone(meta,gpt,soc)

#%% Plot total ecosystem carbon stratified by BGC zone
pgp.GP_TEC_ByBGC_CN(meta,gpt,dBGC)

#%% Plot biomass by BGC zone
pgp.GP_BiomassByBGC_CNV(meta,gpt,dBGC)

#%% Plot average biomass dynammics
pgp.GP_BiomassDynamics_Mean_CN(meta,gpt,dBGC)

#%% Age responses of biomass and net growth
dAC=ugp.CalcStatsByAgeClass(meta,gpt)
pgp.GP_AgeResponsesByRegion(meta,dAC)

#%%
pgp.GP_DeadWoodByBGC_CNV(meta,gpt,dBGC)

#%% Plot net growth by BGC zone (CN and YSM side by side)
# def GP_GrowthNet_ByBGC_CN(meta,dBGC):
# 	d0=copy.deepcopy(dBGC['data'])
# 	lab0=dBGC['code'].copy()
# 	typ='CN'
# 	
# 	# Remove classes with inadequate data
# 	ikp=np.where(d0[typ]['Ctot Net']['N']>=3)[0]
# 	for v in d0[typ].keys():
# 		for k in d0[typ][v].keys():
# 			d0[typ][v][k]=d0[typ][v][k][ikp]
# 	lab2=lab0[ikp]
# 	
# 	# Put in order
# 	ord=np.argsort(d0[typ]['Ctot Net']['mu'])
# 	lab2=np.flip(lab2[ord])
# 	for v in d0[typ].keys():
# 		for k in d0[typ][v].keys():
# 			d0[typ][v][k]=np.flip(d0[typ][v][k][ord])
# 	
# 	wa=np.nansum(d0[typ]['Ctot Net']['mu']*d0[typ]['Area']['sum'])/np.nansum(d0[typ]['Area']['sum'])
# 	lab2=np.append(lab2,'Area\nweighted')
# 	
# 	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(14,7)); cl=np.array([[0.75,0.85,1],[0.17,0.39,0.64],[0.85,1,0.5],[0.45,0.75,1],[0.6,1,0]]); bw=0.4
# 	ax.bar(np.arange(lab2.size)-bw/2,d0[typ]['Ctot Net']['mu'],bw,facecolor=cl[0,:])
# 	for i in range(lab2.size):
# 		ax.errorbar(i-bw/2,d0[typ]['Ctot Net']['mu'][i],yerr=d0[typ]['Ctot Net']['se'][i],color=gp['cla'],fmt='none',capsize=2,lw=0.5)
# 	ax.bar(lab2.size-bw/2,wa,bw,facecolor=cl[0,:],label='CMI/NFI')
# 	ax.errorbar(lab2.size-bw/2,wa,yerr=np.nanmean(d0[typ]['Ctot Net']['se']),color=gp['cla'],fmt='none',capsize=2,lw=0.25)
# 	# for i in range(uZ2.size):
# 	# if d0[typ]['Ctot Net']['mu'][i]>0:
# 	# adj=0.1
# 	# else:
# 	# adj=-0.1
# 	# ax.text(i+0.24,d0[typ]['Ctot Net']['mu'][i]+adj,str(d0[typ]['Ctot Net']['N'][i].astype(int)),color='k',ha='center',va='center',fontsize=5,fontweight='normal')
# 	
# 	d0=copy.deepcopy(dBGC['data'])
# 	typ='YSM'
# 	for v in d0[typ].keys():
# 		for k in d0[typ][v].keys():
# 			d0[typ][v][k]=d0[typ][v][k][ikp]
# 			d0[typ][v][k]=np.flip(d0[typ][v][k][ord])
# 	wa=np.nansum(d0[typ]['Ctot Net']['mu']*d0[typ]['Area']['sum'])/np.nansum(d0[typ]['Area']['sum'])
# 	for i in range(uZ2.size):
# 		if np.isnan(d0[typ]['Ctot Net']['mu'][i])==True:
# 			continue
# 		ax.bar(i+bw/2,d0[typ]['Ctot Net']['mu'][i],bw,facecolor=cl[2,:],label='')
# 		ax.errorbar(i+bw/2,d0[typ]['Ctot Net']['mu'][i],yerr=d0[typ]['Ctot Net']['se'][i],color=gp['cla'],fmt='none',capsize=2,lw=0.5)
# 	ax.bar(uZ2.size+bw/2,wa,bw,facecolor=cl[2,:],label='YSM')
# 	ax.errorbar(uZ2.size+bw/2,wa,yerr=np.nanmean(d0[typ]['Ctot Net']['se']),color=gp['cla'],fmt='none',capsize=2,lw=0.25)
# 	
# 	ax.plot([-1,20],[0,0],'-k',color=gp['cla'],lw=0.5)
# 	plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
# 	ax.set(yticks=np.arange(-5,20,0.5),xticks=np.arange(uZ2.size+1),xticklabels=lab2,ylabel='Net growth (tC ha$^{-1}$ yr$^{-1}$)',xlim=[-0.5,uZ2.size+1-0.5],ylim=[-1.5,5])
# 	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
# 	plt.tight_layout()
# 	gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_NetGrowth_ByGBCZone_CN','png',900)
# 	return

#%% Plot net growth by BGC zone (YSM)
d0=copy.deepcopy(d)
lab0=lab.copy()
typ='YSM'

# Remove classes with inadequate data
ind=np.where(d0[typ]['Ctot Net']['N']>=3)[0]
for v in d0[typ].keys():
    for k in d0[typ][v].keys():
        d0[typ][v][k]=d0[typ][v][k][ind]
uZ2=uZ[ind]
lab2=lab0[ind]

# Put in order
ord=np.argsort(d0[typ]['Ctot Net']['mu'])
uZo=uZ2[ord]
lab2=np.flip(lab2[ord])
for v in d0[typ]:
    for k in d0[typ][v].keys():
        d0[typ][v][k]=np.flip(d0[typ][v][k][ord])

wa=np.sum(d0[typ]['Ctot Net']['mu']*d0[typ]['Area']['sum'])/np.sum(d0[typ]['Area']['sum'])
lab2=np.append(lab2,'Area\nweighted')

plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(14,6)); cl=np.array([[0.47,0.69,0.94],[0.17,0.39,0.64],[1,0.5,0],[0.45,0.75,1],[0.6,1,0]])
ax.bar(np.arange(uZ2.size),d0[typ]['Ctot Net']['mu'],facecolor=cl[0,:],label='')
for i in range(uZ2.size):
    ax.errorbar(i,d0[typ]['Ctot Net']['mu'][i],yerr=d0[typ]['Ctot Net']['se'][i],color=gp['cla'],fmt='none',capsize=2,lw=0.5)
ax.bar(uZ2.size,wa,facecolor=cl[1,:],label='Weighted average')
ax.errorbar(uZ2.size,wa,yerr=np.mean(d0[typ]['Ctot Net']['se']),color=gp['cla'],fmt='none',capsize=2,lw=0.5)
ax.plot([-1,20],[0,0],'-k',color=gp['cla'],lw=0.5)
ax.set(yticks=np.arange(-5,5.5,0.5),xticks=np.arange(uZ2.size+1),xticklabels=lab2,ylabel='Net growth (tC ha$^{-1}$ yr$^{-1}$)',xlim=[-0.5,uZ2.size+1-0.5],ylim=[0,5.25])
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
# for i in range(uZ2.size):
#     if d0[typ]['Ctot Net']['mu'][i]>0:
#         adj=0.1
#     else:
#         adj=-0.1
#     ax.text(i+0.24,d0[typ]['Ctot Net']['mu'][i]+adj,str(d0[typ]['Ctot Net']['N'][i].astype(int)),color='k',ha='center',va='center',fontsize=5,fontweight='normal')
plt.tight_layout()
gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_NetGrowth_ByGBCZone_YSM','png',900)

#%% Age responses of gross growth and mortality
bw=25; bin=np.arange(bw,250+bw,bw)
lw=0.5; ms=3; mec=[0.27,0.49,0.77]; mfc=[1,1,1]; cl=[0.27,0.49,0.77];
plt.close('all'); fig,ax1=plt.subplots(1,2,figsize=gu.cm2inch(15,7))
# Coast
ind=np.where( (gpt['PTF CN']==1) & (gpt['Ecozone BC L1']==meta['LUT']['GP']['Ecozone BC L1']['CWH']) | (gpt['PTF CN']==1) & (gpt['Ecozone BC L1']==meta['LUT']['GP']['Ecozone BC L1']['MH']) )[0]
x=gpt['Age Med t0'][ind]
y=gpt['Ctot G Surv'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
ax1[0].bar(bin,mu,0.9*bw,fc=[0.85,1,0.5],ec='none',label='Growth of survivors')
y=gpt['Ctot G Recr'][ind]
N2,mu2,med2,sig2,se2=gu.discres(x,y,bw,bin)
ax1[0].bar(bin,mu2,0.9*bw,fc=[0.3,0.75,0.],ec='none',label='Growth from recruitment',bottom=mu)
ax1[0].set(ylabel='Biomass flux (tC ha$^-$$^1$ yr$^-$$^1$)',xlabel='Age, years',xticks=np.arange(0,400,bw),yticks=np.arange(-20,30,2),xlim=[0,250+bw],ylim=[-3,8])
ax1[0].yaxis.set_ticks_position('both'); ax1[0].xaxis.set_ticks_position('both'); ax1[0].tick_params(length=meta['Graphics']['gp']['tickl'])
y=gpt['Ctot Mort Nat'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
ax1[0].bar(bin,-mu,0.9*bw,fc=[0.8,0.6,0.4],ec='none',label='Mortality natural')
y=gpt['Ctot Mort Harv'][ind]
N2,mu2,med2,sig2,se2=gu.discres(x,y,bw,bin)
ax1[0].bar(bin,-mu2,0.9*bw,fc=[0.4,0.3,0.2],ec='none',label='Mortality harvest',bottom=-mu)
x=gpt['Age Med t0'][ind]
y=gpt['Ctot Net'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
ax1[0].plot(bin,mu,'-ko',ms=2.5,mfc='w',mec='k',label='Net growth')
ax1[0].plot([0,500],[0,0],'k-')
ax1[0].legend(loc='upper right',frameon=False,facecolor='w',edgecolor='w')
# Interior
#'Ctot G Surv HGT10'
ind=np.where( (gpt['PTF CN']==1) & (gpt['Ecozone BC L1']!=meta['LUT']['GP']['Ecozone BC L1']['CWH']) & (gpt['Ecozone BC L1']!=meta['LUT']['GP']['Ecozone BC L1']['MH']) )[0]
x=gpt['Age Med t0'][ind]
y=gpt['Ctot G Surv'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
ax1[1].bar(bin,mu,0.9*bw,fc=[0.85,1,0.5],ec='none',label='Growth of survivors')
y=gpt['Ctot G Recr'][ind]
N2,mu2,med2,sig2,se2=gu.discres(x,y,bw,bin)
ax1[1].bar(bin,mu2,0.9*bw,fc=[0.3,0.75,0.],ec='none',label='Growth from recruitment',bottom=mu)
ax1[1].set(ylabel='Biomass flux (tC ha$^-$$^1$ yr$^-$$^1$)',xlabel='Age, years',xticks=np.arange(0,400,bw),yticks=np.arange(-20,30,2),xlim=[0,250+bw],ylim=[-3,8])
ax1[1].yaxis.set_ticks_position('both'); ax1[0].xaxis.set_ticks_position('both'); ax1[0].tick_params(length=meta['Graphics']['gp']['tickl'])
y=gpt['Ctot Mort Nat'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
ax1[1].bar(bin,-mu,0.9*bw,fc=[0.8,0.6,0.4],ec='none',label='Mortality natural')
y=gpt['Ctot Mort Harv'][ind]
N2,mu2,med2,sig2,se2=gu.discres(x,y,bw,bin)
ax1[1].bar(bin,-mu2,0.9*bw,fc=[0.4,0.3,0.2],ec='none',label='Mortality harvest',bottom=-mu)
y=gpt['Ctot Net'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
ax1[1].plot(bin,mu,'-ko',ms=2.5,mfc='w',mec='k',label='Net growth')
ax1[1].plot([0,500],[0,0],'k-')
ax1[1].yaxis.set_ticks_position('both'); ax1[1].xaxis.set_ticks_position('both'); ax1[1].tick_params(length=meta['Graphics']['gp']['tickl'])
gu.axletters(ax1,plt,0.04,0.93,FontColor=meta['Graphics']['gp']['cla'],LetterStyle='Caps',FontWeight='Bold')
plt.tight_layout()
#gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_AgeResponseGrossGrowthMortalityByRegion_CN','png',900)

#%% Age responses of gross growth and mortality
bw=25; bin=np.arange(bw,250+bw,bw)
lw=0.5; ms=3; mec=[0.27,0.49,0.77]; mfc=[1,1,1]; cl=[0.27,0.49,0.77];
plt.close('all'); fig,ax1=plt.subplots(1,2,figsize=gu.cm2inch(15,7))
# Coast
ind=np.where( (gpt['PTF YSM']==1) & (gpt['Ecozone BC L1']==meta['LUT']['GP']['Ecozone BC L1']['CWH']) | (gpt['PTF YSM']==1) & (gpt['Ecozone BC L1']==meta['LUT']['GP']['Ecozone BC L1']['MH']) )[0]
x=gpt['Age Med t0'][ind]
y=gpt['Ctot G Surv'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
ax1[0].bar(bin,mu,0.9*bw,fc=[0.85,1,0.5],ec='none',label='Growth of survivors')
y=gpt['Ctot G Recr'][ind]
N2,mu2,med2,sig2,se2=gu.discres(x,y,bw,bin)
ax1[0].bar(bin,mu2,0.9*bw,fc=[0.3,0.75,0.],ec='none',label='Growth from recruitment',bottom=mu)
ax1[0].set(ylabel='Biomass flux (tC ha$^-$$^1$ yr$^-$$^1$)',xlabel='Age, years',xticks=np.arange(0,400,bw),yticks=np.arange(-20,30,2),xlim=[0,250+bw],ylim=[-3,8])
ax1[0].yaxis.set_ticks_position('both'); ax1[0].xaxis.set_ticks_position('both'); ax1[0].tick_params(length=meta['Graphics']['gp']['tickl'])
y=gpt['Ctot Mort Nat'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
ax1[0].bar(bin,-mu,0.9*bw,fc=[0.8,0.6,0.4],ec='none',label='Mortality natural')
y=gpt['Ctot Mort Harv'][ind]
N2,mu2,med2,sig2,se2=gu.discres(x,y,bw,bin)
ax1[0].bar(bin,-mu2,0.9*bw,fc=[0.4,0.3,0.2],ec='none',label='Mortality harvest',bottom=-mu)
x=gpt['Age Med t0'][ind]
y=gpt['Ctot Net'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
ax1[0].plot(bin,mu,'-ko',ms=2.5,mfc='w',mec='k',label='Net growth')
ax1[0].plot([0,500],[0,0],'k-')
ax1[0].legend(loc='upper right',frameon=False,facecolor='w',edgecolor='w')
# Interior
ind=np.where( (gpt['PTF YSM']==1) & (gpt['Ecozone BC L1']!=meta['LUT']['GP']['Ecozone BC L1']['CWH']) & (gpt['Ecozone BC L1']!=meta['LUT']['GP']['Ecozone BC L1']['MH']) )[0]
x=gpt['Age Med t0'][ind]
y=gpt['Ctot G Surv'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
ax1[1].bar(bin,mu,0.9*bw,fc=[0.85,1,0.5],ec='none',label='Growth of survivors')
y=gpt['Ctot G Recr'][ind]
N2,mu2,med2,sig2,se2=gu.discres(x,y,bw,bin)
ax1[1].bar(bin,mu2,0.9*bw,fc=[0.3,0.75,0.],ec='none',label='Growth from recruitment',bottom=mu)
ax1[1].set(ylabel='Biomass flux (tC ha$^-$$^1$ yr$^-$$^1$)',xlabel='Age, years',xticks=np.arange(0,400,bw),yticks=np.arange(-20,30,2),xlim=[0,250+bw],ylim=[-3,8])
ax1[1].yaxis.set_ticks_position('both'); ax1[0].xaxis.set_ticks_position('both'); ax1[0].tick_params(length=meta['Graphics']['gp']['tickl'])
y=gpt['Ctot Mort Nat'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
ax1[1].bar(bin,-mu,0.9*bw,fc=[0.8,0.6,0.4],ec='none',label='Mortality natural')
y=gpt['Ctot Mort Harv'][ind]
N2,mu2,med2,sig2,se2=gu.discres(x,y,bw,bin)
ax1[1].bar(bin,-mu2,0.9*bw,fc=[0.4,0.3,0.2],ec='none',label='Mortality harvest',bottom=-mu)
y=gpt['Ctot Net'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
ax1[1].plot(bin,mu,'-ko',ms=2.5,mfc='w',mec='k',label='Net growth')
ax1[1].plot([0,500],[0,0],'k-')
ax1[1].yaxis.set_ticks_position('both'); ax1[1].xaxis.set_ticks_position('both'); ax1[1].tick_params(length=meta['Graphics']['gp']['tickl'])
gu.axletters(ax1,plt,0.04,0.93,FontColor=meta['Graphics']['gp']['cla'],LetterStyle='Caps',FontWeight='Bold')
plt.tight_layout()
gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_AgeResponseGrossGrowthMortalityByRegion_YSM','png',900)

#%% Age response of survivor growth (how much is associated with trees > 10m)

#ind=np.where( (gpt['Ecozone BC L1']!=meta['LUT']['GP']['Ecozone BC L1']['CWH']) & (gpt['Ecozone BC L1']!=meta['LUT']['GP']['Ecozone BC L1']['MH']) )[0]
ind=np.where( (gpt['Ecozone BC L1']!=-999) )[0]

bw=5; bin=np.arange(bw,100+bw,bw)
iT=np.where(bin<=50)[0]
lw=0.5; ms=3; mec=[0.27,0.49,0.77]; mfc=[1,1,1]; cl=[0.27,0.49,0.77];
plt.close('all'); fig,ax1=plt.subplots(1,figsize=gu.cm2inch(7.8,7))
x=gpt['Age Med t0'][ind]
y=gpt['Ctot G Surv'][ind]-gpt['Ctot G Surv HGT10'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
ax1.bar(bin[iT],mu[iT],0.9*bw,fc=[0.85,1,0.5],ec='none',label='Growth of survivors (tree height < 10m)')
y=gpt['Ctot G Surv HGT10'][ind]
N2,mu2,med2,sig2,se2=gu.discres(x,y,bw,bin)
ax1.bar(bin[iT],mu2[iT],0.9*bw,bottom=mu[iT],fc=[0.3,0.75,0.],ec='none',label='Growth of survivors (tree height > 10m)')
for i in range(bin.size):
    if bin[i]>50:
        continue
    a=mu[i]+mu2[i]
    if np.isnan(a)==False:
        ax1.text(bin[i],0.05+a,str(int(mu2[i]/a*100)) + '%',ha='center',fontsize=7)
ax1.set(ylabel='Growth of survivors (tC ha$^-$$^1$ yr$^-$$^1$)',xlabel='Age, years',xticks=np.arange(0,400,bw),yticks=np.arange(-20,30,0.5),xlim=[0,50+bw/1.5],ylim=[0,4])
ax1.legend(loc='best',frameon=False,facecolor='w',edgecolor='w')
ax1.yaxis.set_ticks_position('both'); ax1.xaxis.set_ticks_position('both'); ax1.tick_params(length=meta['Graphics']['gp']['tickl'])
plt.tight_layout()
gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_AgeResponseGrossGrowthAndAdvancedRegenContribution','png',900)

#%%
def MortalityRelativeVsAgeByBGC_CNY(meta,gpt):
	ord=np.flip(np.argsort(meta['Param']['BE']['BGC Zone Averages']['Area Treed (Mha)']))
	lab=np.array(['' for _ in range(ord.size)],dtype=object)
	plt.close('all'); fig,ax=plt.subplots(3,3,figsize=gu.cm2inch(22,11)); cnt=0
	for i in range(3):
		for j in range(3):
			zone=meta['Param']['BE']['BGC Zone Averages']['Name'][ord[cnt]]
			lab[cnt]=zone
			ind=np.where( (gpt['PTF CNY']==1) & (gpt['Ecozone BC L1']==meta['LUT']['GP']['Ecozone BC L1'][zone]) )[0]
			x=gpt['Age Mean t0'][ind]
			y=gpt['Ctot Mort'][ind]/gpt['Ctot L t0'][ind]*100
			bw=50; bin=np.arange(25,400,bw)
			N,mu,med,sig,se=gu.discres(x,y,bw,bin)
			mu[mu==np.inf]=np.nan
			se[se==np.inf]=np.nan
			ax[i,j].plot(bin,mu,'ko',lw=0.5,ms=3,mec=[0.27,0.44,0.79],mfc=[0.27,0.44,0.79])
			ax[i,j].errorbar(bin,mu,yerr=se,color=[0.27,0.44,0.79],ls='none',lw=1,capsize=2)
			ax[i,j].set(ylabel='Mortality rate (%)',xlabel='Stand age, years',xlim=[0,400],ylim=[0,np.nanmax(mu)+np.nanmax(se)+0.2])
			if (i==0) & (j==0):
				ax[i,j].legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
			ax[i,j].yaxis.set_ticks_position('both'); ax[i,j].xaxis.set_ticks_position('both'); ax[i,j].tick_params(length=meta['Graphics']['gp']['tickl'])
			plt.tight_layout()
			cnt=cnt+1
	gu.axletters(ax,plt,0.025,0.89,FontColor=meta['Graphics']['gp']['cla'],LetterStyle='NoPar',FontWeight='Bold',Labels=lab,LabelSpacer=0.035)
	#if meta['Graphics']['Print Figures']=='On':
	#	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_FullCompareAgeDistByBGC_CN_' + str(iScn+1),'png',900)
	return
MortalityRelativeVsAgeByBGC_CNY(meta,gpt)

#%% Biomass loss by burn severity class analysis

# meta['LUT']['Derived']['burnsev_comp1']
sc=z['bsr_sc']
fy=z['bsr_yr']
tsf=gpt['Year t0']-fy
v='Ctot L t0'
#v='Age Mean t0'
uSC=np.array([1,2,3,4]) #np.unique(sc)
tth=10
mu=np.zeros(uSC.size)
yr1=np.zeros(uSC.size)
yr2=np.zeros(uSC.size)
sig=np.zeros(uSC.size)
for iU in range(uSC.size):
    ind1=np.where( (sc==uSC[iU]) & (tsf>=-tth) & (tsf<0) & (gpt['Year IBM']>=0) )
    ind2=np.where( (sc==uSC[iU]) & (tsf>0) & (tsf<=tth) & (gpt['Year IBM']>=0) )
    mu[iU]=(np.nanmean(gpt[v][ind2])-np.nanmean(gpt[v][ind1]))/np.nanmean(gpt[v][ind1])*100
    yr1[iU]=np.nanmean(gpt['Year t0'][ind1])
    yr2[iU]=np.nanmean(gpt['Year t0'][ind2])
    #sig[iU]=np.mean(gpt[k][ind2])-np.mean(gpt[k][ind1])

u,N=gu.CountByCategories(sc[fy>0],'Percent')
#plt.hist(fy[fy>0])

# Control
ind1=np.where( (np.abs(gpt['Year t0']-np.mean(yr1))<=5) & (z['fire_yr']==0) )
ind2=np.where( (np.abs(gpt['Year t0']-np.mean(yr2))<=5) & (z['fire_yr']==0) )
muC=(np.nanmean(gpt[v][ind2])-np.nanmean(gpt[v][ind1]))/np.nanmean(gpt[v][ind1])*100
print(muC)

plt.close('all'); fig,ax1=plt.subplots(1,figsize=gu.cm2inch(7.8,5.5)); bw=0.4
ax1.bar(np.append(0,uSC)-bw/2,np.append(muC,mu),bw,facecolor=[0.75,0.85,1],label='Without correction')
ax1.bar(uSC+bw/2,mu-muC,bw,facecolor=[0.85,1,0.5],label='Control subtracted')
ax1.set(xticks=np.arange(5),yticks=np.arange(-100,100,10),xticklabels=['Control','Unburned','Low','Moderate','High'],ylabel='Change in biomass (%)',ylim=[-100,0])
ax1.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
ax1.yaxis.set_ticks_position('both'); ax1.xaxis.set_ticks_position('both'); ax1.tick_params(length=meta['Graphics']['gp']['tickl'])
gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_BiomassLossByBurnSeverityClass','png',900)

#%%

fy=z['fire_yr']
tsf=gpt['Year t0']-fy
v='Ctot L t0'
tth=5
mu=np.zeros(1)
sig=np.zeros(1)
ind1=np.where( (tsf>=-tth) & (tsf<0) & (gpt['Year IBM']>0) )
ind2=np.where( (tsf>0) & (tsf<=tth) & (gpt['Year IBM']>0) )
mu=(np.nanmean(gpt[v][ind2])-np.nanmean(gpt[v][ind1]))/np.nanmean(gpt[v][ind1])*100
#sig[iU]=np.mean(gpt[k][ind2])-np.mean(gpt[k][ind1])
print(mu)

plt.hist(fy[fy>0])

#%% Plot Profiles - import data

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
hy=z['harv_yr_comp1']
tsh=gpt['Year t0']-hy
ikp=np.where(gpt['Year Wildfire']<2200)[0]
prfl_h={}
for k in gpt.keys():
    prfl_h[k]={}
    prfl_h[k]['N'],prfl_h[k]['mu'],prfl_h[k]['med'],prfl_h[k]['sig'],prfl_h[k]['se']=gu.discres(tsh[ikp],gpt[k][ikp],bw,bin)

#%% Plot profiles (Age)

it0=np.where(bin<0)[0]; it1=np.where( (bin>0) )[0]
#v='Age VRI t0'
v='Age Med t0'

plt.close('all');
fig,ax=plt.subplots(1,3,figsize=gu.cm2inch(15.5,6)); ms=2
ax[0].plot([0,0],[0,300],'k-',lw=2,color=[0.8,0.8,0.8])
ax[0].plot(bin,prfl_wf[v]['mu'],'ko',ms=ms,mfc='k',mec='k',label='Live')
ax[0].errorbar(bin,prfl_wf[v]['mu'],yerr=prfl_wf[v]['se'],color='k',fmt='none',capsize=1.5,lw=0.5)
mu0=np.nanmean(prfl_wf[v]['mu'][it0])
mu1=np.nanmean(prfl_wf[v]['mu'][it1])
ax[0].plot(bin[it0],mu0*np.ones(it0.size),'k--',color='k',lw=0.5)
ax[0].plot(bin[it1],mu1*np.ones(it1.size),'k-.',color='k',lw=0.5)
#ax[0].text(np.mean(bin[it1]),1.82*mu1,'Change = ' + str(np.round(mu1-mu0,decimals=1)) + '\n(' + str(np.round( (mu1-mu0)/mu0*100 ,decimals=1)) + '%)',ha='center',fontsize=7,style='normal',color='k')
ax[0].set(ylabel='Age, years',xlabel='Years since wildfire',ylim=[0,200],xlim=[-40,80])
ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].tick_params(length=meta['Graphics']['gp']['tickl'])

ax[1].plot([0,0],[0,300],'k-',lw=2,color=[0.8,0.8,0.8])
ax[1].plot(bin,prfl_ibm[v]['mu'],'ko',ms=ms,mfc='k',mec='k',label='Live')
ax[1].errorbar(bin,prfl_ibm[v]['mu'],yerr=prfl_ibm[v]['se'],color='k',fmt='none',capsize=1.5,lw=0.5)
mu0=np.nanmean(prfl_ibm[v]['mu'][it0])
mu1=np.nanmean(prfl_ibm[v]['mu'][it1])
ax[1].plot(bin[it0],mu0*np.ones(it0.size),'k--',color='k',lw=0.5)
ax[1].plot(bin[it1],mu1*np.ones(it1.size),'k-.',color='k',lw=0.5)
#ax[1].text(np.mean(bin[it1]),1.3*mu1,'Change = ' + str(np.round(mu1-mu0,decimals=1)) + '\n(' + str(np.round( (mu1-mu0)/mu0*100 ,decimals=1)) + '%)',ha='center',fontsize=7,style='normal',color='k')
ax[1].set(ylabel='Age, years',xlabel='Years since severe outbreak',ylim=[0,200],xlim=[-40,80])
#ax[0].legend(loc='upper right',frameon=False,facecolor='w',edgecolor='w')
ax[1].yaxis.set_ticks_position('both'); ax[1].xaxis.set_ticks_position('both'); ax[1].tick_params(length=meta['Graphics']['gp']['tickl'])

ax[2].plot([0,0],[0,300],'k-',lw=2,color=[0.8,0.8,0.8])
ax[2].plot(bin,prfl_h[v]['mu'],'ko',ms=ms,mfc='k',mec='k',label='Live')
ax[2].errorbar(bin,prfl_h[v]['mu'],yerr=prfl_h[v]['se'],color='k',fmt='none',capsize=1.5,lw=0.5)
mu0=np.nanmean(prfl_h[v]['mu'][it0])
mu1=np.nanmean(prfl_h[v]['mu'][it1])
ax[2].plot(bin[it0],mu0*np.ones(it0.size),'k--',color='k',lw=0.5)
ax[2].plot(bin[it1],mu1*np.ones(it1.size),'k-.',color='k',lw=0.5)
#ax[2].text(np.mean(bin[it1]),1.3*mu1,'Change = ' + str(np.round(mu1-mu0,decimals=1)) + '\n(' + str(np.round( (mu1-mu0)/mu0*100 ,decimals=1)) + '%)',ha='center',fontsize=7,style='normal',color='k')
ax[2].set(ylabel='Age, years',xlabel='Years since harvest',ylim=[0,200],xlim=[-40,80])
#ax[0].legend(loc='upper right',frameon=False,facecolor='w',edgecolor='w')
ax[2].yaxis.set_ticks_position('both'); ax[2].xaxis.set_ticks_position('both'); ax[2].tick_params(length=meta['Graphics']['gp']['tickl'])
gu.axletters(ax,plt,0.04,0.92,FontColor=meta['Graphics']['gp']['cla'],LetterStyle='NoPar',FontWeight='Bold')#,Labels=['Wildfire','MPB','Wildfire','MPB'])
plt.tight_layout()
gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_Disturbance_Profiles_Age','png',900)

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
ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=meta['Graphics']['gp']['tickl'])

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
ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=meta['Graphics']['gp']['tickl'])

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
ax[0,2].yaxis.set_ticks_position('both'); ax[0,2].xaxis.set_ticks_position('both'); ax[0,2].tick_params(length=meta['Graphics']['gp']['tickl'])

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
ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=meta['Graphics']['gp']['tickl'])

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
ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=meta['Graphics']['gp']['tickl'])

ax[1,2].plot([0,0],[0,200],'k-',lw=2,color=[0.8,0.8,0.8])
ax[1,2].plot(bin,prfl_h[v]['mu'],'ko',ms=ms,mfc='k',mec='k',label='Live')
ax[1,2].errorbar(bin,prfl_h[v]['mu'],yerr=prfl_h[v]['se'],color='k',fmt='none',capsize=1.5,lw=0.5)
mu0=np.nanmean(prfl_h[v]['mu'][it0])
mu1=np.nanmean(prfl_h[v]['mu'][it1])
ax[1,2].plot(bin[it0],mu0*np.ones(it0.size),'k--',color='k',lw=0.5)
ax[1,2].plot(bin[it1],mu1*np.ones(it1.size),'k-.',color='k',lw=0.5)
#ax[1,2].text(np.mean(bin[it1]),0.7*mu1,'Change = ' + str(np.round(mu1-mu0,decimals=1)) + '\n(' + str(np.round( (mu1-mu0)/mu0*100 ,decimals=1)) + '%)',ha='center',fontsize=7,style='normal',color='k')
ax[1,2].set(ylabel='Dead wood (tC ha$^1$)',xlabel='Years since harvest',ylim=[0,40],xlim=[-40,80])
ax[1,2].yaxis.set_ticks_position('both'); ax[1,2].xaxis.set_ticks_position('both'); ax[1,2].tick_params(length=meta['Graphics']['gp']['tickl'])

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
ax[2,0].yaxis.set_ticks_position('both'); ax[2,0].xaxis.set_ticks_position('both'); ax[2,0].tick_params(length=meta['Graphics']['gp']['tickl'])

ax[2,1].plot([0,0],[-200,200],'k-',lw=2,color=[0.8,0.8,0.8])
ax[2,1].plot(bin,prfl_ibm[v]['mu'],'ko',ms=ms,mfc='k',mec='k',label='Live')
ax[2,1].errorbar(bin,prfl_ibm[v]['mu'],yerr=prfl_ibm[v]['se'],color='k',fmt='none',capsize=1.5,lw=0.5)
mu0=np.nanmean(prfl_ibm[v]['mu'][it0])
mu1=np.nanmean(prfl_ibm[v]['mu'][it1])
ax[2,1].plot(bin[it0],mu0*np.ones(it0.size),'k--',color='k',lw=0.5)
ax[2,1].plot(bin[it1],mu1*np.ones(it1.size),'k-.',color='k',lw=0.5)
#ax[2,1].text(np.mean(bin[it1]),1.55*mu1,'Change = ' + str(np.round(mu1-mu0,decimals=1)) + '\n(' + str(np.round( (mu1-mu0)/mu0*100 ,decimals=1)) + '%)',ha='center',fontsize=7,style='normal',color='k')
ax[2,1].set(ylabel='Net growth (tC ha$^{-1}$ yr$^{-1}$)',xlabel='Years since severe outbreak',ylim=[-2,5],xlim=[-40,80])
ax[2,1].yaxis.set_ticks_position('both'); ax[2,1].xaxis.set_ticks_position('both'); ax[2,1].tick_params(length=meta['Graphics']['gp']['tickl'])

ax[2,2].plot([0,0],[-200,200],'k-',lw=2,color=[0.8,0.8,0.8])
ax[2,2].plot(bin,prfl_h[v]['mu'],'ko',ms=ms,mfc='k',mec='k',label='Live')
ax[2,2].errorbar(bin,prfl_h[v]['mu'],yerr=prfl_h[v]['se'],color='k',fmt='none',capsize=1.5,lw=0.5)
mu0=np.nanmean(prfl_h[v]['mu'][it0])
mu1=np.nanmean(prfl_h[v]['mu'][it1])
ax[2,2].plot(bin[it0],mu0*np.ones(it0.size),'k--',color='k',lw=0.5)
ax[2,2].plot(bin[it1],mu1*np.ones(it1.size),'k-.',color='k',lw=0.5)
#ax[2,2].text(np.mean(bin[it1]),1.55*mu1,'Change = ' + str(np.round(mu1-mu0,decimals=1)) + '\n(' + str(np.round( (mu1-mu0)/mu0*100 ,decimals=1)) + '%)',ha='center',fontsize=7,style='normal',color='k')
ax[2,2].set(ylabel='Net growth (tC ha$^{-1}$ yr$^{-1}$)',xlabel='Years since harvest',ylim=[-2,5],xlim=[-40,80])
ax[2,2].yaxis.set_ticks_position('both'); ax[2,2].xaxis.set_ticks_position('both'); ax[2,2].tick_params(length=meta['Graphics']['gp']['tickl'])

gu.axletters(ax,plt,0.04,0.88,FontColor=meta['Graphics']['gp']['cla'],LetterStyle='NoPar',FontWeight='Bold')#,Labels=['Wildfire','MPB','Wildfire','MPB'])
plt.tight_layout()
gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\QA_Disturbance_Profiles','png',900)

#%% Mortality associations
gpg.GP_MortalityAssociations(meta,gpt)


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
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_Initial Tree Density By GBC Zone','png',900)

#%% Biomass by tree density class

d={}
for k in gpt.keys():
    d[k]=np.zeros(3)
    cnt=0
    for k2 in meta['LUT']['Derived']['tdc'].keys():
        ind=np.where( (gpt['PTF CN']==1) & (gpt['Tree Density Class (VRI)']==meta['LUT']['Derived']['tdc'][k2]) )[0]
        d[k][cnt]=np.nanmean(gpt[k][ind])
        cnt=cnt+1
d['Ctot L t0']
d['N L t0']
d['Vws L t0']
d['Ctot Net']/d['N L t0']*1000
d['Ctot Mort Harv']

plt.hist(gpt['Vws L t0'],bins=np.arange(0,1700,100))

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

#%% Differences between growth of forests that have been harvested or not
# Ability of harvest to cause site bias in productivity of young stands

bw=25; bin=np.arange(bw,150+bw,bw)
brw=23; lw=0.5; ms=3; cl1=np.array([0.27,0.49,0.77]); cl2=np.array([0.85,1,0.5])
plt.close('all'); fig,ax1=plt.subplots(1,figsize=gu.cm2inch(7.8,6))

ind=np.where( (gpt['PTF CNY']==1) & (z['harv_yr_comp1']==0) )[0]
x=gpt['Age Med t0'][ind]
y=gpt['Ctot Net'][ind]
N,mu1,med,sig,se=gu.discres(x,y,bw,bin)
ax1.bar(bin-5,mu1,10,fc=cl1,ec='none',label='No harvest history',zorder=-1)
for i in range(bin.size):
    ax1.errorbar(bin[i]-5,mu1[i],yerr=se[i],color=0.5*cl1,fmt='none',capsize=2,lw=0.5)

ind=np.where( (gpt['PTF CNY']==1) & (z['harv_yr_comp1']!=0) & (z['harv_yr_comp1']<1995) )[0]
x=gpt['Age Med t0'][ind]
y=gpt['Ctot Net'][ind]
N,mu2,med,sig,se=gu.discres(x,y,bw,bin)
ax1.bar(bin+5,mu2,10,fc=cl2,ec='none',label='Harvested',zorder=-1)
for i in range(bin.size):
    ax1.errorbar(bin[i]+5,mu2[i],yerr=se[i],color=0.5*cl2,fmt='none',capsize=2,lw=0.5)
mu1a=np.nanmean(mu1[0:2])
mu2a=np.nanmean(mu2[0:2])
print(mu1a)
print(mu2a)
print((mu2a-mu1a)/mu1a*100)
ax1.plot([0,500],[0,0],'-k',lw=lw)
ax1.set(xlabel='Age, years',xticks=np.arange(0,400,bw),ylabel='Net growth (tC ha$^-$$^1$ yr$^-$$^1$)',xlim=[0,175],ylim=[-2,4])
ax1.yaxis.set_ticks_position('both'); ax1.xaxis.set_ticks_position('both'); ax1.tick_params(length=meta['Graphics']['gp']['tickl'])
ax1.legend(loc='upper right',frameon=False,facecolor=None,edgecolor='w')
plt.tight_layout()
gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_AgeResponseHarvVsNoHarv','png',900)

#%% Deciduous analysis

u=np.unique(gpt['Ecozone BC L1'][gpt['Ecozone BC L1']>0])
lab=np.array(['' for _ in range(u.size)],dtype=object)
for iU in range(u.size):
    lab[iU]=ugp.lut_id2cd(meta['LUT']['GP']['Ecozone BC L1'],u[iU])
   
d={}
d['SS']={'sum':np.zeros(u.size)}
d['SS Net']={'sum':np.zeros(u.size)}
d['Pct with D']={'sum':np.zeros(u.size)}
for iU in range(u.size):
    ind0=np.where( (gpt['PTF CN']==1) & (gpt['Ecozone BC L1']==u[iU]) )[0]
    ind1=np.where( (gpt['PTF CN']==1) & (gpt['Ecozone BC L1']==u[iU]) & (gpt['Deciduous L %N t0']>0) )[0]
    d['Pct with D']['sum'][iU]=ind1.size/ind0.size*100
    d['SS']['sum'][iU]=ind0.size
    ind2=np.where( (gpt['PTF CN']==1) & (gpt['Ctot G Surv']>0) & (gpt['Ecozone BC L1']==u[iU]) )[0]
    d['SS Net']['sum'][iU]=ind2.size

vL=['Ctot L t0','Ctot Net','Age Mean t0']
for v in vL:
    d[v]={}
    d[v]['Cmu']=np.zeros(u.size)
    d[v]['Cse']=np.zeros(u.size)
    d[v]['Dmu']=np.zeros(u.size)
    d[v]['Dse']=np.zeros(u.size)

#ind=np.where( (gpt['PTF CN']==1) & (gpt['Ecozone BC L1']==meta['LUT']['GP']['Ecozone BC L1']['MH']) )[0]
#gpt['Deciduous L %N t0'][ind]

for v in vL:
    cnt=0
    for iU in range(u.size):
        ind=np.where( (gpt['PTF CN']==1) & (gpt['Ecozone BC L1']==u[iU]) & (gpt['Deciduous L %N t0']<5) )[0]
        d[v]['Cmu'][cnt]=np.nanmean(gpt[v][ind])
        d[v]['Cse'][cnt]=np.nanstd(gpt[v][ind])/np.sqrt(ind.size)    
        ind=np.where( (gpt['PTF CN']==1) & (gpt['Ecozone BC L1']==u[iU]) & (gpt['Deciduous L %N t0']>=5) )[0]
        d[v]['Dmu'][cnt]=np.nanmean(gpt[v][ind])
        d[v]['Dse'][cnt]=np.nanstd(gpt[v][ind])/np.sqrt(ind.size)
        cnt=cnt+1

# Plot proportion with deciduous
d0=copy.deepcopy(d)
lab0=lab.copy()

# Remove zones with too few data
ind=np.where(d0['SS']['sum']>10)[0]
u0=u[ind]
for k in d0.keys():
    for v in d0[k].keys():
        d0[k][v]=d0[k][v][ind]
lab0=lab[ind]
     
ord=np.argsort(d0['Pct with D']['sum'])
for k in d0.keys():
    for v in d0[k].keys():
        d0[k][v]=np.flip(d0[k][v][ord])
lab0=np.flip(lab0[ord])

cl=np.array([ [0.3,0.45,0.76] ])

plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,5))
ax.bar(np.arange(u0.size),d0['Pct with D']['sum'],facecolor=cl[0,:],label='')
ax.set(xticks=np.arange(u0.size),xticklabels=lab0,ylabel='Proportion of stands with >5% \nbroadleaf deciduous trees (%)',yticks=np.arange(0,110,10),xlim=[-0.5,u0.size-0.5],ylim=[0,100])
ax.legend(loc='upper right',frameon=False,facecolor=[1,1,1],labelspacing=0.25)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
plt.tight_layout()
gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_ProportionWithDeciduousTrees_CN','png',900)

#%% Deciduousness on harvested vs. non-harvested stands

u=np.unique(gpt['Ecozone BC L1'][gpt['Ecozone BC L1']>0])
lab=np.array(['' for _ in range(u.size)],dtype=object)
for iU in range(u.size):
	lab[iU]=ugp.lut_id2cd(meta['LUT']['GP']['Ecozone BC L1'],u[iU])

d={}
d['SS']={'sum':np.zeros(u.size)}
d['Pct H']={'mu':np.zeros(u.size)}
d['Pct NH']={'mu':np.zeros(u.size)}
for iU in range(u.size):
	ind=np.where( (gpt['PTF CNY']==1) & (gpt['Ecozone BC L1']==u[iU]) & (gpt['Harvest Mask']==0) )[0]
	d['Pct NH']['mu'][iU]=np.nanmean(gpt['Deciduous L %N t0'][ind])
	d['SS']['sum'][iU]=ind.size
	ind=np.where( (gpt['PTF CNY']==1) & (gpt['Ecozone BC L1']==u[iU]) & (gpt['Harvest Mask']==1) )[0]
	d['Pct H']['mu'][iU]=np.nanmean(gpt['Deciduous L %N t0'][ind])

# Plot proportion with deciduous
d0=copy.deepcopy(d)
lab0=lab.copy()

# Remove zones with too few data
ind=np.where(d0['SS']['sum']>5)[0]
u0=u[ind]
for k in d0.keys():
	for v in d0[k].keys():
		d0[k][v]=d0[k][v][ind]
lab0=lab[ind]

ord=np.argsort(d0['Pct NH']['mu'])
for k in d0.keys():
	for v in d0[k].keys():
		d0[k][v]=np.flip(d0[k][v][ord])
lab0=np.flip(lab0[ord])

cl=np.array([ [0.3,0.45,0.76],[0.85,1,0.5] ])
bw=0.3

plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,5))
#ax.bar(np.arange(u0.size),d0['Pct NH']['mu'],facecolor=cl[0,:],label='')
ax.bar(np.arange(u0.size)-bw/2,d0['Pct NH']['mu'],bw,facecolor=cl[0,:],label='Stands with no recorded harvest')
ax.bar(np.arange(u0.size)+bw/2,d0['Pct H']['mu'],bw,facecolor=cl[1,:],label='Stands that have been harvested')
ax.set(xticks=np.arange(u0.size),xticklabels=lab0,ylabel='Mean proportion of broadleaf\ndeciduous trees (%)',yticks=np.arange(0,110,5),xlim=[-0.5,u0.size-0.5],ylim=[0,50])
ax.legend(loc='upper right',frameon=False,facecolor=[1,1,1],labelspacing=0.25)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
plt.tight_layout()
#gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_DecidBroadleafProportionByHarvestHistoryAndBGC_CN','png',900)

#%% Plot deciduous vs. conifer age
d0=copy.deepcopy(d)
lab0=lab.copy()

# Remove zones with too few data
ind=np.where(d0['SS']['sum']>10)[0]
u0=u[ind]
for k in d0.keys():
    for v in d0[k].keys():
        d0[k][v]=d0[k][v][ind]
lab0=lab[ind]
   
ord=np.argsort(d0['Age Mean t0']['Cmu'])
for k in d0.keys():
    for v in d0[k].keys():
        d0[k][v]=np.flip(d0[k][v][ord])
lab0=np.flip(lab0[ord])

cl=np.array([ [0.3,0.45,0.76],[0.85,1,0.5] ])
bw=0.3

plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,6))
ax.bar(np.arange(u0.size)-bw/2,d0['Age Mean t0']['Cmu'],bw,facecolor=cl[0,:],label='')
ax.bar(np.arange(u0.size)+bw/2,d0['Age Mean t0']['Dmu'],bw,facecolor=cl[1,:],label='')
ax.errorbar(np.arange(u0.size)-bw/2,d0['Age Mean t0']['Cmu'],yerr=d0['Age Mean t0']['Cse'],color=gp['cla'],fmt='none',capsize=2,lw=0.5)
ax.errorbar(np.arange(u0.size)+bw/2,d0['Age Mean t0']['Dmu'],yerr=d0['Age Mean t0']['Dse'],color=gp['cla'],fmt='none',capsize=2,lw=0.5)
ax.set(xticks=np.arange(u0.size),xticklabels=lab0,ylabel='Stand age (years)',yticks=np.arange(0,500,50),xlim=[-0.5,u0.size-0.5],ylim=[0,300])
ax.legend(loc='upper right',frameon=False,facecolor=[1,1,1],labelspacing=0.25)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
plt.tight_layout()
gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_CompareDecidConifAge_CN','png',900)

#%% Plot deciduous vs. conifer tree biomass
d0=copy.deepcopy(d)
lab0=lab.copy()

# Remove zones with too few data
ind=np.where(d0['SS']['sum']>10)[0]
u0=u[ind]
for k in d0.keys():
    for v in d0[k].keys():
        d0[k][v]=d0[k][v][ind]
lab0=lab[ind]

v0='Ctot L t0'
ord=np.argsort(d0[v0]['Cmu'])
for k in d0.keys():
    for v in d0[k].keys():
        d0[k][v]=np.flip(d0[k][v][ord])
lab0=np.flip(lab0[ord])

cl=np.array([ [0.3,0.45,0.76],[0.85,1,0.5] ])
bw=0.3

plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,6))
ax.bar(np.arange(u0.size)-bw/2,d0[v0]['Cmu'],bw,facecolor=cl[0,:],label='> 95% coniferous')
ax.bar(np.arange(u0.size)+bw/2,d0[v0]['Dmu'],bw,facecolor=cl[1,:],label='> 5% broadleaf deciduous')
ax.errorbar(np.arange(u0.size)-bw/2,d0[v0]['Cmu'],yerr=d0[v0]['Cse'],color=gp['cla'],fmt='none',capsize=2,lw=0.5)
ax.errorbar(np.arange(u0.size)+bw/2,d0[v0]['Dmu'],yerr=d0[v0]['Dse'],color=gp['cla'],fmt='none',capsize=2,lw=0.5)
ax.set(xticks=np.arange(u0.size),xticklabels=lab0,ylabel='Tree biomass (tC/ha)',yticks=np.arange(0,300,20),xlim=[-0.5,u0.size-0.5],ylim=[0,200])
ax.legend(loc='upper right',frameon=False,facecolor=[1,1,1],labelspacing=0.25)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
plt.tight_layout()
gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_CompareDecidConifBiomass_CN','png',900)

#%% Plot deciduous vs. conifer tree biomass production
d0=copy.deepcopy(d)
lab0=lab.copy()

# Remove zones with too few data
ind=np.where(d0['SS Net']['sum']>10)[0]
u0=u[ind]
for k in d0.keys():
    for v in d0[k].keys():
        d0[k][v]=d0[k][v][ind]
lab0=lab[ind]

v0='Ctot Net'
ord=np.argsort(d0[v0]['Cmu'])
for k in d0.keys():
    for v in d0[k].keys():
        d0[k][v]=np.flip(d0[k][v][ord])
lab0=np.flip(lab0[ord])

cl=np.array([ [0.3,0.45,0.76],[0.85,1,0.5] ])
bw=0.3

plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,6))
ax.plot([-10,20],[0,0],'k-',lw=0.5)
ax.bar(np.arange(u0.size)-bw/2,d0[v0]['Cmu'],bw,facecolor=cl[0,:],label='>95% coniferous')
ax.bar(np.arange(u0.size)+bw/2,d0[v0]['Dmu'],bw,facecolor=cl[1,:],label='>5% broadleaf deciduous')
ax.errorbar(np.arange(u0.size)-bw/2,d0[v0]['Cmu'],yerr=d0[v0]['Cse'],color=gp['cla'],fmt='none',capsize=2,lw=0.5)
ax.errorbar(np.arange(u0.size)+bw/2,d0[v0]['Dmu'],yerr=d0[v0]['Dse'],color=gp['cla'],fmt='none',capsize=2,lw=0.5)
ax.set(xticks=np.arange(u0.size),xticklabels=lab0,ylabel='Net biomass production (tC/ha/yr)',yticks=np.arange(-10,300,1),xlim=[-0.5,u0.size-0.5],ylim=[-3,5])
ax.legend(loc='upper right',frameon=False,facecolor=[1,1,1],labelspacing=0.25)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
plt.tight_layout()
gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_CompareDecidConifNetBiomassProd_CN','png',900)

#%% Deciduous analysis (harvest footprint)

ind0=np.where( (gpt['PTF CN']==1) & (gpt['Harvest Mask']==1) )[0]
ind1=np.where( (gpt['PTF CN']==1) & (gpt['Harvest Mask']==1) & (gpt['Deciduous L %N t0']>10) )[0]
ind2=np.where( (gpt['PTF CN']==1) & (gpt['Harvest Mask']==0) )[0]
ind3=np.where( (gpt['PTF CN']==1) & (gpt['Harvest Mask']==0) & (gpt['Deciduous L %N t0']>10) )[0]
ind1.size/ind0.size*100
ind3.size/ind2.size*100

ind0=np.where( (gpt['PTF CN']==1) & (gpt['Harvest Mask']==1) )[0]
ind1=np.where( (gpt['PTF CN']==1) & (gpt['Harvest Mask']==1) & (gpt['Deciduous L %N t0']>10) )[0]
mu1=np.nanmean(gpt['Ctot L t0'][ind0])
mu2=np.nanmean(gpt['Ctot L t0'][ind1])
se1=np.nanstd(gpt['Ctot L t0'][ind0])/np.sqrt(ind0.size)
se2=np.nanstd(gpt['Ctot L t0'][ind1])/np.sqrt(ind1.size)

ind0=np.where( (gpt['PTF CN']==1) & (gpt['Harvest Mask']==0) )[0]
ind1=np.where( (gpt['PTF CN']==1) & (gpt['Harvest Mask']==0) & (gpt['Deciduous L %N t0']>10) )[0]
mu3=np.nanmean(gpt['Ctot L t0'][ind0])
mu4=np.nanmean(gpt['Ctot L t0'][ind1])
se3=np.nanstd(gpt['Ctot L t0'][ind0])/np.sqrt(ind0.size)
se4=np.nanstd(gpt['Ctot L t0'][ind1])/np.sqrt(ind1.size)

cl=np.array([ [0.3,0.45,0.76],[0.85,1,0.5] ])
bw=0.3

plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(10,6))
ax.bar(1-bw/1.85,mu1,bw,facecolor=cl[0,:],label='100% Coniferous')
ax.bar(1+bw/1.85,mu2,bw,facecolor=cl[1,:],label='>0% Deciduous')
ax.bar(1.7-bw/1.85,mu3,bw,facecolor=cl[0,:],label='')
ax.bar(1.7+bw/1.85,mu4,bw,facecolor=cl[1,:],label='')
ax.text(1+bw/1.85,mu1+21,'+' + str(int((mu2-mu1)/mu1*100)) + '(%)',fontsize=10,ha='center')
ax.text(1.7+bw/1.85,mu4+7,'+' + str(int((mu4-mu3)/mu3*100)) + '(%)',fontsize=10,ha='center')
ax.errorbar(1-bw/1.85,mu1,yerr=se1,color=gp['cla'],fmt='none',capsize=2,lw=0.5)
ax.errorbar(1+bw/1.85,mu2,yerr=se2,color=gp['cla'],fmt='none',capsize=2,lw=0.5)
ax.errorbar(1.7-bw/1.85,mu3,yerr=se3,color=gp['cla'],fmt='none',capsize=2,lw=0.5)
ax.errorbar(1.7+bw/1.85,mu4,yerr=se4,color=gp['cla'],fmt='none',capsize=2,lw=0.5)
ax.set(xticks=[1,1.7],xticklabels=['Areas that have\nbeen harvested','Areas that have not\n been harvested'],ylabel='Tree biomass (tC/ha)',yticks=np.arange(-10,300,10),xlim=[0.6,2.1],ylim=[0,100])
ax.legend(loc='upper left',frameon=False,facecolor=[1,1,1],labelspacing=0.25)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
plt.tight_layout()
gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_CompareDecidConifNetBiomassHarvestFoot_CN','png',900)


#%%

