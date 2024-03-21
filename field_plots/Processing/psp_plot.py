#%% Import modules
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import copy
import pandas as pd
import fcgadgets.macgyver.util_general as gu
import fcgadgets.bc1ha.bc1ha_util as u1ha
import fcexplore.field_plots.Processing.psp_util as ugp

#%% Sample frequency per BGC zone
def GP_SamplingFrequencyByBGC(meta,gpt,soc):
	u=np.unique(gpt['Ecozone BC L1']); u=u[(u>0) & (u!=10)]
	plt.close('all'); fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(22,9)); barw=0.35;
	# VRI
	d={'lab':np.array(['' for _ in range(u.size)],dtype=object),'N First':np.zeros(u.size),'N Rem':np.zeros(u.size),'Area':np.zeros(u.size)}
	for iU in range(u.size):
		d['lab'][iU]=ugp.lut_id2cd(meta['LUT']['GP']['Ecozone BC L1'],u[iU])
		#ind0=np.where(meta['Param']['BE']['ByBGCZ']['Name']==d['lab'][iU])[0]
		ind1=np.where( (gpt['PTF VRI']==1) & (gpt['Ecozone BC L1']==u[iU]) & (gpt['Year t0']>0) )
		ind2=np.where( (gpt['PTF VRI']==1) & (gpt['Ecozone BC L1']==u[iU]) & (gpt['Year t0']>0) & (gpt['Year t1']>=0) )
		d['Area'][iU]=meta['Param']['BE']['ByBGCZ'][d['lab'][iU]]['Area Treed (Mha)']
		d['N First'][iU]=ind1[0].size
		d['N Rem'][iU]=ind2[0].size
	ax[0,0].bar(np.arange(u.size)-barw/2,d['N First']/(d['Area']),barw,fc=[0.75,0.85,1],ec=None,label='Single measurement')
	ax[0,0].bar(np.arange(u.size)+barw/2,d['N Rem']/(d['Area']),barw,fc=[0.85,1,0.5],ec=None,label='With remeasurement')
	for iU in range(u.size):
		ax[0,0].text(iU-barw/2,d['N First'][iU]/(d['Area'][iU])+1.5,str(int(d['N First'][iU])),color=0.5*np.array([0.7,0.8,1]),fontsize=5,ha='center')
		ax[0,0].text(iU+barw/2,d['N Rem'][iU]/(d['Area'][iU])+1.5,str(int(d['N Rem'][iU])),color=0.5*np.array([0.8,1,0.7]),fontsize=5,ha='center')
	ax[0,0].set(xticks=np.arange(0,u.size),xticklabels=d['lab'],ylabel='Sampling frequency\n (plots Mha$^{-1}$)',xlim=[-0.5,u.size-1+0.5],ylim=[0,250])
	ax[0,0].set_xticklabels(d['lab'],rotation=90)
	ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=meta['Graphics']['gp']['tickl'])
	#plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
	# CMI/NFI
	d={'lab':np.array(['' for _ in range(u.size)],dtype=object),'N First':np.zeros(u.size),'N Rem':np.zeros(u.size),'Area':np.zeros(u.size)}
	for iU in range(u.size):
		d['lab'][iU]=ugp.lut_id2cd(meta['LUT']['GP']['Ecozone BC L1'],u[iU])
		#ind0=np.where(meta['Param']['BE']['ByBGCZ']['Name']==d['lab'][iU])[0]
		ind1=np.where( (gpt['PTF CN']==1) & (gpt['Ecozone BC L1']==u[iU]) & (gpt['Year t0']>0) )
		ind2=np.where( (gpt['PTF CN']==1) & (gpt['Ecozone BC L1']==u[iU]) & (gpt['Year t0']>0) & (gpt['Year t1']>=0) )
		d['Area'][iU]=meta['Param']['BE']['ByBGCZ'][d['lab'][iU]]['Area Treed (Mha)']
		d['N First'][iU]=ind1[0].size
		d['N Rem'][iU]=ind2[0].size
	ax[0,1].bar(np.arange(u.size)-barw/2,d['N First']/(d['Area']),barw,fc=[0.75,0.85,1],ec=None,label='Single measurement')
	ax[0,1].bar(np.arange(u.size)+barw/2,d['N Rem']/(d['Area']),barw,fc=[0.85,1,0.5],ec=None,label='With remeasurement')
	for iU in range(u.size):
		ax[0,1].text(iU-barw/2,d['N First'][iU]/(d['Area'][iU])+1.5,str(int(d['N First'][iU])),color=0.5*np.array([0.7,0.8,1]),fontsize=5,ha='center')
		ax[0,1].text(iU+barw/2,d['N Rem'][iU]/(d['Area'][iU])+1.5,str(int(d['N Rem'][iU])),color=0.5*np.array([0.8,1,0.7]),fontsize=5,ha='center')
	ax[0,1].set(xticks=np.arange(0,u.size),xticklabels=d['lab'],ylabel='Sampling frequency\n (plots Mha$^{-1}$)',xlim=[-0.5,u.size-1+0.5],ylim=[0,80])
	ax[0,1].set_xticklabels(d['lab'],rotation=90)
	ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=meta['Graphics']['gp']['tickl'])
	ax[0,1].legend(loc='upper right',frameon=False,facecolor=[1,1,1],labelspacing=0.25)
	# YSM
	d={'lab':np.array(['' for _ in range(u.size)],dtype=object),'N First':np.zeros(u.size),'N Rem':np.zeros(u.size),'Area':np.zeros(u.size)}
	for iU in range(u.size):
		d['lab'][iU]=ugp.lut_id2cd(meta['LUT']['GP']['Ecozone BC L1'],u[iU])
		#ind0=np.where(meta['Param']['BE']['ByBGCZ']['Name']==d['lab'][iU])[0]
		ind1=np.where( (gpt['PTF YSM']==1) & (gpt['Ecozone BC L1']==u[iU]) & (gpt['Year t0']>0) )
		ind2=np.where( (gpt['PTF YSM']==1) & (gpt['Ecozone BC L1']==u[iU]) & (gpt['Year t0']>0) & (gpt['Year t1']>=0) )
		d['Area'][iU]=meta['Param']['BE']['ByBGCZ'][d['lab'][iU]]['Area Treed (Mha)']
		d['N First'][iU]=ind1[0].size
		d['N Rem'][iU]=ind2[0].size
	ax[1,0].bar(np.arange(u.size)-barw/2,d['N First']/(d['Area']),barw,fc=[0.75,0.85,1],ec=None,label='Single measurement')
	ax[1,0].bar(np.arange(u.size)+barw/2,d['N Rem']/(d['Area']),barw,fc=[0.85,1,0.5],ec=None,label='With remeasurement')
	for iU in range(u.size):
		ax[1,0].text(iU-barw/2,d['N First'][iU]/(d['Area'][iU])+1.5,str(int(d['N First'][iU])),color=0.5*np.array([0.7,0.8,1]),fontsize=5,ha='center')
		ax[1,0].text(iU+barw/2,d['N Rem'][iU]/(d['Area'][iU])+1.5,str(int(d['N Rem'][iU])),color=0.5*np.array([0.8,1,0.7]),fontsize=5,ha='center')
	ax[1,0].set(xticks=np.arange(0,u.size),xticklabels=d['lab'],ylabel='Sampling frequency\n (plots Mha$^{-1}$)',xlim=[-0.5,u.size-1+0.5],ylim=[0,80])
	ax[1,0].set_xticklabels(d['lab'],rotation=90)
	ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=meta['Graphics']['gp']['tickl'])
	# Soils
	d={'lab':np.array(['' for _ in range(u.size)],dtype=object),'N First':np.zeros(u.size),'N Rem':np.zeros(u.size),'Area':np.zeros(u.size)}
	for iU in range(u.size):
		d['lab'][iU]=ugp.lut_id2cd(meta['LUT']['GP']['Ecozone BC L1'],u[iU])
		#ind0=np.where(meta['Param']['BE']['ByBGCZ']['Name']==d['lab'][iU])[0]
		ind1=np.where( (soc['bgcz']==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'][d['lab'][iU]]) & (soc['TOT_C_THA']>0) )
		d['Area'][iU]=meta['Param']['BE']['ByBGCZ'][d['lab'][iU]]['Area Treed (Mha)']
		d['N First'][iU]=ind1[0].size
		d['N Rem'][iU]=0
	ax[1,1].bar(np.arange(u.size)-barw/2,d['N First']/(d['Area']),barw,fc=[0.75,0.85,1],ec=None,label='Single measurement')
	#ax[1,1].bar(np.arange(u.size)+barw/2,d['N Rem']/(d['Area']),barw,fc=[0.85,1,1.5],ec=None,label='With remeasurement')
	for iU in range(u.size):
		ax[1,1].text(iU-barw/2,d['N First'][iU]/(d['Area'][iU])+3,str(int(d['N First'][iU])),color=0.5*np.array([0.7,0.8,1]),fontsize=5,ha='center')
		#ax[1,1].text(iU+barw/2,d['N Rem'][iU]/(d['Area'][iU])+1.5,str(int(d['N Rem'][iU])),color=0.5*np.array([0.8,1,1.7]),fontsize=5,ha='center')
	ax[1,1].set(xticks=np.arange(0,u.size),xticklabels=d['lab'],ylabel='Sampling frequency\n (plots Mha$^{-1}$)',xlim=[-0.5,u.size-1+0.5],ylim=[0,850])
	ax[1,1].set_xticklabels(d['lab'],rotation=90)
	ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=meta['Graphics']['gp']['tickl'])
	plt.tight_layout()
	gu.axletters(ax,plt,0.02,0.89,FontColor=meta['Graphics']['gp']['cla'],LetterStyle='NoPar',FontWeight='Bold')
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\GP_SamplingFrequencyByBGC','png',900)

	return

#%%
def GP_Biomass_vs_VolumeWholeStem(meta,gpt):
	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(7.8,7))
	ax.plot([0,2000],[0,2000],'-',lw=2,color=[0.8,0.8,0.8])
	ind=np.where( (gpt['PTF CN']==1) & (gpt['Vws L t0']>0) & (gpt['Ctot L t0']>0) )[0]
	ax.plot(gpt['Vws L t0'][ind],gpt['Ctot L t0'][ind],'b.',ms=4,mec='w',mfc=[0.27,0.44,0.79],markeredgewidth=0.25)
	rs,txt=gu.GetRegStats(gpt['Vws L t0'][ind],gpt['Ctot L t0'][ind])
	ax.plot(rs['xhat Line'],rs['yhat Line'],'-k')
	ax.set(ylabel='Tree biomass (tC ha$^{-1}$)',xlabel='Whole-stem volume (m$^{3}$ ha$^{-1}$)',xlim=[0,1700],ylim=[0,1700])
	ax.text(925,480,txt,fontsize=8)
	ax.text(1200,1200,'1:1',fontsize=8,ha='center')
	#plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_Biomass_vs_VolumeWholeStem','png',900)
	return

#%%
def GP_QA_Age_vs_Age(meta,gpt):
	x=gpt['Age VRI t0']
	y=gpt['Age Mean t0']
	ikp=np.where( (x>=0) & (x<=2000) & (y>=0) & (y<=2000) )
	x=x[ikp]
	y=y[ikp]
	rs,txt=gu.GetRegStats(x,y)
	bw=20; bin=np.arange(10,420,bw); N,mu,med,sig,se=gu.discres(x,y,bw,bin)
	
	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(7.8,7.8))
	plt.plot([0,2000],[0,2000],'k-',lw=2,color=[0.8,0.8,0.8])
	plt.plot(x,y,'k.',mfc=[0.25,0.5,1],mec=[0.5,0.75,1],ms=3)
	plt.plot(bin,mu,'ks',mew=0.5,mec='w',mfc='k')
	plt.plot(rs['xhat Line'],rs['yhat Line'],'k-')
	ax.set(xlabel='Age from VRI (years)',ylabel='Age from cores (years)',xlim=[0,1000],ylim=[0,1000])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	plt.tight_layout()
	ax.text(900,100,txt,ha='right',fontsize=7)
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_Age_VRI_vs_PlotCores_CNVY','png',900)
	return

#%%
def GP_AgeClassDistributionAll(meta,gpt):
	x=np.arange(0,401,1)
	y=z['age_vri']
	y=y[y>=0]
	kde=stats.gaussian_kde(y)
	p1=kde(x)

	lw=1.25
	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(14,6))
	#plt.plot(x,p1/np.sum(p1)*100,'b-',lw=lw,color=[0.27,0.49,0.77],label='VRI R1 Comp Polygons')
	ind=np.where( (gpt['Plot Type']==meta['LUT']['GP']['Plot Type']['VRI']) & (gpt['Age Mean t0']>=0) )[0]
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
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_AgeDistribution','png',900)
	return

#%% Age class distribution (by BGC Zone)
def GP_AgeDistByBGC_CN(meta,gpt):
	x=np.arange(0,501,1); yt=np.arange(0,2,0.1)
	ord=np.flip(np.argsort(meta['Param']['Raw']['ByBGCZ']['Area Treed (Mha)']))
	lab=np.array(['' for _ in range(ord.size)],dtype=object)
	plt.close('all'); fig,ax=plt.subplots(3,3,figsize=gu.cm2inch(22,11)); cnt=0
	for i in range(3):
		for j in range(3):
			zone=meta['Param']['Raw']['ByBGCZ']['Name'][ord[cnt]]
			#ind=np.where( (gpt['Ecozone BC L1']==meta['LUT']['GP']['Ecozone BC L1'][zone]) & (gpt['Plot Type']==meta['LUT']['GP']['Plot Type']['VRI']) & (gpt['Age Mean t0']>=0) )[0]
			#kde=stats.gaussian_kde(gpt['Age Mean t0'][ind])
			#p=kde(x); y1=p/np.sum(p)*100
			#ax[i,j].plot(x,y1,'k-',lw=0.75,color=[0.5,0,1],label='VRI ground plots')
	
			ind=np.where( (gpt['Ecozone BC L1']==meta['LUT']['GP']['Ecozone BC L1'][zone]) & (gpt['PTF CN']==1) & (gpt['Age Mean t0']>=0) )[0]
			kde=stats.gaussian_kde(gpt['Age Mean t0'][ind])
			p=kde(x); y2=p/np.sum(p)*100
			ax[i,j].fill_between(x,y2,color=[0.7,0.8,0.95])
			ax[i,j].plot(x,y2,'b-',lw=0.75,color=[0.27,0.44,0.79],label='CMI + NFI network\ntree cores')
			#ym=np.maximum(np.max(y1),np.max(y2)); indYT=np.where(yt>ym)[0]
			ym=np.max(y2); indYT=np.where(yt>ym)[0]
			mu=np.mean(gpt['Age Mean t0'][ind])
			ax[i,j].plot([mu,mu],[0,20],'b--',color=[0.27,0.44,0.79],label='Mean')
			ax[i,j].set(ylabel='Frequency (%)',xlabel='Stand age, years',xlim=[0,500],ylim=[0,yt[indYT[1]]])
			if (i==0) & (j==0):
				ax[i,j].legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
			ax[i,j].yaxis.set_ticks_position('both'); ax[i,j].xaxis.set_ticks_position('both'); ax[i,j].tick_params(length=meta['Graphics']['gp']['tickl'])
			plt.tight_layout()
			lab[cnt]=zone
			cnt=cnt+1
	gu.axletters(ax,plt,0.025,0.89,FontColor=meta['Graphics']['gp']['cla'],LetterStyle='NoPar',FontWeight='Bold',Labels=lab,LabelSpacer=0.035)
	#if meta['Graphics']['Print Figures']=='On':
	#	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\QA_FullCompareAgeDistByBGC_CN_' + str(iScn+1),'png',900)
	return

#%%
def GP_ExportTableStatsByPlotType(meta,gpt,**kwargs):
	d={}
	ind=np.where( (gpt['PTF VRI']==1) )[0]
	for k in gpt.keys():
		d[k]=np.round(np.nanmean(gpt[k][ind]),decimals=2)
	ind=np.where( (gpt['PTF YSM']==1) )[0]
	for k in gpt.keys():
		d[k]=np.append(d[k],np.round(np.nanmean(gpt[k][ind]),decimals=2))
	ind=np.where( (gpt['PTF CN']==1) )[0]
	for k in gpt.keys():
		d[k]=np.append(d[k],np.round(np.nanmean(gpt[k][ind]),decimals=2))
	df=pd.DataFrame(d,index=[0,1,2])
	if kwargs['export']==True:
		df.to_excel(meta['Paths']['GP']['DB'] + '\\Processed\\SummarySL_ByPlotType.xlsx')
	return df

#%%
def GP_BiomassByBGC_CNV(meta,gpt,dBGC):
	d0=copy.deepcopy(dBGC['data'])
	id0=copy.deepcopy(dBGC['id'])
	lab0=copy.deepcopy(dBGC['code'])
	
	typ='CN'
	d0[typ]['Ctot L t0 2']={}
	d0[typ]['Ctot L t0 2']['mu']=d0[typ]['Cbk L t0']['mu']+d0[typ]['Cbr L t0']['mu']+d0[typ]['Cf L t0']['mu']+d0[typ]['Cr L t0']['mu']+d0[typ]['Csw L t0']['mu']
	d0[typ]['Ctot L t0 2']['se']=d0[typ]['Cbk L t0']['se']+d0[typ]['Cbr L t0']['se']+d0[typ]['Cf L t0']['se']+d0[typ]['Cr L t0']['se']+d0[typ]['Csw L t0']['se']
	
	ord=np.argsort(d0[typ]['Ctot L t0 2']['mu'])
	lab0=np.flip(lab0[ord])
	for v in d0[typ].keys():
		for k in d0[typ][v].keys():
			d0[typ][v][k]=np.flip(d0[typ][v][k][ord])
	lab2=np.append(lab0,'Area\nweighted')
	cl=np.array([[0.15,0.05,0.05],[0.45,0.3,0.3],[0.85,0.75,0.65],[0.6,0.8,0.4],[0.95,0.9,0.8]])
	#cl=np.array([[0.65,0.85,0.05],[0.75,0.5,0.95],[0.6,0.85,1],[0,0.6,0],[0.85,0.65,0.15]])
	bw=0.35; spc=0.6*bw;
	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(22,7));
	ax.bar(np.arange(id0.size)-spc,d0[typ]['Cr L t0']['mu'],bw,facecolor=cl[0,:])
	ax.bar(np.arange(id0.size)-spc,d0[typ]['Cbk L t0']['mu'],bw,bottom=d0[typ]['Cr L t0']['mu'],facecolor=cl[1,:])
	ax.bar(np.arange(id0.size)-spc,d0[typ]['Cbr L t0']['mu'],bw,bottom=d0[typ]['Cr L t0']['mu']+d0[typ]['Cbk L t0']['mu'],facecolor=cl[2,:])
	ax.bar(np.arange(id0.size)-spc,d0[typ]['Cf L t0']['mu'],bw,bottom=d0[typ]['Cr L t0']['mu']+d0[typ]['Cbk L t0']['mu']+d0[typ]['Cbr L t0']['mu'],facecolor=cl[3,:])
	ax.bar(np.arange(id0.size)-spc,d0[typ]['Csw L t0']['mu'],bw,bottom=d0[typ]['Cr L t0']['mu']+d0[typ]['Cbk L t0']['mu']+d0[typ]['Cbr L t0']['mu']+d0[typ]['Cf L t0']['mu'],facecolor=cl[4,:])
	for i in range(lab0.size):
		ax.errorbar(i-spc,d0[typ]['Ctot L t0 2']['mu'][i],yerr=d0[typ]['Ctot L t0 2']['se'][i],ecolor='k',fmt='none',capsize=2,lw=0.5)
	
	dW={}
	for k in d0[typ].keys():
		if k=='Area':
			continue
		ikp=np.where(np.isnan(d0[typ][k]['mu'])==False)[0]
		dW[k]=np.sum(d0[typ][k]['mu'][ikp]*d0[typ]['Area']['sum'][ikp])/np.sum(d0[typ]['Area']['sum'][ikp])
	totW=dW['Csw L t0']+dW['Cbk L t0']+dW['Cbr L t0']+dW['Cf L t0']+dW['Cr L t0']
	ax.bar(id0.size-spc,dW['Cr L t0'],bw,facecolor=cl[0,:],label='Roots (' + str(int(dW['Cr L t0']/totW*100)) + '%)')
	ax.bar(id0.size-spc,dW['Cbk L t0'],bw,bottom=dW['Cr L t0'],facecolor=cl[1,:],label='Bark (' + str(int(dW['Cbk L t0']/totW*100)) + '%)')
	ax.bar(id0.size-spc,dW['Cbr L t0'],bw,bottom=dW['Cr L t0']+dW['Cbk L t0'],facecolor=cl[2,:],label='Branches (' + str(int(dW['Cbr L t0']/totW*100)) + '%)')
	ax.bar(id0.size-spc,dW['Cf L t0'],bw,bottom=dW['Cr L t0']+dW['Cbk L t0']+dW['Cbr L t0'],facecolor=cl[3,:],label='Foliage (' + str(int(dW['Cf L t0']/totW*100)) + '%)')
	ax.bar(id0.size-spc,dW['Csw L t0'],bw,bottom=dW['Cr L t0']+dW['Cbk L t0']+dW['Cbr L t0']+dW['Cf L t0'],facecolor=cl[4,:],label='Stemwood (' + str(int(dW['Csw L t0']/totW*100)) + '%)')
	
	typ='VRI'
	d0[typ]['Ctot L t0 2']={}
	d0[typ]['Ctot L t0 2']['mu']=d0[typ]['Cbk L t0']['mu']+d0[typ]['Cbr L t0']['mu']+d0[typ]['Cf L t0']['mu']+d0[typ]['Cr L t0']['mu']+d0[typ]['Csw L t0']['mu']
	d0[typ]['Ctot L t0 2']['se']=d0[typ]['Cbk L t0']['se']+d0[typ]['Cbr L t0']['se']+d0[typ]['Cf L t0']['se']+d0[typ]['Cr L t0']['se']+d0[typ]['Csw L t0']['se']
	ord=np.argsort(d0[typ]['Ctot L t0 2']['mu'])

	#lab=np.flip(lab0[ord])
	for v in d0[typ].keys():
		for k in d0[typ][v].keys():
			d0[typ][v][k]=np.flip(d0[typ][v][k][ord])
	#lab2=np.append(lab,'Area\nweighted')
	ax.bar(np.arange(id0.size)+spc,d0[typ]['Cr L t0']['mu'],bw,facecolor=cl[0,:])
	ax.bar(np.arange(id0.size)+spc,d0[typ]['Cbk L t0']['mu'],bw,bottom=d0[typ]['Cr L t0']['mu'],facecolor=cl[1,:])
	ax.bar(np.arange(id0.size)+spc,d0[typ]['Cbr L t0']['mu'],bw,bottom=d0[typ]['Cr L t0']['mu']+d0[typ]['Cbk L t0']['mu'],facecolor=cl[2,:])
	ax.bar(np.arange(id0.size)+spc,d0[typ]['Cf L t0']['mu'],bw,bottom=d0[typ]['Cr L t0']['mu']+d0[typ]['Cbk L t0']['mu']+d0[typ]['Cbr L t0']['mu'],facecolor=cl[3,:])
	ax.bar(np.arange(id0.size)+spc,d0[typ]['Csw L t0']['mu'],bw,bottom=d0[typ]['Cr L t0']['mu']+d0[typ]['Cbk L t0']['mu']+d0[typ]['Cbr L t0']['mu']+d0[typ]['Cf L t0']['mu'],facecolor=cl[4,:])
	for i in range(id0.size):
		ax.errorbar(i+spc,d0[typ]['Ctot L t0 2']['mu'][i],yerr=d0[typ]['Ctot L t0 2']['se'][i],ecolor='k',fmt='none',capsize=2,lw=0.5)
		dW={}
		for k in d0[typ].keys():
			if k=='Area':
				continue
			ikp=np.where(np.isnan(d0[typ][k]['mu'])==False)[0]
			dW[k]=np.sum(d0[typ][k]['mu'][ikp]*d0[typ]['Area']['sum'][ikp])/np.sum(d0[typ]['Area']['sum'][ikp])
	totW=dW['Csw L t0']+dW['Cbk L t0']+dW['Cbr L t0']+dW['Cf L t0']+dW['Cr L t0']
	ax.bar(id0.size+spc,dW['Cr L t0'],bw,facecolor=cl[0,:])
	ax.bar(id0.size+spc,dW['Cbk L t0'],bw,bottom=dW['Cr L t0'],facecolor=cl[1,:])
	ax.bar(id0.size+spc,dW['Cbr L t0'],bw,bottom=dW['Cr L t0']+dW['Cbk L t0'],facecolor=cl[2,:])
	ax.bar(id0.size+spc,dW['Cf L t0'],bw,bottom=dW['Cr L t0']+dW['Cbk L t0']+dW['Cbr L t0'],facecolor=cl[3,:])
	ax.bar(id0.size+spc,dW['Csw L t0'],bw,bottom=dW['Cr L t0']+dW['Cbk L t0']+dW['Cbr L t0']+dW['Cf L t0'],facecolor=cl[4,:])
	
	ax.set(yticks=np.arange(0,550,20),xticks=np.arange(id0.size+1),xticklabels=lab2,ylabel='Biomass (tC ha$^{-1}$)',xlim=[-0.5,id0.size+1-0.5],ylim=[0,180])
	plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	
	#for iU in range(u.size):
	#	ax.text(i,6,str(d[typ]['Csw L t0']['N'][iU].astype(int)),color='k',ha='center',fontsize=7,fontweight='normal')
	gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_BiomassByBGC_CNV','png',900)
	return

#%%
def GP_DeadWoodByBGC_CNV(meta,gpt,dBGC):
	d0=copy.deepcopy(dBGC['data'])
	id0=copy.deepcopy(dBGC['id'])
	lab0=copy.deepcopy(dBGC['code'])
	
	typ='CN'
	ord=np.argsort(d0[typ]['Ctot D t0']['mu'])
	lab0=np.flip(lab0[ord])
	for v in d0[typ].keys():
		for k in d0[typ][v].keys():
			d0[typ][v][k]=np.flip(np.nan_to_num(d0[typ][v][k][ord]))
	lab2=np.append(lab0,'Area\nweighted')
	cl=np.array([[0.7,0.8,0.95],[0.15,0.05,0.05]])

	bw=0.6; cler=[0.7,0.7,0.7]
	spc=0.0*bw; #spc=0.6*bw;
	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(22,5));
	ax.bar(np.arange(id0.size)-spc,d0[typ]['Ctot D t0']['mu']-d0[typ]['Ctot D Fallen t0']['mu'],bw,facecolor=cl[0,:])
	ax.bar(np.arange(id0.size)-spc,d0[typ]['Ctot D Fallen t0']['mu'],bw,bottom=d0[typ]['Ctot D t0']['mu']-d0[typ]['Ctot D Fallen t0']['mu'],facecolor=cl[1,:])
	for i in range(lab0.size):
		ax.errorbar(i-spc,d0[typ]['Ctot D t0']['mu'][i],yerr=d0[typ]['Ctot D t0']['se'][i],color='k',fmt='none',capsize=2,lw=0.5)
	
	dW={}
	for k in d0[typ].keys():
		if k=='Area':
			continue
		ikp=np.where(np.isnan(d0[typ][k]['mu'])==False)[0]
		dW[k]=np.sum(d0[typ][k]['mu'][ikp]*d0[typ]['Area']['sum'][ikp])/np.sum(d0[typ]['Area']['sum'][ikp])
	#totW=dW['Ctot D t0']
	ax.bar(id0.size-spc,dW['Ctot D t0']-dW['Ctot D Fallen t0'],bw,facecolor=cl[0,:],label='Standing (100%)')
	ax.bar(id0.size-spc,dW['Ctot D Fallen t0'],bw,bottom=dW['Ctot D t0']-dW['Ctot D Fallen t0'],facecolor=cl[1,:],label='Fallen (0%)')

# 	typ='VRI'
# 	for v in d0[typ].keys():
# 		for k in d0[typ][v].keys():
# 			d0[typ][v][k]=np.flip(np.nan_to_num(d0[typ][v][k][ord]))
# 	ax.bar(np.arange(id0.size)+spc,d0[typ]['Ctot D t0']['mu']-d0[typ]['Ctot D Fallen t0']['mu'],bw,facecolor=cl[0,:])
# 	ax.bar(np.arange(id0.size)+spc,d0[typ]['Ctot D Fallen t0']['mu'],bw,bottom=d0[typ]['Ctot D t0']['mu']-d0[typ]['Ctot D Fallen t0']['mu'],facecolor=cl[1,:])
# 	for i in range(id0.size):
# 		ax.errorbar(i+spc,d0[typ]['Ctot D t0']['mu'][i],yerr=d0[typ]['Ctot D t0']['se'][i],color=cler,fmt='none',capsize=2,lw=0.5)
# 		dW={}
# 		for k in d0[typ].keys():
# 			if k=='Area':
# 				continue
# 			ikp=np.where(np.isnan(d0[typ][k]['mu'])==False)[0]
# 			dW[k]=np.sum(d0[typ][k]['mu'][ikp]*d0[typ]['Area']['sum'][ikp])/np.sum(d0[typ]['Area']['sum'][ikp])
# 	totW=dW['Ctot D t0']
# 	ax.bar(id0.size+spc,dW['Ctot D t0'],bw,facecolor=cl[0,:])
# 	ax.bar(id0.size+spc,dW['Ctot D Fallen t0'],bw,bottom=dW['Ctot D t0']-dW['Ctot D Fallen t0'],facecolor=cl[1,:])
	ax.set(yticks=np.arange(0,550,10),xticks=np.arange(id0.size+1),xticklabels=lab2,ylabel='Dead wood (tC ha$^{-1}$)',xlim=[-0.5,id0.size+1-0.5],ylim=[0,40])
	plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_DeadWoodByBGC_CNV','png',900)
	return

#%% Plot total ecosystem carbon stratified by BGC zone
def GP_TEC_ByBGC_CN(meta,gpt,dBGC):
	d0=copy.deepcopy(dBGC['data'])
	id0=copy.deepcopy(dBGC['id'])
	lab0=copy.deepcopy(dBGC['code'])
	
	# Plot Total Ecosystem Carbon (mean/ha)
	ord=np.argsort(d0['Soil']['TEC']['mu'])
	for k in d0.keys():
		for v in d0[k].keys():
			for typ in d0[k][v].keys():
				d0[k][v][typ]=np.flip(d0[k][v][typ][ord])
	lab0=np.flip(lab0[ord])
	id0=np.flip(id0[ord])
	
	cl=np.array([ [0.45,0.3,0.3],[0.15,0.05,0.05],[0.85,0.75,0.65],[0.8,0,0],[0.6,1,0],[1,0.5,0],[0.45,1,1],[0.3,0.7,0] ])
	plt.close('all'); fig,ax=plt.subplots(2,1,figsize=gu.cm2inch(22,10))
	ax[0].bar(np.arange(id0.size),d0['Soil']['MIN_C_THA']['mu'],facecolor=cl[0,:],label='Mineral soil horizon (100 cm depth)')
	ax[0].bar(np.arange(id0.size),d0['Soil']['ORG_C_THA']['mu'],bottom=d0['Soil']['MIN_C_THA']['mu'],facecolor=cl[1,:],label='Organic soil horizon')
	ax[0].bar(np.arange(id0.size),d0['CN']['Ctot D t0']['mu'],bottom=d0['Soil']['ORG_C_THA']['mu']+d0['Soil']['MIN_C_THA']['mu'],facecolor=cl[2,:],label='Standing dead wood')
	y=d0['Soil']['MIN_C_THA']['mu']+d0['Soil']['ORG_C_THA']['mu']+d0['CN']['Ctot D t0']['mu']
	ax[0].bar(np.arange(id0.size),d0['CN']['Cr L t0']['mu'],bottom=y,facecolor=cl[3,:],label='Roots')
	ax[0].bar(np.arange(id0.size),d0['CN']['Cbr L t0']['mu'],bottom=y+d0['CN']['Cr L t0']['mu'],facecolor=cl[4,:],label='Bark')
	ax[0].bar(np.arange(id0.size),d0['CN']['Cbk L t0']['mu'],bottom=y+d0['CN']['Cr L t0']['mu']+d0['CN']['Cbr L t0']['mu'],facecolor=cl[5,:],label='Branches')
	ax[0].bar(np.arange(id0.size),d0['CN']['Cf L t0']['mu'],bottom=y+d0['CN']['Cr L t0']['mu']+d0['CN']['Cbr L t0']['mu']+d0['CN']['Cbk L t0']['mu'],facecolor=cl[6,:],label='Foliage')
	ax[0].bar(np.arange(id0.size),d0['CN']['Csw L t0']['mu'],bottom=y+d0['CN']['Cr L t0']['mu']+d0['CN']['Cbr L t0']['mu']+d0['CN']['Cbk L t0']['mu']+d0['CN']['Cf L t0']['mu'],facecolor=cl[7,:],label='Stemwood')
	#ax.errorbar(np.arange(id0.size),ds['mu SOC tot'],yerr=ds['se'],color=gp['cla'],fmt='none',capsize=2)
	#for i in range(id0.size):
	#ax.text(i,10,str(ds['N'][i].astype(int)),color=gp['cla'],ha='center',fontsize=8)
	ax[0].set(xticks=np.arange(id0.size),xticklabels=lab0,ylabel='Average stock (tC per hectare)',xlim=[-0.5,id0.size-0.5],ylim=[0,575])
	ax[0].legend(loc='upper right',frameon=False,facecolor=[1,1,1],labelspacing=0.25)
	ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].tick_params(length=meta['Graphics']['gp']['tickl'])
	#gu.PrintFig(meta['Paths']['Figs'] + '\\TEC Per-hectare Mean','png',900)
	
	# Plot Total Ecosystem Carbon (total)
	ord=np.argsort(d0['Soil']['TEC']['sum'])
	for k in d0.keys():
		for v in d0[k].keys():
			for typ in d0[k][v].keys():
				d0[k][v][typ]=np.flip(d0[k][v][typ][ord])
	lab0=np.flip(lab0[ord])
	#tot=np.sum(d0['Soil']['TEC']['sum'])
	ax[1].bar(np.arange(id0.size),d0['Soil']['MIN_C_THA']['sum'],facecolor=cl[0,:],label='Mineral soil horizon (100 cm depth)')
	ax[1].bar(np.arange(id0.size),d0['Soil']['ORG_C_THA']['sum'],bottom=d0['Soil']['MIN_C_THA']['sum'],facecolor=cl[1,:],label='Organic soil horizon')
	ax[1].bar(np.arange(id0.size),d0['CN']['Ctot D t0']['sum'],bottom=d0['Soil']['ORG_C_THA']['sum']+d0['Soil']['MIN_C_THA']['sum'],facecolor=cl[2,:],label='Standing + fallen dead wood')
	y=d0['Soil']['MIN_C_THA']['sum']+d0['Soil']['ORG_C_THA']['sum']+d0['CN']['Ctot D t0']['sum']
	ax[1].bar(np.arange(id0.size),d0['CN']['Cr L t0']['sum'],bottom=y,facecolor=cl[3,:],label='Roots')
	ax[1].bar(np.arange(id0.size),d0['CN']['Cbr L t0']['sum'],bottom=y+d0['CN']['Cr L t0']['sum'],facecolor=cl[4,:],label='Bark')
	ax[1].bar(np.arange(id0.size),d0['CN']['Cbk L t0']['sum'],bottom=y+d0['CN']['Cr L t0']['sum']+d0['CN']['Cbr L t0']['sum'],facecolor=cl[5,:],label='Branches')
	ax[1].bar(np.arange(id0.size),d0['CN']['Cf L t0']['sum'],bottom=y+d0['CN']['Cr L t0']['sum']+d0['CN']['Cbr L t0']['sum']+d0['CN']['Cbk L t0']['sum'],facecolor=cl[6,:],label='Foliage')
	ax[1].bar(np.arange(id0.size),d0['CN']['Csw L t0']['sum'],bottom=y+d0['CN']['Cr L t0']['sum']+d0['CN']['Cbr L t0']['sum']+d0['CN']['Cbk L t0']['sum']+d0['CN']['Cf L t0']['sum'],facecolor=cl[7,:],label='Stemwood')
	#ax[1].text(2,3.75,'Total = ' + str(np.round(tot,decimals=1)) + ' (billion tonnes of C)',fontsize=10)
	ax[1].set(xticks=np.arange(id0.size),xticklabels=lab0,ylabel='Total stock (billion tonnes C)',xlim=[-0.5,id0.size-0.5],ylim=[0,5])
	#plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
	ax[1].yaxis.set_ticks_position('both'); ax[1].xaxis.set_ticks_position('both'); ax[1].tick_params(length=meta['Graphics']['gp']['tickl'])
	plt.tight_layout()
	gu.axletters(ax,plt,0.0135,0.91,FontColor=meta['Graphics']['gp']['cla'],LetterStyle='NoPar',FontWeight='Bold')
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_TEC_ByBGC_CN','png',900)
	return

#%%
def GP_BiomassDynamics_Mean_CN(meta,gpt,dBGC):
	d0=copy.deepcopy(dBGC['data'])
	typ='CN'
	
	mu={}
	se={}
	for v in d0[typ].keys():
		ikp=np.where(d0[typ][v]['N']>=0)
		mu[v]=np.nansum(d0[typ][v]['mu'][ikp]*d0[typ]['Area']['sum'][ikp])/np.sum(d0[typ]['Area']['sum'][ikp])
		se[v]=np.nansum(d0[typ][v]['se'][ikp]*d0[typ]['Area']['sum'][ikp])/np.sum(d0[typ]['Area']['sum'][ikp])
	
	lab=['Survivor\ngrowth','Recruitment\ngrowth','Natural\nmortality','Harvest\nmortality','Net\ngrowth','Litterfall','Net primary\nproduction'] #
	clf=[0.7,0.8,0.95]; cle=[0,0,0]; barw=0.6
	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(12,6));
	yt=np.arange(-2,5,0.5); yl=np.array([-2,3.5])
	ax.plot([0,10],[0,0],'-k',color=meta['Graphics']['gp']['cla'],lw=0.5)
	ax.bar(1,mu['Ctot G Surv'],barw,facecolor=clf,label='Growth survivors')
	ax.bar(2,mu['Ctot G Recr'],barw,facecolor=clf,label='Growth recruitment')
	ax.bar(3,-mu['Ctot Mort Nat'],barw,facecolor=clf,label='Natural\nmortality')
	ax.bar(4,-mu['Ctot Mort Harv'],barw,facecolor=clf,label='Harvest')
	ax.bar(5,mu['Ctot Net'],barw,facecolor=clf,label='Net')
	ax.bar(6,-mu['Ctot Litf'],barw,facecolor=clf,label='Litterfall')
	ax.bar(7,mu['Ctot G Surv']+mu['Ctot G Recr']+mu['Ctot Litf'],barw,facecolor=clf,label='NPP')

	ax.errorbar(1,mu['Ctot G Surv'],yerr=se['Ctot G Surv'],color=cle,fmt='none',capsize=2,lw=0.5)
	ax.errorbar(2,mu['Ctot G Recr'],yerr=se['Ctot G Recr'],color=cle,fmt='none',capsize=2,lw=0.5)
	ax.errorbar(3,-mu['Ctot Mort Nat'],yerr=se['Ctot Mort Nat'],color=cle,fmt='none',capsize=2,lw=0.5)
	ax.errorbar(4,-mu['Ctot Mort Harv'],yerr=se['Ctot Mort Harv'],color=cle,fmt='none',capsize=2,lw=0.5)
	ax.errorbar(5,mu['Ctot Net'],yerr=se['Ctot Net'],color=cle,fmt='none',capsize=2,lw=0.5)
	ax.errorbar(6,-mu['Ctot Litf'],yerr=se['Ctot Litf'],color=cle,fmt='none',capsize=2,lw=0.5)
	ax.errorbar(7,mu['Ctot G Surv']+mu['Ctot G Recr']+mu['Ctot Litf'],yerr=se['Ctot G Surv']+se['Ctot G Recr']+se['Ctot Litf'],color=cle,fmt='none',capsize=2,lw=0.5)
	ax.set(xticks=np.arange(1,len(lab)+1),xticklabels=lab,ylabel='Carbon balance of trees (tC ha$^{-1}$ yr$^{-1}$)',xlim=[0.5,7.5],yticks=yt,ylim=yl)
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	ax2=ax.twinx()
	ax2.set_ylabel('Carbon balance of trees (MtCO$_{2}$e yr$^{-1}$)')
	A=62000000 # Consistent with LC2 layer from BC LSCS 2023
	yt2=A*yl/1e6*3.667; #yt2.astype('int16')
	ax2.set(yticks=np.arange(-500,1500,100),ylim=[yt2[0],yt2[-1]])
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_BiomassDynamics_Mean_CN','png',900)
	return

#%%
def GP_AgeResponsesByRegion(meta,dAC):
	ptf='PTF CN'
	lw=0.5; ms=3; bw=17;
	plt.close('all'); fig,ax1=plt.subplots(1,2,figsize=gu.cm2inch(22,6))
	reg='Coast'
	ax1[0].errorbar(dAC['bin'],dAC['data'][ptf][reg]['Ctot L t0']['mu'],yerr=dAC['data'][ptf][reg]['Ctot L t0']['se'],color='k',fmt='none',capsize=1.5,lw=0.25)
	ax1[0].plot(dAC['bin'],dAC['data'][ptf][reg]['Ctot L t0']['mu'],'-ko',ms=ms,lw=lw,mew=lw,color='k',mfc='w',mec='k',label='Biomass',zorder=1)
	ax1[0].set(ylabel='Biomass (tC ha$^-$$^1$)',xlabel='Age, years',xticks=np.arange(0,400,dAC['bw']),yticks=np.arange(0,300,20),xlim=[0,250+bw],ylim=[0,240])
	ax1[0].yaxis.set_ticks_position('both'); ax1[0].xaxis.set_ticks_position('both'); ax1[0].tick_params(length=meta['Graphics']['gp']['tickl'])
	ax2=ax1[0].twinx()

	ax2.bar(dAC['bin'],dAC['data'][ptf][reg]['Ctot Net']['mu'],0.9*bw,fc=[0.7,0.8,0.95],ec='none',label='Net growth',zorder=-1)
	ax2.plot([0,500],[0,0],'-k',lw=lw)
	ax2.set(ylabel='Net growth (tC ha$^-$$^1$ yr$^-$$^1$)',xlabel='',ylim=[-1.75,5.5])
	ax1[0].set_zorder(ax2.get_zorder()+1)
	ax1[0].patch.set_visible(False)
	ax2.tick_params(length=meta['Graphics']['gp']['tickl'])

	reg='Interior'
	ax1[1].errorbar(dAC['bin'],dAC['data'][ptf][reg]['Ctot L t0']['mu'],yerr=dAC['data'][ptf][reg]['Ctot L t0']['se'],color='k',fmt='none',capsize=1.5,lw=0.25)
	ax1[1].plot(dAC['bin'],dAC['data'][ptf][reg]['Ctot L t0']['mu'],'-ko',ms=ms,lw=lw,mew=lw,color='k',mfc='w',mec='k',label='CMI/NFI',zorder=1)
	ax1[1].set(ylabel='Biomass (tC ha$^-$$^1$)',xlabel='Age, years',xticks=np.arange(0,400,dAC['bw']),yticks=np.arange(0,300,20),xlim=[0,250+bw],ylim=[0,240])
	ax1[1].yaxis.set_ticks_position('both'); ax1[1].xaxis.set_ticks_position('both'); ax1[1].tick_params(length=meta['Graphics']['gp']['tickl'])
	ax3=ax1[1].twinx()
	ax3.bar(dAC['bin'],dAC['data'][ptf][reg]['Ctot Net']['mu'],0.9*bw,fc=[0.7,0.8,0.95],ec='none',label='Net growth CMI/NFI',zorder=-1)
	ax3.plot([0,500],[0,0],'-k',lw=lw)
	ax3.set(ylabel='Net growth (tC ha$^-$$^1$ yr$^-$$^1$)',xlabel='',ylim=[-1.75,5.5])
	ax1[1].set_zorder(ax3.get_zorder()+1)
	ax1[1].patch.set_visible(False)
	ax3.tick_params(length=meta['Graphics']['gp']['tickl'])
	#ax1[1].legend(loc='upper right',frameon=False,facecolor=None,edgecolor='w')
	gu.axletters(ax1,plt,0.03,0.93,FontColor=meta['Graphics']['gp']['cla'],LetterStyle='NoPar',FontWeight='Bold',Labels=['Coast','Interior'],LabelSpacer=0.025)
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_AgeResponsesByRegion','png',900)
	return

#%%
def GP_MortalityAssociations(meta,gpt):
	# Estimate harvest mortality
	gpt['Ctot Mort+Lost Harv']=np.zeros(gpt['PTF CN'].size)
	ind=np.where( (gpt['Occ_Harv']==1) )[0]
	gpt['Ctot Mort+Lost Harv'][ind]=gpt['Ctot Mort+Lost'][ind]
	
	gpt['N Mort+Lost Harv']=np.zeros(gpt['PTF CN'].size)
	ind=np.where( (gpt['Occ_Harv']==1) )[0]
	gpt['N Mort+Lost Harv'][ind]=gpt['N Mort+Lost'][ind]
	
	daL=['No damage','Beetles','Defoliators','Diseases','Plant competition','Animal browsing','Fire','Frost, snow, ice, hail','Water stress','Wind','Silviculture','Defects and breakage','Flooding, lightning, slides','Unknown','Harv']
	daL_lab=['No\ndamage\ndetected','Beetles','Defoliators','Diseases','Plant\ncompetition','Animal\nbrowsing','Fire','Frost,\nsnow, ice\n and hail','Water\nstress','Wind', 'Silvi-\nculture','Defects\nand\nbreakage','Flooding\nlightning\nslides','Damage\nof uknown\ncause','Harvest']
	
	# Make sure that all mortality is accounted for by damage agents, assign unaccounted mortality to no damage
	ind=np.where( (gpt['PTF CN']==1) & (gpt['Ctot Mort+Lost']>0) & (gpt['Ctot Mort+Lost No damage']==0) & (gpt['Ctot Mort+Lost Beetles']==0) & (gpt['Ctot Mort+Lost Defoliators']==0) & \
	              (gpt['Ctot Mort+Lost Diseases']==0) & (gpt['Ctot Mort+Lost Plant competition']==0) & (gpt['Ctot Mort+Lost Animal browsing']==0) & \
	              (gpt['Ctot Mort+Lost Fire']==0) & (gpt['Ctot Mort+Lost Frost, snow, ice, hail']==0) & (gpt['Ctot Mort+Lost Water stress']==0) & \
	              (gpt['Ctot Mort+Lost Silviculture']==0) & (gpt['Ctot Mort+Lost Defects and breakage']==0) & (gpt['Ctot Mort+Lost Flooding, lightning, slides']==0) & \
	              (gpt['Ctot Mort+Lost Unknown']==0) & (gpt['Ctot Mort+Lost Harv']==0) )[0]
	gpt['Ctot Mort+Lost No damage'][ind]=gpt['Ctot Mort+Lost'][ind]
	
	# Make sure that all mortality is accounted for by damage agents, assign unaccounted mortality to no damage
	ind=np.where( (gpt['PTF CN']==1) & (gpt['N Mort+Lost']>0) & (gpt['N Mort+Lost No damage']==0) & (gpt['N Mort+Lost Beetles']==0) & (gpt['N Mort+Lost Defoliators']==0) & \
	              (gpt['N Mort+Lost Diseases']==0) & (gpt['N Mort+Lost Plant competition']==0) & (gpt['N Mort+Lost Animal browsing']==0) & \
	              (gpt['N Mort+Lost Fire']==0) & (gpt['N Mort+Lost Frost, snow, ice, hail']==0) & (gpt['N Mort+Lost Water stress']==0) & \
	              (gpt['N Mort+Lost Silviculture']==0) & (gpt['N Mort+Lost Defects and breakage']==0) & (gpt['N Mort+Lost Flooding, lightning, slides']==0) & \
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
	#print(np.sum(d['Sum'])/d['Sum All'])
	
	mup=d['Mean']/np.sum(d['Mean'])*100
	mup2=d['N Mean']/np.sum(d['N Mean'])*100
	#mup=d['Sum']#/np.sum(d['Mean'])*100
	ord=np.flip(np.argsort(mup))
	mup=mup[ord]
	mup2=mup2[ord]
	lab=np.array(daL_lab)[ord]
	
	plt.close('all');
	fig,ax=plt.subplots(1,figsize=gu.cm2inch(18,5.5)); barw=0.35
	ax.bar(np.arange(len(daL))-barw/1.85,mup,barw,fc=[0.7,0.8,0.95],ec=None,label='Biomass loss due to mortality')
	ax.bar(np.arange(len(daL))+barw/1.85,mup2,barw,fc=[0.85,1,0.5],ec=None,label='Number of trees')
	ax.set(xticks=np.arange(0,len(lab)),xticklabels=lab,ylabel='Frequency (%)',xlim=[-0.5,len(lab)-1+0.5])# ,ylim=[0,35]
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	plt.legend(loc='upper right',frameon=False,facecolor=[1,1,1],labelspacing=0.25)
	plt.tight_layout()
	if meta['Graphics']['Print Figures']=='On':
		gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_MortalityAssociations_CN','png',900)
	return

#%%