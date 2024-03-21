
#%% Import modules

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import gdal
import pyproj
import geopandas as gpd
from scipy.optimize import curve_fit
import scipy.io as spio
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats, linalg
from subprocess import call
from fcgadgets.macgyver import utilities_general as gu
import warnings
warnings.filterwarnings('ignore')

#%% Import data

PathProject=r'C:\Users\rhember\Documents\Data\EP703'
PathFigures=r'G:\My Drive\Figures\EP703'

# Import site information
dfS=pd.read_excel(PathProject + '\\EP703_SiteInfo.xlsx',sheet_name='Sheet1')

# Unique sites
uSite=dfS['ID_Instal'].unique()

# Import stand-level data
sobs=gu.ipickle(PathProject + '\\Processed\\EP703_SL.pkl')
#sobs=pd.read_excel(PathProject + '\\Processed\\EP703_SL.xlsx')

# Import tree-level data
tobs=gu.ipickle(PathProject + '\\Processed\\EP703_TL.pkl')
#tobs=pd.read_excel(PathProject + '\\Processed\\EP703_TOBS.xlsx')

#%% Add derived variables, adjustments

# Create an adjusted value of TSF_t0 that will ensure time intervals are not double
# counting the year of measurement
tobs['TSF_t0_adj']=tobs['TSF_t0']+1

# Add an indicator of leading species to tree level structure 
tobs['SpcLead']=np.zeros(tobs['ID_Site'].size)
uSite=np.unique(sobs['ID_Site'])
Nhw=0
for i in range(uSite.size):
	indS=np.where(sobs['ID_Site']==uSite[i])[0]
	indT=np.where(tobs['ID_Site']==uSite[i])[0]
	if sobs['Spc1_ID_t0'][indS[0]]=='FD':
		tobs['SpcLead'][indT]=1
	else:
		tobs['SpcLead'][indT]=2
		Nhw=Nhw+1

tobs['RGR BA']=(np.log(tobs['BA_t1'])-np.log(tobs['BA_t0']))/tobs['DT']

#%% Individual plots


uCC=np.array([1,2,3,4])

DA_Pm=np.zeros(uCC.size)
DA_Bsw_G=np.zeros(uCC.size)
DA_BA_G=np.zeros(uCC.size)
DA_BA_RGR=np.zeros(uCC.size)
DA_H_obs_G=np.zeros(uCC.size)

DR_Pm=np.zeros(uCC.size)
DR_Bsw_G=np.zeros(uCC.size)
DR_BA_G=np.zeros(uCC.size)
DR_BA_RGR=np.zeros(uCC.size)
DR_H_obs_G=np.zeros(uCC.size)
for iCC in range(uCC.size):
	iC=np.where( (tobs['N_Dose']==0) & (tobs['TSF_t0']>=0) & (tobs['TSF_t1']<=18) & (tobs['CrownClass_t0']==uCC[iCC]) & (tobs['ID_Spc_t0']==2) )[0]
	iF=np.where( (tobs['N_Dose']==225) & (tobs['TSF_t0']>=0) & (tobs['TSF_t1']<=18) & (tobs['CrownClass_t0']==uCC[iCC]) & (tobs['ID_Spc_t0']==2) )[0]
	Pm_C=np.nansum(tobs['Mort'][iC])/iC.size*100
	Pm_F=np.nansum(tobs['Mort'][iF])/iF.size*100
	Bsw_G_C=np.nanmean(tobs['Bsw_G'][iC])
	Bsw_G_F=np.nanmean(tobs['Bsw_G'][iF])
	BA_G_C=np.nanmean(tobs['BA_G'][iC])
	BA_G_F=np.nanmean(tobs['BA_G'][iF])
	BA_RGR_C=np.nanmean(tobs['RGR BA'][iC])
	BA_RGR_F=np.nanmean(tobs['RGR BA'][iF])
	H_obs_G_C=np.nanmean(tobs['H_obs_G'][iC])
	H_obs_G_F=np.nanmean(tobs['H_obs_G'][iF])
	
	DA_Pm[iCC]=Pm_F-Pm_C
	#G[iCC]=np.mean(tobs['RGR BA'][ind])
	DA_Bsw_G[iCC]=Bsw_G_F-Bsw_G_C
	DA_BA_G[iCC]=BA_G_F-BA_G_C
	DA_BA_RGR[iCC]=BA_RGR_F-BA_RGR_C
	DA_H_obs_G[iCC]=H_obs_G_F-H_obs_G_C
	
	DR_Pm[iCC]=(Pm_F-Pm_C)/Pm_C*100
	#G[iCC]=np.mean(tobs['RGR BA'][ind])
	DR_Bsw_G[iCC]=(Bsw_G_F-Bsw_G_C)/Bsw_G_C*100
	DR_BA_G[iCC]=(BA_G_F-BA_G_C)/BA_G_C*100
	DR_BA_RGR[iCC]=(BA_RGR_F-BA_RGR_C)/BA_RGR_C*100
	DR_H_obs_G[iCC]=(H_obs_G_F-H_obs_G_C)/H_obs_G_C*100

plt.close('all')
plt.figure(); plt.bar(uCC,DR_Pm,0.85)
plt.figure(); plt.bar(uCC,DR_BA_G,0.85)
#plt.figure(); plt.bar(uCC,DR_BA_RGR,0.85)
plt.figure(); plt.bar(uCC,DR_H_obs_G,0.85)

#%%



id_Inst=3
ind=np.where( (tobs['ID_Site']==id_Inst) & (tobs['N_Dose']==0) )[0]
uP=np.unique(tobs['ID_Plot'][ind])

uCC=np.array([1,2,3,4])

M=np.zeros(uCC.size)
G=np.zeros(uCC.size)
for iCC in range(uCC.size):
	ind=np.where( (tobs['ID_Site']==id_Inst) & (tobs['ID_Plot']==uP[1]) & (tobs['TSF_t0']==0) & (tobs['CrownClass_t0']==uCC[iCC]) )[0]
	M[iCC]=np.sum(tobs['Mort'][ind])
	#G[iCC]=np.mean(tobs['RGR BA'][ind])
	G[iCC]=np.mean(tobs['Bsw_G'][ind])

plt.close('all')
plt.figure(); plt.bar(uCC,M,0.85)
plt.figure(); plt.bar(uCC,G,0.85)

#%%

ind=np.where( (tobs['ID_Site']==id_Inst) & (tobs['ID_Plot']==uP[1]) & (tobs['TSF_t0']==0) )[0]
tobs['Mort'][ind]
tobs['Bsw_G'][ind]
plt.close('all')
plt.hist(tobs['Bsw_t0'][ind])
plt.hist(tobs['H_obs_G'][ind])

#%% Set graphics parameters
		
params={'font.sans-serif':'Calibri',
		'font.size':7,
		'axes.edgecolor':'black',
		'axes.labelsize':7,
		'axes.labelcolor':'black',
		'axes.titlesize':7,
		'axes.linewidth':0.5,		
		'lines.linewidth':0.5,
		'text.color':'black',
		'xtick.color':'black',		
		'xtick.labelsize':7,
		'xtick.major.width':0.5,
		'xtick.major.size':3,
		'xtick.direction':'in',
		'ytick.color':'black',
		'ytick.labelsize':7,
		'ytick.major.width':0.5,
		'ytick.major.size':3,
		'ytick.direction':'in',
		'legend.fontsize':7,		
		'savefig.dpi':300,
		'savefig.transparent':True}
plt.rcParams.update(params)

#%% Histograms of individual-tree data

plt.close('all'); ms=2; Alpha=0.14
fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(15.5,6.25))

ind0=np.where( (tobs['Mort']==0) & (tobs['Rec']==0) )[0]
ind1=np.where( (tobs['Mort']==0) & (tobs['Rec']==0) & (tobs['BA_G']<0) )[0]
ax[0].hist(tobs['BA_G'][ind0],np.arange(-20,100,5))
ax[0].set(position=[0.085,0.12,0.41,0.84],xlabel='Basal area growth (m2/yr)',ylabel='Frequency')

ind0=np.where( (tobs['Mort']==0) & (tobs['Rec']==0) )[0]
ind1=np.where( (tobs['Mort']==0) & (tobs['Rec']==0) & (tobs['H_gf_G']<0) )[0]
ind1.size/np.sum(ind0.size+ind1.size)
ax[1].hist(tobs['H_gf_G'][ind0],np.arange(-1,2,0.05))
ax[1].set(position=[0.58,0.12,0.41,0.84],xlabel='Height growth (m/yr)',ylabel='Frequency')

gu.axletters(ax,plt,0.03,0.92)
plt.savefig(PathFigures + '\\HistogramsIndTree.png',format='png',dpi=900)


# Height growth
plt.close('all'); ms=2; Alpha=0.14
fig,ax=plt.subplots(2,1,figsize=gu.cm2inch(10,12))
ind0=np.where( (tobs['Mort']==0) & (tobs['Rec']==0) )[0]
ax[0].hist(tobs['H_obs_G'][ind0],np.arange(-1,2,0.05))
ax[0].set(xlabel='Height growth, obs (m/yr)',ylabel='Frequency')
ax[1].hist(tobs['H_gf_G'][ind0],np.arange(-1,2,0.05))
ax[1].set(xlabel='Height growth, gf (m/yr)',ylabel='Frequency')


#%% Plot timeline of height measurements

uSpc=[1,2]
uDose=[0,225,450]
uSite=np.unique(tobs['ID_Site'])
uTSF_t0=np.unique(tobs['TSF_t0'])
Data={}
for iSpc in range(len(uSpc)):
	Data[uSpc[iSpc]]={}
	for iDose in range(len(uDose)):
		Data0=np.nan*np.ones((7000,uTSF_t0.size))
		cnt=0
		for iSite in range(uSite.size):
			indSite=np.where( (tobs['ID_Site']==uSite[iSite]) & (tobs['ID_Spc_t0']==uSpc[iSpc]) & (tobs['N_Dose']==uDose[iDose]) )[0]
			uPlot=np.unique(tobs['ID_Plot'][indSite])
			for iPlot in range(uPlot.size):
				indPlot=np.where( (tobs['ID_Site']==uSite[iSite]) & (tobs['ID_Spc_t0']==uSpc[iSpc]) & (tobs['N_Dose']==uDose[iDose]) & (tobs['ID_Plot']==uPlot[iPlot]) & (tobs['H_obs_t0']>0) )[0]
				uTree=np.unique(tobs['ID_Tree'][indPlot])
				for iTree in range(uTree.size):
					indTree=np.where( (tobs['ID_Site']==uSite[iSite]) & (tobs['ID_Spc_t0']==uSpc[iSpc]) & (tobs['N_Dose']==uDose[iDose]) & (tobs['ID_Plot']==uPlot[iPlot]) & (tobs['ID_Tree']==uTree[iTree]) & (tobs['H_obs_t0']>0) )[0]
					c,ia,ib=np.intersect1d(tobs['TSF_t0'][indTree],uTSF_t0,return_indices=True)
					x=tobs['TSF_t0'][indTree]
					y=cnt*np.ones(x.size)
					Data0[cnt,ib]=1
					cnt=cnt+1
		# Remove excess nans
		ikp=np.where(np.nansum(Data0,axis=1)>0)[0]
		Data0=Data0[ikp,:]		
		
		# Sort by start and then duration
		Data0S=np.zeros(Data0.shape)
		cnt=0
		for i in range(Data0.shape[1]):
			ind=np.where(Data0[:,i]>0)[0]
			D=np.zeros(ind.size)
			for j in range(ind.size):		
				t=uTSF_t0[Data0[ind[j]]==1]
				D[j]=np.max(t)-np.min(t)
			idx=np.argsort(-D)
			for j in range(ind.size):		
				Data0S[cnt,:]=Data0[ind[idx[j]],:]
				Data0[ind[idx[j]],:]=0
				cnt=cnt+1
		
		# Add sorted values to data structure
		Data[uSpc[iSpc]][uDose[iDose]]=Data0S
 
# Plot
tsf_max=18
plt.close('all'); 
fig,ax=plt.subplots(3,2,figsize=gu.cm2inch(15.5,14))
Spc=1; Dose=0; z=Data[Spc][Dose].copy()
ind=np.where( (np.nansum(z[:,0:7],axis=1)>0) )[0]; z=z[ind,:]
for i in range(z.shape[0]):
	ind=np.where(z[i,:]>0)[0]
	st=uTSF_t0[ind[0]]
	wdth=uTSF_t0[ind[-1]]-uTSF_t0[ind[0]]
	adj=np.maximum(0,st+wdth-tsf_max)
	ax[0,0].barh(i,wdth-adj,1,left=st,facecolor=[0.8,0.8,0.8],clip_on=True)
ax[0,0].set(position=[0.07,0.71,0.42,0.27],xlim=[-1.5,19],xticks=[-1,0,3,6,9,12,15,18,22],
	   xlabel='',ylabel='Trees with height')
Spc=2; Dose=0; z=Data[Spc][Dose].copy()
ind=np.where( (np.nansum(z[:,0:7],axis=1)>0) )[0]; z=z[ind,:]
for i in range(z.shape[0]):
	ind=np.where(z[i,:]>0)[0]
	st=uTSF_t0[ind[0]]
	wdth=uTSF_t0[ind[-1]]-uTSF_t0[ind[0]]
	adj=np.maximum(0,st+wdth-tsf_max)
	ax[0,1].barh(i,wdth-adj,1,left=st,facecolor=[0.8,0.8,0.8])
ax[0,1].set(position=[0.57,0.71,0.42,0.27],xlim=[-1.5,19],xticks=[-1,0,3,6,9,12,15,18,22],
	   xlabel='',ylabel='Trees with height')
Spc=1; Dose=225; z=Data[Spc][Dose].copy()
ind=np.where( (np.nansum(z[:,0:7],axis=1)>0) )[0]; z=z[ind,:]
for i in range(z.shape[0]):
	ind=np.where(z[i,:]>0)[0]
	st=uTSF_t0[ind[0]]
	wdth=uTSF_t0[ind[-1]]-uTSF_t0[ind[0]]
	adj=np.maximum(0,st+wdth-tsf_max)
	ax[1,0].barh(i,wdth-adj,1,left=st,facecolor=[0.8,0.8,0.8])
ax[1,0].set(position=[0.07,0.39,0.42,0.27],xlim=[-1.5,19],xticks=[-1,0,3,6,9,12,15,18,22],
	   xlabel='',ylabel='Trees with height')
Spc=2; Dose=225; z=Data[Spc][Dose].copy()
ind=np.where( (np.nansum(z[:,0:7],axis=1)>0) )[0]; z=z[ind,:]
for i in range(z.shape[0]):
	ind=np.where(z[i,:]>0)[0]
	st=uTSF_t0[ind[0]]
	wdth=uTSF_t0[ind[-1]]-uTSF_t0[ind[0]]
	adj=np.maximum(0,st+wdth-tsf_max)
	ax[1,1].barh(i,wdth-adj,1,left=st,facecolor=[0.8,0.8,0.8])
ax[1,1].set(position=[0.57,0.39,0.42,0.27],xlim=[-1.5,19],xticks=[-1,0,3,6,9,12,15,18,22],
	   xlabel='',ylabel='Trees with height')
Spc=1; Dose=450; z=Data[Spc][Dose].copy()
ind=np.where( (np.nansum(z[:,0:7],axis=1)>0) )[0]; z=z[ind,:]
for i in range(z.shape[0]):
	ind=np.where(z[i,:]>0)[0]
	st=uTSF_t0[ind[0]]
	wdth=uTSF_t0[ind[-1]]-uTSF_t0[ind[0]]
	adj=np.maximum(0,st+wdth-tsf_max)
	ax[2,0].barh(i,wdth-adj,1,left=st,facecolor=[0.8,0.8,0.8])
ax[2,0].set(position=[0.07,0.07,0.42,0.27],xlim=[-1.5,19],xticks=[-1,0,3,6,9,12,15,18,22],
	   xlabel='Time since application, years',ylabel='Trees with height')
Spc=2; Dose=450; z=Data[Spc][Dose].copy()
ind=np.where( (np.nansum(z[:,0:7],axis=1)>0) )[0]; z=z[ind,:]
for i in range(z.shape[0]):
	ind=np.where(z[i,:]>0)[0]
	st=uTSF_t0[ind[0]]
	wdth=uTSF_t0[ind[-1]]-uTSF_t0[ind[0]]
	adj=np.maximum(0,st+wdth-tsf_max)
	ax[2,1].barh(i,wdth-adj,1,left=st,facecolor=[0.8,0.8,0.8])
ax[2,1].set(position=[0.57,0.07,0.42,0.27],xlim=[-1.5,19],xticks=[-1,0,3,6,9,12,15,18,22],
	   xlabel='Time since application, years',ylabel='Trees with height')

gu.axletters(ax,plt,0.03,0.9,['Douglas-fir, Control','Western hemlock, Control','Douglas-fir, 225', 
							  'Western hemlock, 225','Douglas-fir, 450','Western hemlock, 450'],0.05)
#plt.savefig(PathFigures + '\\TL_HeightVsTime.png',format='png',dpi=900)
#plt.savefig(PathFigures + '\\TL_HeightVsTime.eps',format='eps')
gu.PrintFig(PathFigures + '\\TL_HeightVsTime','png',900)

#%% Plot timeline of basal area measurements

uSpc=[1,2]
uDose=[0,225,450]
uSite=np.unique(tobs['ID_Site'])
uTSF_t0=np.unique(tobs['TSF_t0'])
Data={}
for iSpc in range(len(uSpc)):
	Data[uSpc[iSpc]]={}
	for iDose in range(len(uDose)):
		Data0=np.nan*np.ones((20000,uTSF_t0.size))
		cnt=0
		for iSite in range(uSite.size):
			indSite=np.where( (tobs['ID_Site']==uSite[iSite]) & (tobs['ID_Spc_t0']==uSpc[iSpc]) & (tobs['N_Dose']==uDose[iDose]) )[0]
			uPlot=np.unique(tobs['ID_Plot'][indSite])
			for iPlot in range(uPlot.size):
				indPlot=np.where( (tobs['ID_Site']==uSite[iSite]) & (tobs['ID_Spc_t0']==uSpc[iSpc]) & (tobs['N_Dose']==uDose[iDose]) & (tobs['ID_Plot']==uPlot[iPlot]) & (tobs['BA_t0']>0) )[0]
				uTree=np.unique(tobs['ID_Tree'][indPlot])
				for iTree in range(uTree.size):
					indTree=np.where( (tobs['ID_Site']==uSite[iSite]) & (tobs['ID_Spc_t0']==uSpc[iSpc]) & (tobs['N_Dose']==uDose[iDose]) & (tobs['ID_Plot']==uPlot[iPlot]) & (tobs['ID_Tree']==uTree[iTree]) & (tobs['BA_t0']>0) )[0]
					c,ia,ib=np.intersect1d(tobs['TSF_t0'][indTree],uTSF_t0,return_indices=True)
					x=tobs['TSF_t0'][indTree]
					y=cnt*np.ones(x.size)
					Data0[cnt,ib]=1
					cnt=cnt+1
		# Remove excess nans
		ikp=np.where(np.nansum(Data0,axis=1)>0)[0]
		Data0=Data0[ikp,:]		
		
		# Sort by start and then duration
		Data0S=np.zeros(Data0.shape)
		cnt=0
		for i in range(Data0.shape[1]):
			ind=np.where(Data0[:,i]>0)[0]
			D=np.zeros(ind.size)
			for j in range(ind.size):		
				t=uTSF_t0[Data0[ind[j]]==1]
				D[j]=np.max(t)-np.min(t)
			idx=np.argsort(-D)
			for j in range(ind.size):		
				Data0S[cnt,:]=Data0[ind[idx[j]],:]
				Data0[ind[idx[j]],:]=0
				cnt=cnt+1
		
		# Add sorted values to data structure
		Data[uSpc[iSpc]][uDose[iDose]]=Data0S
 
# Plot
tsf_max=18
plt.close('all'); 
fig,ax=plt.subplots(3,2,figsize=gu.cm2inch(15.5,14))
Spc=1; Dose=0; z=Data[Spc][Dose].copy()
ind=np.where( (np.nansum(z[:,0:7],axis=1)>0) )[0]; z=z[ind,:]
for i in range(z.shape[0]):
	ind=np.where(z[i,:]>0)[0]
	st=uTSF_t0[ind[0]]
	wdth=uTSF_t0[ind[-1]]-uTSF_t0[ind[0]]
	adj=np.maximum(0,st+wdth-tsf_max)
	ax[0,0].barh(i,wdth-adj,1,left=st,facecolor=[0.8,0.8,0.8],clip_on=True)
ax[0,0].set(position=[0.07,0.71,0.42,0.27],xlim=[-1.5,19],xticks=[-1,0,3,6,9,12,15,18,22],
	   xlabel='',ylabel='Trees with basal area')
Spc=2; Dose=0; z=Data[Spc][Dose].copy()
ind=np.where( (np.nansum(z[:,0:7],axis=1)>0) )[0]; z=z[ind,:]
for i in range(z.shape[0]):
	ind=np.where(z[i,:]>0)[0]
	st=uTSF_t0[ind[0]]
	wdth=uTSF_t0[ind[-1]]-uTSF_t0[ind[0]]
	adj=np.maximum(0,st+wdth-tsf_max)
	ax[0,1].barh(i,wdth-adj,1,left=st,facecolor=[0.8,0.8,0.8])
ax[0,1].set(position=[0.57,0.71,0.42,0.27],xlim=[-1.5,19],xticks=[-1,0,3,6,9,12,15,18,22],
	   xlabel='',ylabel='Trees with basal area')
Spc=1; Dose=225; z=Data[Spc][Dose].copy()
ind=np.where( (np.nansum(z[:,0:7],axis=1)>0) )[0]; z=z[ind,:]
for i in range(z.shape[0]):
	ind=np.where(z[i,:]>0)[0]
	st=uTSF_t0[ind[0]]
	wdth=uTSF_t0[ind[-1]]-uTSF_t0[ind[0]]
	adj=np.maximum(0,st+wdth-tsf_max)
	ax[1,0].barh(i,wdth-adj,1,left=st,facecolor=[0.8,0.8,0.8])
ax[1,0].set(position=[0.07,0.39,0.42,0.27],xlim=[-1.5,19],xticks=[-1,0,3,6,9,12,15,18,22],
	   xlabel='',ylabel='Trees with basal area')
Spc=2; Dose=225; z=Data[Spc][Dose].copy()
ind=np.where( (np.nansum(z[:,0:7],axis=1)>0) )[0]; z=z[ind,:]
for i in range(z.shape[0]):
	ind=np.where(z[i,:]>0)[0]
	st=uTSF_t0[ind[0]]
	wdth=uTSF_t0[ind[-1]]-uTSF_t0[ind[0]]
	adj=np.maximum(0,st+wdth-tsf_max)
	ax[1,1].barh(i,wdth-adj,1,left=st,facecolor=[0.8,0.8,0.8])
ax[1,1].set(position=[0.57,0.39,0.42,0.27],xlim=[-1.5,19],xticks=[-1,0,3,6,9,12,15,18,22],
	   xlabel='',ylabel='Trees with basal area')
Spc=1; Dose=450; z=Data[Spc][Dose].copy()
ind=np.where( (np.nansum(z[:,0:7],axis=1)>0) )[0]; z=z[ind,:]
for i in range(z.shape[0]):
	ind=np.where(z[i,:]>0)[0]
	st=uTSF_t0[ind[0]]
	wdth=uTSF_t0[ind[-1]]-uTSF_t0[ind[0]]
	adj=np.maximum(0,st+wdth-tsf_max)
	ax[2,0].barh(i,wdth-adj,1,left=st,facecolor=[0.8,0.8,0.8])
ax[2,0].set(position=[0.07,0.07,0.42,0.27],xlim=[-1.5,19],xticks=[-1,0,3,6,9,12,15,18,22],
	   xlabel='Time since application, years',ylabel='Trees basal area')
Spc=2; Dose=450; z=Data[Spc][Dose].copy()
ind=np.where( (np.nansum(z[:,0:7],axis=1)>0) )[0]; z=z[ind,:]
for i in range(z.shape[0]):
	ind=np.where(z[i,:]>0)[0]
	st=uTSF_t0[ind[0]]
	wdth=uTSF_t0[ind[-1]]-uTSF_t0[ind[0]]
	adj=np.maximum(0,st+wdth-tsf_max)
	ax[2,1].barh(i,wdth-adj,1,left=st,facecolor=[0.8,0.8,0.8])
ax[2,1].set(position=[0.57,0.07,0.42,0.27],xlim=[-1.5,19],xticks=[-1,0,3,6,9,12,15,18,22],
	   xlabel='Time since application, years',ylabel='Trees basal area')

gu.axletters(ax,plt,0.03,0.9,['FDC, C Plots','HW, C Plots','FDC, 225N Plots','HW, 225N Plots','FDC, 450N Plots','HW, 450N Plots'],0.05)
#plt.savefig(PathFigures + '\\TL_HeightVsTime.png',format='png',dpi=900)
#plt.savefig(PathFigures + '\\TL_HeightVsTime.eps',format='eps')
gu.PrintFig(PathFigures + '\\TL_BasalAreaVsTime','emf',dpi=600)


#%% Response of survivor growth to N addition

# Some sites tracked surviving trees between TSF=-1 and TSF=6.
# After looking at it a bunch of different ways, I decided to extract the stats
# for those specific trees. 

uSpc=[1,2]
uDose=[225,450,675]
uTSF=[3,6,9,12,18]
vars=['ID_Site','Year_t0','Year_t1','TSF_t0','TSF_t1',
	  'SA_t0_C','SN_t0_C','N_H_C','BA_t0_C','BA_G_C','H_obs_t0_C','H_obs_G_C','H_mod_G_C','Bsw_obs_t0_C','Bsw_obs_G_C','Bsw_mod_G_C',
	  'SA_t0_F','SN_t0_F','N_H_F','BA_t0_F','BA_G_F','H_obs_t0_F','H_obs_G_F','H_mod_G_F','Bsw_obs_t0_F','Bsw_obs_G_F','Bsw_mod_G_F']

rs={}
for iSpc in range(len(uSpc)):
	rs[uSpc[iSpc]]={}
	for iDose in range(len(uDose)):
		rs[uSpc[iSpc]][uDose[iDose]]={}
		
		for iTSF in range(len(uTSF)):
			rs[uSpc[iSpc]][uDose[iDose]][uTSF[iTSF]]={}
			
			uSite=np.unique(tobs['ID_Site'][tobs['SpcLead']==uSpc[iSpc]])
			d={}
			for v in vars:
				d[v]=np.nan*np.ones(uSite.size)		
		
			for iSite in range(uSite.size):	
				# Index to first measurement
				# Use TSF=-1 for all but site 29 for western hemlock or TSF=0
				iC0=np.where((tobs['ID_Spc_t0']==uSpc[iSpc]) & (tobs['ID_Site']==uSite[iSite]) & (tobs['TSF_t0']==-1) & (tobs['N_Dose']==0) & (tobs['Mort']==0) & (tobs['Rec']==0) )[0]
				iF0=np.where((tobs['ID_Spc_t0']==uSpc[iSpc]) & (tobs['ID_Site']==uSite[iSite]) & (tobs['TSF_t0']==-1) & (tobs['N_Dose']==uDose[iDose]) & (tobs['Mort']==0) & (tobs['Rec']==0) )[0]
				iC0b=np.where((tobs['ID_Spc_t0']==uSpc[iSpc]) & (tobs['ID_Site']==uSite[iSite]) & (tobs['TSF_t0']==0) & (tobs['N_Dose']==0) & (tobs['Mort']==0) & (tobs['Rec']==0) )[0]
				iF0b=np.where((tobs['ID_Spc_t0']==uSpc[iSpc]) & (tobs['ID_Site']==uSite[iSite]) & (tobs['TSF_t0']==0) & (tobs['N_Dose']==uDose[iDose]) & (tobs['Mort']==0) & (tobs['Rec']==0) )[0]
				if iC0b.size>iC0.size: 
					iC0=iC0b
					iF0=iF0b
				# Only continue if there are data
				if (iC0.size==0) | (iF0.size==0): 
					continue
				# Index to second measurement
				iC1=np.where((tobs['ID_Spc_t0']==uSpc[iSpc]) & (tobs['ID_Site']==uSite[iSite]) & (tobs['TSF_t0']==uTSF[iTSF]) & (tobs['N_Dose']==0) & (tobs['Mort']==0) & (tobs['Rec']==0) )[0]
				iF1=np.where((tobs['ID_Spc_t0']==uSpc[iSpc]) & (tobs['ID_Site']==uSite[iSite]) & (tobs['TSF_t0']==uTSF[iTSF]) & (tobs['N_Dose']==uDose[iDose]) & (tobs['Mort']==0) & (tobs['Rec']==0) )[0]
				if (iC1.size==0) | (iF1.size==0): 
					continue
				# Time interval
				dt=tobs['Year_t0'][iC1[0]]-tobs['Year_t0'][iC0[0]]
				
				# Index to survivors
				cC,iaC,ibC=np.intersect1d(tobs['ID_Tree'][iC0],tobs['ID_Tree'][iC1],return_indices=True)
				cF,iaF,ibF=np.intersect1d(tobs['ID_Tree'][iF0],tobs['ID_Tree'][iF1],return_indices=True)
	
				# Variables for control
				C_BA_t0=tobs['BA_t0'][iC0[iaC]]
				C_BA_G=(tobs['BA_t0'][iC1[ibC]]-tobs['BA_t0'][iC0[iaC]])/dt
				C_H_obs_t0=tobs['H_obs_t0'][iC0[iaC]]
				C_H_obs_G=(tobs['H_obs_t0'][iC1[ibC]]-tobs['H_obs_t0'][iC0[iaC]])/dt
				C_H_mod_G=(tobs['H_mod_t0'][iC1[ibC]]-tobs['H_mod_t0'][iC0[iaC]])/dt
				C_Bsw_t0=tobs['Bsw_t0'][iC0[iaC]]
				C_Bsw_G=(tobs['Bsw_t0'][iC1[ibC]]-tobs['Bsw_t0'][iC0[iaC]])/dt
				C_Bsw_mod_G=(tobs['Bsw_mod_t0'][iC1[ibC]]-tobs['Bsw_mod_t0'][iC0[iaC]])/dt
				
				# Make sure modelled response is based on the same sample as observed
				ind=np.where( np.isnan(C_H_obs_G)==True )[0]
				C_H_mod_G[ind]=np.nan
				C_Bsw_mod_G[ind]=np.nan
				
				# Remove bad data				
				ind=np.where( (C_BA_t0<0) )[0]
				C_BA_t0[ind]=np.nan
				C_BA_G[ind]=np.nan
				C_H_obs_t0[ind]=np.nan
				C_H_obs_G[ind]=np.nan				
				C_Bsw_t0[ind]=np.nan
				C_Bsw_G[ind]=np.nan
				C_Bsw_mod_G[ind]=np.nan
				
				# Remove modelled values from observed biomass
				C_Bsw_obs_t0=C_Bsw_t0.copy()
				C_Bsw_obs_G=C_Bsw_G.copy()
				ind=np.where(np.isnan(C_H_obs_G)==True)[0]
				C_Bsw_obs_t0[ind]=np.nan
				C_Bsw_obs_G[ind]=np.nan				
				
				# Variables for treatment 
				F_BA_t0=tobs['BA_t0'][iF0[iaF]]
				F_BA_G=(tobs['BA_t0'][iF1[ibF]]-tobs['BA_t0'][iF0[iaF]])/dt
				F_H_obs_t0=tobs['H_obs_t0'][iF0[iaF]]
				F_H_obs_G=(tobs['H_obs_t0'][iF1[ibF]]-tobs['H_obs_t0'][iF0[iaF]])/dt
				F_H_mod_G=(tobs['H_mod_t0'][iF1[ibF]]-tobs['H_mod_t0'][iF0[iaF]])/dt
				F_Bsw_t0=tobs['Bsw_t0'][iF0[iaF]]
				F_Bsw_G=(tobs['Bsw_t0'][iF1[ibF]]-tobs['Bsw_t0'][iF0[iaF]])/dt
				F_Bsw_mod_G=(tobs['Bsw_mod_t0'][iF1[ibF]]-tobs['Bsw_mod_t0'][iF0[iaF]])/dt
				
				# Make sure modelled response is based on the same sample as observed
				ind=np.where( np.isnan(F_H_obs_G)==True )[0]
				F_H_mod_G[ind]=np.nan
				F_Bsw_mod_G[ind]=np.nan
				
				# Remove bad data
				#ind=np.where( (F_BA_t0<0) | (F_BA_G<0) | (F_H_obs_G<0) )[0]
				ind=np.where( (F_BA_t0<0) )[0]
				F_BA_t0[ind]=np.nan
				F_BA_G[ind]=np.nan
				F_H_obs_t0[ind]=np.nan
				F_H_obs_G[ind]=np.nan				
				F_Bsw_t0[ind]=np.nan
				F_Bsw_G[ind]=np.nan
				F_Bsw_mod_G[ind]=np.nan
				
				# Remove modelled values from observed biomass
				F_Bsw_obs_t0=F_Bsw_t0.copy()
				F_Bsw_obs_G=F_Bsw_G.copy()
				ind=np.where(np.isnan(F_H_obs_G)==True)[0]
				F_Bsw_obs_t0[ind]=np.nan
				F_Bsw_obs_G[ind]=np.nan				
	
				d['ID_Site'][iSite]=uSite[iSite]
				d['Year_t0'][iSite]=np.nanmean(tobs['Year_t0'][iC0])
				d['Year_t1'][iSite]=np.nanmean(tobs['Year_t1'][iC1])
				d['TSF_t0'][iSite]=np.nanmean(tobs['TSF_t0'][iC0])
				d['TSF_t1'][iSite]=np.nanmean(tobs['TSF_t1'][iC1])
				d['SA_t0_C'][iSite]=np.nanmean(tobs['SA_t0'][iC0])
				d['SN_t0_C'][iSite]=np.nanmean(tobs['SN_t0'][iC0])
				d['N_H_C'][iSite]=np.sum(np.isnan(C_H_obs_G)==False)
				d['BA_t0_C'][iSite]=np.nanmean(C_BA_t0)
				d['BA_G_C'][iSite]=np.nanmean(C_BA_G)
				d['H_obs_t0_C'][iSite]=np.nanmean(C_H_obs_t0)
				d['H_obs_G_C'][iSite]=100*np.nanmean(C_H_obs_G)
				d['H_mod_G_C'][iSite]=100*np.nanmean(C_H_mod_G)
				d['Bsw_obs_t0_C'][iSite]=np.nanmean(C_Bsw_obs_t0)
				d['Bsw_obs_G_C'][iSite]=np.nanmean(C_Bsw_obs_G)
				d['Bsw_mod_G_C'][iSite]=np.nanmean(C_Bsw_mod_G)
				d['SA_t0_F'][iSite]=np.nanmean(tobs['SA_t0'][iF0])
				d['SN_t0_F'][iSite]=np.nanmean(tobs['SN_t0'][iF0])
				d['N_H_F'][iSite]=np.sum(np.isnan(F_H_obs_G)==False)
				d['BA_t0_F'][iSite]=np.nanmean(F_BA_t0)
				d['BA_G_F'][iSite]=np.nanmean(F_BA_G)
				d['H_obs_t0_F'][iSite]=np.nanmean(F_H_obs_t0)
				d['H_obs_G_F'][iSite]=100*np.nanmean(F_H_obs_G)
				d['H_mod_G_F'][iSite]=100*np.nanmean(F_H_mod_G)
				d['Bsw_obs_t0_F'][iSite]=np.nanmean(F_Bsw_obs_t0)
				d['Bsw_obs_G_F'][iSite]=np.nanmean(F_Bsw_obs_G)
				d['Bsw_mod_G_F'][iSite]=np.nanmean(F_Bsw_mod_G)
			
			d['SA_t0_DA']=np.round(d['SA_t0_F']-d['SA_t0_C'],decimals=2)
			d['SN_t0_DA']=np.round(d['SN_t0_F']-d['SN_t0_C'],decimals=2)
			d['BA_t0_DA']=np.round(d['BA_t0_F']-d['BA_t0_C'],decimals=2)
			d['BA_t0_DR']=np.round((d['BA_t0_F']-d['BA_t0_C'])/d['BA_t0_C']*100,decimals=2)
			d['BA_G_DA']=np.round(d['BA_G_F']-d['BA_G_C'],decimals=2)
			d['BA_G_DR']=np.round((d['BA_G_F']-d['BA_G_C'])/d['BA_G_C']*100,decimals=2)
			d['H_obs_t0_DA']=np.round(d['H_obs_t0_F']-d['H_obs_t0_C'],decimals=2)
			d['H_obs_t0_DR']=np.round((d['H_obs_t0_F']-d['H_obs_t0_C'])/d['H_obs_t0_C']*100,decimals=2)
			d['H_obs_G_DA']=np.round(d['H_obs_G_F']-d['H_obs_G_C'],decimals=2)
			d['H_mod_G_DA']=np.round(d['H_mod_G_F']-d['H_mod_G_C'],decimals=2)
			d['H_obs_G_DR']=np.round((d['H_obs_G_F']-d['H_obs_G_C'])/d['H_obs_G_C']*100,decimals=2)
			d['Bsw_obs_t0_DA']=np.round(d['Bsw_obs_t0_F']-d['Bsw_obs_t0_C'],decimals=2)
			d['Bsw_obs_t0_DR']=np.round((d['Bsw_obs_t0_F']-d['Bsw_obs_t0_C'])/d['Bsw_obs_t0_C']*100,decimals=2)
			d['Bsw_obs_G_DA']=np.round(d['Bsw_obs_G_F']-d['Bsw_obs_G_C'],decimals=2)
			d['Bsw_obs_G_DR']=np.round((d['Bsw_obs_G_F']-d['Bsw_obs_G_C'])/d['Bsw_obs_G_C']*100,decimals=2)
			d['Bsw_mod_G_DA']=np.round(d['Bsw_mod_G_F']-d['Bsw_mod_G_C'],decimals=2)
			d['Bsw_mod_G_DR']=np.round((d['Bsw_mod_G_F']-d['Bsw_mod_G_C'])/d['Bsw_mod_G_C']*100,decimals=2)
			rs[uSpc[iSpc]][uDose[iDose]][uTSF[iTSF]]=d

#%% Create function to extract response
			
def GetStats(Spc,Dose,Var):
	n=np.array([]); mu=np.array([]); se=np.array([])
	for iTSF in range(len(uTSF)):
		try:
			y=rs[Spc][Dose][uTSF[iTSF]][Var]
		except:
			y=np.nan		
		n=np.append(n,np.nansum(rs[Spc][Dose][uTSF[iTSF]]['N_H_C']))
		mu=np.append(mu,np.nanmean(y))	
		se=np.append(se,np.nanstd(y)/np.sqrt(np.sum(~np.isnan(y))))
	return n,mu,se

#%% Plot actual response

plt.close('all'); ms=2.5; Alpha=0.18; cl1=[0.27,0.49,0.77]; cl2=[0.2,0.8,0]; cl3=[0,1,1]
fig,ax=plt.subplots(3,2,figsize=gu.cm2inch(15.5,15))
ax[0,0].plot([0,20],[0,0],'--k')
n,mu,se=GetStats(1,225,'BA_G_DA')
ax[0,0].errorbar(uTSF,mu,yerr=se,capsize=3,linestyle='',ecolor=cl1)
ax[0,0].plot(uTSF,mu,'-o',markersize=ms,color=cl1,label=r'225 kg N ha$^-$$^1$')
n,mu,se=GetStats(1,450,'BA_G_DA')
ax[0,0].errorbar(uTSF,mu,yerr=se,capsize=3,linestyle='',ecolor=cl2)
ax[0,0].plot(uTSF,mu,'-s',markersize=ms,color=cl2,label=r'450 kg N ha$^-$$^1$')
#n,mu,se=GetStats(1,675,'BA_G_DA')
#ax[0,0].errorbar(uTSF,mu,yerr=se,capsize=3,linestyle='',ecolor=cl3)
#ax[0,0].plot(uTSF,mu,'-d',markersize=ms,color=cl3,label=r'675 kgN ha$^-$$^1$')
ax[0,0].set(position=[0.07,0.72,0.42,0.27],xlim=[0,20],xticks=[0,3,6,9,12,15,18],xlabel='Time since N addition, years',ylim=[-1,5],ylabel=r'Basal area increment (cm$^2$ yr$^-$$^1$)')
ax[0,0].legend(loc='upper right')

ax[0,1].plot([0,20],[0,0],'--k')
n,mu,se=GetStats(2,225,'BA_G_DA')
ax[0,1].errorbar(uTSF,mu,yerr=se,capsize=3,linestyle='',ecolor=cl1)
ax[0,1].plot(uTSF,mu,'-o',markersize=ms,color=cl1,label='Gross growth')
n,mu,se=GetStats(2,450,'BA_G_DA')
ax[0,1].errorbar(uTSF,mu,yerr=se,capsize=3,linestyle='',ecolor=cl2)
ax[0,1].plot(uTSF,mu,'-s',markersize=ms,color=cl2,label='Gross growth')
ax[0,1].set(position=[0.57,0.72,0.42,0.27],xlim=[0,20],xticks=[0,3,6,9,12,15,18],xlabel='Time since N addition, years',ylim=[-1,5],ylabel='Basal area increment (cm$^2$ yr$^-$$^1$)')

ax[1,0].plot([0,20],[0,0],'--k')
n,mu,se=GetStats(1,225,'H_obs_G_DA'); 
ax[1,0].errorbar(uTSF,mu,yerr=se,capsize=3,linestyle='',ecolor=cl1)
ax[1,0].plot(uTSF,mu,'-o',markersize=ms,color=cl1,label='Gross growth')
n,mu,se=GetStats(1,450,'H_obs_G_DA')
ax[1,0].errorbar(uTSF,mu,yerr=se,capsize=3,linestyle='',ecolor=cl2)
ax[1,0].plot(uTSF,mu,'-s',markersize=ms,color=cl2,label='Gross growth')
ax[1,0].set(position=[0.07,0.39,0.42,0.27],xlim=[0,20],xticks=[0,3,6,9,12,15,18],xlabel='Time since N addition, years',ylim=[-5,15],yticks=np.arange(-5,20,5),ylabel='Height increment (cm yr$^-$$^1$)')

ax[1,1].plot([0,20],[0,0],'--k')
n,mu,se=GetStats(2,225,'H_obs_G_DA')
ax[1,1].errorbar(uTSF,mu,yerr=se,capsize=3,linestyle='',ecolor=cl1)
ax[1,1].plot(uTSF,mu,'-o',markersize=ms,color=cl1,label='Gross growth')
n,mu,se=GetStats(2,450,'H_obs_G_DA')
ax[1,1].errorbar(uTSF,mu,yerr=se,capsize=3,linestyle='',ecolor=cl2)
ax[1,1].plot(uTSF,mu,'-s',markersize=ms,color=cl2,label='Gross growth')
ax[1,1].set(position=[0.57,0.39,0.42,0.27],xlim=[0,20],xticks=[0,3,6,9,12,15,18],xlabel='Time since N addition, years',ylim=[-5,15],yticks=np.arange(-5,20,5),ylabel='Height increment (cm yr$^-$$^1$)')

ax[2,0].plot([0,20],[0,0],'--k')
n,mu,se=GetStats(1,225,'Bsw_obs_G_DA')
ax[2,0].errorbar(uTSF,mu,yerr=se,capsize=3,linestyle='',ecolor=cl1)
ax[2,0].plot(uTSF,mu,'-o',markersize=ms,color=cl1,label='Gross growth')
n,mu,se=GetStats(1,450,'Bsw_obs_G_DA')
ax[2,0].errorbar(uTSF,mu,yerr=se,capsize=3,linestyle='',ecolor=cl2)
ax[2,0].plot(uTSF,mu,'-s',markersize=ms,color=cl2,label='Gross growth')
ax[2,0].set(position=[0.07,0.06,0.42,0.27],xlim=[0,20],xticks=[0,3,6,9,12,15,18],xlabel='Time since N addition, years',ylim=[-0.5,2],ylabel='Biomass growth (kg C yr$^-$$^1$)')

ax[2,1].plot([0,20],[0,0],'--k')
n,mu,se=GetStats(2,225,'Bsw_obs_G_DA')
ax[2,1].errorbar(uTSF,mu,yerr=se,capsize=3,linestyle='',ecolor=cl1)
ax[2,1].plot(uTSF,mu,'-o',markersize=ms,color=cl1,label='Gross growth')
n,mu,se=GetStats(2,450,'Bsw_obs_G_DA')
ax[2,1].errorbar(uTSF,mu,yerr=se,capsize=3,linestyle='',ecolor=cl2)
ax[2,1].plot(uTSF,mu,'-s',markersize=ms,color=cl2,label='Gross growth')
ax[2,1].set(position=[0.57,0.06,0.42,0.27],xlim=[0,20],xticks=[0,3,6,9,12,15,18],xlabel='Time since N addition, years',ylim=[-0.5,2],ylabel='Biomass growth (kg C yr$^-$$^1$)')

gu.axletters(ax,plt,0.03,0.9,['Coastal Douglas-fir','Western hemlock','','','',''],0.05)
gu.PrintFig(PathFigures + '\\TL_Growth_vs_TSF_DA','emf',600)

#%% Plot relative response 

plt.close('all'); ms=2.5; Alpha=0.18; cl1=[0.27,0.49,0.77]; cl2=[0.2,0.8,0]
fig,ax=plt.subplots(3,2,figsize=gu.cm2inch(15.5,15))
ax[0,0].plot([0,20],[0,0],'--k')
n,mu,se=GetStats(1,225,'BA_G_DR')
ax[0,0].errorbar(uTSF,mu,yerr=se,capsize=3,linestyle='',ecolor=cl1)
ax[0,0].plot(uTSF,mu,'-o',markersize=ms,color=cl1,label='225 kg N ha$^-$$^1$')
n,mu,se=GetStats(1,450,'BA_G_DR')
ax[0,0].errorbar(uTSF,mu,yerr=se,capsize=3,linestyle='',ecolor=cl2)
ax[0,0].plot(uTSF,mu,'-s',markersize=ms,color=cl2,label='450 kg N ha$^-$$^1$')
ax[0,0].set(position=[0.07,0.72,0.42,0.27],xlim=[0,20],xticks=[0,3,6,9,12,15,18],xlabel='Time since N addition, years',ylim=[-10,70],ylabel='Basal area increment (%)')
ax[0,0].legend(loc='upper right')

ax[0,1].plot([0,20],[0,0],'--k')
n,mu,se=GetStats(2,225,'BA_G_DR')
ax[0,1].errorbar(uTSF,mu,yerr=se,capsize=3,linestyle='',ecolor=cl1)
ax[0,1].plot(uTSF,mu,'-o',markersize=ms,color=cl1,label='Gross growth')
n,mu,se=GetStats(2,450,'BA_G_DR')
ax[0,1].errorbar(uTSF,mu,yerr=se,capsize=3,linestyle='',ecolor=cl2)
ax[0,1].plot(uTSF,mu,'-s',markersize=ms,color=cl2,label='Gross growth')
ax[0,1].set(position=[0.57,0.72,0.42,0.27],xlim=[0,20],xticks=[0,3,6,9,12,15,18],xlabel='Time since N addition, years',ylim=[-10,70],ylabel='Basal area increment (%)')

ax[1,0].plot([0,20],[0,0],'--k')
n,mu,se=GetStats(1,225,'H_obs_G_DR'); 
ax[1,0].errorbar(uTSF,mu,yerr=se,capsize=3,linestyle='',ecolor=cl1)
ax[1,0].plot(uTSF,mu,'-o',markersize=ms,color=cl1,label='Gross growth')
n,mu,se=GetStats(1,450,'H_obs_G_DR')
ax[1,0].errorbar(uTSF,mu,yerr=se,capsize=3,linestyle='',ecolor=cl2)
ax[1,0].plot(uTSF,mu,'-s',markersize=ms,color=cl2,label='Gross growth')
ax[1,0].set(position=[0.07,0.39,0.42,0.27],xlim=[0,20],xticks=[0,3,6,9,12,15,18],xlabel='Time since N addition, years',ylim=[-10,70],ylabel='Height increment (%)')

ax[1,1].plot([0,20],[0,0],'--k')
n,mu,se=GetStats(2,225,'H_obs_G_DR')
ax[1,1].errorbar(uTSF,mu,yerr=se,capsize=3,linestyle='',ecolor=cl1)
ax[1,1].plot(uTSF,mu,'-o',markersize=ms,color=cl1,label='Gross growth')
n,mu,se=GetStats(2,450,'H_obs_G_DR')
ax[1,1].errorbar(uTSF,mu,yerr=se,capsize=3,linestyle='',ecolor=cl2)
ax[1,1].plot(uTSF,mu,'-s',markersize=ms,color=cl2,label='Gross growth')
ax[1,1].set(position=[0.57,0.39,0.42,0.27],xlim=[0,20],xticks=[0,3,6,9,12,15,18],xlabel='Time since N addition, years',ylim=[-10,70],ylabel='Height increment (%)')

ax[2,0].plot([0,20],[0,0],'--k')
n,mu,se=GetStats(1,225,'Bsw_obs_G_DR')
ax[2,0].errorbar(uTSF,mu,yerr=se,capsize=3,linestyle='',ecolor=cl1)
ax[2,0].plot(uTSF,mu,'-o',markersize=ms,color=cl1,label='Gross growth')
n,mu,se=GetStats(1,450,'Bsw_obs_G_DR')
ax[2,0].errorbar(uTSF,mu,yerr=se,capsize=3,linestyle='',ecolor=cl2)
ax[2,0].plot(uTSF,mu,'-s',markersize=ms,color=cl2,label='Gross growth')
ax[2,0].set(position=[0.07,0.06,0.42,0.27],xlim=[0,20],xticks=[0,3,6,9,12,15,18],xlabel='Time since N addition, years',ylim=[-10,70],ylabel='Biomass growth (%)')

ax[2,1].plot([0,20],[0,0],'--k')
n,mu,se=GetStats(2,225,'Bsw_obs_G_DR')
ax[2,1].errorbar(uTSF,mu,yerr=se,capsize=3,linestyle='',ecolor=cl1)
ax[2,1].plot(uTSF,mu,'-o',markersize=ms,color=cl1,label='Gross growth')
n,mu,se=GetStats(2,450,'Bsw_obs_G_DR')
ax[2,1].errorbar(uTSF,mu,yerr=se,capsize=3,linestyle='',ecolor=cl2)
ax[2,1].plot(uTSF,mu,'-s',markersize=ms,color=cl2,label='Gross growth')
ax[2,1].set(position=[0.57,0.06,0.42,0.27],xlim=[0,20],xticks=[0,3,6,9,12,15,18],xlabel='Time since N addition, years',ylim=[-10,70],ylabel='Biomass growth (%)')

gu.axletters(ax,plt,0.03,0.9,['Coastal Douglas-fir','Western hemlock','','','',''],0.05)
#plt.savefig(PathProject + '\\Figures\\DA_vs_TSF_' + spc + '.png',format='png',dpi=900)
#gu.PrintFig(PathFigures + '\\TL_Growth_vs_TSF_DR')

#%% Plot modelled vs. obs actual response

Dose=225

plt.close('all'); ms=2.5; Alpha=0.18; cl1=[0.27,0.49,0.77]; cl2=[0.2,0.8,0]
fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(15,10))
ax[0,0].plot([0,20],[0,0],'--k')
n,mu,se=GetStats(1,Dose,'H_obs_G_DA')
ax[0,0].errorbar(uTSF,mu,yerr=se,capsize=3,linestyle='',ecolor=cl1)
ax[0,0].plot(uTSF,mu,'-o',markersize=ms,color=cl1,label='Observed')
n,mu,se=GetStats(1,Dose,'H_mod_G_DA')
ax[0,0].errorbar(uTSF,mu,yerr=se,capsize=3,linestyle='',ecolor=cl2)
ax[0,0].plot(uTSF,mu,'-s',markersize=ms,color=cl2,label='Modelled')
ax[0,0].set(position=[0.07,0.57,0.42,0.4],xlim=[0,20],xticks=[0,3,6,9,12,15,18],xlabel='Time since N addition, years',ylim=[-5,15],yticks=np.arange(-5,20,5),ylabel='Height increment (cm yr$^-$$^1$)')
ax[0,0].legend(loc='upper right')

ax[0,1].plot([0,20],[0,0],'--k')
n,mu,se=GetStats(2,Dose,'H_obs_G_DA')
ax[0,1].errorbar(uTSF,mu,yerr=se,capsize=3,linestyle='',ecolor=cl1)
ax[0,1].plot(uTSF,mu,'-o',markersize=ms,color=cl1,label='Observed')
n,mu,se=GetStats(2,Dose,'H_mod_G_DA')
ax[0,1].errorbar(uTSF,mu,yerr=se,capsize=3,linestyle='',ecolor=cl2)
ax[0,1].plot(uTSF,mu,'-s',markersize=ms,color=cl2,label='Modelled')
ax[0,1].set(position=[0.57,0.57,0.42,0.4],xlim=[0,20],xticks=[0,3,6,9,12,15,18],xlabel='Time since N addition, years',ylim=[-5,15],yticks=np.arange(-5,20,5),ylabel='Height increment (cm yr$^-$$^1$)')

ax[1,0].plot([0,20],[0,0],'--k')
n,mu,se=GetStats(1,Dose,'Bsw_obs_G_DA')
ax[1,0].errorbar(uTSF,mu,yerr=se,capsize=3,linestyle='',ecolor=cl1)
ax[1,0].plot(uTSF,mu,'-o',markersize=ms,color=cl1,label='Observed')
n,mu,se=GetStats(1,Dose,'Bsw_mod_G_DA')
ax[1,0].errorbar(uTSF,mu,yerr=se,capsize=3,linestyle='',ecolor=cl2)
ax[1,0].plot(uTSF,mu,'-s',markersize=ms,color=cl2,label='Modelled')
ax[1,0].set(position=[0.07,0.08,0.42,0.4],xlim=[0,20],xticks=[0,3,6,9,12,15,18],xlabel='Time since N addition, years',ylim=[-1,2],ylabel='Biomass growth (kg C yr$^-$$^1$)')

ax[1,1].plot([0,20],[0,0],'--k')
n,mu,se=GetStats(2,Dose,'Bsw_obs_G_DA')
ax[1,1].errorbar(uTSF,mu,yerr=se,capsize=3,linestyle='',ecolor=cl1)
ax[1,1].plot(uTSF,mu,'-o',markersize=ms,color=cl1,label='Observed')
n,mu,se=GetStats(2,Dose,'Bsw_mod_G_DA')
ax[1,1].errorbar(uTSF,mu,yerr=se,capsize=3,linestyle='',ecolor=cl2)
ax[1,1].plot(uTSF,mu,'-s',markersize=ms,color=cl2,label='Modelled')
ax[1,1].set(position=[0.57,0.08,0.42,0.4],xlim=[0,20],xticks=[0,3,6,9,12,15,18],xlabel='Time since N addition, years',ylim=[-1,2],ylabel='Biomass growth (kg C yr$^-$$^1$)')

gu.axletters(ax,plt,0.03,0.9,['Coastal Douglas-fir','Western hemlock','',''],0.05)
#plt.savefig(PathFigures + '\\TL_Growth_ObsVsMod.png',format='png',dpi=900)
#gu.PrintFig(PathFigures + '\\TL_Growth_ObsVsMod','emf',300)


#%% Average differences

Dose=225
n,mu1,se=GetStats(1,Dose,'H_obs_G_DA')
n,mu2,se=GetStats(1,Dose,'H_mod_G_DA')
da1=np.mean(mu1)-np.mean(mu2)

n,mu1,se=GetStats(2,Dose,'Bsw_obs_G_DA')
n,mu2,se=GetStats(2,Dose,'Bsw_mod_G_DA')
da2=np.mean(mu1)-np.mean(mu2)

# Sample size
N=[]
n,muC,se=GetStats(1,225,'H_obs_t0_C'); N.append(n[0])
n,muC,se=GetStats(1,225,'H_obs_t0_F'); N.append(n[0])
n,muC,se=GetStats(1,450,'H_obs_t0_C'); N.append(n[0])
n,muC,se=GetStats(1,450,'H_obs_t0_F'); N.append(n[0])
n,muC,se=GetStats(2,225,'H_obs_t0_C'); N.append(n[0])
n,muC,se=GetStats(2,225,'H_obs_t0_F'); N.append(n[0])
n,muC,se=GetStats(2,450,'H_obs_t0_C'); N.append(n[0])
n,muC,se=GetStats(2,450,'H_obs_t0_F'); N.append(n[0])

   
# Compare mean size at beginning
spc=1; dose=225;
n,muC,se=GetStats(spc,dose,'SN_t0_C')
n,muF,se=GetStats(spc,dose,'SN_t0_F')
print([muC[0],muF[0]])
n,muC,se=GetStats(spc,dose,'H_obs_t0_C')
n,muF,se=GetStats(spc,dose,'H_obs_t0_F')
print([muC[0],muF[0]])
n,muC,se=GetStats(spc,dose,'Bsw_obs_t0_C')
n,muF,se=GetStats(spc,dose,'Bsw_obs_t0_F')
print([muC[0],muF[0]])

np.mean((muF[0]-muC[0])/muC[0]*100)




#%% Response of survivor growth to N addition (by crown class)

uSpc=[1,2]
uDose=[225,450]
uCC=[1,2,3,4]
uTSF=[3,6,9,12,18]
vars=['ID_Site','Year_t0','Year_t1','TSF_t0','TSF_t1',
	  'SA_t0_C','SN_t0_C','N_H_C','BA_t0_C','BA_G_C','H_obs_t0_C','H_obs_G_C','H_mod_G_C','Bsw_obs_t0_C','Bsw_obs_G_C','Bsw_mod_G_C',
	  'SA_t0_F','SN_t0_F','N_H_F','BA_t0_F','BA_G_F','H_obs_t0_F','H_obs_G_F','H_mod_G_F','Bsw_obs_t0_F','Bsw_obs_G_F','Bsw_mod_G_F']

rs={}
for iSpc in range(len(uSpc)):
	rs[uSpc[iSpc]]={}
	for iDose in range(len(uDose)):
		rs[uSpc[iSpc]][uDose[iDose]]={}		
		for iCC in range(len(uCC)):
			rs[uSpc[iSpc]][uDose[iDose]][uCC[iCC]]={}		
			for iTSF in range(len(uTSF)):
				rs[uSpc[iSpc]][uDose[iDose]][uCC[iCC]][uTSF[iTSF]]={}
			
				uSite=np.unique(tobs['ID_Site'][tobs['SpcLead']==uSpc[iSpc]])
				d={}
				for v in vars:
					d[v]=np.nan*np.ones(uSite.size)		
		
				for iSite in range(uSite.size):	
					# Use TSF=-1 for all but site 29 for western hemlock or TSF=0
					iC0=np.where((tobs['ID_Spc_t0']==uSpc[iSpc]) & (tobs['ID_Site']==uSite[iSite]) & (tobs['CrownClass_t0']==uCC[iCC]) & (tobs['TSF_t0']==-1) & (tobs['N_Dose']==0) & (tobs['Mort']==0) & (tobs['Rec']==0) )[0]
					iF0=np.where((tobs['ID_Spc_t0']==uSpc[iSpc]) & (tobs['ID_Site']==uSite[iSite]) & (tobs['CrownClass_t0']==uCC[iCC]) & (tobs['TSF_t0']==-1) & (tobs['N_Dose']==uDose[iDose]) & (tobs['Mort']==0) & (tobs['Rec']==0) )[0]
					iC0b=np.where((tobs['ID_Spc_t0']==uSpc[iSpc]) & (tobs['ID_Site']==uSite[iSite]) & (tobs['CrownClass_t0']==uCC[iCC]) & (tobs['TSF_t0']==0) & (tobs['N_Dose']==0) & (tobs['Mort']==0) & (tobs['Rec']==0) )[0]
					iF0b=np.where((tobs['ID_Spc_t0']==uSpc[iSpc]) & (tobs['ID_Site']==uSite[iSite]) & (tobs['CrownClass_t0']==uCC[iCC]) & (tobs['TSF_t0']==0) & (tobs['N_Dose']==uDose[iDose]) & (tobs['Mort']==0) & (tobs['Rec']==0) )[0]
					if iC0b.size>iC0.size: 
						iC0=iC0b
						iF0=iF0b
					# Only continue if there are data
					if (iC0.size==0) | (iF0.size==0): 
						continue
					iC1=np.where((tobs['ID_Spc_t0']==uSpc[iSpc]) & (tobs['ID_Site']==uSite[iSite]) & (tobs['CrownClass_t0']==uCC[iCC]) & (tobs['TSF_t0']==uTSF[iTSF]) & (tobs['N_Dose']==0) & (tobs['Mort']==0) & (tobs['Rec']==0) )[0]
					iF1=np.where((tobs['ID_Spc_t0']==uSpc[iSpc]) & (tobs['ID_Site']==uSite[iSite]) & (tobs['CrownClass_t0']==uCC[iCC]) & (tobs['TSF_t0']==uTSF[iTSF]) & (tobs['N_Dose']==uDose[iDose]) & (tobs['Mort']==0) & (tobs['Rec']==0) )[0]
					if (iC1.size==0) | (iF1.size==0): 
						continue
	
					dt=tobs['Year_t0'][iC1[0]]-tobs['Year_t0'][iC0[0]]
					cC,iaC,ibC=np.intersect1d(tobs['ID_Tree'][iC0],tobs['ID_Tree'][iC1],return_indices=True)
					cF,iaF,ibF=np.intersect1d(tobs['ID_Tree'][iF0],tobs['ID_Tree'][iF1],return_indices=True)
	
					C_BA_t0=tobs['BA_t0'][iC0[iaC]]
					C_BA_G=(tobs['BA_t0'][iC1[ibC]]-tobs['BA_t0'][iC0[iaC]])/dt
					C_H_obs_t0=tobs['H_obs_t0'][iC0[iaC]]
					C_H_obs_G=(tobs['H_obs_t0'][iC1[ibC]]-tobs['H_obs_t0'][iC0[iaC]])/dt
					C_H_mod_G=(tobs['H_mod_t0'][iC1[ibC]]-tobs['H_mod_t0'][iC0[iaC]])/dt
					C_Bsw_t0=tobs['Bsw_t0'][iC0[iaC]]
					C_Bsw_G=(tobs['Bsw_t0'][iC1[ibC]]-tobs['Bsw_t0'][iC0[iaC]])/dt
					C_Bsw_mod_G=(tobs['Bsw_mod_t0'][iC1[ibC]]-tobs['Bsw_mod_t0'][iC0[iaC]])/dt
				
					ind=np.where( np.isnan(C_H_obs_G)==True )[0]
					C_H_mod_G[ind]=np.nan
					C_Bsw_mod_G[ind]=np.nan
				
					ind=np.where( (C_BA_t0<0) )[0]
					C_BA_t0[ind]=np.nan
					C_BA_G[ind]=np.nan
					C_H_obs_t0[ind]=np.nan
					C_H_obs_G[ind]=np.nan				
					C_Bsw_t0[ind]=np.nan
					C_Bsw_G[ind]=np.nan
					C_Bsw_mod_G[ind]=np.nan
				
					C_Bsw_obs_t0=C_Bsw_t0
					C_Bsw_obs_G=C_Bsw_G	
					ind=np.where(np.isnan(C_H_obs_G)==True)[0]
					C_Bsw_obs_t0[ind]=np.nan
					C_Bsw_obs_G[ind]=np.nan
	
					F_BA_t0=tobs['BA_t0'][iF0[iaF]]
					F_BA_G=(tobs['BA_t0'][iF1[ibF]]-tobs['BA_t0'][iF0[iaF]])/dt
					F_H_obs_t0=tobs['H_obs_t0'][iF0[iaF]]
					F_H_obs_G=(tobs['H_obs_t0'][iF1[ibF]]-tobs['H_obs_t0'][iF0[iaF]])/dt
					F_H_mod_G=(tobs['H_mod_t0'][iF1[ibF]]-tobs['H_mod_t0'][iF0[iaF]])/dt
					F_Bsw_t0=tobs['Bsw_t0'][iF0[iaF]]
					F_Bsw_G=(tobs['Bsw_t0'][iF1[ibF]]-tobs['Bsw_t0'][iF0[iaF]])/dt
					F_Bsw_mod_G=(tobs['Bsw_mod_t0'][iF1[ibF]]-tobs['Bsw_mod_t0'][iF0[iaF]])/dt
				
					ind=np.where( np.isnan(F_H_obs_G)==True )[0]
					F_H_mod_G[ind]=np.nan
					F_Bsw_mod_G[ind]=np.nan
				
					ind=np.where( (F_BA_t0<0) )[0]
					F_BA_t0[ind]=np.nan
					F_BA_G[ind]=np.nan
					F_H_obs_t0[ind]=np.nan
					F_H_obs_G[ind]=np.nan				
					F_Bsw_t0[ind]=np.nan
					F_Bsw_G[ind]=np.nan
					F_Bsw_mod_G[ind]=np.nan
					
					F_Bsw_obs_t0=F_Bsw_t0
					F_Bsw_obs_G=F_Bsw_G
					ind=np.where(np.isnan(F_H_obs_G)==True)[0]
					F_Bsw_obs_t0[ind]=np.nan
					F_Bsw_obs_G[ind]=np.nan
	
					d['ID_Site'][iSite]=uSite[iSite]
					d['Year_t0'][iSite]=np.nanmean(tobs['Year_t0'][iC0])
					d['Year_t1'][iSite]=np.nanmean(tobs['Year_t1'][iC1])
					d['TSF_t0'][iSite]=np.nanmean(tobs['TSF_t0'][iC0])
					d['TSF_t1'][iSite]=np.nanmean(tobs['TSF_t1'][iC1])
					d['SA_t0_C'][iSite]=np.nanmean(tobs['SA_t0'][iC0])
					d['SN_t0_C'][iSite]=np.nanmean(tobs['SN_t0'][iC0])
					d['N_H_C'][iSite]=np.sum(np.isnan(C_H_obs_G)==False)
					d['BA_t0_C'][iSite]=np.nanmean(C_BA_t0)
					d['BA_G_C'][iSite]=np.nanmean(C_BA_G)
					d['H_obs_t0_C'][iSite]=np.nanmean(C_H_obs_t0)
					d['H_obs_G_C'][iSite]=100*np.nanmean(C_H_obs_G)
					d['H_mod_G_C'][iSite]=100*np.nanmean(C_H_mod_G)
					d['Bsw_obs_t0_C'][iSite]=np.nanmean(C_Bsw_obs_t0)
					d['Bsw_obs_G_C'][iSite]=np.nanmean(C_Bsw_obs_G)
					d['Bsw_mod_G_C'][iSite]=np.nanmean(C_Bsw_obs_G)
					d['SA_t0_F'][iSite]=np.nanmean(tobs['SA_t0'][iF0])
					d['SN_t0_F'][iSite]=np.nanmean(tobs['SN_t0'][iF0])
					d['N_H_F'][iSite]=np.sum(np.isnan(F_H_obs_G)==False)
					d['BA_t0_F'][iSite]=np.nanmean(F_BA_t0)
					d['BA_G_F'][iSite]=np.nanmean(F_BA_G)
					d['H_obs_t0_F'][iSite]=np.nanmean(F_H_obs_t0)
					d['H_obs_G_F'][iSite]=100*np.nanmean(F_H_obs_G)
					d['H_mod_G_F'][iSite]=100*np.nanmean(F_H_mod_G)
					d['Bsw_obs_t0_F'][iSite]=np.nanmean(F_Bsw_obs_t0)
					d['Bsw_obs_G_F'][iSite]=np.nanmean(F_Bsw_obs_G)
					d['Bsw_mod_G_F'][iSite]=np.nanmean(F_Bsw_obs_G)
			
				d['SA_t0_DA']=np.round(d['SA_t0_F']-d['SA_t0_C'],decimals=2)
				d['SN_t0_DA']=np.round(d['SN_t0_F']-d['SN_t0_C'],decimals=2)
				d['BA_t0_DA']=np.round(d['BA_t0_F']-d['BA_t0_C'],decimals=2)
				d['BA_t0_DR']=np.round((d['BA_t0_F']-d['BA_t0_C'])/d['BA_t0_C']*100,decimals=2)
				d['BA_G_DA']=np.round(d['BA_G_F']-d['BA_G_C'],decimals=2)
				d['BA_G_DR']=np.round((d['BA_G_F']-d['BA_G_C'])/d['BA_G_C']*100,decimals=2)
				d['H_obs_t0_DA']=np.round(d['H_obs_t0_F']-d['H_obs_t0_C'],decimals=2)
				d['H_obs_t0_DR']=np.round((d['H_obs_t0_F']-d['H_obs_t0_C'])/d['H_obs_t0_C']*100,decimals=2)
				d['H_obs_G_DA']=np.round(d['H_obs_G_F']-d['H_obs_G_C'],decimals=2)
				d['H_mod_G_DA']=np.round(d['H_mod_G_F']-d['H_mod_G_C'],decimals=2)
				d['H_obs_G_DR']=np.round((d['H_obs_G_F']-d['H_obs_G_C'])/d['H_obs_G_C']*100,decimals=2)
				d['Bsw_obs_t0_DA']=np.round(d['Bsw_obs_t0_F']-d['Bsw_obs_t0_C'],decimals=2)
				d['Bsw_obs_t0_DR']=np.round((d['Bsw_obs_t0_F']-d['Bsw_obs_t0_C'])/d['Bsw_obs_t0_C']*100,decimals=2)
				d['Bsw_obs_G_DA']=np.round(d['Bsw_obs_G_F']-d['Bsw_obs_G_C'],decimals=2)
				d['Bsw_obs_G_DR']=np.round((d['Bsw_obs_G_F']-d['Bsw_obs_G_C'])/d['Bsw_obs_G_C']*100,decimals=2)
				d['Bsw_mod_G_DA']=np.round(d['Bsw_mod_G_F']-d['Bsw_mod_G_C'],decimals=2)
				d['Bsw_mod_G_DR']=np.round((d['Bsw_mod_G_F']-d['Bsw_mod_G_C'])/d['Bsw_mod_G_C']*100,decimals=2)
				rs[uSpc[iSpc]][uDose[iDose]][uCC[iCC]][uTSF[iTSF]]=d

#%% Fnction to get response
				
def GetStats(Spc,Dose,CC,Var):
	n=np.array([]); mu=np.array([]); se=np.array([])
	for iTSF in range(len(uTSF)):
		try:
			y=rs[Spc][Dose][CC][uTSF[iTSF]][Var]
		except:
			y=np.nan		
		n=np.append(n,np.nansum(rs[Spc][Dose][CC][uTSF[iTSF]]['N_H_C']))
		mu=np.append(mu,np.nanmean(y))	
		se=np.append(se,np.nanstd(y)/np.sqrt(np.sum(~np.isnan(y))))
	return n,mu,se

#%% Plot response of survivor growth by crown class

plt.close('all'); ms=2.5; Alpha=0.18; cl1=[0.27,0.49,0.77]; cl2=[0.2,0.8,0]
fig,ax=plt.subplots(4,2,figsize=gu.cm2inch(15.5,15))
ax[0,0].plot([0,20],[0,0],'--k')
n,mu,se=GetStats(1,225,1,'Bsw_obs_G_DA')
ax[0,0].errorbar(uTSF,mu,yerr=se,capsize=3,linestyle='',ecolor=cl1)
ax[0,0].plot(uTSF,mu,'-o',markersize=ms,color=cl1,label='225 kgN / ha')
n,mu,se=GetStats(1,450,1,'Bsw_obs_G_DA')
ax[0,0].errorbar(uTSF,mu,yerr=se,capsize=3,linestyle='',ecolor=cl2)
ax[0,0].plot(uTSF,mu,'-s',markersize=ms,color=cl2,label='450 kgN / ha')
ax[0,0].set(position=[0.07,0.77,0.42,0.21],xlim=[0,20],xticks=[0,3,6,9,12,15,18],xlabel='',ylim=[-1,5],ylabel='Biomass growth (kgC/yr)')
ax[0,0].legend(loc='upper right')

ax[0,1].plot([0,20],[0,0],'--k')
n,mu,se=GetStats(2,225,1,'Bsw_obs_G_DA')
ax[0,1].errorbar(uTSF,mu,yerr=se,capsize=3,linestyle='',ecolor=cl1)
ax[0,1].plot(uTSF,mu,'-o',markersize=ms,color=cl1,label='225 kgN / ha')
n,mu,se=GetStats(2,450,1,'Bsw_obs_G_DA')
ax[0,1].errorbar(uTSF,mu,yerr=se,capsize=3,linestyle='',ecolor=cl2)
ax[0,1].plot(uTSF,mu,'-s',markersize=ms,color=cl2,label='225 kgN / ha')
ax[0,1].set(position=[0.57,0.77,0.42,0.21],xlim=[0,20],xticks=[0,3,6,9,12,15,18],xlabel='',ylim=[-1,5],ylabel='Biomass growth (kgC/yr)')

ax[1,0].plot([0,20],[0,0],'--k')
n,mu,se=GetStats(1,225,2,'Bsw_obs_G_DA')
ax[1,0].errorbar(uTSF,mu,yerr=se,capsize=3,linestyle='',ecolor=cl1)
ax[1,0].plot(uTSF,mu,'-o',markersize=ms,color=cl1,label='225 kgN / ha')
n,mu,se=GetStats(1,450,2,'Bsw_obs_G_DA')
ax[1,0].errorbar(uTSF,mu,yerr=se,capsize=3,linestyle='',ecolor=cl2)
ax[1,0].plot(uTSF,mu,'-s',markersize=ms,color=cl2,label='225 kgN / ha')
ax[1,0].set(position=[0.07,0.53,0.42,0.21],xlim=[0,20],xticks=[0,3,6,9,12,15,18],xlabel='',ylim=[-1,2],ylabel='Biomass growth (kgC/yr)')

ax[1,1].plot([0,20],[0,0],'--k')
n,mu,se=GetStats(2,225,2,'Bsw_obs_G_DA')
ax[1,1].errorbar(uTSF,mu,yerr=se,capsize=3,linestyle='',ecolor=cl1)
ax[1,1].plot(uTSF,mu,'-o',markersize=ms,color=cl1,label='225 kgN / ha')
n,mu,se=GetStats(2,450,2,'Bsw_obs_G_DA')
ax[1,1].errorbar(uTSF,mu,yerr=se,capsize=3,linestyle='',ecolor=cl2)
ax[1,1].plot(uTSF,mu,'-s',markersize=ms,color=cl2,label='225 kgN / ha')
ax[1,1].set(position=[0.57,0.53,0.42,0.21],xlim=[0,20],xticks=[0,3,6,9,12,15,18],xlabel='',ylim=[-1,2],ylabel='Biomass growth (kgC/yr)')

ax[2,0].plot([0,20],[0,0],'--k')
n,mu,se=GetStats(1,225,3,'Bsw_obs_G_DA')
ax[2,0].errorbar(uTSF,mu,yerr=se,capsize=3,linestyle='',ecolor=cl1)
ax[2,0].plot(uTSF,mu,'-o',markersize=ms,color=cl1,label='225 kgN / ha')
n,mu,se=GetStats(1,450,3,'Bsw_obs_G_DA')
ax[2,0].errorbar(uTSF,mu,yerr=se,capsize=3,linestyle='',ecolor=cl2)
ax[2,0].plot(uTSF,mu,'-s',markersize=ms,color=cl2,label='225 kgN / ha')
ax[2,0].set(position=[0.07,0.29,0.42,0.21],xlim=[0,20],xticks=[0,3,6,9,12,15,18],xlabel='',ylim=[-1,2],ylabel='Biomass growth (kgC/yr)')

ax[2,1].plot([0,20],[0,0],'--k')
n,mu,se=GetStats(2,225,3,'Bsw_obs_G_DA')
ax[2,1].errorbar(uTSF,mu,yerr=se,capsize=3,linestyle='',ecolor=cl1)
ax[2,1].plot(uTSF,mu,'-o',markersize=ms,color=cl1,label='225 kgN / ha')
n,mu,se=GetStats(2,450,3,'Bsw_obs_G_DA')
ax[2,1].errorbar(uTSF,mu,yerr=se,capsize=3,linestyle='',ecolor=cl2)
ax[2,1].plot(uTSF,mu,'-s',markersize=ms,color=cl2,label='225 kgN / ha')
ax[2,1].set(position=[0.57,0.29,0.42,0.21],xlim=[0,20],xticks=[0,3,6,9,12,15,18],xlabel='',ylim=[-1,2],ylabel='Biomass growth (kgC/yr)')

ax[3,0].plot([0,20],[0,0],'--k')
n,mu,se=GetStats(1,225,4,'Bsw_obs_G_DA')
ax[3,0].errorbar(uTSF,mu,yerr=se,capsize=3,linestyle='',ecolor=cl1)
ax[3,0].plot(uTSF,mu,'-o',markersize=ms,color=cl1,label='225 kgN / ha')
n,mu,se=GetStats(1,450,4,'Bsw_obs_G_DA')
ax[3,0].errorbar(uTSF,mu,yerr=se,capsize=3,linestyle='',ecolor=cl2)
ax[3,0].plot(uTSF,mu,'-s',markersize=ms,color=cl2,label='225 kgN / ha')
ax[3,0].set(position=[0.07,0.05,0.42,0.21],xlim=[0,20],xticks=[0,3,6,9,12,15,18],xlabel='Time since N addition, years',ylim=[-1,1],ylabel='Biomass growth (kgC/yr)')

ax[3,1].plot([0,20],[0,0],'--k')
n,mu,se=GetStats(2,225,4,'Bsw_obs_G_DA')
ax[3,1].errorbar(uTSF,mu,yerr=se,capsize=3,linestyle='',ecolor=cl1)
ax[3,1].plot(uTSF,mu,'-o',markersize=ms,color=cl1,label='225 kgN / ha')
n,mu,se=GetStats(2,450,4,'Bsw_obs_G_DA')
ax[3,1].errorbar(uTSF,mu,yerr=se,capsize=3,linestyle='',ecolor=cl2)
ax[3,1].plot(uTSF,mu,'-s',markersize=ms,color=cl2,label='225 kgN / ha')
ax[3,1].set(position=[0.57,0.05,0.42,0.21],xlim=[0,20],xticks=[0,3,6,9,12,15,18],xlabel='Time since N addition, years',ylim=[-1,1],ylabel='Biomass growth (kgC/yr)')

gu.axletters(ax,plt,0.03,0.9,['Dominant','Dominant','Co-dominant','Co-dominant','Intermediate','Intermediate','Suppressed','Suppressed'],0.05)

#plt.SaveFig(PathFigures + '\\TL_Growth_vs_TSF_DA_ByCC_.png',format='png',dpi=900)
 


#%% MODELLING

Dose=225
Spc=2
#uSite=np.unique(tobs.loc[tobs['SpcLead']==Spc,'ID_Site'])
uSite=np.unique(tobs['ID_Site'])

d={}
d['BA_t0_C']=np.array([])
d['BA_G_C']=np.array([])
d['BA_t0_F']=np.array([])
d['BA_G_F']=np.array([])
d['H_obs_t0_C']=np.array([])
d['H_obs_G_C']=np.array([])
d['H_obs_t0_F']=np.array([])
d['H_obs_G_F']=np.array([])
d['Bsw_obs_t0_C']=np.array([])
d['Bsw_obs_G_C']=np.array([])
d['Bsw_obs_t0_F']=np.array([])
d['Bsw_obs_G_F']=np.array([])

for iSite in range(uSite.size):
	
	# Use TSF=-1 for all but site 29 for western hemlock or TSF=0
	iC0=np.where((tobs['ID_Spc_t0']==Spc) & (tobs['ID_Site']==uSite[iSite]) & (tobs['TSF_t0']==-1) & (tobs['N_Dose']==0) & (tobs['Mort']==0) & (tobs['Rec']==0) )[0]
	iF0=np.where((tobs['ID_Spc_t0']==Spc) & (tobs['ID_Site']==uSite[iSite]) & (tobs['TSF_t0']==-1) & (tobs['N_Dose']==Dose) & (tobs['Mort']==0) & (tobs['Rec']==0) )[0]
	iC0b=np.where((tobs['ID_Spc_t0']==Spc) & (tobs['ID_Site']==uSite[iSite]) & (tobs['TSF_t0']==0) & (tobs['N_Dose']==0) & (tobs['Mort']==0) & (tobs['Rec']==0) )[0]
	iF0b=np.where((tobs['ID_Spc_t0']==Spc) & (tobs['ID_Site']==uSite[iSite]) & (tobs['TSF_t0']==0) & (tobs['N_Dose']==Dose) & (tobs['Mort']==0) & (tobs['Rec']==0) )[0]
	if iC0b.size>iC0.size: 
		iC0=iC0b
		iF0=iF0b
	# Only continue if there are data
	if (iC0.size==0) | (iF0.size==0): 
		continue
	iC1=np.where((tobs['ID_Spc_t0']==Spc) & (tobs['ID_Site']==uSite[iSite]) & (tobs['TSF_t0']==6) & (tobs['N_Dose']==0) & (tobs['Mort']==0) & (tobs['Rec']==0) )[0]
	iF1=np.where((tobs['ID_Spc_t0']==Spc) & (tobs['ID_Site']==uSite[iSite]) & (tobs['TSF_t0']==6) & (tobs['N_Dose']==Dose) & (tobs['Mort']==0) & (tobs['Rec']==0) )[0]
	if (iC1.size==0) | (iF1.size==0): 
		continue
	
	dt=tobs.loc[iC1[0],'Year_t0']-tobs.loc[iC0[0],'Year_t0']
	
	cC,iaC,ibC=np.intersect1d(tobs.loc[iC0,'ID_Tree'],tobs.loc[iC1,'ID_Tree'],return_indices=True)
	cF,iaF,ibF=np.intersect1d(tobs.loc[iF0,'ID_Tree'],tobs.loc[iF1,'ID_Tree'],return_indices=True)
	
	C_BA_t0=tobs.loc[iC0[iaC],'BA_t0'].values
	C_BA_G=(tobs.loc[iC1[ibC],'BA_t0'].values-tobs.loc[iC0[iaC],'BA_t0'].values)/dt
	C_H_obs_t0=tobs.loc[iC0[iaC],'H_obs_t0'].values
	C_H_obs_G=(tobs.loc[iC1[ibC],'H_obs_t0'].values-tobs.loc[iC0[iaC],'H_obs_t0'].values)/dt
	C_Bsw_t0=tobs.loc[iC0[iaC],'Bsw_t0'].values
	C_Bsw_G=(tobs.loc[iC1[ibC],'Bsw_t0'].values-tobs.loc[iC0[iaC],'Bsw_t0'].values)/dt
	
	ind=np.where( (C_BA_t0<100) | (C_BA_G<0) | (C_H_obs_G<0) )[0]
	C_BA_t0[ind]=np.nan
	C_BA_G[ind]=np.nan
	C_H_obs_t0[ind]=np.nan
	C_H_obs_G[ind]=np.nan
	C_Bsw_t0[ind]=np.nan
	C_Bsw_G[ind]=np.nan
	
	C_Bsw_obs_t0=C_Bsw_t0
	C_Bsw_obs_G=C_Bsw_G	
	ind=np.where(np.isnan(C_H_obs_G)==True)[0]
	C_Bsw_obs_t0[ind]=np.nan
	C_Bsw_obs_G[ind]=np.nan
	
	F_BA_t0=tobs.loc[iF0[iaF],'BA_t0'].values
	F_BA_G=(tobs.loc[iF1[ibF],'BA_t0'].values-tobs.loc[iF0[iaF],'BA_t0'].values)/dt
	F_H_obs_t0=tobs.loc[iF0[iaF],'H_obs_t0'].values
	F_H_obs_G=(tobs.loc[iF1[ibF],'H_obs_t0'].values-tobs.loc[iF0[iaF],'H_obs_t0'].values)/dt
	F_Bsw_t0=tobs.loc[iF0[iaF],'Bsw_t0'].values
	F_Bsw_G=(tobs.loc[iF1[ibF],'Bsw_t0'].values-tobs.loc[iF0[iaF],'Bsw_t0'].values)/dt
	
	ind=np.where( (F_BA_t0<100) | (F_BA_G<0) | (F_H_obs_G<0) )[0]
	F_BA_t0[ind]=np.nan
	F_BA_G[ind]=np.nan
	F_H_obs_t0[ind]=np.nan
	F_H_obs_G[ind]=np.nan
	F_Bsw_t0[ind]=np.nan
	F_Bsw_G[ind]=np.nan
	
	F_Bsw_obs_t0=F_Bsw_t0
	F_Bsw_obs_G=F_Bsw_G
	ind=np.where(np.isnan(F_H_obs_G)==True)[0]
	F_Bsw_obs_t0[ind]=np.nan
	F_Bsw_obs_G[ind]=np.nan
	
	d['BA_t0_C']=np.append(d['BA_t0_C'],C_BA_t0)
	d['BA_G_C']=np.append(d['BA_G_C'],C_BA_G)
	d['BA_t0_F']=np.append(d['BA_t0_F'],F_BA_t0)
	d['BA_G_F']=np.append(d['BA_G_F'],F_BA_G)
	d['H_obs_t0_C']=np.append(d['H_obs_t0_C'],C_H_obs_t0)
	d['H_obs_G_C']=np.append(d['H_obs_G_C'],C_H_obs_G)
	d['H_obs_t0_F']=np.append(d['H_obs_t0_F'],F_H_obs_t0)
	d['H_obs_G_F']=np.append(d['H_obs_G_F'],F_H_obs_G)
	d['Bsw_obs_t0_C']=np.append(d['Bsw_obs_t0_C'],C_Bsw_obs_t0)
	d['Bsw_obs_G_C']=np.append(d['Bsw_obs_G_C'],C_Bsw_obs_G)
	d['Bsw_obs_t0_F']=np.append(d['Bsw_obs_t0_F'],F_Bsw_obs_t0)
	d['Bsw_obs_G_F']=np.append(d['Bsw_obs_G_F'],F_Bsw_obs_G)


plt.close('all')
plt.plot(d['Bsw_obs_t0_C'],d['Bsw_obs_G_C'],'o')
plt.plot(d['Bsw_obs_t0_F'],d['Bsw_obs_G_F'],'s')

plt.close('all')
plt.plot(np.log(d['Bsw_obs_t0_C']),np.log(d['Bsw_obs_G_C']),'b.')
plt.plot(np.log(d['Bsw_obs_t0_F']),np.log(d['Bsw_obs_G_F']),'r.')








def func(x,a,b):
	return a*x**b
from scipy.optimize import curve_fit
xhat=np.arange(0,900,1)

iFit=np.where((np.isnan(d['Bsw_obs_t0_C']+d['Bsw_obs_G_C'])==False))[0]
x=d['Bsw_obs_t0_C'][iFit].astype(float); 
y=d['Bsw_obs_G_C'][iFit].astype(float)		
poptG,pcovG=curve_fit(func,x,y,[5,0.25])
yhatC=func(xhat,poptG[0],poptG[1])

iFit=np.where((np.isnan(d['Bsw_obs_t0_F']+d['Bsw_obs_G_F'])==False))[0]
x=d['Bsw_obs_t0_F'][iFit].astype(float); 
y=d['Bsw_obs_G_F'][iFit].astype(float)		
poptG,pcovG=curve_fit(func,x,y,[5,0.25])
yhatF=func(xhat,poptG[0],poptG[1])

plt.close('all')
fig,ax=plt.subplots(1,2)
ax[0].plot(d['Bsw_obs_t0_C'],d['Bsw_obs_G_C'],'b.')
ax[0].plot(d['Bsw_obs_t0_F'],d['Bsw_obs_G_F'],'r.')
ax[0].plot(xhat,yhatC,'b-')
ax[0].plot(xhat,yhatF,'r--')

ax[1].plot(xhat,yhatF-yhatC,'-',linewidth=2)





def func(x,a,b,c):
	return a+b*x+c*np.log(x)
from scipy.optimize import curve_fit
xhat=np.arange(0,900,1)

iFit=np.where( (np.isnan(d['Bsw_obs_t0_C']+d['Bsw_obs_G_C'])==False) & (d['Bsw_obs_G_C']>0) )[0]
x=d['Bsw_obs_t0_C'][iFit].astype(float); 
y=np.log(d['Bsw_obs_G_C'][iFit].astype(float))
poptG,pcovG=curve_fit(func,x,y,[0.5,0.05,0.2])
yhatC=func(xhat,poptG[0],poptG[1],poptG[2])

iFit=np.where( (np.isnan(d['Bsw_obs_t0_F']+d['Bsw_obs_G_F'])==False) & (d['Bsw_obs_G_F']>0) )[0]
x=d['Bsw_obs_t0_F'][iFit].astype(float); 
y=np.log(d['Bsw_obs_G_F'][iFit].astype(float))
poptG,pcovG=curve_fit(func,x,y,[0.5,0.05,0.2])
yhatF=func(xhat,poptG[0],poptG[1],poptG[2])

plt.close('all')
plt.plot(xhat,yhatC,'-')
plt.plot(xhat,yhatF,'--')

fig,ax=plt.subplots(1)
plt.plot(xhat,yhatF-yhatC,'--',linewidth=2)



def discres(x,y,bw,bin):
	N=np.nan*np.ones(bin.size)
	mu=np.nan*np.ones(bin.size)
	sig=np.nan*np.ones(bin.size)
	se=np.nan*np.ones(bin.size)
	for i in range(bin.size):
		ind=np.where(np.abs(x-bin[i])<=bw/2)[0]
		N[i]=ind.size
		if ind.size==0:
			continue
		mu[i]=np.nanmean(y[ind])
		sig[i]=np.nanstd(y[ind])
		se[i]=np.nanstd(y[ind])/np.sqrt(ind.size)	
	return N,mu,sig,se






#%% EVALUATE HEIGHT FILLING

uSite=np.unique(tobs['ID_Site'])
iSite=0
ind0=np.where( (tobs['ID_Site']==uSite[iSite]) )[0]

uPlot=np.unique(tobs.loc[ind0,'ID_Plot'])
iPlot=1
ind1=np.where( (tobs['ID_Site']==uSite[iSite]) & (tobs['ID_Plot']==uPlot[iPlot]) )[0]

uTree=np.unique(tobs.loc[ind1,'ID_Tree'])

plt.close('all')
fig,ax=plt.subplots(1)
for iTree in range(0,uTree.size+1,10):	
	ind2=np.where( (tobs['ID_Site']==uSite[iSite]) & (tobs['ID_Plot']==uPlot[iPlot]) & (tobs['ID_Tree']==uTree[iTree]) )[0]
	ba=np.append(tobs.loc[ind2,'BA_t0'],tobs.loc[ind2[-1],'BA_t1'])
	h_gf=np.append(tobs.loc[ind2,'H_gf_t0'],tobs.loc[ind2[-1],'H_gf_t1'])
	h_obs=np.append(tobs.loc[ind2,'H_obs_t0'],tobs.loc[ind2[-1],'H_obs_t1'])
	r=np.random.random((1,3)).flatten()
	ind=np.where(h_gf>0)[0]
	ax.plot(ba[ind],h_gf[ind],'.-',color=r)
	ax.plot(ba[ind],h_obs[ind],'o',color=r,markeredgecolor='w',markersize=8)
	ax.set(xlim=[0,2000],ylim=[0,40])
	ax.grid(True)
	
	ind3=np.where( (tobs['ID_Site']==uSite[iSite]) & (tobs['ID_Plot']==uPlot[iPlot]) & (tobs['ID_Tree']==uTree[iTree]) & (tobs['H_obs_t0']>0) )[0]
	r=np.random.random((1,3)).flatten()
	plt.plot(tobs.loc[ind2,'BA_t0'],tobs.loc[ind2,'H_mod_t0'],'.-',color=r)
	plt.plot(tobs.loc[ind3,'BA_t0'],tobs.loc[ind3,'H_obs_t0'],'o',color=r,markeredgecolor='w',markersize=8)
	#plt.plot(tobs.loc[ind2,'SA_t0'],tobs.loc[ind2,'H_gf_t0'],'-')
	

#%% HEIGHT MEASUREMENTS BY CROWN CLASS AND TIME HORIZON

ID_Spc=tobs['ID_Spc_t0'].values
Mort=tobs['Mort'].values
Rec=tobs['Rec'].values
TSF_t0_adj=tobs['TSF_t0_adj'].values
CC=tobs['CrownClass_t0'].values
Y=tobs['H_obs_G'].values

uCC=np.unique(CC)
uCC=np.array(np.arange(1,5))
Nobs=np.nan*np.ones(uCC.size)
Ntot=np.nan*np.ones(uCC.size)
Prcnt=np.nan*np.ones(uCC.size)
for i in range(uCC.size):
	ind1=np.where( (ID_Spc>0) & (CC==uCC[i]) & (Mort==0) & (Rec==0) & (TSF_t0_adj>=1) & (TSF_t0_adj<=9) & (Y>-2) )[0]
	ind2=np.where( (ID_Spc>0) & (CC==uCC[i]) & (Mort==0) & (Rec==0) & (TSF_t0_adj>=1) & (TSF_t0_adj<=9) )[0]	
	Ntot[i]=ind2.size
	Nobs[i]=ind1.size
	Prcnt[i]=ind1.size/ind2.size*100	
	
Prcnt

ind=np.where( (ID_Spc==1) & (CC==2) & (Mort==0) & (Rec==0) )[0]
plt.hist(Y[ind],np.arange(-10,10,0.1))




'''============================================================================

# Tree-level Models

============================================================================'''




bw=10; bin0=np.arange(10,100,bw)
yC=np.nan*np.ones((bin0.size,uSite.size))
yF=np.nan*np.ones((bin0.size,uSite.size))
for iSite in range(uSite.size):
	iC=np.where( (tobs['ID_Site']==uSite[iSite]) & (tobs['N_Dose']==0) & (tobs['ID_Spc_t0']==1) & (tobs['TSF_t0']>=0) & (tobs['TSF_t0']<=9) & (tobs['Mort']==0) & (tobs['Rec']==0) & (tobs['Bsw_G']>-10) & (tobs['H_gf_G']>-50) & (np.isnan(tobs['H_gf_t0']+tobs['Bsw_t0']+tobs['SBswLT_t0']+tobs['TSF_t0']+tobs['N_Dose'])==False) )[0]
	iF=np.where( (tobs['ID_Site']==uSite[iSite]) & (tobs['N_Dose']==225) & (tobs['ID_Spc_t0']==1) & (tobs['TSF_t0']>=0) & (tobs['TSF_t0']<=9) & (tobs['Mort']==0) & (tobs['Rec']==0) & (tobs['Bsw_G']>-10) & (tobs['H_gf_G']>-50) & (np.isnan(tobs['H_gf_t0']+tobs['Bsw_t0']+tobs['SBswLT_t0']+tobs['TSF_t0']+tobs['N_Dose'])==False) )[0]	
	if iF.size==0: continue
	for i in range(bin0.size):
		ind2=np.where( (np.abs(tobs.loc[iC,'SA_t0']-bin0[i])<=bw/2) )[0]
		yC[i,iSite]=np.mean(tobs.loc[iC[ind2],'Bsw_G'])
		ind2=np.where( (np.abs(tobs.loc[iF,'SA_t0']-bin0[i])<=bw/2) )[0]
		yF[i,iSite]=np.mean(tobs.loc[iF[ind2],'Bsw_G'])
	

plt.close('all')
plt.plot(bin0,np.nanmean(yC,axis=1),'-o')
plt.plot(bin0,np.nanmean(yF,axis=1),'-s')

plt.close('all')
plt.plot(bin0,np.nanmean(yF,axis=1)-np.nanmean(yC,axis=1),'-o')










bw=25; bin0=np.arange(25,500,bw)
yC=np.nan*np.ones((bin0.size,uSite.size))
yF=np.nan*np.ones((bin0.size,uSite.size))
for iSite in range(uSite.size):
	iF=np.where( (tobs['ID_Site']==uSite[iSite]) & (tobs['N_Dose']==225) & (tobs['ID_Spc_t0']==1) & (tobs['TSF_t0']>=0) & (tobs['TSF_t0']<=6) & (tobs['Mort']==0) & (tobs['Rec']==0) & (tobs['Bsw_G']>-10) & (tobs['H_gf_G']>-50) & (np.isnan(tobs['H_gf_t0']+tobs['Bsw_t0']+tobs['SBswLT_t0']+tobs['TSF_t0']+tobs['N_Dose'])==False) )[0]
	iC=np.where( (tobs['ID_Site']==uSite[iSite]) & (tobs['N_Dose']==0) & (tobs['ID_Spc_t0']==1) & (tobs['TSF_t0']>=0) & (tobs['TSF_t0']<=6) & (tobs['Mort']==0) & (tobs['Rec']==0) & (tobs['Bsw_G']>-10) & (tobs['H_gf_G']>-50) & (np.isnan(tobs['H_gf_t0']+tobs['Bsw_t0']+tobs['SBswLT_t0']+tobs['TSF_t0']+tobs['N_Dose'])==False) )[0]
	if iF.size==0: continue
	for i in range(bin0.size):
		ind2=np.where( (np.abs(tobs.loc[iF,'SBswLT_t0']-bin0[i])<=bw/2) )[0]
		yF[i,iSite]=np.mean(tobs.loc[iF[ind2],'Bsw_G'])
		ind2=np.where( (np.abs(tobs.loc[iC,'SBswLT_t0']-bin0[i])<=bw/2) )[0]
		yC[i,iSite]=np.mean(tobs.loc[iC[ind2],'Bsw_G'])
	

plt.close('all')
plt.plot(bin0,np.nanmean(yC,axis=1),'-o')
plt.plot(bin0,np.nanmean(yF,axis=1),'-s')

plt.close('all')
plt.plot(bin0,np.nanmean(yF,axis=1)-np.nanmean(yC,axis=1),'-o')


bw=20; bin0=np.arange(10,500,bw)
yC=np.nan*np.ones((bin0.size,uSite.size))
yF=np.nan*np.ones((bin0.size,uSite.size))
for iSite in range(uSite.size):
	iF=np.where( (tobs['ID_Site']==uSite[iSite]) & (tobs['N_Dose']==225) & (tobs['ID_Spc_t0']==2) & (tobs['TSF_t0']>=0) & (tobs['TSF_t0']<=6) & (tobs['Mort']==0) & (tobs['Rec']==0) & (tobs['Bsw_G']>-10) & (tobs['H_gf_G']>-50) & (np.isnan(tobs['H_gf_t0']+tobs['Bsw_t0']+tobs['SBswLT_t0']+tobs['TSF_t0']+tobs['N_Dose'])==False) )[0]
	iC=np.where( (tobs['ID_Site']==uSite[iSite]) & (tobs['N_Dose']==0) & (tobs['ID_Spc_t0']==2) & (tobs['TSF_t0']>=0) & (tobs['TSF_t0']<=6) & (tobs['Mort']==0) & (tobs['Rec']==0) & (tobs['Bsw_G']>-10) & (tobs['H_gf_G']>-50) & (np.isnan(tobs['H_gf_t0']+tobs['Bsw_t0']+tobs['SBswLT_t0']+tobs['TSF_t0']+tobs['N_Dose'])==False) )[0]
	if iF.size==0: continue
	for i in range(bin0.size):
		ind2=np.where( (np.abs(tobs.loc[iF,'SBswLT_t0']-bin0[i])<=bw/2) )[0]
		yF[i,iSite]=np.mean(tobs.loc[iF[ind2],'TopBroke_t1'])
		ind2=np.where( (np.abs(tobs.loc[iC,'SBswLT_t0']-bin0[i])<=bw/2) )[0]
		yC[i,iSite]=np.mean(tobs.loc[iC[ind2],'TopBroke_t1'])
	

plt.close('all')
plt.plot(bin0,np.nanmean(yC,axis=1),'-o')
plt.plot(bin0,np.nanmean(yF,axis=1),'-s')

plt.close('all')
plt.plot(bin0,np.nanmean(yF,axis=1)-np.nanmean(yC,axis=1),'-o')




# Filter
ind=np.where((tobs['ID_Spc_t0']==1) & \
			 (tobs['TSF_t0']>=0) & \
			 (tobs['TSF_t0']<=10) & \
			 (tobs['N_Dose']==0) & \
			 (tobs['Mort']==0) & \
			 (tobs['Rec']==0) & \
			 (tobs['BA_G']>0) & \
			 (tobs['H_G']>0) & \
			 (np.isnan(tobs['H_t0']+tobs['SA_t0']+tobs['SN_t0']+tobs['Bsw_t0']+tobs['SBswLT_t0']+tobs['TSF_t0']+tobs['N_Dose'])==False))[0]
y=tobs.loc[ind,'Bsw_G'].astype(float)
X=np.zeros((ind.shape[0],1))
X[:,0]=tobs.loc[ind,'TSF_t1'].astype(float)
md=sm.OLS(y,X).fit()
md.summary()


ind=np.where((tobs['ID_Spc_t0']==1) & \
			 (tobs['TSF_t0']>=0) & \
			 (tobs['TSF_t0']<=10) & \
			 (tobs['N_Dose']==225) & \
			 (tobs['Mort']==0) & \
			 (tobs['Rec']==0) & \
			 (tobs['BA_G']>0) & \
			 (tobs['H_G']>0) & \
			 (np.isnan(tobs['H_t0']+tobs['SA_t0']+tobs['SN_t0']+tobs['Bsw_t0']+tobs['SBswLT_t0']+tobs['TSF_t0']+tobs['N_Dose'])==False))[0]
y=tobs.loc[ind,'BA_t1'].astype(float)
X=np.zeros((ind.shape[0],1))
X[:,0]=tobs.loc[ind,'TSF_t1'].astype(float)
md=sm.OLS(y,X).fit()
md.summary()







# Filter
ind=np.where((tobs['ID_Spc_t0']==1) & \
			 (tobs['TSF_t0']>=0) & \
			 (tobs['TSF_t0']<=20) & \
			 (tobs['N_Dose']<300) & \
			 (tobs['Mort']==0) & \
			 (tobs['Rec']==0) & \
			 (tobs['Bsw_G']>0) & \
			 (tobs['H_obs_G']>0) & \
			 (tobs['DA_t1'].isna()) & \
			 (tobs['TopBroke_t1'].isna()) & \
			 (tobs['TopDead_t1'].isna()) & \
			 (np.isnan(tobs['Bsw_G']+tobs['SA_t0']+tobs['SN_t0']+tobs['Bsw_t0']+tobs['SBswLT_t0']+tobs['TSF_t0']+tobs['N_Dose'])==False))[0]
#(tobs['H_obs_G']!=999) & \

# Dependent variable
#y=tobs.loc[ind,'H_obs_G'].astype(float)
y=tobs.loc[ind,'Bsw_G'].astype(float)
#y=np.log(tobs.loc[ind,'Bsw_G'].astype(float))

# Predictor variables
X=np.zeros((ind.shape[0],5))
X[:,0]=tobs.loc[ind,'Bsw_t0'].astype(float)
X[:,1]=np.log(tobs.loc[ind,'Bsw_t0'].astype(float))
X[:,2]=tobs.loc[ind,'SBswLT_t0'].astype(float)
X[:,3]=tobs.loc[ind,'SA_t0'].astype(float)
X[:,4]=tobs.loc[ind,'N_Dose'].astype(float)

# Standardize predictor variables
mu=np.zeros(X.shape[1])
sig=np.zeros(X.shape[1])
Xz=X.copy()
for i in range(X.shape[1]):
	mu[i]=np.mean(X[:,i])
	sig[i]=np.std(X[:,i])
	Xz[:,i]=(X[:,i]-mu[i])/sig[i]
md=sm.OLS(y,Xz).fit()
md.summary()

X=np.array([mu[0],np.log(mu[0]),mu[1],mu[2],0]); X=(X-mu)/sig; 
yhat=md.get_prediction(X).summary_frame(alpha=0.05); y0=yhat['mean'].values

X=np.array([mu[0],np.log(mu[0]),mu[1],mu[2],225]); X=(X-mu)/sig;
yhat=md.get_prediction(X).summary_frame(alpha=0.05); y1=yhat['mean'].values

da=y1-y0
dr=(y1-y0)/y0*100
print(y0)
print(y1)
print(da)
print(dr)


uTSF=tobs['TSF_t0'].unique()
n=np.zeros(uTSF.size)
G=np.zeros(uTSF.size)
for i in range(uTSF.size):
	n[i]=np.where(tobs['TSF_t0']==uTSF[i])[0].size
	#G[i]=np.where(tobs['Bsw_G']==uTSF[i])[0].size
plt.plot(uTSF,n,'o')
plt.plot(uTSF,G,'o')




# Filter
ind=np.where((tobs['ID_Spc_t0']==1) & \
			 (tobs['TSF_t0']>=0) & \
			 (tobs['TSF_t0']<=12) & \
			 (tobs['N_Dose']<300) & \
			 (tobs['Mort']==0) & (tobs['Rec']==0) & \
			 (tobs['BA_G']>0) & \
			 (tobs['H_G']>0) & \
			 (tobs['H_obs_G']>0) & \
			 (tobs['Bsw_G']>0) & (tobs['Bsw_G']<100) & \
			 (tobs['DA_t1'].isna()) & \
			 (tobs['TopBroke_t1'].isna()) & \
			 (tobs['TopDead_t1'].isna()) & \
			 (np.isnan(tobs['Bsw_G']+tobs['SA_t0']+tobs['SN_t0']+tobs['Bsw_t0']+tobs['SBswLT_t0']+tobs['TSF_t0']+tobs['N_Dose'])==False))[0]

# Dependent variable
y=tobs.loc[ind,'Bsw_G'].astype(float)
#y=np.log(tobs.loc[ind,'Bsw_G'].astype(float))

# Predictor variables
X=np.zeros((ind.shape[0],6))
X[:,0]=tobs.loc[ind,'Bsw_t0'].astype(float)
X[:,1]=np.log(tobs.loc[ind,'Bsw_t0'].astype(float))
X[:,2]=tobs.loc[ind,'SBswLT_t0'].astype(float)
X[:,3]=tobs.loc[ind,'SA_t0'].astype(float)
X[:,4]=tobs.loc[ind,'N_Dose'].astype(float)
X[:,5]=tobs.loc[ind,'TSF_t0'].astype(float)
#X[:,6]=tobs.loc[ind,'TSF_t0'].astype(float)**2

# Standardize predictor variables
mu=np.zeros(X.shape[1])
sig=np.zeros(X.shape[1])
Xz=X.copy()
for i in range(X.shape[1]):
	mu[i]=np.mean(X[:,i])
	sig[i]=np.std(X[:,i])
	Xz[:,i]=(X[:,i]-mu[i])/sig[i]

# Add interactions
#Xz=np.c_[Xz,Xz[:,4]*Xz[:,5],Xz[:,4]*Xz[:,6]]
Xz=np.c_[Xz,Xz[:,4]*Xz[:,5]]

# Fit model
md=sm.OLS(y,Xz).fit()
md.summary()


# Response to N addition across range of TSF
tsf=np.arange(1,30,0.5)
da=np.zeros(tsf.size)
dr=np.zeros(tsf.size)
for i in range(tsf.size):
	X=np.array([mu[0],np.log(mu[0]),mu[1],mu[2],0,tsf[i]]); X=(X-mu)/sig; X=np.append(X,[X[4]*X[5]])
	yhat=md.get_prediction(X).summary_frame(alpha=0.05); y0=yhat['mean'].values
	X=np.array([mu[0],np.log(mu[0]),mu[1],mu[2],225,tsf[i]]); X=(X-mu)/sig; X=np.append(X,[X[4]*X[5]])
	yhat=md.get_prediction(X).summary_frame(alpha=0.05); y1=yhat['mean'].values
	da[i]=y1-y0
	dr[i]=(y1-y0)/y0*100
plt.close('all'); plt.plot(tsf,da,'o');plt.grid()
#plt.plot(tsf,dr,'o')


X=np.array([mu[0],np.log(mu[0]),mu[1],mu[2],0,5]); X=(X-mu)/sig; X=np.append(X,X[4]*X[5])
yhat=md.get_prediction(X).summary_frame(alpha=0.05); y0=yhat['mean'].values

X=np.array([mu[0],np.log(mu[0]),mu[1],mu[2],225,5]); X=(X-mu)/sig; X=np.append(X,X[4]*X[5])
yhat=md.get_prediction(X).summary_frame(alpha=0.05); y1=yhat['mean'].values

da=y1-y0
dr=(y1-y0)/y0*100
print(y0)
print(y1)
print(da)
print(dr)


















'''============================================================================
TIME SINCE FERTILIZATION RESPONSE (INDIVIDUAL TREE - SURVIVORS ONLY)
The problem with this is that hemlock started at TSF_t0=-1 instead of TSF_t0==0
and did not measure heights at TSF_t0=3. 
============================================================================'''

#uDose=np.array([225,450,675])
uDose=np.array([225])

# Define a continuous annual TSF variable
tsf=np.arange(1,41,1)

uTSF=np.unique(np.column_stack((tobs['TSF_t0_adj'],tobs['TSF_t1'])).astype(float),axis=0)
uTSF=uTSF[1:,:]
#uTSF=uTSF[0:4,:]

uSpc=[1,2]

# Add variables
indO=np.where( (np.isnan(tobs['H_obs_G'])==False) & (tobs['H_obs_G']!=-999) )[0]
indM=np.where( (np.isnan(tobs['H_obs_G'])==True) | (tobs['H_obs_G']==-999) )[0]

tobs['BA_RGR']=np.nan*np.ones(tobs['Bsw_G'].size)
tobs['H_obs_RGR']=np.nan*np.ones(tobs['Bsw_G'].size)
tobs['Bsw_RGR']=np.nan*np.ones(tobs['Bsw_G'].size)
tobs['Bsw_obs_t0']=np.nan*np.ones(tobs['Bsw_G'].size)
tobs['Bsw_obs_t1']=np.nan*np.ones(tobs['Bsw_G'].size)
tobs['Bsw_obs_G']=np.nan*np.ones(tobs['Bsw_G'].size)
tobs['Bsw_obs_RGR']=np.nan*np.ones(tobs['Bsw_G'].size)
tobs['H_mod_RGR']=np.nan*np.ones(tobs['Bsw_G'].size)
tobs['Bsw_mod_t0']=np.nan*np.ones(tobs['Bsw_G'].size)
tobs['Bsw_mod_t1']=np.nan*np.ones(tobs['Bsw_G'].size)
tobs['Bsw_mod_G']=np.nan*np.ones(tobs['Bsw_G'].size)
tobs['Bsw_mod_RGR']=np.nan*np.ones(tobs['Bsw_G'].size)

tobs['BA_RGR']=(np.log(tobs['BA_t1'])-np.log(tobs['BA_t0']))/tobs['DT']
tobs['Bsw_RGR']=(np.log(tobs['Bsw_t1'])-np.log(tobs['Bsw_t0']))/tobs['DT']

tobs.loc[indO,'Bsw_obs_t0']=tobs.loc[indO,'Bsw_t0']
tobs.loc[indO,'Bsw_obs_t1']=tobs.loc[indO,'Bsw_t1']
tobs.loc[indO,'Bsw_obs_G']=tobs.loc[indO,'Bsw_G']
tobs.loc[indO,'H_obs_RGR']=(np.log(tobs.loc[indO,'H_obs_t1'])-np.log(tobs.loc[indO,'H_obs_t0']))/tobs.loc[indO,'DT']
tobs.loc[indO,'Bsw_obs_RGR']=(np.log(tobs.loc[indO,'Bsw_obs_t1'])-np.log(tobs.loc[indO,'Bsw_obs_t0']))/tobs.loc[indO,'DT']
tobs.loc[indM,'Bsw_mod_RGR']=(np.log(tobs.loc[indO,'Bsw_mod_t1'])-np.log(tobs.loc[indO,'Bsw_mod_t0']))/tobs.loc[indO,'DT']
tobs.loc[indM,'Bsw_mod_t0']=tobs.loc[indM,'Bsw_t0']
tobs.loc[indM,'Bsw_mod_t1']=tobs.loc[indM,'Bsw_t1']
tobs.loc[indM,'Bsw_mod_G']=tobs.loc[indM,'Bsw_G']
tobs.loc[indM,'Bsw_mod_RGR']=(np.log(tobs.loc[indM,'Bsw_mod_t1'])-np.log(tobs.loc[indM,'Bsw_mod_t0']))/tobs.loc[indM,'DT']

# Dependent variables
uVar=['BA_t0','BA_t1','BA_G',
	  'H_obs_t0','H_obs_t1','H_obs_G',
	  'H_mod_G','H_gf_G',
	  'Bsw_obs_t0','Bsw_obs_t1','Bsw_G','Bsw_obs_G']

#uVar=['H_obs_t0','H_mod_t0']

ID_Site=tobs['ID_Site'].values
ID_Spc=tobs['ID_Spc_t0'].values
TSF_t0_adj=tobs['TSF_t0_adj'].values.astype(float)
TSF_t1=tobs['TSF_t1'].values.astype(float)
N_Dose=tobs['N_Dose'].values.astype(float)
Mort=tobs['Mort'].values
Rec=tobs['Rec'].values
CC=tobs['CrownClass_t0'].values

# Calculate responses
rsT={}
for iSpc in range(len(uSpc)):
	rsT[uSpc[iSpc]]={}
	for iVar in range(len(uVar)):
		Y=tobs[uVar[iVar]].values
		rsT[uSpc[iSpc]][uVar[iVar]]={}
		for iDose in range(uDose.size):
			print(iSpc,iVar,iDose)
			rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]={}
			# Calculate the mean of control plots and treatment plots at each installation
			N_C=np.nan*np.ones((tsf.size,uSite.size))
			N_F=np.nan*np.ones((tsf.size,uSite.size))
			Y_C_mu=np.nan*np.ones((tsf.size,uSite.size))
			Y_F_mu=np.nan*np.ones((tsf.size,uSite.size))
			for iSite in range(uSite.size):					   
				for iTSF in range(uTSF.shape[0]):					
					iC=np.where( (ID_Site==uSite[iSite]) & (ID_Spc==uSpc[iSpc]) & (Mort==0) & (Rec==0) & (N_Dose==0) & (TSF_t0_adj>=uTSF[iTSF,0]) & (TSF_t1<=uTSF[iTSF,1]) )[0]
					iF=np.where( (ID_Site==uSite[iSite]) & (ID_Spc==uSpc[iSpc]) & (Mort==0) & (Rec==0) & (N_Dose==uDose[iDose]) & (TSF_t0_adj>=uTSF[iTSF,0]) & (TSF_t1<=uTSF[iTSF,1]) )[0]
					if (iC.size==0) | (iF.size==0): continue	
					it=np.where( (tsf>=uTSF[iTSF,0]) & (tsf<=uTSF[iTSF,1]) )[0]
					N_C[it,iSite]=iC.size
					N_F[it,iSite]=iF.size
					Y_C_mu[it,iSite]=np.nanmean(Y[iC])
					Y_F_mu[it,iSite]=np.nanmean(Y[iF])
			
			# Exclude plots if they have no site-paired control/treatment
			ind=np.where((np.isnan(Y_C_mu+Y_F_mu)==True))
			Y_C_mu[ind]=np.nan
			Y_F_mu[ind]=np.nan
			N_C[ind]=0
			N_F[ind]=0
			
			# Calculate difference variables for each installation
			Y_DA=Y_F_mu-Y_C_mu
			Y_DR=(Y_F_mu-Y_C_mu)/Y_C_mu*100
			
			# Summarize
			SampleSize=np.zeros(Y_C_mu.shape[0])
			for i in range(SampleSize.size):
				SampleSize[i]=np.where(np.isnan(Y_C_mu[i,:]+Y_F_mu[i,:])==False)[0].size
			rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]['N_Inst']=SampleSize
			rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]['N_Tree_C']=np.nansum(N_C,axis=1)
			rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]['N_Tree_F']=np.nansum(N_F,axis=1)
			
			rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]['Y_C_mu']=np.nanmean(Y_C_mu,axis=1)
			rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]['Y_F_mu']=np.nanmean(Y_F_mu,axis=1)
			rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]['Y_C_med']=np.nanmedian(Y_C_mu,axis=1)
			rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]['Y_F_med']=np.nanmedian(Y_F_mu,axis=1)
			
			rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]['DA_med']=np.nanmedian(Y_DA,axis=1)
			rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]['DA_mu']=np.nanmean(Y_DA,axis=1)
			rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]['DA_se']=np.nanstd(Y_DA,axis=1)/np.sqrt(SampleSize)
			
			# Relative respons based on installation averages
			#rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]['DR_med']=np.nanmedian(DR,axis=1)
			rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]['DR_FromMeans']=(np.nanmean(Y_F_mu,axis=1)-np.nanmean(Y_C_mu,axis=1))/np.nanmean(Y_C_mu,axis=1)*100
			#rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]['DR_se']=np.nanstd(DR,axis=1)/np.sqrt(SampleSize)
			
			# Relative response based on averaging indiviudal installation relative responses
			rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]['DR_med']=np.nanmedian(Y_DR,axis=1)
			rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]['DR_mu']=np.nanmean(Y_DR,axis=1)
			rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]['DR_se']=np.nanstd(Y_DR,axis=1)/np.sqrt(SampleSize)
			
 
# Table statistics for manuscript

vd='Bsw_obs_G'
spc=2;
dose=225; 
vr='DA_mu'
mu=[]
it=np.where( (tsf>=1) & (tsf<=3) )[0]; mu.append(np.nanmean(rsT[spc][vd][dose][vr][it]))
it=np.where( (tsf>=1) & (tsf<=6) )[0]; mu.append(np.nanmean(rsT[spc][vd][dose][vr][it]))
it=np.where( (tsf>=1) & (tsf<=9) )[0]; mu.append(np.nanmean(rsT[spc][vd][dose][vr][it]))
it=np.where( (tsf>=1) & (tsf<=12) )[0]; mu.append(np.nanmean(rsT[spc][vd][dose][vr][it]))
it=np.where( (tsf>=1) & (tsf<=20) )[0]; mu.append(np.nanmean(rsT[spc][vd][dose][vr][it]))
mu




'''============================================================================
RESPONSE BY CROWN CLASS (INDIVIDUAL TREE - SURVIVORS ONLY)
The problem with this is that hemlock started at TSF_t0=-1 instead of TSF_t0==0
and did not measure heights at TSF_t0=3. 
============================================================================'''

uDose=np.array([225])

# Define a continuous annual TSF variable
tsf=np.arange(1,41,1)

uTSF=np.unique(np.column_stack((tobs['TSF_t0_adj'],tobs['TSF_t1'])).astype(float),axis=0)
uTSF=uTSF[1:10,:]

uSpc=[1,2]

uCC=[1,2,3,4]

# Dependent variables
uVar=['BA_G','H_obs_t0','H_obs_G','H_mod_G','H_gf_G']

# Calculate responses
rsT={}
for iSpc in range(len(uSpc)):
	rsT[uSpc[iSpc]]={}
	for iVar in range(len(uVar)):
		Y=tobs[uVar[iVar]].values
		rsT[uSpc[iSpc]][uVar[iVar]]={}
		for iDose in range(uDose.size):
			rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]={}
			for iCC in range(len(uCC)):
				print(iSpc,iVar,iDose)
				rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]][uCC[iCC]]={}
				# Calculate the mean of control plots and treatment plots at each installation
				N_C=np.nan*np.ones((tsf.size,uSite.size))
				N_F=np.nan*np.ones((tsf.size,uSite.size))
				Y_C_mu=np.nan*np.ones((tsf.size,uSite.size))
				Y_F_mu=np.nan*np.ones((tsf.size,uSite.size))
				for iSite in range(uSite.size):					   
					for iTSF in range(uTSF.shape[0]):					
						iC=np.where( (CC==uCC[iCC]) & (ID_Site==uSite[iSite]) & (ID_Spc==uSpc[iSpc]) & (Mort==0) & (Rec==0) & (N_Dose==0) & (TSF_t0_adj>=uTSF[iTSF,0]) & (TSF_t1<=uTSF[iTSF,1]) )[0]
						iF=np.where( (CC==uCC[iCC]) & (ID_Site==uSite[iSite]) & (ID_Spc==uSpc[iSpc]) & (Mort==0) & (Rec==0) & (N_Dose==uDose[iDose]) & (TSF_t0_adj>=uTSF[iTSF,0]) & (TSF_t1<=uTSF[iTSF,1]) )[0]
						if (iC.size==0) | (iF.size==0): continue	
						it=np.where( (tsf>=uTSF[iTSF,0]) & (tsf<=uTSF[iTSF,1]) )[0]
						N_C[it,iSite]=iC.size
						N_F[it,iSite]=iF.size
						Y_C_mu[it,iSite]=np.nanmean(Y[iC])
						Y_F_mu[it,iSite]=np.nanmean(Y[iF])
			
				# Exclude plots if they have no site-paired control/treatment
				ind=np.where((np.isnan(Y_C_mu+Y_F_mu)==True))
				Y_C_mu[ind]=np.nan
				Y_F_mu[ind]=np.nan
				N_C[ind]=0
				N_F[ind]=0
				
				# Calculate difference variables for each installation
				Y_DA=Y_F_mu-Y_C_mu
				Y_DR=(Y_F_mu-Y_C_mu)/Y_C_mu*100
				
				# Summarize
				SampleSize=np.zeros(Y_C_mu.shape[0])
				for i in range(SampleSize.size):
					SampleSize[i]=np.where(np.isnan(Y_C_mu[i,:]+Y_F_mu[i,:])==False)[0].size
				rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]][uCC[iCC]]['N_Inst']=SampleSize
				rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]][uCC[iCC]]['N_Tree_C']=np.nansum(N_C,axis=1)
				rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]][uCC[iCC]]['N_Tree_F']=np.nansum(N_F,axis=1)
			
				rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]][uCC[iCC]]['Y_C_mu']=np.nanmean(Y_C_mu,axis=1)
				rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]][uCC[iCC]]['Y_F_mu']=np.nanmean(Y_F_mu,axis=1)
				rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]][uCC[iCC]]['Y_C_med']=np.nanmedian(Y_C_mu,axis=1)
				rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]][uCC[iCC]]['Y_F_med']=np.nanmedian(Y_F_mu,axis=1)
			
				rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]][uCC[iCC]]['DA_med']=np.nanmedian(Y_DA,axis=1)
				rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]][uCC[iCC]]['DA_mu']=np.nanmean(Y_DA,axis=1)
				rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]][uCC[iCC]]['DA_se']=np.nanstd(Y_DA,axis=1)/np.sqrt(SampleSize)
			
				# Relative respons based on installation averages
				#rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]['DR_med']=np.nanmedian(DR,axis=1)
				rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]][uCC[iCC]]['DR_FromMeans']=(np.nanmean(Y_F_mu,axis=1)-np.nanmean(Y_C_mu,axis=1))/np.nanmean(Y_C_mu,axis=1)*100
				#rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]['DR_se']=np.nanstd(DR,axis=1)/np.sqrt(SampleSize)
			
				# Relative response based on averaging indiviudal installation relative responses
				rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]][uCC[iCC]]['DR_med']=np.nanmedian(Y_DR,axis=1)
				rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]][uCC[iCC]]['DR_mu']=np.nanmean(Y_DR,axis=1)
				rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]][uCC[iCC]]['DR_se']=np.nanstd(Y_DR,axis=1)/np.sqrt(SampleSize)
				

# Table statistics for manuscript

vd='H_obs_G'
spc=2;
dose=225; 
cc=1
vr='DR_FromMeans'
mu=[]
it=np.where( (tsf>=1) & (tsf<=3) )[0]; mu.append(np.nanmean(rsT[spc][vd][dose][cc][vr][it]))
it=np.where( (tsf>=1) & (tsf<=6) )[0]; mu.append(np.nanmean(rsT[spc][vd][dose][cc][vr][it]))
it=np.where( (tsf>=1) & (tsf<=9) )[0]; mu.append(np.nanmean(rsT[spc][vd][dose][cc][vr][it]))
it=np.where( (tsf>=1) & (tsf<=12) )[0]; mu.append(np.nanmean(rsT[spc][vd][dose][cc][vr][it]))
it=np.where( (tsf>=1) & (tsf<=20) )[0]; mu.append(np.nanmean(rsT[spc][vd][dose][cc][vr][it]))
mu









'''============================================================================
QUALITY ASSURANCE
============================================================================'''

#------------------------------------------------------------------------------
# Save site summaries to spreadsheet for inspection
#------------------------------------------------------------------------------

cols=['ID_Site','SpcLead','Year_t0','TSF_t0','TSF_t1','C1','C2','C3','C4','C5','C6','Cmu','Cse','F1','F2','F3','F4','F5','F6','Fmu','Fse','Dmu']

uSite=np.unique(tobs['ID_Site'])

sy='H_obs_t0'
dose=225

df0=pd.DataFrame(data=np.nan*np.ones((360,len(cols))),columns=cols)
cnt=0
for i in range(uSite.size):
	ind=np.where((tobs['ID_Site']==uSite[i]))[0]
	if ind.size==0: continue
	uTSF=np.unique(np.column_stack((tobs.loc[ind,'TSF_t0'],tobs.loc[ind,'TSF_t1'])).astype(float),axis=0)
	for j in range(uTSF.shape[0]):
		df0.loc[cnt,'ID_Site']=uSite[i]		
		df0.loc[cnt,'TSF_t0']=uTSF[j,0]
		df0.loc[cnt,'TSF_t1']=uTSF[j,1]
		iC=np.where((tobs['ID_Spc_t0']==2) & (tobs['ID_Site']==uSite[i]) & (tobs['TSF_t0']==uTSF[j,0]) & (tobs['N_Dose']==0) & (tobs['Mort']==0) & (tobs['Rec']==0) )[0]
		iF=np.where((tobs['ID_Spc_t0']==2) & (tobs['ID_Site']==uSite[i]) & (tobs['TSF_t0']==uTSF[j,0]) & (tobs['N_Dose']==dose) & (tobs['Mort']==0) & (tobs['Rec']==0) )[0]
		if (iC.size==0) | (iF.size==0): continue
		df0.loc[cnt,'Year_t0']=tobs.loc[iC[0],'Year_t0']		
		df0.loc[cnt,'SpcLead']=tobs.loc[iC[0],'SpcLead']
		for k in range(iC.size):
			if k==0: df0.loc[cnt,'C1']=tobs.loc[iC[k],sy]
			elif k==1: df0.loc[cnt,'C2']=tobs.loc[iC[k],sy]
			elif k==2: df0.loc[cnt,'C3']=tobs.loc[iC[k],sy]
			elif k==3: df0.loc[cnt,'C4']=tobs.loc[iC[k],sy]
			elif k==4: df0.loc[cnt,'C5']=tobs.loc[iC[k],sy]
			elif k==5: df0.loc[cnt,'C6']=tobs.loc[iC[k],sy]				
		df0.loc[cnt,'Cmu']=np.nanmean(tobs.loc[iC,sy])
		df0.loc[cnt,'Cse']=np.nanstd(tobs.loc[iC,sy])/np.sqrt(iC.size)
		for k in range(iF.size):
			if k==0: df0.loc[cnt,'F1']=tobs.loc[iF[k],sy]
			elif k==1: df0.loc[cnt,'F2']=tobs.loc[iF[k],sy]
			elif k==2: df0.loc[cnt,'F3']=tobs.loc[iF[k],sy]
			elif k==3: df0.loc[cnt,'F4']=tobs.loc[iF[k],sy]
			elif k==4: df0.loc[cnt,'F5']=tobs.loc[iF[k],sy]
			elif k==5: df0.loc[cnt,'F6']=tobs.loc[iF[k],sy]
		df0.loc[cnt,'Fmu']=np.nanmean(tobs.loc[iF,sy])
		df0.loc[cnt,'Fse']=np.nanstd(tobs.loc[iF,sy])/np.sqrt(iF.size)
		df0.loc[cnt,'Dmu']=df0.loc[cnt,'Fmu']-df0.loc[cnt,'Cmu']		
		cnt=cnt+1
	
df0.to_excel(PathProject + '\\Processed\\TL_Summary_' + sy + '.xlsx',index=False)


