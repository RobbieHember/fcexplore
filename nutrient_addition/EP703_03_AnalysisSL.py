
'''============================================================================

ANALYSIS OF EP703 IN SUPPORT OF ESTIMATING THE GHG BENEFIT FROM COASTAL 
FERTILIZATION

Steps:
    1) Import data
    2) Quality assurance
    3) Descriptive statistics
    4) N response vs. time since N addition

Inputs:
    1) EP703 stand-level dataframe from EP703_01_ProcessRawData.py

Outputs:
    1) Graphics for article
    2) Spreadsheets for inspection and article tables
    
============================================================================'''

#%% Import Modules

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
from fcgadgets.utilities import utilities_general as gu

#%% Set figure properties

fs=7
params={'font.sans-serif':'Calibri',
        'font.size':fs,
        'axes.edgecolor':'black',
        'axes.labelsize':fs,
        'axes.labelcolor':'black',
        'axes.titlesize':fs,
        'axes.linewidth':0.5,        
        'lines.linewidth':0.5,
        'text.color':'black',
        'xtick.color':'black',        
        'xtick.labelsize':fs,
        'xtick.major.width':0.5,
        'xtick.major.size':3,
        'xtick.direction':'in',
        'ytick.color':'black',
        'ytick.labelsize':fs,
        'ytick.major.width':0.5,
        'ytick.major.size':3,
        'ytick.direction':'in',
        'legend.fontsize':fs,        
        'savefig.dpi':300,
        'savefig.transparent':True}
plt.rcParams.update(params)

#%% Import paths

PathProject=r'C:\Users\rhember\Documents\Data\EP703'
PathFigures=r'C:\Users\rhember\OneDrive - Government of BC\Figures\EP703'
PathManuscriptSI=r'G:\My Drive\Manuscripts\CCIPB_FertilizationGHGBenefit_EP703\SI'

#%% Import data

# Import site information
dS=gu.ReadExcel(PathProject + '\\EP703_SiteInfo.xlsx')

# Unique sites
uSite=np.unique(dS['ID_Instal'])

# Import stand-level data
sobs=gu.ipickle(PathProject + '\\Processed\\EP703_SL.pkl')

#%% Derived variables and adjustments

# Create an adjusted value of TSF_t0 that will ensure time intervals are not double
# counting the year of measurement
sobs['TSF_t0_adj']=sobs['TSF_t0']+1

sobs['Bsw_Net_RGR']=(np.log(sobs['Bsw_t1'])-np.log(sobs['Bsw_t0']))/sobs['DT']

# Initial biomass
u=np.unique(np.column_stack((sobs['ID_Site'],sobs['ID_Plot'])).astype(float),axis=0)
sobs['Bsw_Init']=np.nan*np.ones(sobs['Bsw_t0'].size)
for i in range(u.shape[0]):
    ind=np.where( (sobs['ID_Site']==u[i,0]) & (sobs['ID_Plot']==u[i,1]) )[0]
    atsf=np.abs(sobs['TSF_t0'][ind])
    if ind.size!=0:
        ind2=np.where( atsf==np.min(atsf) )[0]
        sobs['Bsw_Init'][ind]=sobs['Bsw_t0'][ind[ind2[0]]]

# Relative change in stand density (%/yr)
sobs['dN_rel']=sobs['N_Net']/sobs['N_t0']*100

sobs['RGR']=(np.log(sobs['Bsw_t1'])-np.log(sobs['Bsw_t0']))/sobs['DT']

sobs['Bsw_G_ave']=sobs['Bsw_G']/sobs['N_t1']*1000


#%% Descriptive statistics

# Number of installations and C-F plot pairs (FD leading)
ind=np.where( (sobs['Spc1_ID_t0']=='FD') & (sobs['N_Dose']==0) & (sobs['TSF_t0']==0) )[0]
np.unique(sobs['ID_Site'][ind]).size
np.unique(np.column_stack((sobs['ID_Site'][ind],sobs['ID_Plot'][ind])),axis=0).size

# Number of installations and C-F plot pairs (HW leading)
ind=np.where( (sobs['Spc1_ID_t0']=='HW') & (sobs['N_Dose']==0) & (sobs['TSF_t0']==0) )[0]
np.unique(sobs['ID_Site'][ind]).size
np.unique(np.column_stack((sobs['ID_Site'][ind],sobs['ID_Plot'][ind])),axis=0).size

# Biomass by pool
ind=np.where( (sobs['Spc1_ID_t0']=='FD') & (sobs['N_Dose']==0) )[0]
np.nanmean(sobs['Bsw_t0'][ind])
np.nanmean(sobs['Bbk_t0'][ind])
np.nanmean(sobs['Bbr_t0'][ind])
np.nanmean(sobs['Bf_t0'][ind])

# Get initial stand attributes at plot/treatment level
uDose=np.unique(sobs['N_Dose'])

d={}
vs=['Spc1_ID_t0','Spc2_ID_t0','Spc1_Pct_t0','Spc2_Pct_t0','N_Dose','SI','Age_t0',
    'N_t0','Dam_t0','H_obs_t0','H_gf_t0','Bsw_t0']
for i in range(len(vs)):
    d[vs[i]]=np.array([])

for iSite in range(uSite.size):
    ind=np.where(sobs['ID_Site']==uSite[iSite])[0]
    for iDose in range(uDose.size):
        ind=np.where( (sobs['ID_Site']==uSite[iSite]) & (sobs['N_Dose']==uDose[iDose]) )[0]
        ind=np.where( (sobs['ID_Site']==uSite[iSite]) & (sobs['N_Dose']==uDose[iDose]) & (sobs['TSF_t0_adj']==0) )[0]
        if ind.size==0: continue
        for k in d.keys():
            try:
                d[k]=np.append(d[k],np.nanmean(sobs[k][ind]))
            except:
                d[k]=np.append(d[k],sobs[k][ind[0]])

df=pd.DataFrame(columns=['Variable','FD_C_mean','FD_C_sig','FD_225_mean','FD_225_sig',
                         'FD_450_mean','FD_450_sig','FD_675_mean','FD_675_sig',
                         'HW_C_mean','HW_C_sig','HW_225_mean','HW_225_sig',
                         'HW_450_mean','HW_450_sig','HW_675_mean','HW_675_sig'])

iFDC=np.where( (d['Spc1_ID_t0']=='FD') & (d['N_Dose']==0) )[0]
iFD225=np.where( (d['Spc1_ID_t0']=='FD') & (d['N_Dose']==225) )[0]
iFD450=np.where( (d['Spc1_ID_t0']=='FD') & (d['N_Dose']==450) )[0]
iFD675=np.where( (d['Spc1_ID_t0']=='FD') & (d['N_Dose']==675) )[0]
iHWC=np.where( (d['Spc1_ID_t0']=='HW') & (d['N_Dose']==0) )[0]
iHW225=np.where( (d['Spc1_ID_t0']=='HW') & (d['N_Dose']==225) )[0]
iHW450=np.where( (d['Spc1_ID_t0']=='HW') & (d['N_Dose']==450) )[0]
iHW675=np.where( (d['Spc1_ID_t0']=='HW') & (d['N_Dose']==675) )[0]
for k in d.keys():
    if (k=='Spc1_ID_t0') | (k=='Spc2_ID_t0'): continue
    df=df.append({'Variable':k,
               'FD_C_mean':np.nanmean(d[k][iFDC]),'FD_C_sig':np.nanstd(d[k][iFDC]),
               'FD_225_mean':np.nanmean(d[k][iFD225]),'FD_225_sig':np.nanstd(d[k][iFD225]),
               'FD_450_mean':np.nanmean(d[k][iFD450]),'FD_450_sig':np.nanstd(d[k][iFD450]),
               'FD_675_mean':np.nanmean(d[k][iFD675]),'FD_675_sig':np.nanstd(d[k][iFD675]),
               'HW_C_mean':np.nanmean(d[k][iHWC]),'HW_C_sig':np.nanstd(d[k][iHWC]),
               'HW_225_mean':np.nanmean(d[k][iHW225]),'HW_225_sig':np.nanstd(d[k][iHW225]),
               'HW_450_mean':np.nanmean(d[k][iHW450]),'HW_450_sig':np.nanstd(d[k][iHW450]),
               'HW_675_mean':np.nanmean(d[k][iHW675]),'HW_675_sig':np.nanstd(d[k][iHW675])},ignore_index=True)

#df.to_excel(PathManuscriptSI + '\\InitialStandAttributes.xlsx')


#%% Difference between growth rate of control and fert plots prior to fert

uTSF=np.unique(np.column_stack((sobs['TSF_t0'],sobs['TSF_t1'])).astype(float),axis=0)

dose=450
spc='FD'

iC0=np.where( (sobs['Spc1_ID_t0']==spc) & (sobs['TSF_t0']==-1) & (sobs['N_Dose']==0) )[0]
iF0=np.where( (sobs['Spc1_ID_t0']==spc) & (sobs['TSF_t0']==-1) & (sobs['N_Dose']==dose) )[0]
id_t0_C=np.unique(np.column_stack((sobs['ID_Site'][iC0],sobs['ID_Plot'][iC0])).astype(float),axis=0)
id_t0_F=np.unique(np.column_stack((sobs['ID_Site'][iF0],sobs['ID_Plot'][iF0])).astype(float),axis=0)

Bsw_t0_C=np.nan*np.ones(id_t0_C.shape[0])
for i in range(id_t0_C.shape[0]):
    ind=np.where( (sobs['ID_Site']==id_t0_C[i,0]) & (sobs['ID_Plot']==id_t0_C[i,1]) & (sobs['TSF_t0']==-1) & (sobs['N_Dose']==0) )[0]
    Bsw_t0_C[i]=np.nanmean(sobs['Bsw_t0'][ind])
Bsw_t0_F=np.nan*np.ones(id_t0_F.shape[0])
for i in range(id_t0_F.shape[0]):
    ind=np.where( (sobs['ID_Site']==id_t0_F[i,0]) & (sobs['ID_Plot']==id_t0_F[i,1]) & (sobs['TSF_t0']==-1) & (sobs['N_Dose']>0) )[0]
    Bsw_t0_F[i]=np.nanmean(sobs['Bsw_t0'][ind])

Bsw_t1_C=np.nan*np.ones(id_t0_C.shape[0])
for i in range(id_t0_C.shape[0]):
    ind=np.where( (sobs['ID_Site']==id_t0_C[i,0]) & (sobs['ID_Plot']==id_t0_C[i,1]) & (sobs['TSF_t0']==0) & (sobs['N_Dose']==0) )[0]
    Bsw_t1_C[i]=np.nanmean(sobs['Bsw_t0'][ind])

Bsw_t1_F=np.nan*np.ones(id_t0_F.shape[0])
for i in range(id_t0_F.shape[0]):
    ind=np.where( (sobs['ID_Site']==id_t0_F[i,0]) & (sobs['ID_Plot']==id_t0_F[i,1]) & (sobs['TSF_t0']==0) & (sobs['N_Dose']>0) )[0]
    Bsw_t1_F[i]=np.nanmean(sobs['Bsw_t0'][ind])

Gnet_C=Bsw_t1_C-Bsw_t0_C
Gnet_F=Bsw_t1_F-Bsw_t0_F

ikp=np.where(Gnet_C>-5)[0]
stats_C=np.array([np.nanmean(Gnet_C[ikp]),np.nanstd(Gnet_C[ikp])/np.sqrt(ikp.size)])
ikp=np.where(Gnet_F>-5)[0]
stats_F=np.array([np.nanmean(Gnet_F[ikp]),np.nanstd(Gnet_F[ikp])/np.sqrt(ikp.size)])
print(stats_C)
print(stats_F)

stats_C[0]/stats_F[0]


#%% Nitrogen response of biomass fluxes vs. time since fertilization

#uDose=np.unique(sobs['N_Dose'])
#uDose=uDose[1:]
uDose=np.array([225,450,675])

# Define a continuous annual TSF variable
tsf=np.arange(1,41,1)

uTSF=np.unique(np.column_stack((sobs['TSF_t0_adj'],sobs['TSF_t1'])).astype(float),axis=0)
uTSF=uTSF[1:,:]

uSpc=['FD','HW']

uVar=['Bsw_G','Bsw_M','Bsw_Net']

rsT={}
for iSpc in range(len(uSpc)):
    rsT[uSpc[iSpc]]={}
    for iVar in range(len(uVar)):
        Y=sobs[uVar[iVar]]
        rsT[uSpc[iSpc]][uVar[iVar]]={}
        for iDose in range(uDose.size):
            print(iSpc,iVar,iDose)
            rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]={}            
            Y_C_mu=np.nan*np.ones((tsf.size,uSite.size))
            Y_F_mu=np.nan*np.ones((tsf.size,uSite.size))
            for iSite in range(uSite.size):
                for iTSF in range(uTSF.shape[0]):
                    iC=np.where( (sobs['ID_Site']==uSite[iSite]) & (sobs['Spc1_ID_t0']==uSpc[iSpc]) & (sobs['TSF_t0_adj']>=uTSF[iTSF,0]) & (sobs['TSF_t1']<=uTSF[iTSF,1]) & (sobs['N_Dose']==0) )[0]
                    iF=np.where( (sobs['ID_Site']==uSite[iSite]) & (sobs['Spc1_ID_t0']==uSpc[iSpc]) & (sobs['TSF_t0_adj']>=uTSF[iTSF,0]) & (sobs['TSF_t1']<=uTSF[iTSF,1]) & (sobs['N_Dose']==uDose[iDose]) )[0]
                    if (iC.size==0) | (iF.size==0): continue    
                    it=np.where( (tsf>=uTSF[iTSF,0]) & (tsf<=uTSF[iTSF,1]) )[0]
                    Y_C_mu[it,iSite]=np.nanmean(sobs[uVar[iVar]][iC])
                    Y_F_mu[it,iSite]=np.nanmean(sobs[uVar[iVar]][iF])
            
            # Exclude plots if they have no site-paired control/treatment
            ind=np.where((np.isnan(Y_C_mu+Y_F_mu)==True))
            Y_C_mu[ind]=np.nan
            Y_F_mu[ind]=np.nan
            
            # Calculate differences
            DA=Y_F_mu-Y_C_mu
            DR=DA/Y_C_mu*100
            
            # Summarize
            SampleSize=np.zeros(Y_C_mu.shape[0])
            for i in range(SampleSize.size):
                SampleSize[i]=np.where(np.isnan(Y_C_mu[i,:]+Y_F_mu[i,:])==False)[0].size
            rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]['N_Inst']=SampleSize
            rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]['yC_mu']=np.nanmean(Y_C_mu,axis=1)
            rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]['yF_mu']=np.nanmean(Y_F_mu,axis=1)
            rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]['DA_med']=np.nanmedian(DA,axis=1)
            rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]['DA_mu']=np.nanmean(DA,axis=1)
            rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]['DA_se']=np.nanstd(DA,axis=1)/np.sqrt(SampleSize)
            rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]['DR_med']=np.nanmedian(DR,axis=1)
            rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]['DR_mu']=np.nanmean(DR,axis=1)
            rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]['DR_se']=np.nanstd(DR,axis=1)/np.sqrt(SampleSize)
            rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]['DR_FromMeans']=(np.nanmean(Y_F_mu,axis=1)-np.nanmean(Y_C_mu,axis=1))/np.nanmean(Y_C_mu,axis=1)*100

#plt.plot(rsT['FD']['Bsw_Net'][225]['DA_mu'],'ob-')
  
#------------------------------------------------------------------------------
# Plot
#------------------------------------------------------------------------------
            
spc='FD'

plt.close('all'); ms=2; Alpha=0.14
fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(15.5,10.5))
ax[0,0].plot([-5,42],[0,0],'k-')
ax[0,0].fill_between(tsf,rsT[spc]['Bsw_G'][225]['DA_mu']-rsT[spc]['Bsw_G'][225]['DA_se'],
  rsT[spc]['Bsw_G'][225]['DA_mu']+rsT[spc]['Bsw_G'][225]['DA_se'],color=[0.27,0.49,0.79],alpha=Alpha,linewidth=0)
ax[0,0].plot(tsf,rsT[spc]['Bsw_G'][225]['DA_mu'],'-o',color=[0.27,0.49,0.79],markersize=ms,label='225 kg N ha$^-$$^1$') 
ax[0,0].fill_between(tsf,rsT[spc]['Bsw_G'][450]['DA_mu']-rsT[spc]['Bsw_G'][450]['DA_se'],
  rsT[spc]['Bsw_G'][450]['DA_mu']+rsT[spc]['Bsw_G'][450]['DA_se'],color=[0,0.7,0],alpha=Alpha,linewidth=0)
ax[0,0].plot(tsf,rsT[spc]['Bsw_G'][450]['DA_mu'],'-s',color=[0,0.7,0],markersize=ms,label='450 kg N ha$^-$$^1$')
ax[0,0].fill_between(tsf,rsT[spc]['Bsw_G'][675]['DA_mu']-rsT[spc]['Bsw_G'][675]['DA_se'],
  rsT[spc]['Bsw_G'][675]['DA_mu']+rsT[spc]['Bsw_G'][675]['DA_se'],color=[0.5,0,1],alpha=Alpha,linewidth=0)
ax[0,0].plot(tsf,rsT[spc]['Bsw_G'][675]['DA_mu'],'-d',color=[0.5,0,1],markersize=ms,label='675 kg N ha$^-$$^1$')
ax[0,0].set(position=[0.07,0.57,0.42,0.41],xlim=[-1,41],xlabel='Time since N addition, years',ylim=[-1.8,2.2],ylabel='Gross growth response (Mg C ha$^-$$^1$ yr$^-$$^1$)')
ax[0,0].legend(loc='lower left')
ax[1,0].plot([-5,42],[0,0],'k-')
ax[1,0].fill_between(tsf,np.cumsum(rsT[spc]['Bsw_G'][225]['DA_mu']-rsT[spc]['Bsw_G'][225]['DA_se']),
  np.cumsum(rsT[spc]['Bsw_G'][225]['DA_mu']+rsT[spc]['Bsw_G'][225]['DA_se']),color=[0.27,0.49,0.79],alpha=Alpha,linewidth=0)
ax[1,0].plot(tsf,np.cumsum(rsT[spc]['Bsw_G'][225]['DA_mu']),'-o',color=[0.27,0.49,0.79],markersize=ms) 
ax[1,0].fill_between(tsf,np.cumsum(rsT[spc]['Bsw_G'][450]['DA_mu']-rsT[spc]['Bsw_G'][450]['DA_se']),
  np.cumsum(rsT[spc]['Bsw_G'][450]['DA_mu']+rsT[spc]['Bsw_G'][450]['DA_se']),color=[0,0.7,0],alpha=Alpha,linewidth=0)
ax[1,0].plot(tsf,np.cumsum(rsT[spc]['Bsw_G'][450]['DA_mu']),'-s',color=[0,0.7,0],markersize=ms)
ax[1,0].fill_between(tsf,np.cumsum(rsT[spc]['Bsw_G'][675]['DA_mu']-rsT[spc]['Bsw_G'][675]['DA_se']),
  np.cumsum(rsT[spc]['Bsw_G'][675]['DA_mu']+rsT[spc]['Bsw_G'][675]['DA_se']),color=[0.5,0,1],alpha=Alpha,linewidth=0)
ax[1,0].plot(tsf,np.cumsum(rsT[spc]['Bsw_G'][675]['DA_mu']),'-d',color=[0.5,0,1],markersize=ms)
ax[1,0].set(position=[0.07,0.07,0.42,0.41],xlim=[-1,41],xlabel='Time since N addition, years',ylim=[-15,15],ylabel='Gross growth response (Mg C ha$^-$$^1$ yr$^-$$^1$)')
ax[0,1].plot([-5,42],[0,0],'k-')
ax[0,1].fill_between(tsf,rsT[spc]['Bsw_M'][225]['DA_mu']-rsT[spc]['Bsw_M'][225]['DA_se'],
  rsT[spc]['Bsw_M'][225]['DA_mu']+rsT[spc]['Bsw_M'][225]['DA_se'],color=[0.27,0.49,0.79],alpha=Alpha,linewidth=0)
ax[0,1].plot(tsf,rsT[spc]['Bsw_M'][225]['DA_mu'],'-o',color=[0.27,0.49,0.79],markersize=ms) 
ax[0,1].fill_between(tsf,rsT[spc]['Bsw_M'][450]['DA_mu']-rsT[spc]['Bsw_M'][450]['DA_se'],
  rsT[spc]['Bsw_M'][450]['DA_mu']+rsT[spc]['Bsw_M'][450]['DA_se'],color=[0,0.7,0],alpha=Alpha,linewidth=0)
ax[0,1].plot(tsf,rsT[spc]['Bsw_M'][450]['DA_mu'],'-s',color=[0,0.7,0],markersize=ms)
ax[0,1].fill_between(tsf,rsT[spc]['Bsw_M'][675]['DA_mu']-rsT[spc]['Bsw_M'][675]['DA_se'],
  rsT[spc]['Bsw_M'][675]['DA_mu']+rsT[spc]['Bsw_M'][675]['DA_se'],color=[0.5,0,1],alpha=Alpha,linewidth=0)
ax[0,1].plot(tsf,rsT[spc]['Bsw_M'][675]['DA_mu'],'-d',color=[0.5,0,1],markersize=ms)
ax[0,1].set(position=[0.57,0.57,0.42,0.41],xlim=[-1,41],xlabel='Time since N addition, years',ylim=[-1.8,2.2],ylabel='Mortality response (Mg C ha$^-$$^1$ yr$^-$$^1$)')
ax[1,1].plot([-5,42],[0,0],'k-')
ax[1,1].fill_between(tsf,np.cumsum(rsT[spc]['Bsw_M'][225]['DA_mu']-rsT[spc]['Bsw_M'][225]['DA_se']),
  np.cumsum(rsT[spc]['Bsw_M'][225]['DA_mu']+rsT[spc]['Bsw_M'][225]['DA_se']),color=[0.27,0.49,0.79],alpha=Alpha,linewidth=0)
ax[1,1].plot(tsf,np.cumsum(rsT[spc]['Bsw_M'][225]['DA_mu']),'-o',color=[0.27,0.49,0.79],markersize=ms) 
ax[1,1].fill_between(tsf,np.cumsum(rsT[spc]['Bsw_M'][450]['DA_mu']-rsT[spc]['Bsw_M'][450]['DA_se']),
  np.cumsum(rsT[spc]['Bsw_M'][450]['DA_mu']+rsT[spc]['Bsw_M'][450]['DA_se']),color=[0,0.7,0],alpha=Alpha,linewidth=0)
ax[1,1].plot(tsf,np.cumsum(rsT[spc]['Bsw_M'][450]['DA_mu']),'-s',color=[0,0.7,0],markersize=ms)
ax[1,1].fill_between(tsf,np.cumsum(rsT[spc]['Bsw_M'][675]['DA_mu']-rsT[spc]['Bsw_M'][675]['DA_se']),
  np.cumsum(rsT[spc]['Bsw_M'][675]['DA_mu']+rsT[spc]['Bsw_M'][675]['DA_se']),color=[0.5,0,1],alpha=Alpha,linewidth=0)
ax[1,1].plot(tsf,np.cumsum(rsT[spc]['Bsw_M'][675]['DA_mu']),'-d',color=[0.5,0,1],markersize=ms)
ax[1,1].set(position=[0.57,0.07,0.42,0.41],xlim=[-1,41],xlabel='Time since N addition, years',ylim=[-15,15],ylabel='Mortality response (Mg C ha$^-$$^1$ yr$^-$$^1$)')
gu.axletters(ax,plt,0.03,0.92,['Annual','Annual','Cumulative','Cumulative'],0.05)

#gu.PrintFig(PathFigures + '\\DA_vs_TSF_' + spc,'png',900)

#------------------------------------------------------------------------------
# Table statistics for manuscript
#------------------------------------------------------------------------------

spc='FD';
dose=225;
vd='Bsw_Net'
vr='yF_mu'
mu=[]
#it=np.where( (tsf>=1) & (tsf<=3) )[0];mu.append(np.mean(rsT[spc][vd][dose][vr][it]))
#it=np.where( (tsf>=1) & (tsf<=6) )[0];mu.append(np.mean(rsT[spc][vd][dose][vr][it]))
#it=np.where( (tsf>=1) & (tsf<=9) )[0];mu.append(np.mean(rsT[spc][vd][dose][vr][it]))
#it=np.where( (tsf>=1) & (tsf<=12) )[0];mu.append(np.mean(rsT[spc][vd][dose][vr][it]))
it=np.where( (tsf>=1) & (tsf<=15) )[0];mu.append(np.mean(rsT[spc][vd][dose][vr][it]))
#it=np.where( (tsf>=1) & (tsf<=18) )[0];mu.append(np.mean(rsT[spc][vd][dose][vr][it]))
#it=np.where( (tsf>=1) & (tsf<=20) )[0];mu.append(np.mean(rsT[spc][vd][dose][vr][it]))
#it=np.where( (tsf>=1) & (tsf<=30) )[0];mu.append(np.mean(rsT[spc][vd][dose][vr][it]))
#it=np.where( (tsf>=1) & (tsf<=40) )[0];mu.append(np.mean(rsT[spc][vd][dose][vr][it]))
mu


#%% Nitrogen response of biomass fluxes vs. calendar year

uDose=np.array([225,450,675])

# Define a continuous annual TSF variable
tsf=np.arange(1971,2012,1)

uTSF=np.unique(np.column_stack((sobs['Year_t0'],sobs['Year_t1'])).astype(float),axis=0)
uTSF=uTSF[1:,:]

uSpc=['FD','HW']

uVar=['Bsw_G','Bsw_M','Bsw_Net']

rsT={}
for iSpc in range(len(uSpc)):
    rsT[uSpc[iSpc]]={}
    for iVar in range(len(uVar)):
        Y=sobs[uVar[iVar]]
        rsT[uSpc[iSpc]][uVar[iVar]]={}
        for iDose in range(uDose.size):
            print(iSpc,iVar,iDose)
            rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]={}            
            Y_C_mu=np.nan*np.ones((tsf.size,uSite.size))
            Y_F_mu=np.nan*np.ones((tsf.size,uSite.size))
            for iSite in range(uSite.size):
                for iTSF in range(uTSF.shape[0]):
                    iC=np.where( (sobs['ID_Site']==uSite[iSite]) & (sobs['Spc1_ID_t0']==uSpc[iSpc]) & (sobs['Year_t0']>=uTSF[iTSF,0]) & (sobs['Year_t1']<=uTSF[iTSF,1]) & (sobs['N_Dose']==0) )[0]
                    iF=np.where( (sobs['ID_Site']==uSite[iSite]) & (sobs['Spc1_ID_t0']==uSpc[iSpc]) & (sobs['Year_t0']>=uTSF[iTSF,0]) & (sobs['Year_t1']<=uTSF[iTSF,1]) & (sobs['N_Dose']==uDose[iDose]) )[0]
                    if (iC.size==0) | (iF.size==0): continue    
                    it=np.where( (tsf>=uTSF[iTSF,0]) & (tsf<=uTSF[iTSF,1]) )[0]
                    Y_C_mu[it,iSite]=np.nanmean(sobs[uVar[iVar]][iC])
                    Y_F_mu[it,iSite]=np.nanmean(sobs[uVar[iVar]][iF])
            
            # Exclude plots if they have no site-paired control/treatment
            ind=np.where((np.isnan(Y_C_mu+Y_F_mu)==True))
            Y_C_mu[ind]=np.nan
            Y_F_mu[ind]=np.nan
            
            # Calculate differences
            DA=Y_F_mu-Y_C_mu
            DR=DA/Y_C_mu*100
            
            # Summarize
            SampleSize=np.zeros(Y_C_mu.shape[0])
            for i in range(SampleSize.size):
                SampleSize[i]=np.where(np.isnan(Y_C_mu[i,:]+Y_F_mu[i,:])==False)[0].size
            rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]['N_Inst']=SampleSize
            rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]['yC_mu']=np.nanmean(Y_C_mu,axis=1)
            rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]['yF_mu']=np.nanmean(Y_F_mu,axis=1)
            rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]['DA_med']=np.nanmedian(DA,axis=1)
            rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]['DA_mu']=np.nanmean(DA,axis=1)
            rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]['DA_se']=np.nanstd(DA,axis=1)/np.sqrt(SampleSize)
            rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]['DR_med']=np.nanmedian(DR,axis=1)
            rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]['DR_mu']=np.nanmean(DR,axis=1)
            rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]['DR_se']=np.nanstd(DR,axis=1)/np.sqrt(SampleSize)
            rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]['DR_FromMeans']=(np.nanmean(Y_F_mu,axis=1)-np.nanmean(Y_C_mu,axis=1))/np.nanmean(Y_C_mu,axis=1)*100

#plt.plot(rsT['FD']['Bsw_Net'][225]['DA_mu'],'ob-')
  
#------------------------------------------------------------------------------
# Plot
#------------------------------------------------------------------------------
            
spc='FD'

plt.close('all'); ms=2; Alpha=0.14
fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(15.5,10.5))
ax[0,0].plot([1965,2015],[0,0],'k-')
ax[0,0].fill_between(tsf,rsT[spc]['Bsw_G'][225]['DA_mu']-rsT[spc]['Bsw_G'][225]['DA_se'],
  rsT[spc]['Bsw_G'][225]['DA_mu']+rsT[spc]['Bsw_G'][225]['DA_se'],color=[0.27,0.49,0.79],alpha=Alpha,linewidth=0)
ax[0,0].plot(tsf,rsT[spc]['Bsw_G'][225]['DA_mu'],'-o',color=[0.27,0.49,0.79],markersize=ms,label='225 kg N ha$^-$$^1$') 
ax[0,0].fill_between(tsf,rsT[spc]['Bsw_G'][450]['DA_mu']-rsT[spc]['Bsw_G'][450]['DA_se'],
  rsT[spc]['Bsw_G'][450]['DA_mu']+rsT[spc]['Bsw_G'][450]['DA_se'],color=[0,0.7,0],alpha=Alpha,linewidth=0)
ax[0,0].plot(tsf,rsT[spc]['Bsw_G'][450]['DA_mu'],'-s',color=[0,0.7,0],markersize=ms,label='450 kg N ha$^-$$^1$')
#ax[0,0].fill_between(tsf,rsT[spc]['Bsw_G'][675]['DA_mu']-rsT[spc]['Bsw_G'][675]['DA_se'],
#  rsT[spc]['Bsw_G'][675]['DA_mu']+rsT[spc]['Bsw_G'][675]['DA_se'],color=[0.5,0,1],alpha=Alpha,linewidth=0)
#ax[0,0].plot(tsf,rsT[spc]['Bsw_G'][675]['DA_mu'],'-d',color=[0.5,0,1],markersize=ms,label='675 kg N ha$^-$$^1$')
ax[0,0].set(position=[0.07,0.57,0.42,0.41],xlim=[1970,2012],xlabel='Time since N addition, years',ylim=[-1.8,2.2],ylabel='Gross growth response (Mg C ha$^-$$^1$ yr$^-$$^1$)')
ax[0,0].legend(loc='lower left')

ax[1,0].plot([1965,2015],[0,0],'k-')
ax[1,0].fill_between(tsf,np.cumsum(rsT[spc]['Bsw_G'][225]['DA_mu']-rsT[spc]['Bsw_G'][225]['DA_se']),
  np.cumsum(rsT[spc]['Bsw_G'][225]['DA_mu']+rsT[spc]['Bsw_G'][225]['DA_se']),color=[0.27,0.49,0.79],alpha=Alpha,linewidth=0)
ax[1,0].plot(tsf,np.cumsum(rsT[spc]['Bsw_G'][225]['DA_mu']),'-o',color=[0.27,0.49,0.79],markersize=ms) 
ax[1,0].fill_between(tsf,np.cumsum(rsT[spc]['Bsw_G'][450]['DA_mu']-rsT[spc]['Bsw_G'][450]['DA_se']),
  np.cumsum(rsT[spc]['Bsw_G'][450]['DA_mu']+rsT[spc]['Bsw_G'][450]['DA_se']),color=[0,0.7,0],alpha=Alpha,linewidth=0)
ax[1,0].plot(tsf,np.cumsum(rsT[spc]['Bsw_G'][450]['DA_mu']),'-s',color=[0,0.7,0],markersize=ms)
ax[1,0].fill_between(tsf,np.cumsum(rsT[spc]['Bsw_G'][675]['DA_mu']-rsT[spc]['Bsw_G'][675]['DA_se']),
  np.cumsum(rsT[spc]['Bsw_G'][675]['DA_mu']+rsT[spc]['Bsw_G'][675]['DA_se']),color=[0.5,0,1],alpha=Alpha,linewidth=0)
ax[1,0].plot(tsf,np.cumsum(rsT[spc]['Bsw_G'][675]['DA_mu']),'-d',color=[0.5,0,1],markersize=ms)
ax[1,0].set(position=[0.07,0.07,0.42,0.41],xlim=[1970,2012],xlabel='Time since N addition, years',ylim=[-15,15],ylabel='Gross growth response (Mg C ha$^-$$^1$ yr$^-$$^1$)')
ax[0,1].plot([1965,2015],[0,0],'k-')
ax[0,1].fill_between(tsf,rsT[spc]['Bsw_M'][225]['DA_mu']-rsT[spc]['Bsw_M'][225]['DA_se'],
  rsT[spc]['Bsw_M'][225]['DA_mu']+rsT[spc]['Bsw_M'][225]['DA_se'],color=[0.27,0.49,0.79],alpha=Alpha,linewidth=0)
ax[0,1].plot(tsf,rsT[spc]['Bsw_M'][225]['DA_mu'],'-o',color=[0.27,0.49,0.79],markersize=ms) 
#ax[0,1].fill_between(tsf,rsT[spc]['Bsw_M'][450]['DA_mu']-rsT[spc]['Bsw_M'][450]['DA_se'],
#  rsT[spc]['Bsw_M'][450]['DA_mu']+rsT[spc]['Bsw_M'][450]['DA_se'],color=[0,0.7,0],alpha=Alpha,linewidth=0)
#ax[0,1].plot(tsf,rsT[spc]['Bsw_M'][450]['DA_mu'],'-s',color=[0,0.7,0],markersize=ms)
#ax[0,1].fill_between(tsf,rsT[spc]['Bsw_M'][675]['DA_mu']-rsT[spc]['Bsw_M'][675]['DA_se'],
#  rsT[spc]['Bsw_M'][675]['DA_mu']+rsT[spc]['Bsw_M'][675]['DA_se'],color=[0.5,0,1],alpha=Alpha,linewidth=0)
#ax[0,1].plot(tsf,rsT[spc]['Bsw_M'][675]['DA_mu'],'-d',color=[0.5,0,1],markersize=ms)
ax[0,1].set(position=[0.57,0.57,0.42,0.41],xlim=[1970,2012],xlabel='Time since N addition, years',ylim=[-1.8,2.2],ylabel='Mortality response (Mg C ha$^-$$^1$ yr$^-$$^1$)')
ax[1,1].plot([1965,2015],[0,0],'k-')
ax[1,1].fill_between(tsf,np.cumsum(rsT[spc]['Bsw_M'][225]['DA_mu']-rsT[spc]['Bsw_M'][225]['DA_se']),
  np.cumsum(rsT[spc]['Bsw_M'][225]['DA_mu']+rsT[spc]['Bsw_M'][225]['DA_se']),color=[0.27,0.49,0.79],alpha=Alpha,linewidth=0)
ax[1,1].plot(tsf,np.cumsum(rsT[spc]['Bsw_M'][225]['DA_mu']),'-o',color=[0.27,0.49,0.79],markersize=ms) 
#ax[1,1].fill_between(tsf,np.cumsum(rsT[spc]['Bsw_M'][450]['DA_mu']-rsT[spc]['Bsw_M'][450]['DA_se']),
#  np.cumsum(rsT[spc]['Bsw_M'][450]['DA_mu']+rsT[spc]['Bsw_M'][450]['DA_se']),color=[0,0.7,0],alpha=Alpha,linewidth=0)
#ax[1,1].plot(tsf,np.cumsum(rsT[spc]['Bsw_M'][450]['DA_mu']),'-s',color=[0,0.7,0],markersize=ms)
#ax[1,1].fill_between(tsf,np.cumsum(rsT[spc]['Bsw_M'][675]['DA_mu']-rsT[spc]['Bsw_M'][675]['DA_se']),
#  np.cumsum(rsT[spc]['Bsw_M'][675]['DA_mu']+rsT[spc]['Bsw_M'][675]['DA_se']),color=[0.5,0,1],alpha=Alpha,linewidth=0)
#ax[1,1].plot(tsf,np.cumsum(rsT[spc]['Bsw_M'][675]['DA_mu']),'-d',color=[0.5,0,1],markersize=ms)
ax[1,1].set(position=[0.57,0.07,0.42,0.41],xlim=[1970,2012],xlabel='Time since N addition, years',ylim=[-15,15],ylabel='Mortality response (Mg C ha$^-$$^1$ yr$^-$$^1$)')
gu.axletters(ax,plt,0.03,0.92,['Annual','Annual','Cumulative','Cumulative'],0.05)

#gu.PrintFig(PathFigures + '\\DA_vs_TSF_' + spc,'png',900)







#%% Nitrogen response vs. time since fert (site by site)

uDose=np.array([225,450,675])

# Define a continuous annual TSF variable
tsf=np.arange(-1,41,1)

uTSF=np.unique(np.column_stack((sobs['TSF_t0_adj'],sobs['TSF_t1'])).astype(float),axis=0)
#uTSF=uTSF[1:,:]

uVar=['Bsw_G','Bsw_M','Bsw_Net','N_t0','Bsw_t0','H_gf_t0','BA_t0']

rsT={}
for iVar in range(len(uVar)):
    Y=sobs[uVar[iVar]]
    rsT[uVar[iVar]]={}
    for iDose in range(uDose.size):
        rsT[uVar[iVar]][uDose[iDose]]={}            
        Y_C_mu=np.nan*np.ones((tsf.size,uSite.size))
        Y_F_mu=np.nan*np.ones((tsf.size,uSite.size))
        for iSite in range(uSite.size):
            for iTSF in range(uTSF.shape[0]):
                iC=np.where( (sobs['Spc1_ID_t0']!='X') & (sobs['ID_Site']==uSite[iSite]) & (sobs['TSF_t0_adj']>=uTSF[iTSF,0]) & (sobs['TSF_t1']<=uTSF[iTSF,1]) & (sobs['N_Dose']==0) )[0]
                iF=np.where( (sobs['Spc1_ID_t0']!='X') & (sobs['ID_Site']==uSite[iSite]) & (sobs['TSF_t0_adj']>=uTSF[iTSF,0]) & (sobs['TSF_t1']<=uTSF[iTSF,1]) & (sobs['N_Dose']==uDose[iDose]) )[0]
                if (iC.size==0) | (iF.size==0): continue    
                it=np.where( (tsf>=uTSF[iTSF,0]) & (tsf<=uTSF[iTSF,1]) )[0]
                Y_C_mu[it,iSite]=np.nanmean(sobs[uVar[iVar]][iC])
                Y_F_mu[it,iSite]=np.nanmean(sobs[uVar[iVar]][iF])
            
        # Exclude plots if they have no site-paired control/treatment
        ind=np.where((np.isnan(Y_C_mu+Y_F_mu)==True))
        Y_C_mu[ind]=np.nan
        Y_F_mu[ind]=np.nan
            
        # Calculate differences
        DA=Y_F_mu-Y_C_mu
        DR=DA/Y_C_mu*100
            
        # Summarize
        SampleSize=np.zeros(Y_C_mu.shape[0])
        for i in range(SampleSize.size):
            SampleSize[i]=np.where(np.isnan(Y_C_mu[i,:]+Y_F_mu[i,:])==False)[0].size
        rsT[uVar[iVar]][uDose[iDose]]['N_Inst']=SampleSize
        rsT[uVar[iVar]][uDose[iDose]]['yC_mu']=Y_C_mu
        rsT[uVar[iVar]][uDose[iDose]]['yF_mu']=Y_F_mu
        rsT[uVar[iVar]][uDose[iDose]]['DA_mu']=DA
        rsT[uVar[iVar]][uDose[iDose]]['DR_mu']=DR
   
# Get leading species
LeadSpc=np.zeros(uSite.size)
for i in range(uSite.size):
    ind=np.where(sobs['ID_Site']==uSite[i])[0]
    if ind.size==0:
        continue
    if sobs['Spc1_ID_t0'][ind[0]]=='FD':
        LeadSpc[i]=1
    else:
        LeadSpc[i]=2

#------------------------------------------------------------------------------     
# Plot each installation for tree coring site selection
#------------------------------------------------------------------------------

Dose=225
for i in range(len(LeadSpc)):
    #if LeadSpc[i]!=1:
    #    continue
    plt.close('all')
    fig,ax=plt.subplots(3,3,figsize=gu.cm2inch(18,14))
    ax[0,0].plot(tsf,rsT['N_t0'][Dose]['yC_mu'][:,i],'-b.')
    ax[0,0].plot(tsf,rsT['N_t0'][Dose]['yF_mu'][:,i],'-g.')
    ax[0,0].set(ylabel='Stand density (stems/ha)')
    ax[0,1].plot(tsf,rsT['H_gf_t0'][Dose]['yC_mu'][:,i],'-b.')
    ax[0,1].plot(tsf,rsT['H_gf_t0'][Dose]['yF_mu'][:,i],'-g.')
    ax[0,1].set(ylabel='Height (m)')
    ax[0,2].plot(tsf,rsT['BA_t0'][Dose]['yC_mu'][:,i],'-b.')
    ax[0,2].plot(tsf,rsT['BA_t0'][Dose]['yF_mu'][:,i],'-g.')
    ax[0,2].set(ylabel='Basal area (m2/ha)')
    ind=np.where(np.isnan(rsT['Bsw_t0'][Dose]['yF_mu'][:,i])==False)[0]
    try:
        d=rsT['Bsw_t0'][Dose]['yF_mu'][ind[0],i]-rsT['Bsw_t0'][Dose]['yC_mu'][ind[0],i]
    except:
        d=0
    d=0
    ax[1,0].plot(tsf,rsT['Bsw_t0'][Dose]['yC_mu'][:,i],'-b.')
    ax[1,0].plot(tsf,rsT['Bsw_t0'][Dose]['yF_mu'][:,i]-d,'-g.')
    ax[1,0].set(ylabel='Stemwood biomass (MgC/ha)')
    
    ax[1,1].plot(tsf,rsT['Bsw_G'][Dose]['yC_mu'][:,i],'-b.')
    ax[1,1].plot(tsf,rsT['Bsw_G'][Dose]['yF_mu'][:,i],'-g.')
    ax[1,1].set(ylabel='Stemwood gross growth (MgC/ha/yr)')
    ax[1,2].plot(tsf,rsT['Bsw_Net'][Dose]['yC_mu'][:,i],'-b.')
    ax[1,2].plot(tsf,rsT['Bsw_Net'][Dose]['yF_mu'][:,i],'-g.')
    ax[1,2].set(ylabel='Stemwood net growth (MgC/ha/yr)')
    
    ax[2,0].plot(tsf,0*np.ones(tsf.size),'-k')
    ax[2,0].plot(tsf,rsT['Bsw_G'][Dose]['yF_mu'][:,i]-rsT['Bsw_G'][Dose]['yC_mu'][:,i],'-b.')
    ax[2,0].set(ylabel='Delta gross growth (MgC/ha/yr)')
    
    ax[2,1].plot(tsf,0*np.ones(tsf.size),'-k')
    ax[2,1].plot(tsf,rsT['Bsw_M'][Dose]['yF_mu'][:,i]-rsT['Bsw_M'][Dose]['yC_mu'][:,i],'-b.')
    ax[2,1].set(ylabel='Delta mortality (MgC/ha/yr)')
    
    ax[2,2].plot(tsf,0*np.ones(tsf.size),'-k')
    ax[2,2].plot(tsf,rsT['Bsw_Net'][Dose]['yF_mu'][:,i]-rsT['Bsw_Net'][Dose]['yC_mu'][:,i],'-b.')
    ax[2,2].set(ylabel='Delta net growth (MgC/ha/yr)')
    plt.tight_layout()
    fig.savefig(PathFigures + '\\ForCoringProject\\TS_By_Inst_' + str(uSite[i]) + '.png',format='png',dpi=300)

#------------------------------------------------------------------------------
# Get coordinates for tree coring site selectoin
#------------------------------------------------------------------------------

a=[1,3,4,5,10,11,12,14,16,25,38,41,42,47,71,72]
ll=np.zeros((len(a),2))
for i in range(len(a)):
    ind=np.where(sobs['ID_Site']==a[i])[0]
    ll[i,0]=sobs['Lat'][ind[0]]
    ll[i,1]=sobs['Lon'][ind[0]]


#%% Look at relationship between initial biomass and delta growth

Dose=450
iTSF=np.where( (uTSF[:,0]>0) & (uTSF[:,0]<20) )[0]
cols=['dB','dG','dNet','rdB','rdG','rdNet']
df=pd.DataFrame(data=np.nan*np.ones((len(LeadSpc),len(cols))),columns=cols)
for i in range(len(LeadSpc)):
    df.loc[i,'dB']=np.nanmean(rsT['Bsw_t0'][Dose]['yF_mu'][0:2,i]-rsT['Bsw_t0'][Dose]['yC_mu'][0:2,i])
    df.loc[i,'rdB']=np.nanmean( (rsT['Bsw_t0'][Dose]['yF_mu'][0:2,i]-rsT['Bsw_t0'][Dose]['yC_mu'][0:2,i])/rsT['Bsw_t0'][Dose]['yC_mu'][0:2,i]*100 )
    
    df.loc[i,'dG']=np.nanmean(rsT['Bsw_G'][Dose]['yF_mu'][iTSF,i]-rsT['Bsw_G'][Dose]['yC_mu'][iTSF,i])
    df.loc[i,'dNet']=np.nanmean(rsT['Bsw_Net'][Dose]['yF_mu'][iTSF,i]-rsT['Bsw_Net'][Dose]['yC_mu'][iTSF,i])
    
    df.loc[i,'rdG']=np.nanmean( (rsT['Bsw_G'][Dose]['yF_mu'][iTSF,i]-rsT['Bsw_G'][Dose]['yC_mu'][iTSF,i])/rsT['Bsw_G'][Dose]['yC_mu'][iTSF,i]*100 )
    df.loc[i,'rdNet']=np.nanmean( (rsT['Bsw_Net'][Dose]['yF_mu'][iTSF,i]-rsT['Bsw_Net'][Dose]['yC_mu'][iTSF,i])/rsT['Bsw_Net'][Dose]['yC_mu'][iTSF,i]*100 )


df=df[np.isnan(df['dB'])==False]
df=df.reset_index(drop=True)

plt.close('all')
fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(13,11))
y=df['dG'].values.astype(float); x=df['dB'].values.astype(float);x=sm.tools.tools.add_constant(x)
md=sm.OLS(y,x).fit()
xhat=np.linspace(np.min(x[:,1]),np.max(x[:,1]),10);yhat=md.predict(np.c_[np.ones(xhat.size),xhat])
ax[0,0].plot(xhat,np.zeros(xhat.size),'k-')
ax[0,0].plot(np.zeros(yhat.size),yhat,'k-')
ax[0,0].plot(x[:,1],y,'ko',mec=[1,1,1],mfc=[0,0,0],ms=4)
ax[0,0].plot(xhat,yhat,'r-')
ax[0,0].set(position=[0.075,0.57,0.42,0.4],ylabel='$\Delta$ gross growth (Mg C ha$^-$$^1$ yr$^-$$^1$)',xlabel='$\Delta$ initial biomass (Mg C ha$^-$$^1$)')

y=df['rdG'].values.astype(float); x=df['rdB'].values.astype(float);x=sm.tools.tools.add_constant(x)
md=sm.OLS(y,x).fit()
xhat=np.linspace(np.min(x[:,1]),np.max(x[:,1]),10);yhat=md.predict(np.c_[np.ones(xhat.size),xhat])
ax[0,1].plot(xhat,np.zeros(xhat.size),'k-')
ax[0,1].plot(np.zeros(yhat.size),yhat,'k-')
ax[0,1].plot(x[:,1],y,'ko',mec=[1,1,1],mfc=[0,0,0],ms=4)
ax[0,1].plot(xhat,yhat,'r-')
ax[0,1].set(position=[0.57,0.57,0.42,0.4],ylabel='$\Delta$ gross growth (%)',xlabel='$\Delta$ initial biomass (%)')

y=df['dNet'].values.astype(float); x=df['dB'].values.astype(float);x=sm.tools.tools.add_constant(x)
md=sm.OLS(y,x).fit()
xhat=np.linspace(np.min(x[:,1]),np.max(x[:,1]),10);yhat=md.predict(np.c_[np.ones(xhat.size),xhat])
ax[1,0].plot(xhat,np.zeros(xhat.size),'k-')
ax[1,0].plot(np.zeros(yhat.size),yhat,'k-')
ax[1,0].plot(x[:,1],y,'ko',mec=[1,1,1],mfc=[0,0,0],ms=4)
ax[1,0].plot(xhat,yhat,'r-')
ax[1,0].set(position=[0.075,0.07,0.42,0.4],ylabel='$\Delta$ net growth (Mg C ha$^-$$^1$ yr$^-$$^1$)',xlabel='$\Delta$ initial biomass (Mg C ha$^-$$^1$)')

y=df['rdNet'].values.astype(float); x=df['rdB'].values.astype(float);x=sm.tools.tools.add_constant(x)
md=sm.OLS(y,x).fit()
xhat=np.linspace(np.min(x[:,1]),np.max(x[:,1]),10);yhat=md.predict(np.c_[np.ones(xhat.size),xhat])
ax[1,1].plot(xhat,np.zeros(xhat.size),'k-')
ax[1,1].plot(np.zeros(yhat.size),yhat,'k-')
ax[1,1].plot(x[:,1],y,'ko',mec=[1,1,1],mfc=[0,0,0],ms=4)
ax[1,1].plot(xhat,yhat,'r-')
ax[1,1].set(position=[0.57,0.07,0.42,0.4],ylabel='$\Delta$ net growth (%)',xlabel='$\Delta$ initial biomass (%)')

for i in range(0,2):
    for j in range(0,2):
        ax[i,j].yaxis.set_ticks_position('both'); ax[i,j].xaxis.set_ticks_position('both')
gu.axletters(ax,plt,0.03,0.92)

#gu.PrintFig(PathFigures + '\\SL_DeltaG_Vs_InitB','png',dpi=900)


#%% Nonlinear model

from scipy.optimize import curve_fit

# Function (tree height as a function of basal area)
def func(x,a,b,c):
    y=a+(b-a)*np.exp(-c*x)
    return y
xhat=np.arange(0,20,1)
plt.plot(xhat,func(xhat,0.1,1.7,0.3),'o-')

y=rsT['Bsw_M'][450]['DA_mu']
m,n=y.shape


dose=225
#iSpc=np.where( (LeadSpc==1) & (np.isnan(rsT['Bsw_Net'][dose]['yC_mu'][40,:])==False) )[0]
iSpc=np.where( (LeadSpc>=0) )[0]
plt.close('all')
y=np.nanmean(rsT['Bsw_G'][dose]['DA_mu'][:,iSpc],axis=1)
plt.plot(tsf,y,'-o',linewidth=2)
y=np.nanmean(rsT['Bsw_M'][dose]['DA_mu'][:,iSpc],axis=1)
plt.plot(tsf,y,'-o',linewidth=2)
y=np.nanmean(rsT['Bsw_Net'][dose]['DA_mu'][:,iSpc],axis=1)
plt.plot(tsf,y,'-o',linewidth=2,color=[0,0,0])


dose=225
ikp=np.where( (np.isnan(rsT['Bsw_Net'][dose]['DA_mu'][40,:])==False) )[0]
plt.close('all')
y=np.nanmean(rsT['Bsw_G'][dose]['DA_mu'][:,ikp],axis=1)
plt.plot(y,linewidth=2)
y=np.nanmean(rsT['Bsw_M'][dose]['DA_mu'][:,ikp],axis=1)
plt.plot(y,linewidth=2)
y=np.nanmean(rsT['Bsw_Net'][dose]['DA_mu'][:,ikp],axis=1)
plt.plot(y,linewidth=2,color=[0,0,0])


plt.close('all')
y=np.nan_to_num(np.nanmean(rsT['Bsw_G'][dose]['DA_mu'][:,ikp],axis=1))
plt.plot(np.cumsum(y),linewidth=2)
y=np.nan_to_num(np.nanmean(rsT['Bsw_M'][dose]['DA_mu'][:,ikp],axis=1))
plt.plot(np.cumsum(y),linewidth=2)
y=np.nan_to_num(np.nanmean(rsT['Bsw_Net'][dose]['DA_mu'][:,ikp],axis=1))
plt.plot(np.cumsum(y),linewidth=2,color=[0,0,0])


y=np.nan_to_num(np.nanmean(rsT['Bsw_G'][dose]['DA_mu'][:,iSpc],axis=1))-np.nan_to_num(np.nanmean(rsT['Bsw_M'][dose]['DA_mu'][:,iSpc],axis=1))
plt.plot(np.cumsum(y),'--',linewidth=2,color=[1,0,0])

tsf1=np.tile(np.reshape(tsf,(m,1)),(1,n))
tsf1.shape
y=y.flatten()
tsf1=tsf1.flatten()
ikp=np.where( (np.isnan(y)!=1) )[0]

p,pcov=curve_fit(func,tsf1[ikp],y[ikp],[0.1,1.6,0.2])
print(p)

plt.close('all')
plt.plot(xhat,func(xhat,p[0],p[1],p[2]),'o-')


#%% Nitrogen response of height and basal area increment vs. TSF

uVar=['BA_G','BA_t0','H_obs_G','H_gf_G','H_obs_t0','H_gf_t0','dN_rel']

rsT={}
for iSpc in range(len(uSpc)):
    rsT[uSpc[iSpc]]={}
    for iVar in range(len(uVar)):
        Y=sobs[uVar[iVar]]
        rsT[uSpc[iSpc]][uVar[iVar]]={}
        for iDose in range(uDose.size):
            print(iSpc,iVar,iDose)
            rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]={}            
            Y_C_mu=np.nan*np.ones((tsf.size,uSite.size))
            Y_F_mu=np.nan*np.ones((tsf.size,uSite.size))
            for iSite in range(uSite.size):
                for iTSF in range(uTSF.shape[0]):
                    iC=np.where( (sobs['ID_Site']==uSite[iSite]) & (sobs['Spc1_ID_t0']==uSpc[iSpc]) & (sobs['TSF_t0_adj']>=uTSF[iTSF,0]) & (sobs['TSF_t1']<=uTSF[iTSF,1]) & (sobs['N_Dose']==0) )[0]
                    iF=np.where( (sobs['ID_Site']==uSite[iSite]) & (sobs['Spc1_ID_t0']==uSpc[iSpc]) & (sobs['TSF_t0_adj']>=uTSF[iTSF,0]) & (sobs['TSF_t1']<=uTSF[iTSF,1]) & (sobs['N_Dose']==uDose[iDose]) )[0]
                    if (iC.size==0) | (iF.size==0): 
                        continue    
                    it=np.where( (tsf>=uTSF[iTSF,0]) & (tsf<=uTSF[iTSF,1]) )[0]
                    Y_C_mu[it,iSite]=np.mean(sobs[uVar[iVar]][iC])
                    Y_F_mu[it,iSite]=np.mean(sobs[uVar[iVar]][iF])
            # Exclude plots if they have no site-paired control/treatment
            ind=np.where((np.isnan(Y_C_mu+Y_F_mu)==True))
            Y_C_mu[ind]=np.nan
            Y_F_mu[ind]=np.nan
            # Calculate differences
            DA=Y_F_mu-Y_C_mu
            DR=DA/Y_C_mu*100
            
            # Summarize
            SampleSize=np.zeros(Y_C_mu.shape[0])
            for i in range(SampleSize.size):                
                SampleSize[i]=np.where(np.isnan(Y_C_mu[i,:]+Y_F_mu[i,:])==False)[0].size
            rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]['N_Inst']=SampleSize
            rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]['yC_mu']=np.nanmean(Y_C_mu,axis=1)
            rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]['yF_mu']=np.nanmean(Y_F_mu,axis=1)
            rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]['DA_mu']=np.nanmean(DA,axis=1)
            rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]['DA_se']=np.nanstd(DA,axis=1)/np.sqrt(SampleSize)
            rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]['DR_mu']=np.nanmean(DR,axis=1)
            rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]['DR_se']=np.nanstd(DR,axis=1)/np.sqrt(SampleSize)
            rsT[uSpc[iSpc]][uVar[iVar]][uDose[iDose]]['DR_FromMeans']=(np.nanmean(Y_F_mu,axis=1)-np.nanmean(Y_C_mu,axis=1))/np.nanmean(Y_C_mu,axis=1)*100


# Tabular summary
spc='HW';
dose=225;
vd='dN_rel'
vr='DA_mu'
mu=[]
it=np.where( (tsf>=1) & (tsf<=3) )[0];mu.append(np.mean(rsT[spc][vd][dose][vr][it]))
it=np.where( (tsf>=1) & (tsf<=6) )[0];mu.append(np.mean(rsT[spc][vd][dose][vr][it]))
it=np.where( (tsf>=1) & (tsf<=9) )[0];mu.append(np.mean(rsT[spc][vd][dose][vr][it]))
it=np.where( (tsf>=1) & (tsf<=12) )[0];mu.append(np.mean(rsT[spc][vd][dose][vr][it]))
it=np.where( (tsf>=1) & (tsf<=20) )[0];mu.append(np.mean(rsT[spc][vd][dose][vr][it]))
mu


#%% Nitrogen response stratified by crown class

uDose=np.unique(sobs['N_Dose'])
uDose=uDose[1:]
#uDose=225*np.ones(1)

#uTSF=np.unique(np.column_stack((sobs['TSF_t0'],sobs['TSF_t1'])).astype(float),axis=0)
#uTSF=uTSF[1:,:]

#uTSF=np.unique(np.column_stack((sobs['TSF_t0_adj'],sobs['TSF_t1'])).astype(float),axis=0)
#ind=np.where( (uTSF[:,0]>0) & (uTSF[:,0]<=1000) )[0]
#uTSF=uTSF[ind,:]
uTSF=np.array([3,6,9,12,15,18])

#FullResponse=np.c_[dfR['ID_Site'].values,dfR['DA_mu'].values]

spc=['FD','HW']
vr=['Bsw_G_cc1','Bsw_G_cc2','Bsw_G_cc3','Bsw_G_cc4',
    'Bsw_M_cc1','Bsw_M_cc2','Bsw_M_cc3','Bsw_M_cc4',
    'Bsw_Net_cc1','Bsw_Net_cc2','Bsw_Net_cc3','Bsw_Net_cc4']
rsD={}
for iSpc in range(len(spc)):
    rsD[spc[iSpc]]={}
    for iVar in range(len(vr)):
        rsD[spc[iSpc]][vr[iVar]]={}
        for iTSF in range(uTSF.size):
            rsD[spc[iSpc]][vr[iVar]][iTSF]={}            
            print(iSpc,iVar,iTSF)            
            yCmu0=np.nan*np.ones((uDose.size,uSite.size))
            yFmu0=np.nan*np.ones((uDose.size,uSite.size))
            for iSite in range(uSite.size):
                
                # Exclude responders
                #ind=np.where(FullResponse[:,0]==uSite[iSite])[0]
                #if FullResponse[ind,1]>0.05: continue
            
                for iDose in range(uDose.shape[0]):
                    iC=np.where( (sobs['ID_Site']==uSite[iSite]) & (sobs['Spc1_ID_t0']==spc[iSpc]) & (sobs['TSF_t0_adj']>=1) & (sobs['TSF_t1']<=uTSF[iTSF]) & (sobs['N_Dose']==0) )[0]
                    iF=np.where( (sobs['ID_Site']==uSite[iSite]) & (sobs['Spc1_ID_t0']==spc[iSpc]) & (sobs['TSF_t0_adj']>=1) & (sobs['TSF_t1']<=uTSF[iTSF]) & (sobs['N_Dose']==uDose[iDose]) )[0]
                    if (iC.size==0) | (iF.size==0): continue                        
                    yCmu0[iDose,iSite]=np.mean(sobs[vr[iVar]][iC])
                    yFmu0[iDose,iSite]=np.mean(sobs[vr[iVar]][iF])
            # Exclude plots if they have no site-paired control/treatment
            ind=np.where((np.isnan(yCmu0+yFmu0)==True))
            yCmu0[ind]=np.nan
            yFmu0[ind]=np.nan
            # Calculate differences
            DA=yFmu0-yCmu0
            DR=DA/yCmu0*100
            # Summarize
            SampleSize=np.zeros(yCmu0.shape[0])
            for i in range(SampleSize.size): 
                SampleSize[i]=np.where(np.isnan(yCmu0[i,:]+yFmu0[i,:])==False)[0].size
            rsD[spc[iSpc]][vr[iVar]][iTSF]['Cmu']=np.nanmean(yCmu0,axis=1)
            rsD[spc[iSpc]][vr[iVar]][iTSF]['Fmu']=np.nanmean(yFmu0,axis=1)            
            rsD[spc[iSpc]][vr[iVar]][iTSF]['DA_mu']=np.nanmean(DA,axis=1)
            rsD[spc[iSpc]][vr[iVar]][iTSF]['DA_se']=np.nanstd(DA,axis=1)/np.sqrt(SampleSize)
            rsD[spc[iSpc]][vr[iVar]][iTSF]['DR_mu']=np.nanmean(DR,axis=1)
            rsD[spc[iSpc]][vr[iVar]][iTSF]['DR_se']=np.nanstd(DR,axis=1)/np.sqrt(SampleSize)

#------------------------------------------------------------------------------
# Plot
#------------------------------------------------------------------------------

# Remember this is mean response from TSF=1 on (i.e., 1 to iTSF)
iTSF=np.where(uTSF==15)[0][0]

plt.close('all'); ylim=[-0.3,0.8]; ah=0.26; ms=2.5; Alpha=0.18; bw=0.28; cl1=[0.24,0.49,0.77]; cl2=[0.6,1,0]; cl3=[0.5,0,1]
fig,ax=plt.subplots(3,2,figsize=gu.cm2inch(15.5,13.5))
spc='FD'; dose=225; y=[]; ye=[]
for i in range(4):
    y.append(rsD[spc]['Bsw_G_cc' + str(i+1)][iTSF]['DA_mu'][np.where(uDose==dose)[0]][0])
    ye.append(rsD[spc]['Bsw_G_cc' + str(i+1)][iTSF]['DA_se'][np.where(uDose==dose)[0]][0])
ax[0,0].bar(np.arange(1,5,1)-0.31,y,bw,yerr=ye,facecolor=cl1,capsize=2.5,ecolor=list(0.5*np.array(cl1)))
dose=450; y=[]; ye=[]
for i in range(4):
    y.append(rsD[spc]['Bsw_G_cc' + str(i+1)][iTSF]['DA_mu'][np.where(uDose==dose)[0]][0])
    ye.append(rsD[spc]['Bsw_G_cc' + str(i+1)][iTSF]['DA_se'][np.where(uDose==dose)[0]][0])
ax[0,0].bar(np.arange(1,5,1),y,bw,yerr=ye,facecolor=cl2,capsize=2.5,ecolor=list(0.5*np.array(cl2)))
dose=675; y=[]; ye=[]
for i in range(4):
    y.append(rsD[spc]['Bsw_G_cc' + str(i+1)][iTSF]['DA_mu'][np.where(uDose==dose)[0]][0])
    ye.append(rsD[spc]['Bsw_G_cc' + str(i+1)][iTSF]['DA_se'][np.where(uDose==dose)[0]][0])
ax[0,0].bar(np.arange(1,5,1)+0.31,y,bw,yerr=ye,facecolor=cl3,capsize=2.5,ecolor=list(0.5*np.array(cl3)))
ax[0,0].legend(['225 kg N ha$^-$$^1$','450 kg N ha$^-$$^1$','675 kg N ha$^-$$^1$'],loc=1,frameon=False)
ax[0,0].plot([0.4,4.6],[0,0],'k-')
ax[0,0].set(position=[0.0725,0.715,0.42,ah],xlim=[0.4,4.6],xticks=np.arange(1,5),xticklabels=['Dominant','Co-dominant','Intermediate','Suppressed'],xlabel='Crown class',ylim=ylim,ylabel='Gross growth (Mg C ha$^-$$^1$ yr$^-$$^1$)')

spc='HW'; dose=225; y=[]; ye=[]
for i in range(4):
    y.append(rsD[spc]['Bsw_G_cc' + str(i+1)][iTSF]['DA_mu'][np.where(uDose==dose)[0]][0])
    ye.append(rsD[spc]['Bsw_G_cc' + str(i+1)][iTSF]['DA_se'][np.where(uDose==dose)[0]][0])
ax[0,1].bar(np.arange(1,5,1)-0.31,y,bw,yerr=ye,facecolor=cl1,capsize=2.5,ecolor=list(0.5*np.array(cl1)))
dose=450; y=[]; ye=[]
for i in range(4):
    y.append(rsD[spc]['Bsw_G_cc' + str(i+1)][iTSF]['DA_mu'][np.where(uDose==dose)[0]][0])
    ye.append(rsD[spc]['Bsw_G_cc' + str(i+1)][iTSF]['DA_se'][np.where(uDose==dose)[0]][0])
ax[0,1].bar(np.arange(1,5,1),y,bw,yerr=ye,facecolor=cl2,capsize=2.5,ecolor=list(0.5*np.array(cl2)))
dose=675; y=[]; ye=[]
for i in range(4):
    y.append(rsD[spc]['Bsw_G_cc' + str(i+1)][iTSF]['DA_mu'][np.where(uDose==dose)[0]][0])
    ye.append(rsD[spc]['Bsw_G_cc' + str(i+1)][iTSF]['DA_se'][np.where(uDose==dose)[0]][0])
ax[0,1].bar(np.arange(1,5,1)+0.31,y,bw,yerr=ye,facecolor=cl3,capsize=2.5,ecolor=list(0.5*np.array(cl3)))
ax[0,1].plot([0.4,4.6],[0,0],'k-')
ax[0,1].set(position=[0.57,0.715,0.42,ah],xlim=[0.4,4.6],xticks=np.arange(1,5),xticklabels=['Dominant','Co-dominant','Intermediate','Suppressed'],xlabel='Crown class',ylim=ylim,ylabel='Gross growth (Mg C ha$^-$$^1$ yr$^-$$^1$)')

spc='FD'; dose=225; y=[]; ye=[]
for i in range(4):
    y.append(rsD[spc]['Bsw_M_cc' + str(i+1)][iTSF]['DA_mu'][np.where(uDose==dose)[0]][0])
    ye.append(rsD[spc]['Bsw_M_cc' + str(i+1)][iTSF]['DA_se'][np.where(uDose==dose)[0]][0])
ax[1,0].bar(np.arange(1,5,1)-0.31,y,bw,yerr=ye,facecolor=cl1,capsize=2.5,ecolor=list(0.5*np.array(cl1)))
dose=450; y=[]; ye=[]
for i in range(4):
    y.append(rsD[spc]['Bsw_M_cc' + str(i+1)][iTSF]['DA_mu'][np.where(uDose==dose)[0]][0])
    ye.append(rsD[spc]['Bsw_M_cc' + str(i+1)][iTSF]['DA_se'][np.where(uDose==dose)[0]][0])
ax[1,0].bar(np.arange(1,5,1),y,bw,yerr=ye,facecolor=cl2,capsize=2.5,ecolor=list(0.5*np.array(cl2)))
dose=675; y=[]; ye=[]
for i in range(4):
    y.append(rsD[spc]['Bsw_M_cc' + str(i+1)][iTSF]['DA_mu'][np.where(uDose==dose)[0]][0])
    ye.append(rsD[spc]['Bsw_M_cc' + str(i+1)][iTSF]['DA_se'][np.where(uDose==dose)[0]][0])
ax[1,0].bar(np.arange(1,5,1)+0.31,y,bw,yerr=ye,facecolor=cl3,capsize=2.5,ecolor=list(0.5*np.array(cl3)))
ax[1,0].plot([0.4,4.6],[0,0],'k-')
ax[1,0].set(position=[0.0725,0.39,0.42,ah],xlim=[0.4,4.6],xticks=np.arange(1,5),xticklabels=['Dominant','Co-dominant','Intermediate','Suppressed'],xlabel='Crown class',ylim=ylim,ylabel='Mortality (Mg C ha$^-$$^1$ yr$^-$$^1$)')

spc='HW'; dose=225; y=[]; ye=[]
for i in range(4):
    y.append(rsD[spc]['Bsw_M_cc' + str(i+1)][iTSF]['DA_mu'][np.where(uDose==dose)[0]][0])
    ye.append(rsD[spc]['Bsw_M_cc' + str(i+1)][iTSF]['DA_se'][np.where(uDose==dose)[0]][0])
ax[1,1].bar(np.arange(1,5,1)-0.31,y,bw,yerr=ye,facecolor=cl1,capsize=2.5,ecolor=list(0.5*np.array(cl1)))
dose=450; y=[]; ye=[]
for i in range(4):
    y.append(rsD[spc]['Bsw_M_cc' + str(i+1)][iTSF]['DA_mu'][np.where(uDose==dose)[0]][0])
    ye.append(rsD[spc]['Bsw_M_cc' + str(i+1)][iTSF]['DA_se'][np.where(uDose==dose)[0]][0])
ax[1,1].bar(np.arange(1,5,1),y,bw,yerr=ye,facecolor=cl2,capsize=2.5,ecolor=list(0.5*np.array(cl2)))
dose=675; y=[]; ye=[]
for i in range(4):
    y.append(rsD[spc]['Bsw_M_cc' + str(i+1)][iTSF]['DA_mu'][np.where(uDose==dose)[0]][0])
    ye.append(rsD[spc]['Bsw_M_cc' + str(i+1)][iTSF]['DA_se'][np.where(uDose==dose)[0]][0])
ax[1,1].bar(np.arange(1,5,1)+0.31,y,bw,yerr=ye,facecolor=cl3,capsize=2.5,ecolor=list(0.5*np.array(cl3)))
ax[1,1].plot([0.4,4.6],[0,0],'k-')
ax[1,1].set(position=[0.57,0.39,0.42,ah],xlim=[0.4,4.6],xticks=np.arange(1,5),xticklabels=['Dominant','Co-dominant','Intermediate','Suppressed'],xlabel='Crown class',ylim=ylim,ylabel='Mortality (Mg C ha$^-$$^1$ yr$^-$$^1$)')

spc='FD'; dose=225; y=[]; ye=[]
for i in range(4):
    y.append(rsD[spc]['Bsw_Net_cc' + str(i+1)][iTSF]['DA_mu'][np.where(uDose==dose)[0]][0])
    ye.append(rsD[spc]['Bsw_Net_cc' + str(i+1)][iTSF]['DA_se'][np.where(uDose==dose)[0]][0])
ax[2,0].bar(np.arange(1,5,1)-0.31,y,bw,yerr=ye,facecolor=cl1,capsize=2.5,ecolor=list(0.5*np.array(cl1)))
dose=450; y=[]; ye=[]
for i in range(4):
    y.append(rsD[spc]['Bsw_Net_cc' + str(i+1)][iTSF]['DA_mu'][np.where(uDose==dose)[0]][0])
    ye.append(rsD[spc]['Bsw_Net_cc' + str(i+1)][iTSF]['DA_se'][np.where(uDose==dose)[0]][0])
ax[2,0].bar(np.arange(1,5,1),y,bw,yerr=ye,facecolor=cl2,capsize=2.5,ecolor=list(0.5*np.array(cl2)))
dose=675; y=[]; ye=[]
for i in range(4):
    y.append(rsD[spc]['Bsw_Net_cc' + str(i+1)][iTSF]['DA_mu'][np.where(uDose==dose)[0]][0])
    ye.append(rsD[spc]['Bsw_Net_cc' + str(i+1)][iTSF]['DA_se'][np.where(uDose==dose)[0]][0])
ax[2,0].bar(np.arange(1,5,1)+0.31,y,bw,yerr=ye,facecolor=cl3,capsize=2.5,ecolor=list(0.5*np.array(cl3)))
ax[2,0].plot([0.4,4.6],[0,0],'k-')
ax[2,0].set(position=[0.0725,0.06,0.42,ah],xlim=[0.4,4.6],xticks=np.arange(1,5),xticklabels=['Dominant','Co-dominant','Intermediate','Suppressed'],xlabel='Crown class',ylim=ylim,ylabel='Net growth (Mg C ha$^-$$^1$ yr$^-$$^1$)')

spc='HW'; dose=225; y=[]; ye=[]
for i in range(4):
    y.append(rsD[spc]['Bsw_Net_cc' + str(i+1)][iTSF]['DA_mu'][np.where(uDose==dose)[0]][0])
    ye.append(rsD[spc]['Bsw_Net_cc' + str(i+1)][iTSF]['DA_se'][np.where(uDose==dose)[0]][0])
ax[2,1].bar(np.arange(1,5,1)-0.31,y,bw,yerr=ye,facecolor=cl1,capsize=2.5,ecolor=list(0.5*np.array(cl1)))
dose=450; y=[]; ye=[]
for i in range(4):
    y.append(rsD[spc]['Bsw_Net_cc' + str(i+1)][iTSF]['DA_mu'][np.where(uDose==dose)[0]][0])
    ye.append(rsD[spc]['Bsw_Net_cc' + str(i+1)][iTSF]['DA_se'][np.where(uDose==dose)[0]][0])
ax[2,1].bar(np.arange(1,5,1),y,bw,yerr=ye,facecolor=cl2,capsize=2.5,ecolor=list(0.5*np.array(cl2)))
dose=675; y=[]; ye=[]
for i in range(4):
    y.append(rsD[spc]['Bsw_Net_cc' + str(i+1)][iTSF]['DA_mu'][np.where(uDose==dose)[0]][0])
    ye.append(rsD[spc]['Bsw_Net_cc' + str(i+1)][iTSF]['DA_se'][np.where(uDose==dose)[0]][0])
ax[2,1].bar(np.arange(1,5,1)+0.31,y,bw,yerr=ye,facecolor=cl3,capsize=2.5,ecolor=list(0.5*np.array(cl3)))
ax[2,1].plot([0.4,4.6],[0,0],'k-')
ax[2,1].set(position=[0.57,0.06,0.42,ah],xlim=[0.4,4.6],xticks=np.arange(1,5),xticklabels=['Dominant','Co-dominant','Intermediate','Suppressed'],xlabel='Crown class',ylim=ylim,ylabel='Net growth (Mg C ha$^-$$^1$ yr$^-$$^1$)')

gu.axletters(ax,plt,0.03,0.92,['Douglas-fir','Western hemlock','Douglas-fir','Western hemlock',
                               'Douglas-fir','Western hemlock'],0.06)
    
gu.PrintFig(PathFigures + '\\DA_by_CrownClass_TSF_1to9','png',500)

#%% Nitrogen response stratified by N dose

uDose=np.unique(sobs['N_Dose'])
uDose=uDose[1:]

#uTSF=np.unique(np.column_stack((sobs['TSF_t0'],sobs['TSF_t1'])).astype(float),axis=0)
#uTSF=uTSF[1:,:]

uTSF=np.unique(np.column_stack((sobs['TSF_t0_adj'],sobs['TSF_t1'])).astype(float),axis=0)
ind=np.where( (uTSF[:,0]>0) & (uTSF[:,0]<=20) )[0]
uTSF=uTSF[ind,:]

spc=['FD','HW']
vr=['Bsw_G','Bsw_M','Bsw_Net']
rsD={}
for iSpc in range(len(spc)):
    rsD[spc[iSpc]]={}
    for iVar in range(len(vr)):
        rsD[spc[iSpc]][vr[iVar]]={}
        for iTSF in range(uTSF.shape[0]):
            rsD[spc[iSpc]][vr[iVar]][iTSF]={}            
            print(iSpc,iVar,iTSF)            
            yCmu0=np.nan*np.ones((uDose.size,uSite.size))
            yFmu0=np.nan*np.ones((uDose.size,uSite.size))
            for iSite in range(uSite.size):
                for iDose in range(uDose.shape[0]):
                    #iC=np.where( (sobs['ID_Site']==uSite[iSite]) & (sobs['Spc1_ID_t0']==spc[iSpc]) & (sobs['TSF_t0']>=0) & (sobs['TSF_t1']<=uTSF[iTSF,1]) & (sobs['N_Dose']==0) )[0]
                    #iF=np.where( (sobs['ID_Site']==uSite[iSite]) & (sobs['Spc1_ID_t0']==spc[iSpc]) & (sobs['TSF_t0']>=0) & (sobs['TSF_t1']<=uTSF[iTSF,1]) & (sobs['N_Dose']==uDose[iDose]) )[0]
                    iC=np.where( (sobs['ID_Site']==uSite[iSite]) & (sobs['Spc1_ID_t0']==spc[iSpc]) & (sobs['TSF_t0_adj']>=1) & (sobs['TSF_t1']<=uTSF[iTSF,1]) & (sobs['N_Dose']==0) )[0]
                    iF=np.where( (sobs['ID_Site']==uSite[iSite]) & (sobs['Spc1_ID_t0']==spc[iSpc]) & (sobs['TSF_t0_adj']>=1) & (sobs['TSF_t1']<=uTSF[iTSF,1]) & (sobs['N_Dose']==uDose[iDose]) )[0]
                    if (iC.size==0) | (iF.size==0): continue                        
                    yCmu0[iDose,iSite]=np.mean(sobs[vr[iVar]][iC])
                    yFmu0[iDose,iSite]=np.mean(sobs[vr[iVar]][iF])
            # Exclude plots if they have no site-paired control/treatment
            ind=np.where((np.isnan(yCmu0+yFmu0)==True))
            yCmu0[ind]=np.nan
            yFmu0[ind]=np.nan
            # Calculate differences
            DA=yFmu0-yCmu0
            DR=DA/yCmu0*100
            # Summarize
            SampleSize=np.zeros(yCmu0.shape[0])
            for i in range(SampleSize.size): 
                SampleSize[i]=np.where(np.isnan(yCmu0[i,:]+yFmu0[i,:])==False)[0].size
            rsD[spc[iSpc]][vr[iVar]][iTSF]['Cmu']=np.nanmean(yCmu0,axis=1)
            rsD[spc[iSpc]][vr[iVar]][iTSF]['Fmu']=np.nanmean(yFmu0,axis=1)            
            rsD[spc[iSpc]][vr[iVar]][iTSF]['DA_mu']=np.nanmean(DA,axis=1)
            rsD[spc[iSpc]][vr[iVar]][iTSF]['DA_se']=np.nanstd(DA,axis=1)/np.sqrt(SampleSize)
            rsD[spc[iSpc]][vr[iVar]][iTSF]['DR_mu']=np.nanmean(DR,axis=1)
            rsD[spc[iSpc]][vr[iVar]][iTSF]['DR_se']=np.nanstd(DR,axis=1)/np.sqrt(SampleSize)
        
#------------------------------------------------------------------------------
# Plot
#------------------------------------------------------------------------------

plt.close('all'); ms=2.5; Alpha=0.18; cl1=[0.4,1,0]; cl2=[0.27,0.49,0.77]; cl3=[0.5,0,0]
fig,ax=plt.subplots(2,3,figsize=gu.cm2inch(15.5,11))
spc='FD'
iTSF=np.where(uTSF[:,0]==1)[0][0]
ax[0,0].plot([0,1000],[0,0],'k-')
ax[0,0].errorbar(uDose,rsD[spc]['Bsw_G'][iTSF]['DA_mu'],yerr=rsD[spc]['Bsw_G'][iTSF]['DA_se'],capsize=3,linestyle='',ecolor=cl1)
ax[0,0].plot(uDose,rsD[spc]['Bsw_G'][iTSF]['DA_mu'],'-o',markersize=ms,color=cl1,label='Gross growth')
ax[0,0].errorbar(uDose,rsD[spc]['Bsw_Net'][iTSF]['DA_mu'],yerr=rsD[spc]['Bsw_Net'][0]['DA_se'],capsize=3,linestyle='',ecolor=cl2)
ax[0,0].plot(uDose,rsD[spc]['Bsw_Net'][iTSF]['DA_mu'],'-s',markersize=ms,color=cl2,label='Net growth')
#ax[0,0].errorbar(uDose,-rsD[spc]['Bsw_M'][iTSF]['DA_mu'],yerr=rsD[spc]['Bsw_M'][0]['DA_se'],capsize=3,linestyle='',ecolor=cl3)
#ax[0,0].plot(uDose,-rsD[spc]['Bsw_M'][iTSF]['DA_mu'],'-v',markersize=ms,color=cl3,label='-Mortality')
ax[0,0].set(position=[0.08,0.555,0.27,0.43],xlim=[0,1000],xticks=uDose,ylim=[-0.5,2],ylabel='Response (Mg C ha$^-$$^1$ yr$^-$$^1$)')
ax[0,0].legend(loc='lower left')
iTSF=np.where(uTSF[:,0]==7)[0][0]
ax[0,1].plot([0,1000],[0,0],'k-')
ax[0,1].errorbar(uDose,rsD[spc]['Bsw_G'][iTSF]['DA_mu'],yerr=rsD[spc]['Bsw_G'][iTSF]['DA_se'],capsize=3,linestyle='',ecolor=cl1)
ax[0,1].plot(uDose,rsD[spc]['Bsw_G'][iTSF]['DA_mu'],'-o',markersize=ms,color=cl1)
ax[0,1].errorbar(uDose,rsD[spc]['Bsw_Net'][iTSF]['DA_mu'],yerr=rsD[spc]['Bsw_Net'][iTSF]['DA_se'],capsize=3,linestyle='',ecolor=cl2)
ax[0,1].plot(uDose,rsD[spc]['Bsw_Net'][iTSF]['DA_mu'],'-s',markersize=ms,color=cl2)
ax[0,1].set(position=[0.4,0.555,0.27,0.43],xlim=[0,1000],xticks=uDose,ylim=[-0.5,2])
iTSF=np.where(uTSF[:,0]==10)[0][0]
ax[0,2].plot([0,1000],[0,0],'k-')
ax[0,2].errorbar(uDose,rsD[spc]['Bsw_G'][iTSF]['DA_mu'],yerr=rsD[spc]['Bsw_G'][iTSF]['DA_se'],capsize=3,linestyle='',ecolor=cl1)
ax[0,2].plot(uDose,rsD[spc]['Bsw_G'][iTSF]['DA_mu'],'-o',markersize=ms,color=cl1)
ax[0,2].errorbar(uDose,rsD[spc]['Bsw_Net'][iTSF]['DA_mu'],yerr=rsD[spc]['Bsw_Net'][iTSF]['DA_se'],capsize=3,linestyle='',ecolor=cl2)
ax[0,2].plot(uDose,rsD[spc]['Bsw_Net'][iTSF]['DA_mu'],'-s',markersize=ms,color=cl2)
ax[0,2].set(position=[0.72,0.555,0.27,0.43],xlim=[0,1000],xticks=uDose,ylim=[-0.5,2])
spc='HW'
iTSF=np.where(uTSF[:,0]==1)[0][0]
ax[1,0].plot([0,1000],[0,0],'k-')
ax[1,0].errorbar(uDose,rsD[spc]['Bsw_G'][iTSF]['DA_mu'],yerr=rsD[spc]['Bsw_G'][iTSF]['DA_se'],capsize=3,linestyle='',ecolor=cl1)
ax[1,0].plot(uDose,rsD[spc]['Bsw_G'][iTSF]['DA_mu'],'-o',markersize=ms,color=cl1,label='Gross growth')
ax[1,0].errorbar(uDose,rsD[spc]['Bsw_Net'][iTSF]['DA_mu'],yerr=rsD[spc]['Bsw_Net'][iTSF]['DA_se'],capsize=3,linestyle='',ecolor=cl2)
ax[1,0].plot(uDose,rsD[spc]['Bsw_Net'][iTSF]['DA_mu'],'-s',markersize=ms,color=cl2,label='Net growth')
ax[1,0].set(position=[0.08,0.07,0.27,0.43],xlim=[0,1000],xticks=uDose,xlabel='N application dose (kg N ha$^-$$^1$)',ylim=[-0.5,2],ylabel='Response (Mg C ha$^-$$^1$ yr$^-$$^1$)')
iTSF=np.where(uTSF[:,0]==7)[0][0]
ax[1,1].plot([0,1000],[0,0],'k-')
ax[1,1].errorbar(uDose,rsD[spc]['Bsw_G'][iTSF]['DA_mu'],yerr=rsD[spc]['Bsw_G'][iTSF]['DA_se'],capsize=3,linestyle='',ecolor=cl1)
ax[1,1].plot(uDose,rsD[spc]['Bsw_G'][iTSF]['DA_mu'],'-o',markersize=ms,color=cl1)
ax[1,1].errorbar(uDose,rsD[spc]['Bsw_Net'][iTSF]['DA_mu'],yerr=rsD[spc]['Bsw_Net'][iTSF]['DA_se'],capsize=3,linestyle='',ecolor=cl2)
ax[1,1].plot(uDose,rsD[spc]['Bsw_Net'][iTSF]['DA_mu'],'-s',markersize=ms,color=cl2)
ax[1,1].set(position=[0.4,0.07,0.27,0.43],xlim=[0,1000],xticks=uDose,xlabel='N application dose (kg N ha$^-$$^1$)',ylim=[-0.5,2])
iTSF=np.where(uTSF[:,0]==10)[0][0]
ax[1,2].plot([0,1000],[0,0],'k-')
ax[1,2].errorbar(uDose,rsD[spc]['Bsw_G'][iTSF]['DA_mu'],yerr=rsD[spc]['Bsw_G'][iTSF]['DA_se'],capsize=3,linestyle='',ecolor=cl1)
ax[1,2].plot(uDose,rsD[spc]['Bsw_G'][iTSF]['DA_mu'],'-o',markersize=ms,color=cl1)
ax[1,2].errorbar(uDose,rsD[spc]['Bsw_Net'][iTSF]['DA_mu'],yerr=rsD[spc]['Bsw_Net'][iTSF]['DA_se'],capsize=3,linestyle='',ecolor=cl2)
ax[1,2].plot(uDose,rsD[spc]['Bsw_Net'][iTSF]['DA_mu'],'-s',markersize=ms,color=cl2)
ax[1,2].set(position=[0.72,0.07,0.27,0.43],xlim=[0,1000],xticks=uDose,xlabel='N application dose (kg N ha$^-$$^1$)',ylim=[-0.5,2])
gu.axletters(ax,plt,0.03,0.92,['Douglas-fir (1 to 3 years)','Douglas-fir (1 to 9 years)','Douglas-fir (1 to 12 years)','Western hemlock (1 to 3 years)','Western hemlock (1 to 9 years)','Western hemlock (1 to 12 years)'],0.09)
#gu.PrintFig(PathFigures + '\\DA_vs_Dose.png','emf',500)


#%% Self-thinning relationship

cols=['ID_Site','Spc1_ID_t0','LogN_C','LogN_F','LogS_C','LogS_F']

dose=225
tsf0_max=10
df0=pd.DataFrame(data=np.nan*np.ones((360,len(cols))),columns=cols)
cnt=0
for i in range(uSite.size):
    ind=np.where((sobs['ID_Site']==uSite[i]))[0]
    if ind.size==0: continue

    df0.loc[cnt,'ID_Site']=uSite[i]
    df0.loc[cnt,'Spc1_ID_t0']=sobs['Spc1_ID_t0'][ind[0]]
    
    iC=np.where( (sobs['ID_Site']==uSite[i]) & (sobs['TSF_t0']<=tsf0_max) & (sobs['N_Dose']==0) )[0]
    iF=np.where( (sobs['ID_Site']==uSite[i]) & (sobs['TSF_t0']<=tsf0_max) & (sobs['N_Dose']==dose) )[0]
    df0.loc[cnt,'LogN_C']=np.mean(np.log(sobs['N_t0'][iC]))
    df0.loc[cnt,'LogN_F']=np.mean(np.log(sobs['N_t0'][iF]))
    df0.loc[cnt,'LogS_C']=np.mean(np.log(sobs['Bsw_t0'][iC]*1000/sobs['N_t0'][iC]))
    df0.loc[cnt,'LogS_F']=np.mean(np.log(sobs['Bsw_t0'][iF]*1000/sobs['N_t0'][iF]))
    cnt=cnt+1
ind=np.where(np.isnan(df0['LogN_C']+df0['LogN_F'])==False)[0]
df0=df0.loc[ind]
df0=df0.reset_index(drop=True)

dose=450
tsf0_max=10
df1=pd.DataFrame(data=np.nan*np.ones((360,len(cols))),columns=cols)
cnt=0
for i in range(uSite.size):
    ind=np.where((sobs['ID_Site']==uSite[i]))[0]
    if ind.size==0: continue

    df1.loc[cnt,'ID_Site']=uSite[i]
    df1.loc[cnt,'Spc1_ID_t0']=sobs['Spc1_ID_t0'][ind[0]]
    
    iC=np.where( (sobs['ID_Site']==uSite[i]) & (sobs['TSF_t0']<=tsf0_max) & (sobs['N_Dose']==0) )[0]
    iF=np.where( (sobs['ID_Site']==uSite[i]) & (sobs['TSF_t0']<=tsf0_max) & (sobs['N_Dose']==dose) )[0]
    df1.loc[cnt,'LogN_C']=np.mean(np.log(sobs['N_t0'][iC]))
    df1.loc[cnt,'LogN_F']=np.mean(np.log(sobs['N_t0'][iF]))
    df1.loc[cnt,'LogS_C']=np.mean(np.log(sobs['Bsw_t0'][iC]*1000/sobs['N_t0'][iC]))
    df1.loc[cnt,'LogS_F']=np.mean(np.log(sobs['Bsw_t0'][iF]*1000/sobs['N_t0'][iF]))
    cnt=cnt+1
ind=np.where(np.isnan(df1['LogN_C']+df1['LogN_F'])==False)[0]
df1=df1.loc[ind]
df1=df1.reset_index(drop=True)

# Plot
plt.close('all'); 
fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(15.5,7))

# Douglas-fir
ikp=np.where( (df0['Spc1_ID_t0']=='FD') & (np.abs(df0['LogN_C'])<20) & (np.abs(df0['LogN_F'])<20) & (df0['LogS_C']>1.8) )[0]

xhat=np.arange(2,6.5,0.1)
mdC=sm.RLM(df0.loc[ikp,'LogN_C'].values,np.c_[np.ones(ikp.size),df0.loc[ikp,'LogS_C'].values]).fit()
mdC.summary()
#yhatC=mdC.fittedvalues
yhatC=mdC.predict(np.c_[np.ones(xhat.size),xhat])
#yhatC=mdC.get_prediction(np.c_[np.ones(xhat.size),xhat]).summary_frame(alpha=0.05)
#yhatC=yhatC['mean'].values

mdF=sm.RLM(df0.loc[ikp,'LogN_F'].values,np.c_[np.ones(ikp.size),df0.loc[ikp,'LogS_F'].values]).fit()
mdF.summary()
#yhatF=mdF.fittedvalues
yhatF=mdF.predict(np.c_[np.ones(xhat.size),xhat])

ikp2=np.where( (df1['Spc1_ID_t0']=='FD') & (np.abs(df1['LogN_C'])<20) & (np.abs(df1['LogN_F'])<20) & (df1['LogS_C']>1.8) )[0]
mdF2=sm.RLM(df1.loc[ikp2,'LogN_F'].values,np.c_[np.ones(ikp2.size),df1.loc[ikp2,'LogS_F'].values]).fit()
mdF2.summary()
yhatF2=mdF2.predict(np.c_[np.ones(xhat.size),xhat])

lx=np.arange(2,8,0.01)
lyC=np.zeros(lx.size)
lyF1=np.zeros(lx.size)
lyF2=np.zeros(lx.size)
for i in range(lx.size):
    lyC[i]=np.exp(mdC.predict(np.c_[1,lx[i]]))
    lyF1[i]=np.exp(mdF.predict(np.c_[1,lx[i]]))
    lyF2[i]=np.exp(mdF2.predict(np.c_[1,lx[i]]))
indC=np.where(np.abs(lyC-1000)==np.min(np.abs(lyC-1000)))[0]
indF1=np.where(np.abs(lyF1-1000)==np.min(np.abs(lyF1-1000)))[0]
indF2=np.where(np.abs(lyF2-1000)==np.min(np.abs(lyF2-1000)))[0]
print(np.exp(lx[indC]))
print(np.exp(lx[indF1]))
print(np.exp(lx[indF2]))


#ax[0].plot(df0.loc[ikp,'LogS_C'],df0.loc[ikp,'LogN_C'],'o',color=[0,0,1],markersize=4)
#ax[0].plot(df0.loc[ikp,'LogS_F'],df0.loc[ikp,'LogN_F'],'s',color=[1,0,0],markersize=4)
ax[0].plot([1.2,6.8],[6.908,6.908],'-',color=[0.8,0.8,0.8],linewidth=2.5,label='1000 stems ha$^-$$^1$')
ax[0].plot(xhat,yhatC,'-',color=[0.24,0.49,0.77],linewidth=1,label='Control')
ax[0].plot(xhat,yhatF,'--',color=[0.6,1,0],linewidth=1,label='N225')
ax[0].plot(xhat,yhatF2,'-.',color=[0.5,0,1],linewidth=1,label='N450')
ax[0].set(position=[0.07,0.11,0.42,0.82],ylim=[6,9.5],xlim=[1.2,6.8],xlabel='Log average stemwood biomass',ylabel='Log stand density')
ax[0].legend(loc='upper right')

# Western hemlock
ikp=np.where( (df0['Spc1_ID_t0']=='HW') & (np.abs(df0['LogN_C'])<20) & (np.abs(df0['LogN_F'])<20) & (df0['LogS_C']>1.8) )[0]

xhat=np.arange(1.5,6,0.1)
mdC=sm.OLS(df0.loc[ikp,'LogN_C'].values,np.c_[np.ones(ikp.size),df0.loc[ikp,'LogS_C'].values]).fit()
mdC.summary()
yhatC=mdC.predict(np.c_[np.ones(xhat.size),xhat])

mdF=sm.OLS(df0.loc[ikp,'LogN_F'].values,np.c_[np.ones(ikp.size),df0.loc[ikp,'LogS_F'].values]).fit()
mdF.summary()
yhatF=mdF.predict(np.c_[np.ones(xhat.size),xhat])

ikp2=np.where( (df1['Spc1_ID_t0']=='HW') & (np.abs(df1['LogN_C'])<20) & (np.abs(df1['LogN_F'])<20) & (df1['LogS_C']>1.8) )[0]
mdF2=sm.RLM(df1.loc[ikp2,'LogN_F'].values,np.c_[np.ones(ikp2.size),df1.loc[ikp2,'LogS_F'].values]).fit()
mdF2.summary()
yhatF2=mdF2.predict(np.c_[np.ones(xhat.size),xhat])

lx=np.arange(2,8,0.01)
lyC=np.zeros(lx.size)
lyF1=np.zeros(lx.size)
lyF2=np.zeros(lx.size)
for i in range(lx.size):
    lyC[i]=np.exp(mdC.predict(np.c_[1,lx[i]]))
    lyF1[i]=np.exp(mdF.predict(np.c_[1,lx[i]]))
    lyF2[i]=np.exp(mdF2.predict(np.c_[1,lx[i]]))
indC=np.where(np.abs(lyC-1000)==np.min(np.abs(lyC-1000)))[0]
indF1=np.where(np.abs(lyF1-1000)==np.min(np.abs(lyF1-1000)))[0]
indF2=np.where(np.abs(lyF2-1000)==np.min(np.abs(lyF2-1000)))[0]
print(np.exp(lx[indC]))
print(np.exp(lx[indF1]))
print(np.exp(lx[indF2]))

#ax[1].plot(df0.loc[ikp,'LogS_C'],df0.loc[ikp,'LogN_C'],'o',color=[0,0,1],markersize=4)
#ax[1].plot(df0.loc[ikp,'LogS_F'],df0.loc[ikp,'LogN_F'],'s',color=[1,0,0],markersize=4)
ax[1].plot([1.2,6.8],[6.908,6.908],'-',color=[0.8,0.8,0.8],linewidth=2.5)
ax[1].plot(xhat,yhatC,'-',color=[0.24,0.49,0.77],linewidth=1)
ax[1].plot(xhat,yhatF,'--',color=[0.6,1,0],linewidth=1)
ax[1].plot(xhat,yhatF2,'-.',color=[0.5,0,1],linewidth=1)
ax[1].set(position=[0.57,0.11,0.42,0.82],ylim=[6,9.5],xlim=[1.2,6.8],xlabel='Log average stemwood biomass',ylabel='Log stand density')
gu.axletters(ax,plt,0.03,0.92,['Coastal Douglas-fir','Western hemlock'],0.04)

gu.PrintFig(PathFigures + '\\SelfThinningCurves','emf',500)

# Plot density-size relationships
plt.close('all'); 
fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(15.5,7))
ikp=np.where( (sobs['Spc1_ID_t0']=='FD') & (sobs['N_Dose']==0) )[0]
ax[0].plot(sobs['Bsw_t0'][ikp]/sobs['N_t0'][ikp],sobs['N_t0'][ikp],'.')
ikp=np.where( (sobs['Spc1_ID_t0']=='FD') & (sobs['N_Dose']==225) )[0]
ax[0].plot(sobs['Bsw_t0'][ikp]/sobs['N_t0'][ikp],sobs['N_t0'][ikp],'.')
ax[0].set(yscale='log',xscale='log')
ikp=np.where( (sobs['Spc1_ID_t0']=='HW') & (sobs['N_Dose']==0) )[0]
ax[1].plot(sobs['Bsw_t0'][ikp]/sobs['N_t0'][ikp],sobs['N_t0'][ikp],'.')
ikp=np.where( (sobs['Spc1_ID_t0']=='HW') & (sobs['N_Dose']==225) )[0]
ax[1].plot(sobs['Bsw_t0'][ikp]/sobs['N_t0'][ikp],sobs['N_t0'][ikp],'.')
ax[1].set(yscale='log',xscale='log')



#------------------------------------------------------------------------------
# Mortality from self-thinning relationship
#------------------------------------------------------------------------------

Amax=100
SN0=2200
SB0=1 # Initial biomass (MgC/ha)
SG0=1.5 # Net growth (MgC/ha/yr)
Ba1000=350 # kgDM/tree
Neta=np.array([-1.4,-1.5,-1.6])
rrG=np.arange(1,2.05,0.05)
AveDyingB=0.2
#rrG=np.array([1,1.22])
dSN_rel=np.zeros((rrG.size,Neta.size))
M=np.zeros((rrG.size,Neta.size))
for iNeta in range(Neta.size):
    for k in range(rrG.size):        
        SG=rrG[k]*SG0
        SN=SN0*np.ones(m)
        SB=SB0*np.ones(m)
        Ba=np.zeros(Amax)
        for iA in range(1,Amax):
            SB[iA]=SB[iA-1]+SG; 
            SN[iA]=SN[iA-1]
            Ba[iA]=SB[iA]*1000*2/SN[iA]; 
            Ba_max=Ba1000*(SN[iA]/1000)**Neta[iNeta] # kgC/tree
            if Ba[iA]>Ba_max:
                for j in range(1000):
                    SN[iA]=SN[iA]-1; 
                    Ba[iA]=SB[iA]*1000*2/SN[iA]; 
                    Ba_max=Ba1000*(SN[iA]/1000)**Neta[iNeta] 
                    if Ba[iA]<=Ba_max:
                        break            
        #plt.plot(Ba,SN,'-')
        dSN=-np.append(0,SN[1:]-SN[0:-1])
        dSN_rel[k,iNeta]=np.mean(-np.append(0,SN[1:]-SN[0:-1])/SN*100)
        M[k,iNeta]=AveDyingB*np.mean(dSN*Ba/1000)

cl=[0,0.8,0]
cl=[0.27,0.49,0.77]


plt.close('all'); 
fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(14.5,6.5))
ax[0].plot(rrG,M[:,0]-M[0,0],'k--',color=cl,linewidth=0.75)
ax[0].plot(rrG,M[:,1]-M[0,1],'k-',color=cl,linewidth=1.5)
ax[0].plot(rrG,M[:,2]-M[0,2],'k--',color=cl,linewidth=0.75)
ax[0].set(position=[0.075,0.13,0.42,0.82],xlim=[1,2],ylim=[0,2],xlabel='Growth enhancement (response ratio)',ylabel='Mortality (Mg C ha-1 yr-1)')
ax[0].grid(True)

#plt.close('all'); 
#fig,ax=plt.subplots(1,figsize=gu.cm2inch(8,6.5))
ax[1].plot([1,2],[1,2],'k-',linewidth=3,color=[0.8,0.8,0.8])
ax[1].text(1.64,1.71,'1:1')
ax[1].plot(rrG,M[:,0]/M[0,0],'k--',color=cl,linewidth=0.75)
ax[1].plot(rrG,M[:,1]/M[0,1],'k-',color=cl,linewidth=1.5)
ax[1].plot(rrG,M[:,2]/M[0,2],'k--',color=cl,linewidth=0.75)
ax[1].set(position=[0.565,0.13,0.42,0.82],xlim=[1,2],ylim=[1,2],xlabel='Growth enhancement (response ratio)',ylabel='Mortality (response ratio)')
ax[1].grid(True)
gu.axletters(ax,plt,0.04,0.92)
#plt.savefig(r'G:\My Drive\Figures\Fertilization\MortalityFromSelfThinningRule.png',format='png',dpi=900)



plt.close('all'); 
fig,ax=plt.subplots(1,figsize=gu.cm2inch(8,6.5))
plt.plot(rrG,dSN_rel[:,0]-dSN_rel[0,0],'k--',color=[0.7,0.7,0.7],linewidth=0.75)
plt.plot(rrG,dSN_rel[:,1]-dSN_rel[0,1],'k-',linewidth=1.5)
plt.plot(rrG,dSN_rel[:,2]-dSN_rel[0,2],'k--',color=[0.7,0.7,0.7],linewidth=0.75)
ax.set(position=[0.15,0.13,0.82,0.82],xlim=[1,2],ylim=[0,2],xlabel='Growth enhancement (response ratio)',ylabel='Mortality (Mg C ha-1 yr-1)')
ax.grid(True)


plt.plot(rrG,y/y[0],'o-')



#%% Time series analysis

tv=np.arange(1971,2012,1)
uSpc=['FD','HW']
uDose=[225,450]
uVar=['dN_rel','RGR','Age_t0','N_t0','N_Net','BA_t0','H_obs_t0','Bsw_t0','Bsw_Net','Bsw_G','Bsw_M']
rts={}
for iSpc in range(len(uSpc)):
    rts[uSpc[iSpc]]={}
    for iDose in range(len(uDose)):
        rts[uSpc[iSpc]][uDose[iDose]]={}
        rts[uSpc[iSpc]][uDose[iDose]]['y_c']={}
        rts[uSpc[iSpc]][uDose[iDose]]['y_f']={}
        for iVar in range(len(uVar)):
            rts[uSpc[iSpc]][uDose[iDose]]['y_c'][uVar[iVar]]=np.nan*np.ones((tv.size,sobs['Year_t0'].size))
            rts[uSpc[iSpc]][uDose[iDose]]['y_f'][uVar[iVar]]=np.nan*np.ones((tv.size,sobs['Year_t0'].size))
            cnt_c=0; 
            cnt_f=0;
            for iYr in range(sobs['Year_t0'].size):
                #if (sobs['ID_Site'][iYr]==29) | (sobs['ID_Site'][iYr]==71) | (sobs['ID_Site'][iYr]==77):
                #    continue
                #if sobs['Age_t0'][iYr]>80:
                #    continue
                it=np.where( (tv>=sobs['Year_t0'][iYr]) & (tv<=sobs['Year_t1'][iYr]) )[0]
                if (sobs['Spc1_ID_t0'][iYr]==uSpc[iSpc]) & (sobs['N_Dose'][iYr]==0):
                    #it=np.where( (tv>=sobs['Year_t0'][iYr]) & (tv<=sobs['Year_t1'][iYr]) )[0]                    
                    rts[uSpc[iSpc]][uDose[iDose]]['y_c'][uVar[iVar]][it,cnt_c]=sobs[uVar[iVar]][iYr]; 
                    cnt_c=cnt_c+1
                if (sobs['Spc1_ID_t0'][iYr]==uSpc[iSpc]) & (sobs['N_Dose'][iYr]==uDose[iDose]):
                    #it=np.where( (tv>=sobs['Year_t0'][iYr]) & (tv<=sobs['Year_t1'][iYr]) )[0]
                    rts[uSpc[iSpc]][uDose[iDose]]['y_f'][uVar[iVar]][it,cnt_f]=sobs[uVar[iVar]][iYr]; 
                    cnt_f=cnt_f+1



Spc='HW'
nam='RGR'
Dose=450

plt.close('all')
plt.plot(tv,np.nanmean(rts[Spc][Dose]['y_c'][nam],axis=1),'-o')
plt.plot(tv,np.nanmean(rts[Spc][Dose]['y_f'][nam],axis=1),'-s')

plt.close('all')
plt.plot(tv,np.nanmean(rts[Spc][Dose]['y_f'][nam],axis=1)-np.nanmean(rts[Spc][Dose]['y_c'][nam],axis=1),'-o')



nam='Bsw_Net'
Spc='FD'
plt.close('all')
plt.plot(tv,0*np.ones(tv.size),'k-')
plt.plot(tv,np.nanmean(rts[Spc][225]['y_f'][nam],axis=1)-np.nanmean(rts[Spc][225]['y_c'][nam],axis=1),'-o')
plt.plot(tv,np.nanmean(rts[Spc][450]['y_f'][nam],axis=1)-np.nanmean(rts[Spc][450]['y_c'][nam],axis=1),'-s')


plt.close('all')
plt.plot(tv,np.cumsum(np.nanmean(rts[Spc][225]['y_f'][nam],axis=1)-np.nanmean(rts[Spc][225]['y_c'][nam],axis=1)),'-o')
plt.plot(tv,np.cumsum(np.nanmean(rts[Spc][450]['y_f'][nam],axis=1)-np.nanmean(rts[Spc][450]['y_c'][nam],axis=1)),'-s')


# Damage agents
da=['DA_A_t1', 'DA_AB_t1', 'DA_AD_t1', 'DA_AE_t1', 'DA_AX_t1', 'DA_BT_t1', 'DA_CK_t1', 'DA_D_t1', 
 'DA_DDF_t1', 'DA_DDP_t1', 'DA_DDS_t1', 'DA_DM_t1', 'DA_DN_t1', 'DA_DR_t1', 'DA_DRA_t1', 'DA_DS_t1',
 'DA_DT_t1', 'DA_FCK_t1', 'DA_FK_t1', 'DA_MT_t1', 'DA_N_t1', 'DA_NGC_t1', 'DA_NW_t1', 'DA_NWS_t1', 
 'DA_NX_t1', 'DA_NY_t1', 'DA_PC_t1', 'DA_SC_t1', 'DA_SCB_t1', 'DA_SCM_t1', 'DA_SCT_t1']

muC=np.zeros(len(da))
muF=np.zeros(len(da))
for i in range(len(da)):
    indC=np.where(sobs['N_Dose']==0)[0]
    indF=np.where(sobs['N_Dose']==225)[0]
    muC[i]=np.nanmean(sobs[da[i]][indC])
    muF[i]=np.nanmean(sobs[da[i]][indF])

plt.close('all')
plt.barh(np.arange(0,len(da))-0.25,muC,0.5)
plt.barh(np.arange(0,len(da))+0.25,muF,0.5)



muc=np.zeros(tv.size)
muf=np.zeros(tv.size)
for i in range(len(da)):
    muc=muc+np.nanmean(ytc[da[i]],axis=1)
    muf=muf+np.nanmean(ytf1[da[i]],axis=1)
 
plt.close('all')
plt.plot(tv,muc,'-o')
plt.plot(tv,muf,'-s')

nam='DA_AX_t1'; plt.close('all')
plt.plot(tv,np.nanmean(ytc[nam],axis=1),'-o')
plt.plot(tv,np.nanmean(ytf1[nam],axis=1),'-s')



Gc=np.nan*np.ones((tv.size,5000))
Ac=np.nan*np.ones((tv.size,5000))
Gf1=np.nan*np.ones((tv.size,5000))
cwd=np.nan*np.ones((tv.size,5000))
c_c=0;c_f1=0;
c_cwd=0
for i in range(sobs['Year_t0'].size):
    it=np.where( (tv>=sobs['Year_t0'][i]) & (tv<=sobs['Year_t1'][i]) )[0]
    #cwd[it,c_cwd]=sobs['tmean_mjjas_r'][i]
    
    if (sobs['Spc1_ID_t0'][i]=='FD') & (sobs['N_Dose'][i]==0):
        it=np.where( (tv>=sobs['Year_t0'][i]) & (tv<=sobs['Year_t1'][i]) )[0]
        Gc[it,c_c]=sobs['Bsw_M'][i]
        Ac[it,c_c]=sobs['Age_t0'][i]
        c_c=c_c+1
    if (sobs['Spc1_ID_t0'][i]=='FD') & (sobs['N_Dose'][i]>0) & (sobs['N_Dose'][i]==225):
        it=np.where( (tv>=sobs['Year_t0'][i]) & (tv<=sobs['Year_t1'][i]) )[0]
        Gf1[it,c_f1]=sobs['Bsw_M'][i]    
        c_f1=c_f1+1

dGY=gu.ReadExcel(r'C:\Users\rhember\Documents\fd_yield_tipsy.xlsx')
Amu=np.nanmean(Ac,axis=1)
G_gy=np.zeros(tv.size)
for i in range(Amu.size):
    ind=np.where(dGY['Age']==int(Amu[i]))[0]
    G_gy[i]=dGY['Gnet'][ind]


plt.close('all')    
plt.plot(tv,np.nanmean(Gc,axis=1),'-ko')
plt.plot(tv,np.nanmean(Gf1,axis=1),'-go')

plt.close('all')    
plt.plot(tv,np.nanmean(Gc,axis=1)-gu.movingave(G_gy,5,'centre'),'-ko')
plt.plot(tv,np.nanmean(Gf1,axis=1)-gu.movingave(G_gy,5,'centre'),'-go')

fig,ax=plt.subplots(1)
plt.plot(tv,np.nanmean(Gf1,axis=1)-np.nanmean(Gc,axis=1),'-go')
plt.grid()


plt.plot(tv,zscore(np.nanmean(cwd,axis=1)),'-rs')


fig,ax=plt.subplots(1)
plt.plot(tv,np.nanmean(Gf1-Gc,axis=1),'-go')



fig,ax=plt.subplots(1)
plt.plot(np.cumsum(np.nan_to_num(np.nanmean(Gf1,axis=1)-np.nanmean(Gc,axis=1))),'-go')


#%% A analysis

A=np.arange(1,201,1)
uSpc=['FD','HW']
uDose=[225,450]
uVar=['N_t0','N_Net','BA_t0','H_obs_t0','Bsw_t0','Bsw_Net','Bsw_G','Bsw_M']
rts={}
for iSpc in range(len(uSpc)):
    rts[uSpc[iSpc]]={}
    for iDose in range(len(uDose)):
        rts[uSpc[iSpc]][uDose[iDose]]={}
        rts[uSpc[iSpc]][uDose[iDose]]['y_c']={}
        rts[uSpc[iSpc]][uDose[iDose]]['y_f']={}
        for iVar in range(len(uVar)):
            rts[uSpc[iSpc]][uDose[iDose]]['y_c'][uVar[iVar]]=np.nan*np.ones((A.size,sobs['Year_t0'].size))
            rts[uSpc[iSpc]][uDose[iDose]]['y_f'][uVar[iVar]]=np.nan*np.ones((A.size,sobs['Year_t0'].size))
            cnt_c=0; 
            cnt_f=0;
            for iYr in range(sobs['Year_t0'].size):
                it=np.where( (A>=sobs['Age_t0'][iYr]) & (A<=sobs['Age_t1'][iYr]) )[0]
                if (sobs['Spc1_ID_t0'][iYr]==uSpc[iSpc]) & (sobs['N_Dose'][iYr]==0):
                    #it=np.where( (tv>=sobs['Year_t0'][iYr]) & (tv<=sobs['Year_t1'][iYr]) )[0]                    
                    rts[uSpc[iSpc]][uDose[iDose]]['y_c'][uVar[iVar]][it,cnt_c]=sobs[uVar[iVar]][iYr]; 
                    cnt_c=cnt_c+1
                if (sobs['Spc1_ID_t0'][iYr]==uSpc[iSpc]) & (sobs['N_Dose'][iYr]==uDose[iDose]):
                    #it=np.where( (tv>=sobs['Year_t0'][iYr]) & (tv<=sobs['Year_t1'][iYr]) )[0]
                    rts[uSpc[iSpc]][uDose[iDose]]['y_f'][uVar[iVar]][it,cnt_f]=sobs[uVar[iVar]][iYr]; 
                    cnt_f=cnt_f+1


Spc='HW'
nam='Bsw_G'
Dose=450

plt.close('all')
plt.plot(A,np.nanmean(rts[Spc][Dose]['y_c'][nam],axis=1),'-o')
plt.plot(A,np.nanmean(rts[Spc][Dose]['y_f'][nam],axis=1),'-s')


#%% Plot age responses of stand density and stemwood biomass (for SI)

uV=['N_t0','Dam_t0','H_obs_t0','H_gf_t0','Bsw_t0','H_obs_G']
uS=['FD','HW']
uA=np.unique(sobs['Age_t0'])
uDose=[0,225,450,675]
fA_mu=[None]*len(uS)
fA_sig=[None]*len(uS)
for iS in range(len(uS)):
    fA_mu0=[None]*len(uDose)
    fA_sig0=[None]*len(uDose)
    for iDose in range(len(uDose)):
        fA_mu0[iDose]={}
        fA_sig0[iDose]={}
        for iV in range(len(uV)):
            fA_mu0[iDose][uV[iV]]=np.nan*np.ones(uA.size)
            fA_sig0[iDose][uV[iV]]=np.nan*np.ones(uA.size)
            for j in range(uA.size):
                ind=np.where( (sobs['Age_t0']==uA[j]) & (sobs['Spc1_ID_t0']==uS[iS]) & (sobs['N_Dose']==uDose[iDose]) )[0]
                fA_mu0[iDose][uV[iV]][j]=np.nanmean(sobs[uV[iV]][ind])
                fA_sig0[iDose][uV[iV]][j]=np.nanstd(sobs[uV[iV]][ind])
    fA_mu[iS]=fA_mu0
    fA_sig[iS]=fA_sig0


# Plot
plt.close('all'); ms=2; Alpha=0.14; xlim=[0,125]
fig,ax=plt.subplots(4,2,figsize=gu.cm2inch(15.5,14))
ax[0,0].plot([-5,200],[0,0],'k-')
ax[0,0].fill_between(uA,fA_mu[0][0]['N_t0']-fA_sig[0][0]['N_t0'],
  fA_mu[0][0]['N_t0']+fA_sig[0][0]['N_t0'],color=[0.27,0.49,0.79],alpha=Alpha,linewidth=0)
ax[0,0].plot(uA,fA_mu[0][0]['N_t0'],'-o',color=[0.27,0.49,0.79],markersize=ms) 
ax[0,0].set(position=[0.07,0.78,0.42,0.2],xlim=xlim,xlabel='',ylim=[0,6000],ylabel='Stand density (stems/ha)')

ax[1,0].plot([-5,200],[0,0],'k-')
ax[1,0].fill_between(uA,fA_mu[0][0]['Dam_t0']-fA_sig[0][0]['Dam_t0'],
  fA_mu[0][0]['Dam_t0']+fA_sig[0][0]['Dam_t0'],color=[0.27,0.49,0.79],alpha=Alpha,linewidth=0)
ax[1,0].plot(uA,fA_mu[0][0]['Dam_t0'],'-o',color=[0.27,0.49,0.79],markersize=ms) 
ax[1,0].set(position=[0.07,0.54,0.42,0.2],xlim=xlim,xlabel='',ylim=[0,60],ylabel='Diameter mean (cm)')

ax[2,0].plot([-5,200],[0,0],'k-')
ax[2,0].fill_between(uA,fA_mu[0][0]['H_obs_t0']-fA_sig[0][0]['H_obs_t0'],
  fA_mu[0][0]['H_obs_t0']+fA_sig[0][0]['H_obs_t0'],color=[0.27,0.49,0.79],alpha=Alpha,linewidth=0)
ax[2,0].plot(uA,fA_mu[0][0]['H_obs_t0'],'-o',color=[0.27,0.49,0.79],markersize=ms) 
ax[2,0].set(position=[0.07,0.3,0.42,0.2],xlim=xlim,xlabel='',ylim=[0,50],ylabel='Height mean (m)')

ax[3,0].plot([-5,200],[0,0],'k-')
ax[3,0].fill_between(uA,fA_mu[0][0]['Bsw_t0']-fA_sig[0][0]['Bsw_t0'],
  fA_mu[0][0]['Bsw_t0']+fA_sig[0][0]['Bsw_t0'],color=[0.27,0.49,0.79],alpha=Alpha,linewidth=0)
ax[3,0].plot(uA,fA_mu[0][0]['Bsw_t0'],'-o',color=[0.27,0.49,0.79],markersize=ms) 
ax[3,0].set(position=[0.07,0.06,0.42,0.2],xlim=xlim,xlabel='Stand age, years',ylim=[0,400],ylabel='Stemwood biomas (Mg C/ha)')

ax[0,1].plot([-5,200],[0,0],'k-')
ax[0,1].fill_between(uA,fA_mu[1][0]['N_t0']-fA_sig[1][0]['N_t0'],
  fA_mu[1][0]['N_t0']+fA_sig[1][0]['N_t0'],color=[0.27,0.49,0.79],alpha=Alpha,linewidth=0)
ax[0,1].plot(uA,fA_mu[1][0]['N_t0'],'-o',color=[0.27,0.49,0.79],markersize=ms) 
ax[0,1].set(position=[0.58,0.78,0.42,0.2],xlim=xlim,xlabel='',ylim=[0,6000],ylabel='Stand density (stems/ha)')

ax[1,1].plot([-5,200],[0,0],'k-')
ax[1,1].fill_between(uA,fA_mu[1][0]['Dam_t0']-fA_sig[1][0]['Dam_t0'],
  fA_mu[1][0]['Dam_t0']+fA_sig[1][0]['Dam_t0'],color=[0.27,0.49,0.79],alpha=Alpha,linewidth=0)
ax[1,1].plot(uA,fA_mu[1][0]['Dam_t0'],'-o',color=[0.27,0.49,0.79],markersize=ms) 
ax[1,1].set(position=[0.58,0.54,0.42,0.2],xlim=xlim,xlabel='',ylim=[0,60],ylabel='Diameter mean (cm)')

ax[2,1].plot([-5,200],[0,0],'k-')
ax[2,1].fill_between(uA,fA_mu[1][0]['H_obs_t0']-fA_sig[1][0]['H_obs_t0'],
  fA_mu[1][0]['H_obs_t0']+fA_sig[1][0]['H_obs_t0'],color=[0.27,0.49,0.79],alpha=Alpha,linewidth=0)
ax[2,1].plot(uA,fA_mu[1][0]['H_obs_t0'],'-o',color=[0.27,0.49,0.79],markersize=ms) 
ax[2,1].set(position=[0.58,0.3,0.42,0.2],xlim=xlim,xlabel='',ylim=[0,50],ylabel='Height mean (m)')

ax[3,1].plot([-5,200],[0,0],'k-')
ax[3,1].fill_between(uA,fA_mu[1][0]['Bsw_t0']-fA_sig[1][0]['Bsw_t0'],
  fA_mu[1][0]['Bsw_t0']+fA_sig[1][0]['Bsw_t0'],color=[0.27,0.49,0.79],alpha=Alpha,linewidth=0)
ax[3,1].plot(uA,fA_mu[1][0]['Bsw_t0'],'-o',color=[0.27,0.49,0.79],markersize=ms) 
ax[3,1].set(position=[0.58,0.06,0.42,0.2],xlim=xlim,xlabel='Stand age, years',ylim=[0,400],ylabel='Stemwood biomas (Mg C/ha)')

gu.axletters(ax,plt,0.03,0.85,['Douglas-fir','Western hemlock','','','',''],0.05)
plt.savefig(r'G:\My Drive\Figures\GHGBenefit_EP703\AgeResponse_Yield_StandLevel.png',format='png',dpi=900)


#%% Quality assurance check
# Save site summaries to spreadsheet for inspection.

cols=['ID_Site','Spc1_ID_t0','Year_t0','TSF_t0','TSF_t1', \
      'C1','C2','C3','C4','C5','C6','Cmu','Cse', \
      'F1','F2','F3','F4','F5','F6','Fmu','Fse','DA_mu','DR_mu']

#sy='H_obs_t0'
sy='Bsw_G'
dose=225

df0=pd.DataFrame(data=np.nan*np.ones((360,len(cols))),columns=cols)
cnt=0
for i in range(uSite.size):
    ind=np.where((sobs['ID_Site']==uSite[i]))[0]
    if ind.size==0: continue
    uTSF=np.unique(np.column_stack((sobs['TSF_t0'][ind],sobs['TSF_t1'][ind])).astype(float),axis=0)
    for j in range(uTSF.shape[0]):
        df0.loc[cnt,'ID_Site']=uSite[i]
        df0.loc[cnt,'Spc1_ID_t0']=sobs['Spc1_ID_t0'][ind[0]]
        df0.loc[cnt,'TSF_t0']=uTSF[j,0]
        df0.loc[cnt,'TSF_t1']=uTSF[j,1]
        iC=np.where((sobs['ID_Site']==uSite[i]) & (sobs['TSF_t0']==uTSF[j,0]) & (sobs['N_Dose']==0))[0]
        iF=np.where((sobs['ID_Site']==uSite[i]) & (sobs['TSF_t0']==uTSF[j,0]) & (sobs['N_Dose']==dose))[0]
        df0.loc[cnt,'Year_t0']=sobs['Year_t0'][iC[0]]
        for k in range(iC.size):
            if k==0: df0.loc[cnt,'C1']=sobs[sy][iC[k]]
            elif k==1: df0.loc[cnt,'C2']=sobs[sy][iC[k]]
            elif k==2: df0.loc[cnt,'C3']=sobs[sy][iC[k]]
            elif k==3: df0.loc[cnt,'C4']=sobs[sy][iC[k]]
            elif k==4: df0.loc[cnt,'C5']=sobs[sy][iC[k]]
            elif k==5: df0.loc[cnt,'C6']=sobs[sy][iC[k]]                
        df0.loc[cnt,'Cmu']=np.mean(sobs[sy][iC])
        df0.loc[cnt,'Cse']=np.std(sobs[sy][iC])/np.sqrt(iC.size)
        for k in range(iF.size):
            if k==0: df0.loc[cnt,'F1']=sobs[sy][iF[k]]
            elif k==1: df0.loc[cnt,'F2']=sobs[sy][iF[k]]
            elif k==2: df0.loc[cnt,'F3']=sobs[sy][iF[k]]
            elif k==3: df0.loc[cnt,'F4']=sobs[sy][iF[k]]
            elif k==4: df0.loc[cnt,'F5']=sobs[sy][iF[k]]
            elif k==5: df0.loc[cnt,'F6']=sobs[sy][iF[k]]
        df0.loc[cnt,'Fmu']=np.mean(sobs[sy][iF])
        df0.loc[cnt,'Fse']=np.std(sobs[sy][iF])/np.sqrt(iF.size)
        df0.loc[cnt,'DA_mu']=np.round(df0.loc[cnt,'Fmu']-df0.loc[cnt,'Cmu'],2)
        df0.loc[cnt,'DR_mu']=np.round((df0.loc[cnt,'Fmu']-df0.loc[cnt,'Cmu'])/df0.loc[cnt,'Cmu']*100,0)
        cnt=cnt+1
    
df0.to_excel(PathProject + '\\Processed\\QA_' + sy + '.xlsx',index=False)


#%% Regression analysis

#------------------------------------------------------------------------------
# Build a dataframe for the purpose of conducting regression analysis
#------------------------------------------------------------------------------

cols=['ID_Site','DA_mu','Age_Init','N_Dose','Spc1_ID_t0','SI','N_Init','Bsw_Init','pH_Min_Init','pH_Humus_Init','P_Min_Init','P_Humus_Init','Ca_Min_Init']#,'ndep_ann_r','tmean_ann_r','prcp_ann_r']

yv='Bsw_G'
dose=225
tsf0_max=9

dfR=pd.DataFrame(data=np.nan*np.ones((360,len(cols))),columns=cols)
cnt=0
for i in range(uSite.size):
    ind=np.where((sobs['ID_Site']==uSite[i]))[0]
    if ind.size==0: continue

    dfR.loc[cnt,'ID_Site']=uSite[i]
    
    iC=np.where( (sobs['ID_Site']==uSite[i]) & (sobs['TSF_t0']<=tsf0_max) & (sobs['N_Dose']==0) )[0]
    iF=np.where( (sobs['ID_Site']==uSite[i]) & (sobs['TSF_t0']<=tsf0_max) & (sobs['N_Dose']==dose) )[0]
    if iF.size==0: continue
    dfR.loc[cnt,'DA_mu']=np.mean(sobs[yv][iF])-np.mean(sobs[yv][iC])
    
    dfR.loc[cnt,'Age_Init']=np.min(sobs['Age_t0'][ind])
    dfR.loc[cnt,'N_Dose']=sobs['N_Dose'][iF[0]]
            
    col_adj=4
    
    for j in range(len(cols)-col_adj):
        dfR.loc[cnt,cols[j+col_adj]]=sobs[cols[j+col_adj]][ind[0]]
    cnt=cnt+1
  
dose=450
for i in range(uSite.size):
    ind=np.where((sobs['ID_Site']==uSite[i]))[0]
    if ind.size==0: continue

    dfR.loc[cnt,'ID_Site']=uSite[i]
    
    iC=np.where( (sobs['ID_Site']==uSite[i]) & (sobs['TSF_t0']<=tsf0_max) & (sobs['N_Dose']==0) )[0]
    iF=np.where( (sobs['ID_Site']==uSite[i]) & (sobs['TSF_t0']<=tsf0_max) & (sobs['N_Dose']==dose) )[0]
    if iF.size==0: continue
    dfR.loc[cnt,'DA_mu']=np.mean(sobs[yv][iF])-np.mean(sobs[yv][iC])
    
    dfR.loc[cnt,'Age_Init']=np.min(sobs['Age_t0'][ind])
    dfR.loc[cnt,'N_Dose']=sobs['N_Dose'][iF[0]]
            
    col_adj=4
    
    for j in range(len(cols)-col_adj):
        dfR.loc[cnt,cols[j+col_adj]]=sobs[cols[j+col_adj]][ind[0]]
    cnt=cnt+1

# Drop excess rows
dfR=dfR[np.isnan(dfR['DA_mu'])==False]
dfR=dfR.reset_index(drop=True)

# Drop site (and species)
dfR=dfR.drop(['ID_Site','Spc1_ID_t0'],axis=1)

dfR.shape

#------------------------------------------------------------------------------
# Exclude sites with no response data
#------------------------------------------------------------------------------
    
dfR=dfR.dropna()
dfR=dfR.reset_index(drop=True)

dfR.shape

#------------------------------------------------------------------------------
# Correlation matrix
#------------------------------------------------------------------------------

# Simple bivariate correlation
dfR.corr()

# Partial correlation
cp=gu.PartialCorrelation(dfR.to_numpy())
plt.bar(np.arange(1,cp.shape[0]+1),cp[:,0])

#------------------------------------------------------------------------------
# OLS regression results
#------------------------------------------------------------------------------
 
# Dependent variable
y=dfR['DA_mu'].values.astype(float)

# Predictor variables
xv=['SI','Age_Init','pH_Min_Init','N_Dose']
xv=['Age_Init']
#xv=['N_Dose','pH_Min_Init','P_Min_Init','Ca_Min_Init']
iSI=0;iA=1;iPH=2;iND=3;iT=4;iW=5

X=np.zeros((y.shape[0],len(xv)))
for i in range(len(xv)):
    X[:,i]=dfR[xv[i]].values.astype(float)

# Standardize predictor variables
mu=np.zeros(X.shape[1])
sig=np.zeros(X.shape[1])
Xz=X.copy()
for i in range(X.shape[1]):
    mu[i]=np.mean(X[:,i])
    sig[i]=np.std(X[:,i])
    Xz[:,i]=(X[:,i]-mu[i])/sig[i]

# Add constant
Xz=sm.tools.tools.add_constant(Xz)

# Fit model
md=sm.OLS(y,Xz).fit()
md.summary()


plt.close('all')
# SI
x0=np.min(dfR['SI']); X=np.zeros(len(xv)+1); X[0]=1; X[iSI+1]=(x0-mu[iSI])/sig[iSI]; y0=md.predict(X)
x1=np.max(dfR['SI']); X=np.zeros(len(xv)+1); X[0]=1; X[iSI+1]=(x1-mu[iSI])/sig[iSI]; y1=md.predict(X)
plt.figure(1)
plt.plot([x0,x1],[y0,y1],'-')

# Age
x0=np.min(dfR['Age_Init']); X=np.zeros(len(xv)+1); X[0]=1; X[iA+1]=(x0-mu[iA])/sig[iA]; y0=md.predict(X)
x1=np.max(dfR['Age_Init']); X=np.zeros(len(xv)+1); X[0]=1; X[iA+1]=(x1-mu[iA])/sig[iA]; y1=md.predict(X)
plt.figure(2)
plt.plot([x0,x1],[y0,y1],'-')

# pH
x0=np.min(dfR['pH_Min_Init']); X=np.zeros(len(xv)+1); X[0]=1; X[iPH+1]=(x0-mu[iPH])/sig[iPH]; y0=md.predict(X)
x1=np.max(dfR['pH_Min_Init']); X=np.zeros(len(xv)+1); X[0]=1; X[iPH+1]=(x1-mu[iPH])/sig[iPH]; y1=md.predict(X)
plt.figure(3)
plt.plot([x0,x1],[y0,y1],'-')

# Ndep
x0=np.min(dfR['ndep_ann_r']); X=np.zeros(len(xv)+1); X[0]=1; X[iND+1]=(x0-mu[iND])/sig[iND]; y0=md.predict(X)
x1=np.max(dfR['ndep_ann_r']); X=np.zeros(len(xv)+1); X[0]=1; X[iND+1]=(x1-mu[iND])/sig[iND]; y1=md.predict(X)
plt.figure(4)
plt.plot([x0,x1],[y0,y1],'-')



#------------------------------------------------------------------------------
# Scatterplots
#------------------------------------------------------------------------------


# SI
plt.close('all')
ind=np.where(dfR['Spc1_ID_t0']=='FD')[0]
plt.plot(dfR.loc[ind,'SI'],dfR.loc[ind,'DA_mu'],'o')
ind=np.where(dfR['Spc1_ID_t0']=='HW')[0]
plt.plot(dfR.loc[ind,'SI'],dfR.loc[ind,'DA_mu'],'s')

# Age
plt.close('all')
ind=np.where(dfR['Spc1_ID_t0']=='FD')[0]
plt.plot(dfR.loc[ind,'Age_Init'],dfR.loc[ind,'DA_mu'],'o')
ind=np.where(dfR['Spc1_ID_t0']=='HW')[0]
plt.plot(dfR.loc[ind,'Age_Init'],dfR.loc[ind,'DA_mu'],'s')

# Stand density
plt.close('all')
ind=np.where(dfR['Spc1_ID_t0']=='FD')[0]
plt.plot(dfR.loc[ind,'N_Init'],dfR.loc[ind,'DA_mu'],'o')
ind=np.where(dfR['Spc1_ID_t0']=='HW')[0]
plt.plot(dfR.loc[ind,'N_Init'],dfR.loc[ind,'DA_mu'],'s')

# pH
plt.close('all')
ind=np.where(dfR['Spc1_ID_t0']=='FD')[0]
plt.plot(dfR.loc[ind,'pH_Min_Init'],dfR.loc[ind,'DA_mu'],'o')
ind=np.where(dfR['Spc1_ID_t0']=='HW')[0]
plt.plot(dfR.loc[ind,'pH_Min_Init'],dfR.loc[ind,'DA_mu'],'s')












