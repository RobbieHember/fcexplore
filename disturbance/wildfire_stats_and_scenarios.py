'''

WILDFIRE STATS AND SCENARIOS BY BGC ZONE

See documentation.

'''

#%% Import Python modules

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import geopandas as gpd
import time
import gc as garc
import scipy.stats as stats
import statsmodels.api as sm
from fcgadgets.utilities import utilities_general as gu
from fcgadgets.utilities import utilities_gis as gis
from fcgadgets.utilities import utilities_inventory as invu
from fcgadgets.taz import wildfire_stat_models as wfsm

#%% Path of file to store stats and scenarios

PathData=r'G:\My Drive\Data\Wildfire\Wildfire_Stats_Scenarios_By_BGCZ\Wildfire_Stats_Scenarios_By_BGCZ.pkl'
PathFigures=r'G:\My Drive\Figures\Wildfire\Wildfire_Stats_Sceanrios_By_BGCZ'

#%% Set figure properties

params=gu.Import_GraphicsParameters('spyder_fs7')
plt.rcParams.update(params)

#%% Import BGC zones

zBECZ=gis.OpenGeoTiff(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\VRI\becz.tif')
zBECZf=zBECZ['Data'].flatten()
zBECZs=zBECZ['Data'][0::50,0::50].flatten()
dBECZ=gu.ReadExcel(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\VRI\becz_lut.xlsx')

#%% Total area burned

tv_obs=np.arange(2009,2020,1)
A=np.zeros(tv_obs.size)
for iT in range(tv_obs.size):
    print(tv_obs[iT])
    zFire=gis.OpenGeoTiff(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Disturbances\PROT_HISTORICAL_FIRE_POLYS_SP_' + str(tv_obs[iT]) + '.tif')
    ind=np.where(zFire['Data'].flatten()>0)[0]
    A[iT]=ind.size
    del zFire
    garc.collect()

plt.bar(tv_obs,A/1e6,0.7,facecolor=[0.5,0,0])


#%% Import wildfire data and calculate prob occ

# Initialize variables
tv_obs=np.arange(1920,2020,1)
wfss={}
for iZ in range(dBECZ['ZONE'].size):
    nam=dBECZ['ZONE'][iZ]
    wfss[nam]={}
    wfss[nam]['Ao']=np.zeros(tv_obs.size)
    ind=np.where( (zBECZs==dBECZ['VALUE'][iZ]) )[0]
    wfss[nam]['Oc']=np.zeros((tv_obs.size,ind.size))
    ind=np.where( (zBECZf==dBECZ['VALUE'][iZ]) )[0]
    wfss[nam]['Azone']=ind.size

# Populate variables from rasterized wildfire perimiter data
for iT in range(tv_obs.size):
    print(tv_obs[iT])
    zFire=gis.OpenGeoTiff(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Disturbances\PROT_HISTORICAL_FIRE_POLYS_SP_' + str(tv_obs[iT]) + '.tif')
    zFiref=zFire['Data'].flatten()
    zFires=zFire['Data'][0::50,0::50].flatten()
    for iZ in range(dBECZ['ZONE'].size):
        nam=dBECZ['ZONE'][iZ]
        ind=np.where( (zBECZf==dBECZ['VALUE'][iZ]) & (zFiref>0) )[0]
        wfss[nam]['Ao'][iT]=wfss[nam]['Ao'][iT]+ind.size
        ind=np.where( (zBECZs==dBECZ['VALUE'][iZ]) )[0]
        wfss[nam]['Oc'][iT,:]=zFires[ind]
    del zFire,zFiref,zFires
    garc.collect()

# Annual probability of occurrence
for k in wfss.keys():    
    wfss[k]['Po_obs']=np.sum(wfss[k]['Oc'])/wfss[k]['Oc'].size

# Pareto distribution fits
for zone in wfss.keys():
    Ao=wfss[zone]['Ao']/wfss[zone]['Azone']
    shape,loc,scale=stats.pareto.fit(Ao)
    wfss[zone]['Beta_Pareto']=np.array([shape,loc,scale])

# Find relationship between shape parameter and annual prob occ
n_stand=1000
n_t=2000
modifier=np.arange(0.45,1.5,0.05)
for zone in wfss.keys():    
    shape=np.zeros(modifier.size)
    Po=np.zeros(modifier.size)
    for iMod in range(modifier.size):
        b0=wfss[zone]['Beta_Pareto'].copy()
        b0[0]=modifier[iMod]*b0[0]
        y=wfsm.GenerateDisturbancesFromPareto(n_t,n_stand,b0)
        shape[iMod]=b0[0].copy()
        Po[iMod]=np.sum(y)/y.size*100
    
    plt.close('all')
    fig,ax=plt.subplots(1)
    ax.plot(Po,shape,'-',linewidth=2)
    ax.set(xscale='log',yscale='log',xlabel='Annual probability of occurrence (%/yr)',ylabel='Pareto shape parameter')
    
    # Fit model of prob occ vs
    x=sm.tools.tools.add_constant(np.log(Po))
    y=np.log(shape)
    md=sm.OLS(y,x).fit()
    #md.summary()
    beta=md.params
    wfss[zone]['Pareto_shape_for_Po']=beta
    
    # Check that relationship is good
    #Po_hat=Po.copy()
    #shape_hat=np.exp(beta[1]*np.log(Po_hat)+beta[0])
    #plt.plot(Po_hat,shape_hat,'x')
    

# Remove occurrence data (no need to save)
for zone in wfss.keys():
    del wfss[zone]['Oc']

# Save
gu.opickle(PathData,wfss)

# Plot annual time series
plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,5.5))
ax.bar(tv_obs,wfss['IDF']['Ao']/wfss['IDF']['Azone']*100,1,facecolor=[0.5,0,0])
ax.set(position=[0.06,0.1,0.92,0.86],xlim=[tv_obs[0]-0.5,tv_obs[-1]+1+0.5],xticks=np.arange(tv_obs[0],tv_obs[-1]+2,10),
       ylabel='Probability (% yr$^-$$^1$)')
gu.PrintFig(PathFigures + '\\Wilfire_ts_IDF','png',900)

# Plot annual probability of occurrence
Po_obs=np.zeros(len(wfss.keys())); cnt=0
for k in wfss.keys():
    Po_obs[cnt]=wfss[k]['Po_obs']; cnt=cnt+1
    
plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,5.5))
ax.bar(np.arange(1,Po_obs.size+1),Po_obs*100,facecolor=[0.8,0.8,0.8])
ax.set(position=[0.07,0.1,0.92,0.86],xlim=[0.5,Po_obs.size+.5],xticks=np.arange(1,Po_obs.size+1),
       xticklabels=list(wfss.keys()),ylabel='Probability (% yr$^-$$^1$)')
gu.PrintFig(PathFigures + '\\Wilfire_Occurrence_ByBGCZ','png',900)

#%% Scenario development (deterministic component)

# Import parameters
#wfss=gu.ipickle(PathData)

tv_scn=np.arange(-2000,2201,1)

for zone in wfss.keys():

    Po_obs=wfss[zone]['Po_obs']*100

    wfss[zone]['Po_Det_WF_Scn1']=Po_obs*np.ones(tv_scn.size)
    wfss[zone]['Po_Det_WF_Scn2']=Po_obs*np.ones(tv_scn.size)
    wfss[zone]['Po_Det_WF_Scn3']=Po_obs*np.ones(tv_scn.size)
    wfss[zone]['Po_Det_WF_Scn4']=Po_obs*np.ones(tv_scn.size)
    
    wfss[zone]['Po_Det_WF_Scn1']=1.5*Po_obs
    ind=np.where( (tv_scn<1920-50) )[0]; 
    wfss[zone]['Po_Det_WF_Scn2'][ind]=1.5*Po_obs
    wfss[zone]['Po_Det_WF_Scn3'][ind]=1.5*Po_obs
    wfss[zone]['Po_Det_WF_Scn4'][ind]=1.5*Po_obs

    ind=np.where( (tv_scn>=1920-50) & (tv_scn<1920) )[0];
    wfss[zone]['Po_Det_WF_Scn2'][ind]=np.interp(tv_scn[ind],np.array([1920-50,1920]),np.array([1.5*Po_obs,Po_obs]))
    wfss[zone]['Po_Det_WF_Scn3'][ind]=np.interp(tv_scn[ind],np.array([1920-50,1920]),np.array([1.5*Po_obs,Po_obs]))
    wfss[zone]['Po_Det_WF_Scn4'][ind]=np.interp(tv_scn[ind],np.array([1920-50,1920]),np.array([1.5*Po_obs,Po_obs]))

    ind=np.where( (tv_scn>2020) )[0]; 
    wfss[zone]['Po_Det_WF_Scn2'][ind]=np.interp(tv_scn[ind],np.array([2020,2200]),np.array([Po_obs,1.0*Po_obs]))
    wfss[zone]['Po_Det_WF_Scn3'][ind]=np.interp(tv_scn[ind],np.array([2020,2200]),np.array([Po_obs,1.5*Po_obs]))
    wfss[zone]['Po_Det_WF_Scn4'][ind]=np.interp(tv_scn[ind],np.array([2020,2200]),np.array([Po_obs,2.0*Po_obs]))

# Save
gu.opickle(PathData,wfss)

# Plot deteriministic component of scenarios
zone='IDF'
plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,5.5))
rc=patches.Rectangle((1920,0),100,1,facecolor=[0.85,0.85,0.85])
ax.add_patch(rc)
ax.plot(tv_scn,wfss[zone]['Po_Det_WF_Scn1'],'b-',linewidth=1.5,label='Wildfire occurrence Scn. 1 (pre-industrial baseline)')
ax.plot(tv_scn,wfss[zone]['Po_Det_WF_Scn2'],'g:',linewidth=1.5,label='Wildfire occurrence Scn. 2')
ax.plot(tv_scn,wfss[zone]['Po_Det_WF_Scn3'],'c-.',linewidth=1.5,label='Wildfire occurrence Scn. 3')
ax.plot(tv_scn,wfss[zone]['Po_Det_WF_Scn4'],'r--',linewidth=1.5,label='Wildfire occurrence Scn. 4')
ax.annotate('Modern era',(1920+50,0.65),ha='center')
ax.legend(loc='upper left',frameon=False)
ax.set(position=[0.065,0.17,0.92,0.8],ylim=[0,1],xlim=[1000-0.5,tv_scn[-1]+1+0.5],
       xticks=np.arange(tv_scn[0],tv_scn[-1]+100,100),
       ylabel='Annual probability (% yr$^-$$^1$)',
       xlabel='Time, calendar year')
gu.PrintFig(PathFigures + '\\Wildfire_Scenarios_ts_IDF','png',900)

#%% Apply random component to scenarios

n_stand=2000

for zone in wfss.keys():

    # Initialize annual probability of occurrence (final with deterministic and
    # random components)
    wfss[zone]['Po_WF_Scn1']=np.zeros((tv_scn.size,n_stand))
    
    for iT in range(tv_scn.size):
        # Adjust shape parameter to match specified annual probability of 
        # occurrence from the deterministic component
        b0=wfss[zone]['Beta_Pareto'].copy()       
        b_shape=wfss[zone]['Pareto_shape_for_Po'].copy()
        b0[0]=np.exp(b_shape[1]*np.log(wfss[zone]['Po_Det_WF_Scn1'][iT])+b_shape[0])
        wfss[zone]['Po_WF_Scn1'][iT,:]=wfsm.GenerateDisturbancesFromPareto(1,n_stand,b0)

Po=np.sum(wfss[zone]['Po_WF_Scn1'],axis=1)/n_stand*100
plt.close('all'); 
plt.bar(tv_scn[3000:],Po[3000:],1)
plt.plot(tv_scn[3000:],gu.movingave(Po[3000:],50,'historical'),linewidth=3)


#%% Histogram of burn severity rating in 2017

# Import data
dBSR=gu.ReadExcel(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Disturbances\VEG_BURN_SEVERITY_SP.xlsx')
zBSR=gis.OpenGeoTiff(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Disturbances\VEG_BURN_SEVERITY_SP_2017.tif')

# Check data
#plt.matshow(zBSR.Data)

# Get number of occurrences (exclude unknown class)
u=np.unique(zBSR.Data[zBSR.Data!=5])
N=np.zeros(u.size)
for i in range(u.size):
    ind=np.where(zBSR.Data.flatten()==u[i])[0]
    N[i]=ind.size

# Percent frequency
p=N/np.sum(N)*100
print(np.round(p))

# Define pre-industrial distribution
p_pi=np.array([25,45,25,5])
np.sum(p_pi)

# Plot histogram
plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(10,5.5))
ax.bar(u-0.16,p,0.3,facecolor=[0.8,0.8,0.8],label='Modern era (observed)')
ax.bar(u+0.16,p_pi,0.3,facecolor=[0.25,0.5,1],label='Pre-industrial era')
ax.set(position=[0.1,0.1,0.84,0.86],xticks=np.arange(1,5,1),xticklabels=dBSR['Code'][0:5],
       ylabel='Probability (%)')
ax.legend(loc='upper right',frameon=False)
gu.PrintFig(PathFigures + '\\BurnSeverityDistribution2017','png',900)



#%% Relationship between onset-spread model parameters and annual prob occ

# Import BGC-specific paramters
#wfss=gu.ipickle(r'C:\Users\rhember\Documents\Data\DisturbanceStatsByBGCZ\Wildfire.pkl')

# Import BGC-specific annual probability of occurrence
Pa=gu.ipickle(r'C:\Users\rhember\Documents\Data\DisturbanceStatsByBGCZ\Wildfire_FireballVsAnnProb.pkl')
binO=np.linspace(0.0000001,0.0000004,10)
binS=np.linspace(0.02,0.14,10)

z=np.reshape(Pa[:,2],(binO.size,binS.size))
plt.close('all'); plt.matshow(z)
plt.colorbar()

# Fit model of prob occ as function of prob onset and prob spread
x=np.column_stack((Pa[:,0:2],Pa[:,0]*Pa[:,1]))
x=sm.tools.tools.add_constant(x)
y=Pa[:,2]
md=sm.OLS(y,x).fit()
md.summary()

# Get curves
binOb=np.linspace(0.0000001,0.0000004,200)
binSb=np.linspace(0.02,0.14,200)
wfss['IDF']['OS contour']=np.zeros((4000,2)); cnt_IDF=0
wfss['ESSF']['OS contour']=np.zeros((4000,2)); cnt_ESSF=0
for i in range(binOb.size):
    for j in range(binSb.size):
        x=np.column_stack((1,binOb[i],binSb[j],binOb[i]*binSb[j]))
        y=md.predict(x)
        if np.abs(y-100*wfss['IDF']['Pa'])<0.0001:
            wfss['IDF']['OS contour'][cnt_IDF,:]=[binOb[i],binSb[j]]
            cnt_IDF=cnt_IDF+1
        if np.abs(y-100*wfss['ESSF']['Pa'])<0.0001:
            wfss['ESSF']['OS contour'][cnt_ESSF,:]=[binOb[i],binSb[j]]
            cnt_ESSF=cnt_ESSF+1    
wfss['IDF']['OS contour']=wfss['IDF']['OS contour'][0:cnt_IDF,:]
wfss['ESSF']['OS contour']=wfss['ESSF']['OS contour'][0:cnt_ESSF,:]

plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(7,6.5))
ax.plot([0.0000002,0.0000002],[0.02,0.08],linewidth=4,color=[0.8,0.8,0.8])
ax.plot(wfss['IDF']['OS contour'][:,0],wfss['IDF']['OS contour'][:,1],'g-',label='IDF',linewidth=1.5) 
ax.plot(wfss['ESSF']['OS contour'][:,0],wfss['ESSF']['OS contour'][:,1],'g--',label='ESSF',linewidth=1.5)  
ax.legend(loc='upper right',frameon=False)
ax.set(position=[0.16,0.16,0.8,0.8],xlabel='Probability of onset',ylabel='Probability of spread',
       xlim=[0.0000001,0.0000004],ylim=[0.02,0.08])
#ax.grid(True)
gu.PrintFig(r'G:\My Drive\Figures\DisturbanceStatsByBGCZ\Wilfire_OSpar_vs_Po','png',900)



#%% Charcoal record digitized from Marlon et al. (2012)
# This must already be a moving average of the z-score as it does not exceed abs(1).

dChar=gu.ReadExcel(r'G:\My Drive\Data\DigitizedFigures\Marlonetal2012_Fig2.xlsx')
tvChar=np.arange(500,2000,1)
yChar=np.interp(tvChar,dChar['Time'],dChar['Charcoal influx zscore'])
plt.figure()
plt.plot(tvChar,yChar,'.')

# Adjust so that it is relative mean during observed wildfire record
it_mu=np.where(tvChar>=1920)[0]
yChar_adj=(yChar-np.mean(yChar[it_mu]))*1

plt.plot(tvChar,yChar_adj,'.')
plt.grid('on')

yChar_adj_full=np.ones(tv_full.size)
ind=np.where( (tv_full>=tvChar[0]) & (tv_full<=tvChar[-1]) )[0]
yChar_adj_full[ind]=yChar_adj
ind=np.where( (tv_full<tvChar[0]) )[0]
yChar_adj_full[ind]=0.3
ind=np.where( (tv_full>tvChar[-1]) )[0]
yChar_adj_full[ind]=0
yChar_adj_full_ma=gu.movingave(yChar_adj_full,ivl_ma,'centre')

d=Po_full_ma_zscore-yChar_adj_full_ma


plt.figure()
plt.plot(tv_full,Po_full_ma_zscore,'.')
plt.plot(tv_full,yChar_adj_full_ma,'.')


plt.close('all'); 
plt.plot(Po_full_ma_zscore,'-')

#%%


