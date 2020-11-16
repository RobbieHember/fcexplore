'''

BEETLE STATS AND SCENARIOS BY BGC ZONE

See documentation.

'''

#%% Import Python modules

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pickle
import geopandas as gpd
import openpyxl
import gdal
import time
import gc as garc
import scipy.stats as stats
import statsmodels.api as sm
from numpy import matlib as mb
from shapely.geometry import Point
from fcgadgets.utilities import utilities_general as gu
from fcgadgets.utilities import utilities_gis as gis
from fcgadgets.utilities import utilities_inventory as invu
from fcgadgets.taz import wildfire_stat_models as wfsm
from fcgadgets.taz import general_stat_models as gensm

#%% Path of file to store stats and scenarios

PathData=r'G:\My Drive\Data\Beetle_Stats_Scenarios_By_BGCZ\IBM_Stats_Scenarios_By_BGCZ.pkl'
PathFigures=r'G:\My Drive\Figures\Beetles\Beetle_Stats_Sceanrios_By_BGCZ'

#%% Plotting parameters

params=gu.Import_GraphicsParameters('spyder_fs7')
plt.rcParams.update(params)

flg_prnt=0

#%% Import BGC zones

zBECZ=gis.OpenGeoTiff(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\VRI\becz.tif')
zBECZf=zBECZ['Data'].flatten()
zBECZs=zBECZ['Data'][0::50,0::50].flatten()
dBECZ=gu.ReadExcel(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\VRI\becz_lut.xlsx')

#%% Import beetle data and calculate prob occ

# Initialize variables
tv_obs=np.arange(1950,2020,1)
#ibmss=gu.ipickle(PathData)

ibmss={}
for iZ in range(dBECZ['ZONE'].size):
    nam=dBECZ['ZONE'][iZ]
    ibmss[nam]={}
    ibmss[nam]['Ao']=np.zeros(tv_obs.size)
    ind=np.where( (zBECZs==dBECZ['VALUE'][iZ]) )[0]
    ibmss[nam]['Oc']=np.zeros((tv_obs.size,ind.size))
    ind=np.where( (zBECZf==dBECZ['VALUE'][iZ]) )[0]
    ibmss[nam]['Azone']=ind.size

# Populate variables from rasterized wildfire perimiter data
for iT in range(tv_obs.size):
    print(tv_obs[iT])
    zDist=gis.OpenGeoTiff(r'Z:\!Workgrp\Forest Carbon\Data\BC1ha\Disturbances\PEST_INFESTATION_POLY_IBM_SeverityClass_' + str(tv_obs[iT]) + '.tif')
    zDistf=zDist['Data'].flatten()
    zDists=zDist['Data'][0::50,0::50].flatten()
    for iZ in range(dBECZ['ZONE'].size):
        nam=dBECZ['ZONE'][iZ]
        ind=np.where( (zBECZf==dBECZ['VALUE'][iZ]) & (zDistf>0) )[0]
        ibmss[nam]['Ao'][iT]=ibmss[nam]['Ao'][iT]+ind.size
        ind=np.where( (zBECZs==dBECZ['VALUE'][iZ]) )[0]
        ibmss[nam]['Oc'][iT,:]=zDists[ind]
    del zDist,zDistf,zDists
    garc.collect()

# Annual probability of occurrence
for k in ibmss.keys():   
    ind=np.where(ibmss[k]['Oc'].flatten()>0)[0]
    ibmss[k]['Po_obs']=ind.size/ibmss[k]['Oc'].size

# Pareto distribution fits
for zone in ibmss.keys():
    Ao=ibmss[zone]['Ao']/ibmss[zone]['Azone']
    shape,loc,scale=stats.pareto.fit(Ao)
    ibmss[zone]['Beta_Pareto']=np.array([shape,loc,scale])

# Find relationship between shape parameter and annual prob occ
n_stand=1000
n_t=2000
modifier=np.arange(0.45,1.5,0.05)
for zone in ibmss.keys():    
    shape=np.zeros(modifier.size)
    Po=np.zeros(modifier.size)
    for iMod in range(modifier.size):
        b0=ibmss[zone]['Beta_Pareto'].copy()
        b0[0]=modifier[iMod]*b0[0]
        y=gensm.GenerateDisturbancesFromPareto(n_t,n_stand,b0)
        shape[iMod]=b0[0].copy()
        Po[iMod]=np.sum(y)/y.size*100
    
    if np.sum(Po)>0:
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
        ibmss[zone]['Pareto_shape_for_Po']=beta
    else:
        ibmss[zone]['Pareto_shape_for_Po']=np.array([-999,-999])
    
    # Check that relationship is good
    #Po_hat=Po.copy()
    #shape_hat=np.exp(beta[1]*np.log(Po_hat)+beta[0])
    #plt.plot(Po_hat,shape_hat,'x')

# Remove occurrence data (no need to save)
ibmss_bk=ibmss.copy()
for zone in ibmss.keys():
    del ibmss[zone]['Oc']
    
# Save
gu.opickle(PathData,ibmss)


# Plot PDF
zone='SBPS'
plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,5.5))
x=np.linspace(0,1,100)
b=ibmss[zone]['Beta_Pareto']
b=[0.8,0.001,0.00000001]
y=stats.pareto.pdf(x=x,b=b[0],loc=b[1],scale=b[2])
ax.plot(x,y,'r-',lw=1,label='pareto pdf')
ax.set(yscale='log')

# Plot annual time series
zone='SBPS'
plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,5.5))
ax.bar(tv_obs,ibmss[zone]['Ao']/ibmss[zone]['Azone']*100,1,facecolor=[0.5,0,0])
ax.set(position=[0.06,0.1,0.92,0.86],xlim=[tv_obs[0]-0.5,tv_obs[-1]+1+0.5],xticks=np.arange(tv_obs[0],tv_obs[-1]+2,10),
       ylabel='Probability (% yr$^-$$^1$)')
gu.PrintFig(PathFigures + '\\IBM_ts_IDF','png',900)

# Plot annual probability of occurrence
Po_obs=np.zeros(len(ibmss.keys())); cnt=0
for k in ibmss.keys():
    Po_obs[cnt]=ibmss[k]['Po_obs']; cnt=cnt+1
    
plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,5.5))
ax.bar(np.arange(1,Po_obs.size+1),Po_obs*100,facecolor=[0.8,0.8,0.8])
ax.set(position=[0.07,0.1,0.92,0.86],xlim=[0.5,Po_obs.size+.5],xticks=np.arange(1,Po_obs.size+1),
       xticklabels=list(ibmss.keys()),ylabel='Probability (% yr$^-$$^1$)')
gu.PrintFig(PathFigures + '\\IBM_Occurrence_ByBGCZ','png',900)

#%% Scenario development (deterministic component)

tv_scn=np.arange(-2000,2201,1)
yr_mod_start=1951
rr_pi=0.25

for zone in ibmss.keys():

    Po_obs=ibmss[zone]['Po_obs']*100

    ibmss[zone]['Po_Det_Scn1']=Po_obs*np.ones(tv_scn.size)
    ibmss[zone]['Po_Det_Scn2']=Po_obs*np.ones(tv_scn.size)
    ibmss[zone]['Po_Det_Scn3']=Po_obs*np.ones(tv_scn.size)

    ind=np.where( (tv_scn<yr_mod_start-50) )[0]; 
    ibmss[zone]['Po_Det_Scn1'][ind]=rr_pi*Po_obs
    ibmss[zone]['Po_Det_Scn2'][ind]=rr_pi*Po_obs
    ibmss[zone]['Po_Det_Scn3'][ind]=rr_pi*Po_obs

    ind=np.where( (tv_scn>=yr_mod_start-50) & (tv_scn<yr_mod_start) )[0];
    ibmss[zone]['Po_Det_Scn1'][ind]=np.interp(tv_scn[ind],np.array([yr_mod_start-50,yr_mod_start]),np.array([rr_pi*Po_obs,Po_obs]))
    ibmss[zone]['Po_Det_Scn2'][ind]=np.interp(tv_scn[ind],np.array([yr_mod_start-50,yr_mod_start]),np.array([rr_pi*Po_obs,Po_obs]))
    ibmss[zone]['Po_Det_Scn3'][ind]=np.interp(tv_scn[ind],np.array([yr_mod_start-50,yr_mod_start]),np.array([rr_pi*Po_obs,Po_obs]))

    ind=np.where( (tv_scn>2020) )[0]; 
    ibmss[zone]['Po_Det_Scn1'][ind]=np.interp(tv_scn[ind],np.array([2020,2200]),np.array([Po_obs,1.0*Po_obs]))
    ibmss[zone]['Po_Det_Scn2'][ind]=np.interp(tv_scn[ind],np.array([2020,2200]),np.array([Po_obs,1.5*Po_obs]))
    ibmss[zone]['Po_Det_Scn3'][ind]=np.interp(tv_scn[ind],np.array([2020,2200]),np.array([Po_obs,2.0*Po_obs]))

# Save
gu.opickle(PathData,ibmss)

# Plot deteriministic component of scenarios
zone='SBPB'
plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,5.5))
rc=patches.Rectangle((1951,0),2020-1951,10,facecolor=[0.85,0.85,0.85])
ax.add_patch(rc)
ax.plot(tv_scn,ibmss[zone]['Po_Det_Scn1'],'b-',linewidth=1.5,label='IBM occurrence Scn. 1')
ax.plot(tv_scn,ibmss[zone]['Po_Det_Scn2'],'g--',linewidth=1.5,label='IBM occurrence Scn. 2')
ax.plot(tv_scn,ibmss[zone]['Po_Det_Scn3'],'r-.',linewidth=1.5,label='IBM occurrence Scn. 3')
ax.annotate('Modern era',(1951+35,6.5),ha='center')
ax.legend(loc='upper left',frameon=False)
ax.set(position=[0.065,0.17,0.92,0.8],xlim=[1000-0.5,tv_scn[-1]+1+0.5],
       xticks=np.arange(tv_scn[0],tv_scn[-1],50),
       ylim=[0,10],
       ylabel='Annual probability (% yr$^-$$^1$)',
       xlabel='Time, calendar year')
gu.PrintFig(PathFigures + '\\IBM_Scenarios_ts_IBM_SBPS','png',900)

#%% Manual specification
# Fitted models not working very well - predict huge mortality. Provide alternative option.

for zone in ibmss.keys():
    ibmss[zone]['Beta_Pareto_Alt']=[1.6,0.00002,0.0015]

gu.opickle(PathData,ibmss)

#%% Apply random component to scenarios

n_stand=2000

for zone in ibmss.keys():

    # Initialize annual probability of occurrence (final with deterministic and
    # random components)
    ibmss[zone]['Po_Scn1']=np.zeros((tv_scn.size,n_stand))
    
    for iT in range(tv_scn.size):
        # Adjust shape parameter to match specified annual probability of 
        # occurrence from the deterministic component
        b0=ibmss[zone]['Beta_Pareto'].copy()  
        
        b_shape=ibmss[zone]['Pareto_shape_for_Po'].copy()
        b0[0]=np.exp(b_shape[1]*np.log(ibmss[zone]['Po_Det_Scn1'][iT])+b_shape[0])
        
        # Manual exploration of coefficients:
        b0=[1.6,0.00002,0.0015]
        
        ibmss[zone]['Po_Scn1'][iT,:]=gensm.GenerateDisturbancesFromPareto(1,n_stand,b0)

zone='SBPS'
#zone='ESSF'
Po=np.sum(ibmss[zone]['Po_Scn1'],axis=1)/n_stand*100
plt.close('all'); 
plt.bar(tv_scn[3000:],Po[3000:],1)
plt.plot(tv_scn[3000:],gu.movingave(Po[3000:],50,'historical'),linewidth=3)
print(np.mean(Po[3000:]))

#%% Histogram of burn severity rating in 2017

# Import data
u=np.array(['T','L','M','S','V'])
id=np.array([5,2,3,4,6])

# Get number of occurrences (exclude unknown class)
nams=np.array(list(ibmss.keys()))
p=np.zeros((u.size,nams.size))
cnt=0
for zone in ibmss.keys():
    N=np.zeros(u.size)
    for i in range(u.size):
        ind=np.where(ibmss[zone]['Oc']==id[i])[0]
        N[i]=ind.size
    p[:,cnt]=N/np.sum(N)
    cnt=cnt+1

# Plot histogram
ind=np.where(nams=='SBS')[0]

print(p[:,ind])

plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(12,5.5))
ax.bar(np.arange(1,6,1),p[:,ind].flatten()*100,0.9,facecolor=[0.8,0.8,0.8],label='Modern era (observed)')
#ax.bar(u+0.16,p_pi,0.3,facecolor=[0.25,0.5,1],label='Pre-industrial era')
ax.set(position=[0.1,0.1,0.84,0.86],xticks=np.arange(1,6,1),xticklabels=['Trace','Low','Medium','Severe','Very severe'],
       ylabel='Probability (%)')
#ax.legend(loc='upper right',frameon=False)
gu.PrintFig(PathFigures + '\\IBM_DistributionSBS','png',900)

