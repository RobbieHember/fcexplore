'''

WILDFIRE OCCURRENCE STATS AND SCENARIOS BY BGC ZONE

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
from fcgadgets.macgyver import utilities_general as gu
from fcgadgets.macgyver import utilities_gis as gis
from fcgadgets.macgyver import utilities_inventory as invu
from fcgadgets.taz import aspatial_stat_models as asm

#%% Set figure properties

plt.rcParams.update( gu.Import_GraphicsParameters('spyder_fs6') )

#%% Path of file to store stats and scenarios

PathData=r'C:\Users\rhember\Documents\Data\Taz Datasets\Wildfire Stats and Scenarios\Wildfire_Stats_Scenarios_By_BGCZ.pkl'

PathFigures=r'C:\Users\rhember\OneDrive - Government of BC\Figures\Wildfire\Wildfire_Stats_Sceanrios_By_BGCZ'

#%% Project assumptions

tv_obs=np.arange(1920,2020,1)
tv_scn=np.arange(-2000,2201,1)

#%% Import BGC zones

zBECZ=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\becz.tif')
zBECZf=zBECZ['Data'].flatten()
zBECZs=zBECZ['Data'][0::50,0::50].flatten()
dBECZ=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\becz_lut.xlsx')

#%% Calculate observed historical probability of occurrence from wildfire perimiter
# This takes > 30 min

flg=0
if flg==1:
    
    wfss={}
    for iZ in range(dBECZ['ZONE'].size):
        nam=dBECZ['ZONE'][iZ]
        wfss[nam]={}
        wfss[nam]['Ao']=np.zeros(tv_obs.size)
        ind=np.where( (zBECZs==dBECZ['VALUE'][iZ]) )[0]
        wfss[nam]['Oc']=np.zeros((tv_obs.size,ind.size))
        ind=np.where( (zBECZf==dBECZ['VALUE'][iZ]) )[0]
        wfss[nam]['Azone']=ind.size
    
    for iT in range(tv_obs.size):
        print(tv_obs[iT])
        zFire=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\PROT_HISTORICAL_FIRE_POLYS_SP_' + str(tv_obs[iT]) + '.tif')
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
    
    # Save 
    # gu.opickle(r'C:\Users\rhember\Documents\Data\Wildfire\Taz Calibration\wfss.pkl',wfss)

else:
    wfss=gu.ipickle(r'C:\Users\rhember\Documents\Data\Wildfire\Taz Calibration\wfss.pkl')
    

#%% Calculate annual probability of occurrence 

for zone in wfss.keys():    
    wfss[zone]['Po_obs']=np.sum(wfss[zone]['Oc'])/wfss[zone]['Oc'].size

#%% Bar chart of mean annual probability of occurrence

Po_obs_mu=np.zeros(len(wfss.keys())); cnt=0
for k in wfss.keys():
    Po_obs_mu[cnt]=wfss[k]['Po_obs']; cnt=cnt+1
    
plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(12,5.5))
ax.bar(np.arange(1,Po_obs_mu.size+1),Po_obs_mu*100,facecolor=[0.8,0.8,0.8])
ax.set(position=[0.085,0.1,0.9,0.86],xlim=[0.5,Po_obs_mu.size+.5],xticks=np.arange(1,Po_obs_mu.size+1),
       xticklabels=list(wfss.keys()),ylabel='Probability of occurrence (% yr$^-$$^1$)')
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')
ax.tick_params(length=2)
#gu.PrintFig(PathFigures + '\\Wilfire_Occurrence_ByBGCZ','png',900)

#%% Bar chart of total area burned

A=np.zeros(tv_obs.size); cnt=0
for k in wfss.keys():
    A=A+wfss[k]['Ao']
    
plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(12,5.5))
ax.bar(tv_obs,A/1000,facecolor=[0.8,0.8,0.8])
ax.set(position=[0.085,0.1,0.9,0.86],xlim=[1919,2021],xticks=np.arange(1920,2040,10),ylabel='Area burned (ha yr$^-$$^1$)')
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=1.5)
#gu.PrintFig(PathFigures + '\\Wilfire_Occurrence_ByBGCZ','png',900)

#%% Pareto distribution fits

for zone in wfss.keys():
    Ao=wfss[zone]['Ao']/wfss[zone]['Azone']
    shape,loc,scale=stats.pareto.fit(Ao)
    wfss[zone]['Beta_Pareto']=np.array([shape,loc,scale])

#%% Plot time series of annual time series of observed annual area burned (%/yr)

zone='SBS'
wfss[zone]['Beta_Pareto']

plt.close('all')
fig,ax=plt.subplots(2,1,figsize=gu.cm2inch(15,8))
ax[0].bar(tv_obs,wfss[zone]['Ao']/wfss['IDF']['Azone']*100,1,facecolor=[0.5,0,0])
ax[0].set(position=[0.06,0.58,0.9,0.42],xlim=[tv_obs[0]-0.5,tv_obs[-1]+1+0.5],xticks=np.arange(tv_obs[0],tv_obs[-1]+2,10),ylabel='Probability (% yr$^-$$^1$)')

shape=3; scale=0.0075; loc=-scale;
yhat=stats.pareto.rvs(shape,loc=loc,scale=scale,size=tv_obs.size)
ax[1].bar(tv_obs,yhat*100,1);
ax[1].set(position=[0.06,0.08,0.9,0.42],xlim=[tv_obs[0]-0.5,tv_obs[-1]+1+0.5],xticks=np.arange(tv_obs[0],tv_obs[-1]+2,10),ylabel='Probability (% yr$^-$$^1$)');

#gu.PrintFig(PathFigures + '\\Wilfire_ts_IDF','png',900)

# Plot histogram of annual area burned (%/yr)
n=2000
#Po=stats.pareto.rvs(3.8,-0.009,0.009,n)
Po=stats.pareto.rvs(3,-0.008,0.008,n)
plt.close('all'); 
plt.hist(Po*100,np.arange(0,3,0.1))
print(np.median(Po*100))
print(np.mean(Po*100))

#%% Plot time series of total area burned

A=np.zeros(tv_obs.size); 
for k in wfss.keys():
    A=A+wfss[k]['Ao']
    
plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,6))
ax.bar(tv_obs,A/1000,1,facecolor=[0.29,0.49,0.74])
ax.set(position=[0.08,0.1,0.88,0.84],xlim=[tv_obs[0]-0.5,tv_obs[-1]+1+0.5],xticks=np.arange(tv_obs[0],tv_obs[-1]+2,10),ylabel='Area affected by wildfire (ha yr$^-$$^1$)',xlabel='Time, years')
gu.PrintFig(PathFigures + '\\Wildfire_AreaAffected','png',900)


#%% Look at how mean annual probability of occurrence varies with shape and scale parameters

binShape=np.arange(2,4,0.05)
binScale=np.arange(0,0.02,0.002)
z=np.zeros((binShape.size,binScale.size))
for i in range(binShape.size):
    for j in range(binScale.size):
        Po=stats.pareto.rvs(binShape[i],-binScale[j],binScale[j],5000)
        z[i,j]=np.mean(Po)

plt.close('all'); fig,ax=plt.subplots(1,2)
im=ax[0].matshow(z*100,clim=[0,2],extent=[binScale[0],binScale[-1],binShape[0],binShape[-1]])
ax[0].set(position=[0.1,0.1,0.7,0.8],aspect='auto')
cb=plt.colorbar(im,cax=ax[1])
ax[1].set(position=[0.881,0.1,0.04,0.8])

#%% Generate calibrated Pareto parameters by BGC zone

# Specify shape parameter that will be used in all BGC zones
Shape=3.1

bin=np.arange(0.00005,0.015,0.00001)
for zone in wfss.keys():
    Po_hat_mu=np.zeros(bin.size)
    for iBin in range(bin.size):
        Po=stats.pareto.rvs(Shape,-bin[iBin],bin[iBin],2000)
        Po_hat_mu[iBin]=np.mean(Po)
    d=wfss[zone]['Po_obs']-Po_hat_mu
    iMin=np.where(np.abs(d)==np.min(np.abs(d)))[0]
    Po_hat_mu[iMin]
    Scale=bin[iMin[0]]
    wfss[zone]['Beta_Pareto_Cal']=np.array([Shape,-Scale,Scale])
    #print(Scale)

#%% Relationship between Po and scale paramter at a fixed shape

x=np.zeros(dBECZ['ZONE'].size)
y=np.zeros(dBECZ['ZONE'].size)
c=0
for k in wfss.keys():
    y[c]=wfss[k]['Beta_Pareto_Cal'][2]
    x[c]=wfss[k]['Po_obs']*100
    c=c+1

plt.close('all')
fig,ax=plt.subplots(1)
ax.plot(x,y,'o',linewidth=2)
ax.set(xscale='linear',yscale='linear',xlabel='Annual probability of occurrence (%/yr)',ylabel='Pareto scale parameter')
    
# Fit model of prob occ vs
x1=sm.tools.tools.add_constant(x)
y1=y
md=sm.OLS(y1,x1).fit()
#md.summary()
beta=md.params
xhat=np.linspace(np.min(x),np.max(x),2)
plt.plot(xhat,beta[1]*xhat+beta[0],'r--')

# Add model coefficients
for k in wfss.keys():
    wfss[k]['Pareto_scale_to_match_Po_mu']=beta
    
#%% Scenario development (deterministic component)

# Adjust historical period with variation from Marlon 2012
# Include climate-sensitive future

# Charcoal record digitized from Marlon et al. (2012)
# This must already be a moving average of the z-score as it does not exceed abs(1).
dChar=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\DigitizedFigures\Marlonetal2012_Fig2.xlsx')
tvChar=np.arange(500,2000,1)
yChar=np.interp(tvChar,dChar['Time'],dChar['Charcoal influx zscore'])
#plt.figure()
#plt.plot(tvChar,yChar,'.')

for zone in wfss.keys():
    
    # Observed Po
    Po_obs=wfss[zone]['Po_obs']*100
    
    # Indices
    iCharPre=np.where( (tvChar>=tv_scn[0]) & (tvChar<=1919) )[0]
    iCharCal=np.where( (tvChar>=1919) )[0]
    ind0=np.where( (tvChar>=tvChar[0]) & (tvChar<=1919) )[0]
    ind1=np.where( (tv_scn>=tvChar[0])  & (tv_scn<=1919) )[0]
    indF=np.where(tv_scn>=1920)[0]
    
    # Scenario 1
    ind=np.where( (tv_scn>=500)  & (tv_scn<=1919) )[0]
    #wfss[zone]['Po_Det_WF_Scn1']=np.mean(wfss[zone]['Po_Det_WF_Scn2'][ind])*np.ones(tv_scn.size)
    wfss[zone]['Po_Det_WF_Scn1']=1.0*Po_obs*np.ones(tv_scn.size)
    
    # Scenario 2
    yCharL=0.15*yChar
    d=Po_obs-np.mean(yCharL[iCharCal])
    yCharL=yCharL+d    
    y1=np.mean(yCharL)*np.ones(tv_scn.size)    
    y1[ind1]=yCharL[ind0]    
    y1[indF]=Po_obs
    y1=np.maximum(0,y1)
    
    yCharL=0.25*yChar
    d=Po_obs-np.mean(yCharL[iCharCal])
    yCharL=yCharL+d    
    y2=np.mean(yCharL)*np.ones(tv_scn.size)    
    y2[ind1]=yCharL[ind0]    
    y2[indF]=Po_obs
    y2=np.maximum(0,y2)
    
    yCharL=0.35*yChar
    d=Po_obs-np.mean(yCharL[iCharCal])
    yCharL=yCharL+d    
    y3=np.mean(yCharL)*np.ones(tv_scn.size)    
    y3[ind1]=yCharL[ind0]    
    y3[indF]=Po_obs
    y3=np.maximum(0,y3)
    
    yMin=np.min(np.column_stack((y1,y3)),axis=1)
    yMax=np.max(np.column_stack((y1,y3)),axis=1)
    yMu=y2
    
    flg=0
    if flg==1:
        plt.close('all')
        #plt.plot(tv_scn,y,'b-')
        plt.plot(tv_scn,yMin,'r--')
        plt.plot(tv_scn,yMu,'k-')
        plt.plot(tv_scn,yMax,'r--')
    
    wfss[zone]['Po_Det_WF_Scn2']=yMu
    wfss[zone]['Po_Det_WF_Scn2_Min']=yMin
    wfss[zone]['Po_Det_WF_Scn2_Max']=yMax
    
    # Assume 1.5 and 2.0 x increase by 2100 for Scn3 and 4
    ind=np.where( (tv_scn>2020) )[0]; 
    #wfss[zone]['Po_Det_WF_Scn2'][ind]=np.interp(tv_scn[ind],np.array([2020,2200]),np.array([Po_obs,1.0*Po_obs]))
    #wfss[zone]['Po_Det_WF_Scn3'][ind]=np.interp(tv_scn[ind],np.array([2020,2200]),np.array([Po_obs,2.125*Po_obs]))
    wfss[zone]['Po_Det_WF_Scn2'][ind]=np.interp(tv_scn[ind],np.array([2020,2200]),np.array([Po_obs,3.25*Po_obs]))

    # Scenario 3
    yCharL=0.25*yChar
    d=Po_obs-np.mean(yCharL[iCharCal])
    yCharL=yCharL+d    
    y1=np.mean(yCharL)*np.ones(tv_scn.size)    
    y1[ind1]=yCharL[ind0]    
    y1[indF]=Po_obs
    y1=np.maximum(0,y1)
    
    yCharL=0.5*yChar
    d=Po_obs-np.mean(yCharL[iCharCal])
    yCharL=yCharL+d    
    y2=np.mean(yCharL)*np.ones(tv_scn.size)    
    y2[ind1]=yCharL[ind0]    
    y2[indF]=Po_obs
    y2=np.maximum(0,y2)
    
    yCharL=0.75*yChar
    d=Po_obs-np.mean(yCharL[iCharCal])
    yCharL=yCharL+d    
    y3=np.mean(yCharL)*np.ones(tv_scn.size)    
    y3[ind1]=yCharL[ind0]    
    y3[indF]=Po_obs
    y3=np.maximum(0,y3)
    
    yMin=np.min(np.column_stack((y1,y3)),axis=1)
    yMax=np.max(np.column_stack((y1,y3)),axis=1)
    yMu=y2
    
    wfss[zone]['Po_Det_WF_Scn3']=yMu
    wfss[zone]['Po_Det_WF_Scn3_Min']=yMin
    wfss[zone]['Po_Det_WF_Scn3_Max']=yMax
    
    # Assume 1.5 and 2.0 x increase by 2100 for Scn3 and 4
    ind=np.where( (tv_scn>2020) )[0]; 
    #wfss[zone]['Po_Det_WF_Scn2'][ind]=np.interp(tv_scn[ind],np.array([2020,2200]),np.array([Po_obs,1.0*Po_obs]))
    #wfss[zone]['Po_Det_WF_Scn3'][ind]=np.interp(tv_scn[ind],np.array([2020,2200]),np.array([Po_obs,2.125*Po_obs]))
    wfss[zone]['Po_Det_WF_Scn3'][ind]=np.interp(tv_scn[ind],np.array([2020,2200]),np.array([Po_obs,3.25*Po_obs]))

    # Scenario 4
    wfss[zone]['Po_Det_WF_Scn4']=Po_obs*np.ones(tv_scn.size)
    ind=np.where( (tv_scn<1920-50) )[0]; 
    wfss[zone]['Po_Det_WF_Scn4'][ind]=1.5*Po_obs
    ind=np.where( (tv_scn>=1920-50) & (tv_scn<1920) )[0];
    wfss[zone]['Po_Det_WF_Scn4'][ind]=np.interp(tv_scn[ind],np.array([1920-50,1920]),np.array([1.5*Po_obs,Po_obs]))
    # Assume 1.5 and 2.0 x increase by 2100 for Scn3 and 4
    ind=np.where( (tv_scn>2020) )[0]; 
    wfss[zone]['Po_Det_WF_Scn4'][ind]=np.interp(tv_scn[ind],np.array([2020,2200]),np.array([Po_obs,3.25*Po_obs]))

    # Scenario 5
    wfss[zone]['Po_Det_WF_Scn5']=Po_obs*np.ones(tv_scn.size)
    
    ind=np.where( (tv_scn<1920-50) )[0]; 
    wfss[zone]['Po_Det_WF_Scn5'][ind]=2.5*Po_obs

    ind=np.where( (tv_scn>=1920-50) & (tv_scn<1920) )[0];
    wfss[zone]['Po_Det_WF_Scn5'][ind]=np.interp(tv_scn[ind],np.array([1920-50,1920]),np.array([2.5*Po_obs,Po_obs]))
    
    # Assume 1.5 and 2.0 x increase by 2100 for Scn3 and 4
    ind=np.where( (tv_scn>2020) )[0]; 
    wfss[zone]['Po_Det_WF_Scn5'][ind]=np.interp(tv_scn[ind],np.array([2020,2200]),np.array([Po_obs,3.25*Po_obs]))

#%% Plot deteriministic component of scenarios

zone='SBS'

itH=np.where(tv_scn<1920)[0]
itF=np.where(tv_scn>=1920)[0]

plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(12,5.5)); lw=1
rc=patches.Rectangle((1920,0),100,1.6,facecolor=[0.9,0.9,0.9])
ax.add_patch(rc)

ax.fill_between(tv_scn[itH],wfss[zone]['Po_Det_WF_Scn2_Min'][itH],wfss[zone]['Po_Det_WF_Scn2_Max'][itH],color=[0.9,0.95,1],lw=lw,zorder=1)
ax.fill_between(tv_scn[itH],wfss[zone]['Po_Det_WF_Scn3_Min'][itH],wfss[zone]['Po_Det_WF_Scn3_Max'][itH],color=[0.95,0.85,1],lw=lw,zorder=1)

ax.plot(tv_scn,wfss[zone]['Po_Det_WF_Scn1'],'b-',color=[1,0.5,0],lw=lw,label='Scenario 1')

ax.plot(tv_scn,wfss[zone]['Po_Det_WF_Scn2'],'b--',color=[0.5,0.65,0.8],lw=lw,label='Scenario 2')
#ax.plot(tv_scn[itH],wfss[zone]['Po_Det_WF_Scn2_Min'][itH],'b--',color=[0.5,0.65,0.8],lw=0.5)
#ax.plot(tv_scn[itH],wfss[zone]['Po_Det_WF_Scn2_Max'][itH],'b--',color=[0.5,0.65,0.8],lw=0.5)

ax.plot(tv_scn,wfss[zone]['Po_Det_WF_Scn3'],'b-.',color=[0.5,0.,1],lw=lw,label='Scenario 3')
#ax.plot(tv_scn[itH],wfss[zone]['Po_Det_WF_Scn3_Min'][itH],'b--',color=[0.5,0.,1],lw=0.5)
#ax.plot(tv_scn[itH],wfss[zone]['Po_Det_WF_Scn3_Max'][itH],'b--',color=[0.5,0.,1],lw=0.5)

ax.plot(tv_scn,wfss[zone]['Po_Det_WF_Scn4'],'b:',color=[0.4,0.9,0],lw=lw,label='Scenario 4')
ax.plot(tv_scn,wfss[zone]['Po_Det_WF_Scn5'],'b--',color=[0.5,0,0],lw=lw,label='Scenario 5')

ax.annotate('Modern\nperiod',(1920+50,0.75),ha='center')
ax.legend(loc='upper left',frameon=False)
ax.set(position=[0.09,0.17,0.88,0.8],ylim=[0,1.4],xlim=[1400-0.5,tv_scn[-1]+1+0.5],
       xticks=np.arange(tv_scn[0],tv_scn[-1]+100,100),ylabel='Probability of occurrence (% yr$^-$$^1$)',xlabel='Time, calendar year')
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')
ax.tick_params(length=2)

# gu.PrintFig(PathFigures + '\\Wildfire_Scenarios_ts_IDF_new','png',900)


#%% Save

# Remove occurrence data (no need to save)
flg=0
if flg==1:
    for zone in wfss.keys():
        del wfss[zone]['Oc']

# Save
gu.opickle(PathData,wfss)


#%% Scenario development (deterministic component) (OLD)

# Import parameters
# wfss=gu.ipickle(PathData)

#for zone in wfss.keys():
#
#    Po_obs=wfss[zone]['Po_obs']*100
#
#    wfss[zone]['Po_Det_WF_Scn1']=1.5*Po_obs*np.ones(tv_scn.size)
#    wfss[zone]['Po_Det_WF_Scn2']=Po_obs*np.ones(tv_scn.size)
#    wfss[zone]['Po_Det_WF_Scn3']=Po_obs*np.ones(tv_scn.size)
#    wfss[zone]['Po_Det_WF_Scn4']=Po_obs*np.ones(tv_scn.size)
#    
#    ind=np.where( (tv_scn<1920-50) )[0]; 
#    wfss[zone]['Po_Det_WF_Scn2'][ind]=1.5*Po_obs
#    wfss[zone]['Po_Det_WF_Scn3'][ind]=1.5*Po_obs
#    wfss[zone]['Po_Det_WF_Scn4'][ind]=1.5*Po_obs
#
#    ind=np.where( (tv_scn>=1920-50) & (tv_scn<1920) )[0];
#    wfss[zone]['Po_Det_WF_Scn2'][ind]=np.interp(tv_scn[ind],np.array([1920-50,1920]),np.array([1.5*Po_obs,Po_obs]))
#    wfss[zone]['Po_Det_WF_Scn3'][ind]=np.interp(tv_scn[ind],np.array([1920-50,1920]),np.array([1.5*Po_obs,Po_obs]))
#    wfss[zone]['Po_Det_WF_Scn4'][ind]=np.interp(tv_scn[ind],np.array([1920-50,1920]),np.array([1.5*Po_obs,Po_obs]))
#    
#    # Assume 1.5 and 2.0 x increase by 2100 for Scn3 and 4
#    ind=np.where( (tv_scn>2020) )[0]; 
#    wfss[zone]['Po_Det_WF_Scn2'][ind]=np.interp(tv_scn[ind],np.array([2020,2200]),np.array([Po_obs,1.0*Po_obs]))
#    wfss[zone]['Po_Det_WF_Scn3'][ind]=np.interp(tv_scn[ind],np.array([2020,2200]),np.array([Po_obs,2.125*Po_obs]))
#    wfss[zone]['Po_Det_WF_Scn4'][ind]=np.interp(tv_scn[ind],np.array([2020,2200]),np.array([Po_obs,3.25*Po_obs]))
#
## Check
#ind=np.where(tv_scn==2100)[0]
#ind1=np.where(tv_scn==2000)[0]
#print(wfss[zone]['Po_Det_WF_Scn4'][ind]/wfss[zone]['Po_Det_WF_Scn4'][ind1])




#%% Plot deterministic component of scenarios with other observations

def LookAtCariboo():
    
    # Cariboo reconstruction from dendro scars
    dTR=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\Wildfire\Cariboo Dendro Wildfire Occurrence\Harvey et al 2017 fig 2.xlsx')
    dTR['Year']=np.round(dTR['Year'])
    dTR['Site']=np.round(dTR['Site'])
    tvTR=np.arange(1750,2010)
    uTR=np.unique(dTR['Site'])
    vTR=np.zeros((tvTR.size,uTR.size))
    for i in range(uTR.size):
        ind=np.where(dTR['Site']==uTR[i])[0]
        ic,ia,ib=np.intersect1d(tvTR,dTR['Year'][ind],return_indices=True)
        vTR[ia,i]=1
    
    yTR=np.sum(vTR,axis=1)/vTR.shape[1]
    yTR_ma=gu.movingave(yTR,30,'Historical')
    muTR=np.mean(yTR_ma)
    sigTR=np.std(yTR_ma)
    yTR_ma_z=(yTR_ma-muTR)/sigTR
    
    plt.close('all')
    fig,ax=plt.subplots(1)
    ax.plot(tvTR,yTR_ma_z,lw=1)
    return



#%% Scenario development (deterministic+random component)
# Make sure it is working as epected

#n_stand=2000
#
#for zone in wfss.keys():
#
#    # Initialize annual probability of occurrence (final with deterministic and
#    # random components)
#    Po=np.zeros((tv_scn.size,n_stand))    
#    for iT in range(tv_scn.size):
#        # Adjust scale parameter to match specified annual probability of 
#        # occurrence from the deterministic component
#        b0=wfss[zone]['Beta_Pareto_Cal'].copy()
#        Scale=wfss[zone]['Pareto_scale_to_match_Po_mu'][1]*wfss[zone]['Po_Det_WF_Scn1'][iT]+wfss[zone]['Pareto_scale_to_match_Po_mu'][0]
#        b0[1]=-Scale
#        b0[2]=Scale
#        Po[iT,:]=gensm.GenerateDisturbancesFromPareto(1,n_stand,b0)    
#    # Convert to percent area of occurrence
#    wfss[zone]['PctOcc_DetPlusRand_WF_Scn1']=np.sum(Po,axis=1)/n_stand*100
#    
#    Po=np.zeros((tv_scn.size,n_stand))
#    for iT in range(tv_scn.size):
#        # Adjust shape parameter to match specified annual probability of 
#        # occurrence from the deterministic component
#        b0=wfss[zone]['Beta_Pareto_Cal'].copy()
#        Scale=wfss[zone]['Pareto_scale_to_match_Po_mu'][1]*wfss[zone]['Po_Det_WF_Scn4'][iT]+wfss[zone]['Pareto_scale_to_match_Po_mu'][0]
#        b0[1]=-Scale
#        b0[2]=Scale
#        Po[iT,:]=gensm.GenerateDisturbancesFromPareto(1,n_stand,b0)  
#    # Convert to percent area of occurrence
#    wfss[zone]['PctOcc_DetPlusRand_WF_Scn4']=np.sum(Po,axis=1)/n_stand*100
#
## Plot annual percent occurrence
#zone='IDF'
#t0=1000
#plt.close('all'); 
#plt.bar(tv_scn[t0:],wfss[zone]['PctOcc_DetPlusRand_WF_Scn1'][t0:],1)
#plt.plot(tv_scn[t0:],gu.movingave(wfss[zone]['PctOcc_DetPlusRand_WF_Scn1'][t0:],50,'historical'),linewidth=3)
#plt.bar(tv_scn[t0:],wfss[zone]['PctOcc_DetPlusRand_WF_Scn4'][t0:],1)
#plt.plot(tv_scn[t0:],gu.movingave(wfss[zone]['PctOcc_DetPlusRand_WF_Scn4'][t0:],50,'historical'),linewidth=3)






