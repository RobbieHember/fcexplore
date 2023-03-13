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
import statsmodels.formula.api as smf
import fcgadgets.macgyver.utilities_general as gu
import fcgadgets.macgyver.utilities_gis as gis
import fcgadgets.macgyver.utilities_inventory as invu
import fcgadgets.taz.aspatial_stat_models as asm

#%% Set figure properties

gp=gu.SetGraphics('Manuscript')

#%% Path of file to store stats and scenarios

PathData=r'C:\Users\rhember\Documents\Data\Taz Datasets\Wildfire Stats and Scenarios\Wildfire_Stats_Scenarios_By_BGCZ.pkl'

PathFigures=r'C:\Users\rhember\OneDrive - Government of BC\Figures\Wildfire\Wildfire_Stats_Sceanrios_By_BGCZ'

#%% Project assumptions

tv_obs=np.arange(1920,2020,1)
tv_scn=np.arange(-2000,2201,1)
tv_long=np.arange(1000,2022,1)

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
        wfss[nam]['Area Wildfire']=np.zeros(tv_obs.size)
        ind=np.where( (zBECZs==dBECZ['VALUE'][iZ]) )[0]
        wfss[nam]['Occurrence']=np.zeros( (tv_obs.size,ind.size) )
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
            wfss[nam]['Area Wildfire'][iT]=wfss[nam]['Area Wildfire'][iT]+ind.size
            ind=np.where( (zBECZs==dBECZ['VALUE'][iZ]) )[0]
            wfss[nam]['Occurrence'][iT,:]=zFires[ind]
        del zFire,zFiref,zFires
        garc.collect()

    # Save
    # gu.opickle(r'C:\Users\rhember\Documents\Data\Wildfire\Taz Calibration\wfss.pkl',wfss)

else:
    wfss=gu.ipickle(r'C:\Users\rhember\Documents\Data\Wildfire\Taz Calibration\wfss.pkl')


#%% Calculate annual probability of occurrence

for zone in wfss.keys():
    wfss[zone]['Po_obs']=np.sum(wfss[zone]['Occurrence'])/wfss[zone]['Occurrence'].size

#%% Calculate annual probability of occurrence with logistic regression
# Checked that we get the same answer as above. It is identical

flg=0
if flg==1:
    for zone in wfss.keys():
        d={}
        d['Occurrence']=wfss[zone]['Occurrence'].flatten()
        d['Random']=np.random.normal(loc=0.0,scale=1.0,size=d['Occurrence'].size)
        df=pd.DataFrame(d)
        log_reg=smf.logit("Occurrence ~ Random", data=df).fit()
        log_reg.summary()

        log_reg.params=log_reg.params+2*log_reg.bse
        d2={'Random':np.array([0.0])}
        wfss[zone]['Po_obs2']=log_reg.predict(d2)


    plt.close('all')
    for zone in wfss.keys():
        plt.plot(wfss[zone]['Po_obs'],wfss[zone]['Po_obs2'],'ko')

#%% Adjustments based on model evaluation
# *** Last run turned off. Not effective at improving accuracy of age (implying
# underestimation of disturbance during observation period.) ***
flg=0
if flg==1:
    zone='BWBS'; wfss[zone]['Po_obs']=1.3*wfss[zone]['Po_obs']
    zone='PP'; wfss[zone]['Po_obs']=1.24*wfss[zone]['Po_obs']
    zone='SBPS'; wfss[zone]['Po_obs']=1.18*wfss[zone]['Po_obs']

    zone='CMA'; wfss[zone]['Po_obs']=0.74*wfss[zone]['Po_obs']
    zone='MH'; wfss[zone]['Po_obs']=0.7*wfss[zone]['Po_obs']

#%% Bar chart of mean annual probability of occurrence

Po_obs_mu=np.zeros(len(wfss.keys())); cnt=0
for k in wfss.keys():
    Po_obs_mu[cnt]=wfss[k]['Po_obs']; cnt=cnt+1

plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(12,5))
ax.bar(np.arange(1,Po_obs_mu.size+1),Po_obs_mu*100,facecolor=[0.8,0.8,0.8])
ax.set(position=[0.085,0.125,0.9,0.85],xlim=[0.5,Po_obs_mu.size+.5],xticks=np.arange(1,Po_obs_mu.size+1),
       xticklabels=list(wfss.keys()),xlabel='BGC zone',ylabel='Probability of occurrence (% yr$^-$$^1$)')
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')
ax.tick_params(length=2)
#gu.PrintFig(PathFigures + '\\Wilfire_Occurrence_ByBGCZ','png',900)

#%% Bar chart of total area burned

A_Burned=np.zeros(tv_obs.size); cnt=0
for zone in wfss.keys():
    A_Burned=A_Burned+wfss[zone]['Area Wildfire']

plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(12,5.5))
ax.bar(tv_obs,A_Burned/1000,facecolor=[0.85,0,0])
ax.set(position=[0.09,0.12,0.895,0.86],xlim=[1919,2021],xticks=np.arange(1920,2040,10),ylabel='Area burned (ha yr$^-$$^1$)',xlabel='Time, years')
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=1.5)
#gu.PrintFig(PathFigures + '\\Wildfire_AreaAffected_Total','png',900)

#%% Pareto distribution fits

# The lower the shape parameter, the more extreme events
dS={}
dS['Shape']=np.zeros(16)
dS['Scale']=np.zeros(16)
dS['Mean Po (obs)']=np.zeros(16)
dS['Mean Po (fit)']=np.zeros(16)
dS['Mean Po Delta']=np.zeros(16)
cnt=0
for zone in wfss.keys():
    Ao=wfss[zone]['Area Wildfire']/wfss[zone]['Azone']
    shape,loc,scale=stats.pareto.fit(Ao)
    wfss[zone]['Beta_Pareto']=np.array([shape,loc,scale])
    #print(scale)
    dS['Shape'][cnt]=shape
    dS['Scale'][cnt]=scale
    dS['Mean Po (obs)'][cnt]=100*wfss[zone]['Po_obs']
    dS['Mean Po (fit)'][cnt]=100*np.mean(stats.pareto.rvs(shape,loc,scale,3000))
    dS['Mean Po Delta'][cnt]=dS['Mean Po (obs)'][cnt]-dS['Mean Po (fit)'][cnt]
    cnt=cnt+1

df=pd.DataFrame(dS)

#%% Look at how mean annual probability of occurrence varies with shape and scale parameters

binShape=np.arange(2,4,0.05)
binScale=np.arange(0,0.02,0.002)
z=np.zeros((binShape.size,binScale.size))
for i in range(binShape.size):
    for j in range(binScale.size):
        Po=stats.pareto.rvs(binShape[i],-binScale[j],binScale[j],5000)
        z[i,j]=np.mean(Po)

plt.close('all'); fig,ax=plt.subplots(1,2)
im=ax[0].matshow(z*100,clim=[0.2,0.4],extent=[binScale[0],binScale[-1],binShape[0],binShape[-1]])
ax[0].set(position=[0.1,0.1,0.7,0.8],aspect='auto',yticks=binShape,xticks=binScale,xlabel='Scale',ylabel='Shape')
cb=plt.colorbar(im,cax=ax[1])
ax[1].set(position=[0.881,0.1,0.04,0.8])

#%% Generate calibrated Pareto parameters by BGC zone

# Specify shape parameter that will be used in all BGC zones
#Shape0=2.8
Shape0=3.2

bin=np.arange(0.00005,0.015,0.00001)
for zone in wfss.keys():
    Po_hat_mu=np.zeros(bin.size)
    for iBin in range(bin.size):
        Po=stats.pareto.rvs(Shape0,-bin[iBin],bin[iBin],2000)
        Po_hat_mu[iBin]=np.mean(Po)
    d=wfss[zone]['Po_obs']-Po_hat_mu
    iMin=np.where(np.abs(d)==np.min(np.abs(d)))[0]
    Po_hat_mu[iMin]
    Scale=bin[iMin[0]]
    wfss[zone]['Beta_Pareto_Cal']=np.array([Shape0,-Scale,Scale])
    #print(Scale)

#%% Plot time series of annual time series of observed annual area burned (%/yr)

zone='SBS'

print(wfss[zone]['Beta_Pareto'])
#shape=Shape0; loc=-0.0085; scale=-loc;
#shape=wfss[zone]['Beta_Pareto'][0]; scale=-1*wfss[zone]['Beta_Pareto'][1]; loc=-wfss[zone]['Beta_Pareto'][2]
shape=wfss[zone]['Beta_Pareto_Cal'][0]; loc=wfss[zone]['Beta_Pareto_Cal'][1]; scale=wfss[zone]['Beta_Pareto_Cal'][2]

plt.close('all')
fig,ax=plt.subplots(3,1,figsize=gu.cm2inch(15,8))
ax[0].bar(tv_obs,wfss[zone]['Area Wildfire']/wfss['IDF']['Azone']*100,1,facecolor=[0.5,0,0]);
yhat=stats.pareto.rvs(shape,loc=loc,scale=scale,size=tv_long.size);
it0=np.where(tv_long<1920)[0]
ax[0].bar(tv_long[it0],yhat[it0]*100,1,fc=[1,0.6,0]);
ax[0].set(position=[0.06,0.7,0.93,0.28],xlim=[tv_long[0]-0.5,tv_long[-1]+1+0.5],xticks=np.arange(tv_long[0],tv_obs[-1]+2,100),xticklabels='',ylim=[0,12],ylabel='Area burned (% yr$^-$$^1$)',xlabel='',yticks=np.arange(0,12,2));
ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].tick_params(length=1.5)

yhat=stats.pareto.rvs(shape,loc=loc,scale=scale,size=tv_long.size);
ax[1].bar(tv_long,yhat*100,1,fc=[1,0.6,0]);
ax[1].set(position=[0.06,0.39,0.93,0.28],xlim=[tv_long[0]-0.5,tv_long[-1]+1+0.5],xticks=np.arange(tv_long[0],tv_obs[-1]+2,100),xticklabels='',ylim=[0,12],ylabel='Area burned (% yr$^-$$^1$)',xlabel='',yticks=np.arange(0,12,2));
ax[1].yaxis.set_ticks_position('both'); ax[1].xaxis.set_ticks_position('both'); ax[1].tick_params(length=1.5)

yhat=stats.pareto.rvs(shape,loc=loc,scale=scale,size=tv_long.size);
ax[2].bar(tv_long,yhat*100,1,fc=[1,0.6,0]);
ax[2].set(position=[0.06,0.085,0.93,0.28],xlim=[tv_long[0]-0.5,tv_long[-1]+1+0.5],xticks=np.arange(tv_long[0],tv_obs[-1]+2,100),ylim=[0,12],ylabel='Area burned (% yr$^-$$^1$)',xlabel='Time, years',yticks=np.arange(0,12,2));
ax[2].yaxis.set_ticks_position('both'); ax[2].xaxis.set_ticks_position('both'); ax[2].tick_params(length=1.5)
gu.axletters(ax,plt,0.01,0.85,FontColor=[0,0,0],LetterStyle='Caps',FontWeight='Bold')
# gu.PrintFig(PathFigures + '\\AreaBurnedSimulationTimeSeriesLong_' + zone,'png',900)

#%% Plot histogram of annual area burned (%/yr)

n=10000
#Po=stats.pareto.rvs(3.8,-0.009,0.009,n)
Po=stats.pareto.rvs(shape,loc,scale,n)
plt.close('all');
plt.hist(Po*100,np.arange(0,10,0.2))
print(np.median(Po*100))
print(np.mean(Po*100))

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
    indF2=np.where(tv_scn>2020)[0]

    #--------------------------------------------------------------------------
    # Scenario 1 - constant Po from observations
    #--------------------------------------------------------------------------

    ind=np.where( (tv_scn>=500)  & (tv_scn<=1919) )[0]
    wfss[zone]['Po_Det_WF_Scn1']=1.0*Po_obs*np.ones(tv_scn.size)

    #--------------------------------------------------------------------------
    # Scenario 2 - history informed by fire scars, 2 x observation mean by 2100
    #--------------------------------------------------------------------------

    # Historical

    yCharL=0.15*yChar
    d=Po_obs-np.mean(yCharL[iCharCal])
    yCharL=yCharL+d
    y1=np.mean(yCharL)*np.ones(tv_scn.size)
    y1[ind1]=yCharL[ind0]
    y1[indF]=Po_obs
    y1=np.maximum(0,y1)

    yCharL=0.4*yChar
    d=Po_obs-np.mean(yCharL[iCharCal])
    yCharL=yCharL+d
    y2=np.mean(yCharL)*np.ones(tv_scn.size)
    y2[ind1]=yCharL[ind0]
    y2[indF]=Po_obs
    y2=np.maximum(0,y2)

    yCharL=0.45*yChar
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
    #wfss[zone]['Po_Det_WF_Scn2_Min']=yMin
    #wfss[zone]['Po_Det_WF_Scn2_Max']=yMax

    # Future

    flg=0
    if flg==1:
        # Figure out what factor to use to get a specific increase by 2100
        np.interp(2100,np.array([2020,2200]),np.array([Po_obs,25.75*Po_obs]))/Po_obs

    # Assume 2.0 x increase by 2100
    wfss[zone]['Po_Det_WF_Scn2'][indF2]=np.interp(tv_scn[indF2],np.array([2020,2200]),np.array([Po_obs,3.25*Po_obs]))

    #--------------------------------------------------------------------------
    # Scenario 3 - history informed by fire scars, 4 x observation mean by 2100
    #--------------------------------------------------------------------------

    wfss[zone]['Po_Det_WF_Scn3']=yMu.copy()
    wfss[zone]['Po_Det_WF_Scn3'][indF2]=np.interp(tv_scn[indF2],np.array([2020,2200]),np.array([Po_obs,7.75*Po_obs]))

    #--------------------------------------------------------------------------
    # Scenario 4 - history informed by fire scars, 4 x observation mean by 2100
    #--------------------------------------------------------------------------

    wfss[zone]['Po_Det_WF_Scn4']=yMu.copy()
    wfss[zone]['Po_Det_WF_Scn4'][indF2]=np.interp(tv_scn[indF2],np.array([2020,2200]),np.array([Po_obs,16.75*Po_obs]))

    #--------------------------------------------------------------------------
    # Scenario 5 - history informed by fire scars, 4 x observation mean by 2100
    #--------------------------------------------------------------------------

    wfss[zone]['Po_Det_WF_Scn5']=yMu.copy()
    wfss[zone]['Po_Det_WF_Scn5'][indF2]=np.interp(tv_scn[indF2],np.array([2020,2200]),np.array([Po_obs,25.75*Po_obs]))


#%% Plot deteriministic component of scenarios

zone='SBS'

itH=np.where(tv_scn<1920)[0]
itF=np.where(tv_scn>=1920)[0]

plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(12.5,6.5)); lw=1
rc=patches.Rectangle((1920,0),100,20,facecolor=[0.92,0.92,0.92])
ax.add_patch(rc)

#ax.fill_between(tv_scn[itH],wfss[zone]['Po_Det_WF_Scn2_Min'][itH],wfss[zone]['Po_Det_WF_Scn2_Max'][itH],color=[0.9,0.95,1],lw=lw,zorder=1)
#ax.fill_between(tv_scn[itH],wfss[zone]['Po_Det_WF_Scn3_Min'][itH],wfss[zone]['Po_Det_WF_Scn3_Max'][itH],color=[0.95,0.85,1],lw=lw,zorder=1)

ax.plot(tv_scn,wfss[zone]['Po_Det_WF_Scn1'],'b-',color=[0.27,0.49,0.77],lw=lw,label='Scenario 1: Constant calibrated against observation\nperiod')

ax.plot(tv_scn,wfss[zone]['Po_Det_WF_Scn2'],'b--',color=[1,0.5,0],lw=lw,label='Scenario 2: Pre-observation period calibrated against firescars,\n 2 x observed level by 2100')
#ax.plot(tv_scn[itH],wfss[zone]['Po_Det_WF_Scn2_Min'][itH],'b--',color=[0.5,0.65,0.8],lw=0.5)
#ax.plot(tv_scn[itH],wfss[zone]['Po_Det_WF_Scn2_Max'][itH],'b--',color=[0.5,0.65,0.8],lw=0.5)

ax.plot(tv_scn[indF2],wfss[zone]['Po_Det_WF_Scn3'][indF2],'b-.',color=[1,0,0],lw=lw,label='Scenario 3: Pre-observation period calibrated against firescars,\n 4 x observed level by 2100')
ax.plot(tv_scn[indF2],wfss[zone]['Po_Det_WF_Scn4'][indF2],'b:',color=[0.8,0,0],lw=lw,label='Scenario 4: Pre-observation period calibrated against firescars,\n 8 x observed level by 2100')
ax.plot(tv_scn[indF2],wfss[zone]['Po_Det_WF_Scn5'][indF2],'b--',color=[0.4,0,0],lw=lw,label='Scenario 5: Pre-observation period calibrated against firescars,\n 12 x observed level by 2100')

#ax.annotate('2 x observ. period\nmean in 2100',(2100,0.95),ha='right',color=[1,0,0])
#ax.annotate('4 x observ. period\nmean in 2100',(2100,1.95),ha='right',color=[0.8,0,0])
#ax.annotate('8 x observ. period\nmean in 2100',(2100,3.85),ha='right',color=[0.4,0,0])

ax.annotate('Observation\nperiod',(1920+50,2.0),ha='center',fontweight='bold')
ax.legend(loc='upper left',frameon=False)
ax.set(position=[0.09,0.14,0.88,0.84],ylim=[0,5],xlim=[1400-0.5,tv_scn[-1]+1+0.5],
       xticks=np.arange(tv_scn[0],tv_scn[-1]+100,100),ylabel='Probability of occurrence (% yr$^-$$^1$)',xlabel='Time, calendar year')
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')
ax.tick_params(length=2)

# gu.PrintFig(PathFigures + '\\Wildfire_Scenarios_ts_' + zone,'png',900)


#%% Save

# Remove occurrence data (no need to save)
flg=0
if flg==1:
    for zone in wfss.keys():
        del wfss[zone]['Occurrence']

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






