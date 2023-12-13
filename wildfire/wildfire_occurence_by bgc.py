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
from matplotlib.ticker import ScalarFormatter
import scipy.stats as stats
import statsmodels.api as sm
import statsmodels.formula.api as smf
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.macgyver.util_inventory as uinv
import fcgadgets.bc1ha.bc1ha_util as u1ha
import fcgadgets.taz.aspatial_stat_models as asm
gp=gu.SetGraphics('Manuscript')

#%% Path of file to store stats and scenarios
meta=u1ha.Init()
PathData=meta['Paths']['Model']['Taz Datasets'] + '\\Wildfire Stats and Scenarios\Wildfire_Stats_Scenarios_By_BGCZ.pkl'
PathFigures=r'C:\Users\rhember\OneDrive - Government of BC\Figures\Wildfire\Wildfire_Stats_Sceanrios_By_BGCZ'

#%% Project assumptions
tv_obs=np.arange(1920,2024,1)
tv_scn=np.arange(-2000,2201,1)
tv_long=np.arange(1000,2024,1)
ivl=50

#%% Import raster data
z=u1ha.Import_Raster(meta,[],['bgcz','lc_comp1_2019'])
zZ_Full=z['bgcz']['Data'].flatten()
zLC_Full=z['lc_comp1_2019']['Data'].flatten()
zZ_Sub=z['bgcz']['Data'][0::ivl,0::ivl].flatten()
zLC_Sub=z['lc_comp1_2019']['Data'][0::ivl,0::ivl].flatten()

lut_bgcz=meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']
bgcz_cd=np.array(list(meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'].keys()))
lut_lc=meta['LUT']['Derived']['lc_comp1']

#%% Calculate observed historical probability of occurrence from wildfire perimiter
# This takes 1 min
wfss={}
for iZ in range(bgcz_cd.size):
    nam=bgcz_cd[iZ]
    wfss[nam]={}
    wfss[nam]['Area Wildfire']=np.zeros(tv_obs.size)
    ind=np.where( (zZ_Full==lut_bgcz[nam]) & (zLC_Full==lut_lc['Forest']) )[0]
    wfss[nam]['Area Zone']=ind.size
    ind=np.where( (zZ_Sub==lut_bgcz[nam]) & (zLC_Sub==lut_lc['Forest']) )[0]
    wfss[nam]['Occurrence']=np.zeros((tv_obs.size,ind.size),dtype='int8')

for iY in range(6):
    zFire=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PROT_HISTORICAL_FIRE_POLYS_SP\FIRE_YEAR_' + str(iY+1) + '_Year.tif')['Data']
    idxF=gu.IndicesFromUniqueArrayValues(zFire.flatten())        
    for iZ in range(bgcz_cd.size):
        nam=bgcz_cd[iZ]            
        for yr in idxF.keys():
            if yr==0: continue
            iT=np.where(tv_obs==yr)[0]
            if iT.size==0: continue
            ind=np.where( (zZ_Full[idxF[yr]]==lut_bgcz[nam]) & (zLC_Full[idxF[yr]]==lut_lc['Forest']) )[0]
            wfss[nam]['Area Wildfire'][iT]=wfss[nam]['Area Wildfire'][iT]+ind.size            
        iZoneS=np.where( (zZ_Sub==lut_bgcz[nam]) & (zLC_Sub==lut_lc['Forest']) )[0]
        idxS=gu.IndicesFromUniqueArrayValues(zFire[0::ivl,0::ivl].flatten()[iZoneS])
        for yr in idxS.keys():
            if yr==0: continue
            iT=np.where(tv_obs==yr)[0]
            if iT.size==0: continue                
            wfss[nam]['Occurrence'][iT,idxS[yr]]=1
    del zFire,idxF,idxS
    garc.collect()

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
# *** FOR EXPLORITORY PURPOSES ONLY ***
# *** Last run turned off. Not effective at improving accuracy of age (implying
# underestimation of disturbance during observation period.) ***
flg=0
if flg==1:
    zone='BWBS'; wfss[zone]['Po_obs']=1.3*wfss[zone]['Po_obs']
    zone='PP'; wfss[zone]['Po_obs']=1.24*wfss[zone]['Po_obs']
    zone='SBPS'; wfss[zone]['Po_obs']=1.18*wfss[zone]['Po_obs']

    zone='CMA'; wfss[zone]['Po_obs']=0.74*wfss[zone]['Po_obs']
    zone='MH'; wfss[zone]['Po_obs']=0.7*wfss[zone]['Po_obs']

#%% Calculate total area burned by year
A_Burned=np.zeros(tv_obs.size); cnt=0
for zone in wfss.keys():
    A_Burned=A_Burned+wfss[zone]['Area Wildfire']

# plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(12,5.5))
# ax.bar(tv_obs,A_Burned/1e6,facecolor=[0.75,0.65,0.55])
# ax.set(xticks=np.arange(1920,2040,10),ylabel='Area affected (Million ha yr$^-$$^1$)',xlabel='Time, years',xlim=[1919,2024])
# ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=1.5)
# #gu.PrintFig(PathFigures + '\\Wildfire_AreaAffected_Total','png',900)

#%% Calculate mean annual probability of occurrence by BGC zone
Po_obs_mu=np.zeros(len(wfss.keys())); cnt=0
for k in wfss.keys():
    Po_obs_mu[cnt]=wfss[k]['Po_obs']; cnt=cnt+1

ord=np.flip(np.argsort(Po_obs_mu))
bgcz_cdS=bgcz_cd.copy()[ord]
Po_obs_mu=Po_obs_mu[ord]

# plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(12,5))
# ax.bar(np.arange(1,Po_obs_mu.size+1),Po_obs_mu*100,facecolor=[0.7,0.8,0.95])
# ax.set(xticks=np.arange(1,Po_obs_mu.size+1),
#        xticklabels=bgcz_cdS,xlabel='BGC zone',ylabel='Probability of occurrence (% yr$^-$$^1$)',xlim=[0.5,Po_obs_mu.size+.5],ylim=[0,0.5])
# ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=1.5)
# plt.tight_layout()
# #gu.PrintFig(PathFigures + '\\Wilfire_Occurrence_ByBGCZ','png',900)

#%% Plot Time series and by BGC zone combined
plt.close('all'); fig,ax=plt.subplots(2,1,figsize=gu.cm2inch(14,8))
ax[0].bar(tv_obs,A_Burned/1e6,facecolor=[0.7,0.8,0.95])
ax[0].set(xticks=np.arange(1920,2040,10),ylabel='Area affected (Million ha yr$^-$$^1$)',xlabel='Time, years',xlim=[1919,2024])
ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].tick_params(length=1.5)
ax[1].bar(np.arange(1,Po_obs_mu.size+1),Po_obs_mu*100,facecolor=[0.7,0.8,0.95])
ax[1].set(xticks=np.arange(1,Po_obs_mu.size+1),
       xticklabels=bgcz_cdS,xlabel='BGC zone',ylabel='Probability of occurrence (% yr$^-$$^1$)',xlim=[0.5,Po_obs_mu.size+.5],ylim=[0,0.5])
ax[1].yaxis.set_ticks_position('both'); ax[1].xaxis.set_ticks_position('both'); ax[1].tick_params(length=1.5)
plt.tight_layout()
gu.axletters(ax,plt,0.0175,0.87,FontColor=[0,0,0],LetterStyle='Caps',FontWeight='Bold')
gu.PrintFig(PathFigures + '\\Wildfire_ObservedStats','png',900)

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
    Ao=wfss[zone]['Area Wildfire']/wfss[zone]['Area Zone']
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

it0=np.where(tv_long<1920)[0]

plt.close('all'); fig,ax=plt.subplots(5,1,figsize=gu.cm2inch(15,9))
ax[0].bar(tv_obs,wfss[zone]['Area Wildfire']/wfss[zone]['Area Zone']*100,1,facecolor=[0.5,0,0]);
yhat=stats.pareto.rvs(shape,loc=loc,scale=scale,size=tv_long.size);
ax[0].bar(tv_long[it0],yhat[it0]*100,1,fc=[1,0.6,0]);
ax[0].set(xticks=np.arange(tv_long[0],tv_obs[-1]+2,100),xticklabels='',ylabel='',xlabel='',yticks=np.arange(0,12,2),xlim=[tv_long[0]-0.5,tv_long[-1]+1+0.5],ylim=[0,12]);
ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].tick_params(length=1.5)

ax[1].bar(tv_obs,wfss[zone]['Area Wildfire']/wfss[zone]['Area Zone']*100,1,facecolor=[0.5,0,0]);
yhat=stats.pareto.rvs(shape,loc=loc,scale=scale,size=tv_long.size);
ax[1].bar(tv_long[it0],yhat[it0]*100,1,fc=[1,0.6,0]);
ax[1].set(xticks=np.arange(tv_long[0],tv_obs[-1]+2,100),xticklabels='',ylabel='',xlabel='',yticks=np.arange(0,12,2),xlim=[tv_long[0]-0.5,tv_long[-1]+1+0.5],ylim=[0,12]);
ax[1].yaxis.set_ticks_position('both'); ax[1].xaxis.set_ticks_position('both'); ax[1].tick_params(length=1.5)

ax[2].bar(tv_obs,wfss[zone]['Area Wildfire']/wfss[zone]['Area Zone']*100,1,facecolor=[0.5,0,0]);
yhat=stats.pareto.rvs(shape,loc=loc,scale=scale,size=tv_long.size);
ax[2].bar(tv_long[it0],yhat[it0]*100,1,fc=[1,0.6,0]);
ax[2].set(xticks=np.arange(tv_long[0],tv_obs[-1]+2,100),xticklabels='',ylabel='Area burned (% yr$^-$$^1$)',xlabel='',yticks=np.arange(0,12,2),xlim=[tv_long[0]-0.5,tv_long[-1]+1+0.5],ylim=[0,12])
ax[2].yaxis.set_ticks_position('both'); ax[2].xaxis.set_ticks_position('both'); ax[2].tick_params(length=1.5)

ax[3].bar(tv_obs,wfss[zone]['Area Wildfire']/wfss[zone]['Area Zone']*100,1,facecolor=[0.5,0,0]);
yhat=stats.pareto.rvs(shape,loc=loc,scale=scale,size=tv_long.size);
ax[3].bar(tv_long[it0],yhat[it0]*100,1,fc=[1,0.6,0]);
ax[3].set(xticks=np.arange(tv_long[0],tv_obs[-1]+2,100),xticklabels='',ylabel='',xlabel='',yticks=np.arange(0,12,2),xlim=[tv_long[0]-0.5,tv_long[-1]+1+0.5],ylim=[0,12]);
ax[3].yaxis.set_ticks_position('both'); ax[3].xaxis.set_ticks_position('both'); ax[3].tick_params(length=1.5)

ax[4].bar(tv_obs,wfss[zone]['Area Wildfire']/wfss[zone]['Area Zone']*100,1,facecolor=[0.5,0,0]);
yhat=stats.pareto.rvs(shape,loc=loc,scale=scale,size=tv_long.size);
ax[4].bar(tv_long[it0],yhat[it0]*100,1,fc=[1,0.6,0]);
ax[4].set(xticks=np.arange(tv_long[0],tv_obs[-1]+2,100),ylabel='',xlabel='Time, years',yticks=np.arange(0,12,2),xlim=[tv_long[0]-0.5,tv_long[-1]+1+0.5],ylim=[0,12]);
ax[4].yaxis.set_ticks_position('both'); ax[4].xaxis.set_ticks_position('both'); ax[4].tick_params(length=1.5)
for i in range(5):
    ax[i].set(ylim=[0,6])
plt.tight_layout()
gu.axletters(ax,plt,0.015,0.74,FontColor=[0,0,0],LetterStyle='Caps',FontWeight='Bold')

gu.PrintFig(PathFigures + '\\AreaBurnedSimulationTimeSeriesLong_' + zone,'png',900)

#%% Plot histogram of annual area burned (%/yr)

n=10000
#Po=stats.pareto.rvs(3.8,-0.009,0.009,n)
Po=stats.pareto.rvs(shape,loc,scale,n)
plt.close('all');
plt.hist(Po*100,np.arange(0,10,0.2))
print(np.median(Po*100))
print(np.mean(Po*100))

#%% Relationship between Po and scale paramter at a fixed shape

x=np.zeros(bgcz_cd.size)
y=np.zeros(bgcz_cd.size)
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

#plt.figure(); plt.plot(tvChar,yChar,'-')

for zone in wfss.keys():

    # Observed Po
    Po_obs=wfss[zone]['Po_obs']*100

    # Indices
    iCharPre=np.where( (tvChar>=tv_scn[0]) & (tvChar<=1919) )[0]
    iCharCalPre=np.where( (tvChar>=500)  & (tvChar<=1919) )[0]
    iCharCal=np.where( (tvChar>=1919) )[0]
    ind0=np.where( (tvChar>=tvChar[0]) & (tvChar<=1919) )[0]
    ind1=np.where( (tv_scn>=tvChar[0])  & (tv_scn<=1919) )[0]
    indF=np.where(tv_scn>=1920)[0]
    indF2=np.where(tv_scn>2020)[0]

    # Initialize
    dF={'F0':{},'F1':{},'F2':{},'F3':{}}
    wfss[zone]['Po Det Scenarios']={'H0':dF.copy(),'H1':dF.copy(),'H2':dF.copy(),'H3':dF.copy(),'H4':dF.copy()}

    #--------------------------------------------------------------------------
    # Constant Po from observations
    #--------------------------------------------------------------------------
    
    wfss[zone]['Po Det Scenarios']['H0']['F0']=1.0*Po_obs*np.ones(tv_scn.size)

    dH={}
    dH['H0']=wfss[zone]['Po Det Scenarios']['H0']['F0'].copy()

    #--------------------------------------------------------------------------
    # Historical period informed by fire scars
    #--------------------------------------------------------------------------

    yChar0=yChar.copy()+100
    mu0=np.mean(yChar0[iCharCal])
    mu1=np.mean(yChar0[iCharCalPre])
    yChar0=(yChar0-mu0)/(mu1-mu0)
    #plt.close('all'); plt.plot(yChar0)
    
    ind0=np.where( (tv_scn<tvChar[0]) )[0]
    ind1=np.where( (tv_scn>=tvChar[0])  & (tv_scn<=1919) )[0]
    ind2=np.where( (tvChar>=tvChar[0])  & (tvChar<=1919) )[0]
    ind3=np.where( (tv_scn>1919) )[0]
                  
    dH['H1']=np.zeros(tv_scn.size)
    dH['H1'][ind1]=Po_obs+Po_obs*0.5*yChar0[ind2]
    dH['H1'][ind0]=np.mean(dH['H1'][ind1])
    dH['H1'][ind3]=Po_obs
    #plt.close('all'); plt.plot(dH['H1'])
    
    dH['H2']=np.zeros(tv_scn.size)
    dH['H2'][ind1]=Po_obs+Po_obs*1.0*yChar0[ind2]
    dH['H2'][ind0]=np.mean(dH['H2'][ind1])
    dH['H2'][ind3]=Po_obs
    #plt.close('all'); plt.plot(dH['H2'])
    
    dH['H3']=np.zeros(tv_scn.size)
    dH['H3'][ind1]=Po_obs+Po_obs*1.5*yChar0[ind2]
    dH['H3'][ind0]=np.mean(dH['H3'][ind1])
    dH['H3'][ind3]=Po_obs
    #plt.close('all'); plt.plot(dH['H3'])
    
    dH['H4']=dH['H2'].copy()
    ind4=np.where(tv_scn>1888)[0]
    dH['H4'][ind4]=np.mean(dH['H2'][ind1])
    wfss[zone]['Po Det Scenarios']['H4']['F0']=dH['H4']
    
    #mu=[np.mean(dH['H0'][ind]),np.mean(dH['H1'][ind]),np.mean(dH['H2'][ind]),np.mean(dH['H3'][ind])]
    #print(mu/mu[0]*100)
    
    #dH['H0']=dH['H0']*(np.mean(dH['H2'][ind])/np.mean(dH['H0'][ind]))

    flg=0
    if flg==1:
        plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(14.5,6.5)); lw=1
        rc=patches.Rectangle((1920,0),103,20,facecolor=[0.92,0.92,0.92])
        ax.add_patch(rc)
        iT1=np.where(tv_scn<1920)[0]; iT2=np.where( (tv_scn>=1920) & (tv_scn<=2024) )[0]
        mu1=np.mean(dH['H1'][iT1]); mu2=np.mean(dH['H1'][iT2]); pLo=(mu2-mu1)/mu1*100
        mu1=np.mean(dH['H2'][iT1]); mu2=np.mean(dH['H2'][iT2]); pMu=(mu2-mu1)/mu1*100
        mu1=np.mean(dH['H3'][iT1]); mu2=np.mean(dH['H3'][iT2]); pHi=(mu2-mu1)/mu1*100
        print(pLo)
        print(pMu)
        print(pHi)
        plt.plot(tv_scn,Po_obs*np.ones(tv_scn.size),'k-',lw=0.5,label='H0: Constant observed')
        plt.plot(tv_scn,dH['H1'],'b-',lw=lw,color=[0.27,0.49,0.77],label='H1: PI period = 1.5 x modern period')
        plt.plot(tv_scn,dH['H2'],'g--',lw=lw,color=[1,0.45,0],label='H2: PI period = 2.0 x modern period')
        plt.plot(tv_scn,dH['H3'],'r-.',lw=lw,color=[0.85,0,0],label='H3: PI period = 2.5 x modern period')
        plt.plot(tv_scn,dH['H4'],'r:',lw=lw,color=[0.5,0,0.95],label='H4: Constant H2')
        ax.annotate('Observation\nperiod',(1920+50,0.8),ha='center',fontweight='bold',color=[0.5,0.5,0.5])
        ax.legend(loc='upper left',frameon=False)
        ax.set(xticks=np.arange(tv_scn[0],tv_scn[-1]+100,100),ylabel='Probability of occurrence (% yr$^-$$^1$)',xlabel='Time, calendar year',
               yticks=np.arange(0,3,0.2),ylim=[0,1.2],xlim=[500-0.5,2050+0.5])        
        ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=1.5)
        plt.tight_layout()
        gu.PrintFig(PathFigures + '\\FirescarCalibration_' + zone,'png',900)

    for k in wfss[zone]['Po Det Scenarios']['H1'].keys():
        if k=='F0': continue
        wfss[zone]['Po Det Scenarios']['H1'][k]=dH['H1'].copy()
        wfss[zone]['Po Det Scenarios']['H2'][k]=dH['H2'].copy()
        wfss[zone]['Po Det Scenarios']['H3'][k]=dH['H3'].copy()
        wfss[zone]['Po Det Scenarios']['H4'][k]=dH['H4'].copy()

    #--------------------------------------------------------------------------
    # Future
    #--------------------------------------------------------------------------
    flg=0
    if flg==1:
        # Figure out what factor to use to get a specific increase by 2100
        np.interp(2100,np.array([2020,2200]),np.array([Po_obs,5.5*Po_obs]))/Po_obs

    # 3 x increase by 2100
    for k in wfss[zone]['Po Det Scenarios'].keys():
        if k=='H0': continue
        wfss[zone]['Po Det Scenarios'][k]['F1']=dH[k].copy()
        wfss[zone]['Po Det Scenarios'][k]['F1'][indF2]=np.interp(tv_scn[indF2],np.array([2020,2200]),np.array([Po_obs,5.5*Po_obs]))

    # 5 x observation mean by 2100
    for k in wfss[zone]['Po Det Scenarios'].keys():
        if k=='H0': continue
        wfss[zone]['Po Det Scenarios'][k]['F2']=dH[k].copy()        
        wfss[zone]['Po Det Scenarios'][k]['F2'][indF2]=np.interp(tv_scn[indF2],np.array([2020,2200]),np.array([Po_obs,10*Po_obs]))

    # 20 x observation mean by 2100
    for k in wfss[zone]['Po Det Scenarios'].keys():
        if k=='H0': continue
        wfss[zone]['Po Det Scenarios'][k]['F3']=dH[k].copy()
        wfss[zone]['Po Det Scenarios'][k]['F3'][indF2]=np.interp(tv_scn[indF2],np.array([2020,2200]),np.array([Po_obs,43.75*Po_obs]))

#%% Save

# Remove occurrence data (no need to save)
flg=0
if flg==1:
    for zone in wfss.keys():
        del wfss[zone]['Occurrence']

# Save
gu.opickle(PathData,wfss)

#%% Plot deteriministic component of scenarios
ufcs.Plot_WildfireScenarios(meta)

def Plot_WildfireScenarios(meta):
    PathData=meta['Paths']['Model']['Taz Datasets'] + '\\Wildfire Stats and Scenarios\Wildfire_Stats_Scenarios_By_BGCZ.pkl'
    wfss=gu.ipickle(PathData,wfss)
    zone='SBS'
    
    itH=np.where(tv_scn<1920)[0]
    itF=np.where(tv_scn>=1920)[0]
    
    plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(14.5,6.5)); lw=1
    rc=patches.Rectangle((1920,0),103,20,facecolor=[0.92,0.92,0.92])
    ax.add_patch(rc)
    # Observations
    ax.plot(tv_obs,wfss[zone]['Area Wildfire']/wfss[zone]['Area Zone']*100,'k.',ms=2,lw=0.25,label='Observations')
    # Scenario 4
    ax.plot(tv_scn,wfss[zone]['Po Det Scenarios']['H4']['F0'],'b-',color=[0.47,0.69,0.97],lw=2,label='Wildfire Occurrence Scn-H4F0: Constant PI average at H2 level')
    # Scenario 1
    ax.plot(tv_scn,wfss[zone]['Po Det Scenarios']['H1']['F1'],'b--',color=[1,0.5,0],lw=lw,label='Wildfire Occurrence Scn-H1F1: Pre-observation period calibrated against firescars,\n 3 x observed level by 2100')
    # Scenario 2
    ax.plot(tv_scn,wfss[zone]['Po Det Scenarios']['H2']['F2'],'b-.',color=[1,0,0],lw=lw,label='Wildfire Occurrence Scn-H2F2: Pre-observation period calibrated against firescars,\n 5 x observed level by 2100')
    # Scenario 3
    ax.plot(tv_scn,wfss[zone]['Po Det Scenarios']['H3']['F3'],'b:',color=[0.45,0,0],lw=lw,label='Wildfire Occurrence Scn-H3F3: Pre-observation period calibrated against firescars,\n 20 x observed level by 2100')
    
    ax.annotate('Observation\nperiod',(1920+50,2.25),ha='center',fontweight='bold',color=[0.5,0.5,0.5])
    ax.legend(loc='upper left',frameon=False)
    ax.set(xticks=np.arange(tv_scn[0],tv_scn[-1]+100,100),ylabel='Probability of occurrence (% yr$^-$$^1$)',xlabel='Time, calendar year',
           yscale='linear',yticks=np.arange(0,6,1),ylim=[0,5],xlim=[500-0.5,tv_scn[-1]+1+0.5])
    #ax.get_yaxis().set_minor_formatter(ScalarFormatter())
    #ax.set_ylim([0,5])
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=1.5)
    
    plt.tight_layout()
    gu.PrintFig(PathFigures + '\\Wildfire_Scenarios_ts_' + zone,'png',900)
    return


