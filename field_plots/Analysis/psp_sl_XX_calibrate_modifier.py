"""
PSP - TIME-DEPENDENT MODELLING
"""

#%% Import modules
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import copy
from scipy import stats
import seaborn as sns
import statsmodels.formula.api as smf
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.bc1ha.bc1ha_util as u1ha
import fcexplore.field_plots.Processing.psp_util as ugp
import fcgadgets.macgyver.util_mcmc as umc
gp=gu.SetGraphics('Manuscript')

#%% Import data
meta=u1ha.Init()
meta,gpt=ugp.ImportGroundPlotData(meta,type='Stand')

#%% Filter

flg=2

# Interior
if flg==1:
    gy=gu.ReadExcel(r'C:\Data\TIPSY\for_calibration_interior.xlsx')
    ikp=np.where( (np.isin(gpt['Ecozone BC L1'],[meta['LUT']['GP']['Ecozone BC L1']['SBS'],meta['LUT']['GP']['Ecozone BC L1']['IDF'],meta['LUT']['GP']['Ecozone BC L1']['SBPS'],meta['LUT']['GP']['Ecozone BC L1']['BWBS']])==True) & \
                 (gpt['PTF CN']==1) & (gpt['Year t0']>0) & (gpt['Year t1']>0) & \
                 (gpt['Age Mean t0']>=0) & (np.isnan(gpt['Ctot Net'])==False) )[0]
    print(ikp.size)

# Coast
if flg==2:
    gy=gu.ReadExcel(r'C:\Data\TIPSY\for_calibration_coast.xlsx')
    ikp=np.where( (np.isin(gpt['Ecozone BC L1'],[meta['LUT']['GP']['Ecozone BC L1']['CWH'],meta['LUT']['GP']['Ecozone BC L1']['MH'],meta['LUT']['GP']['Ecozone BC L1']['ICH'],meta['LUT']['GP']['Ecozone BC L1']['CDF']])==True) & \
                 (gpt['PTF CN']==1) & (gpt['Year t0']>0) & (gpt['Year t1']>0) & \
                 (gpt['Age Mean t0']>=0) & (np.isnan(gpt['Ctot Net'])==False) )[0]
    print(ikp.size)

d={}
d['Delta t']=gpt['Delta t'][ikp]
d['Year t0']=gpt['Year t0'][ikp]
d['Year t1']=gpt['Year t1'][ikp]
d['Age Mean t0']=gpt['Age Mean t0'][ikp]
d['Csw L t0']=gpt['Csw L t0'][ikp]
d['Csw Net']=gpt['Csw Net'][ikp]

# Ensure age isn't older than maximum of GY age curve
d['Age Mean t0']=np.minimum(d['Age Mean t0'],np.max(gy['Age Germination']))

# # Stats
# t_mu=np.mean(d['Year t0'])
# t_sig=np.std(d['Year t0'])
# A_mu=np.mean(d['Age Mean t0'])
# A_sig=np.std(d['Age Mean t0'])

#%% Import GY model

# Index for each value in ground plot DB
d['iA0']=np.zeros(ikp.size,dtype='int16')
d['iAI']=[None]*ikp.size
for i in range(ikp.size):
    dt=d['Delta t'][i]
    d['iA0'][i]=np.where(gy['Age Germination']==int(d['Age Mean t0'][i]))[0]
    d['iAI'][i]=np.where( (gy['Age Germination']>=int(d['Age Mean t0'][i])) & (gy['Age Germination']<=int(d['Age Mean t0'][i]+dt)) )[0]    
    #d['Mod Csw L t0 Init'][i]=gy['Csw'][iA0]
    #d['Mod Csw Net Init'][i]=np.mean(gy['Net Growth C'][iAI])
#d['Mod Csw L t0 Init']=np.nan*np.ones(ikp.size)
#d['Mod Csw Net Init']=np.nan*np.ones(ikp.size)

# Matrix version of GY model
gym={}
gym['Age']=np.tile(gy['Age Germination'],(ikp.size,1)).T
gym['Csw']=np.tile(gy['Csw'],(ikp.size,1)).T
gym['Net Growth C']=np.tile(gy['Net Growth C'],(ikp.size,1)).T
gym['Time']=np.zeros(gym['Age'].shape)
for i in range(ikp.size):
    iA=d['iA0'][i]
    yr0=d['Year t0'][i]
    yr=np.arange(yr0-iA,yr0-iA+gym['Age'].shape[0],1)
    gym['Time'][:,i]=yr

#%% Model
def model(theta):
    global ikp
    bS=theta[0]
    bA=-1*theta[1]
    bT=theta[2]
    bAT=-1*theta[3]
      
    gym0=copy.deepcopy(gym)  
    fG=(bS+bA*(gym0['Age']/75)+bT*(gym0['Time']-1850)+bAT*np.minimum(1.5,(gym0['Time']-1850)*(gym0['Age']/75)))-2
    gym0['Net Growth C']=fG+gym0['Net Growth C']
    gym0['Csw']=np.cumsum(gym0['Net Growth C'],axis=0)
    #plt.plot(gym0['Csw'][:,0])
    #plt.plot(gym0['Net Growth C'][:,0])
    
    yG=np.zeros(ikp.size)
    yB=np.zeros(ikp.size)
    for i in range(ikp.size):
        yG[i]=np.mean(gym0['Net Growth C'][ d['iAI'][i],i ])
        yB[i]=gym0['Csw'][ d['iA0'][i],i ]
    #print(np.mean(yG))
    #print(np.mean(yB))
    
    # Find an appropriate adjustment factor for biomass so that the likelihood
    # is similar to that of net growth
    #sfBiomass=1.0/100.0
    #sfBiomass=np.mean(d['Csw Net'])/np.mean(d['Csw L t0'])
    sfBiomass=0.014
    #sfBiomass=np.nanmean(yG)/np.nanmean(yB)
    #np.nanmean(gpt['Ctot L t0'])*sfBiomass
     
    # Biomass error    
    yO=d['Csw L t0']
    yM=yB
    EB,txt=gu.GetRegStats(yO,yM)    
    y=sfBiomass*yO
    yhat=sfBiomass*yM
    
    # Net growth error
    yO=d['Csw Net']
    yM=yG
    EG,txt=gu.GetRegStats(yO,yM)    
    y=np.append(y,yO)    
    yhat=np.append(yhat,yM)
    
    E=np.array([EB['MEC'],EG['MEC']])
    
    return yhat,y,E

#%% Run
global chain,ll,ar,E
#b=np.array([0,0.0035,0.008,0.0055,1])
thetaInit=np.array([0.23,0.0035,0.008,0.0055,1.25]).reshape(1,-1)
thetaMin= np.array([0.0 ,0.0   ,0.0  ,0.0   ,0.01])
thetaMax= np.array([4.0 ,0.2   ,0.2  ,0.2   ,3.0])

#thetaScale=np.array([0.15,0.003,0.006,0.004,0.35])
thetaScale=np.abs(2.0*thetaInit[0,:])
SearchWidth=np.abs(0.2*thetaInit[0,:]) # As search width goes up, acceptance rate goes down
N_Iter=100000
chain=np.zeros((N_Iter,thetaInit.size))
chain[0]=thetaInit
ll=np.zeros(N_Iter)
ar=np.zeros(N_Iter)
E=np.zeros((N_Iter,2))
iStart=0
chain,ll,ar,E=umc.MetropolisHastings2b(thetaInit,thetaScale,thetaMin,thetaMax,SearchWidth,[],model,N_Iter,chain,ll,ar,E,iStart)

# Plotting:
iStartPostBurn=int(0.5*N_Iter)
plt.close('all'); 
fig,ax=plt.subplots(2,1)
ax[0].plot(ll,'-o') # Plot log likelihood
ax[1].plot(gu.BlockMean(ar,int(N_Iter/30)),'-') # Plot acceptance rate

# Plot parameter posteriors
fig,ax=plt.subplots(1,5)
for i in range(5):
    sns.distplot(chain[iStartPostBurn:,i],ax=ax[i])
    ax[i].set_xlabel(f"param[{i}]")
fig.tight_layout()

# Plot chains
fig,ax=plt.subplots(5,1)
for i in range(5):
    ax[i].plot(chain[:,i])
    ax[i].plot(gu.BlockMean(np.arange(0,N_Iter,1),1000),gu.BlockMean(chain[:,i],1000),lw=1.25,color=[1,1,0])
fig.tight_layout()

fig,ax=plt.subplots(2,1)
ax[0].plot(E[:,0])
ax[1].plot(E[:,1])
fig.tight_layout()

fig,ax=plt.subplots(1)
ax.plot(E[:,0],E[:,1],'.')

#%% Save
flg=0
if flg==1:
    d={'chain':chain,'ll':ll,'ar':ar,'E':E}
    gu.opickle(r'C:\Data\MCMC_Chain.pkl',d)



#%%
#plt.close('all'); 
b=np.mean(chain[iStartPostBurn:,:],axis=0)
#b=np.array([0,0.0035,0.008,0.0055,1])
fig,ax=plt.subplots(1)
plt.plot(gy['Age Germination'],gy['Net Growth C'],'k-')
fG=(b[0]-b[1]*(gy['Age Germination']/75)+b[2]*(1900-1850)-b[3]*((1900-1850)*np.minimum(1.5,(gy['Age Germination']/75))))-2
plt.plot(gy['Age Germination'],fG+gy['Net Growth C'],'b-')
fG=(b[0]-b[1]*(gy['Age Germination']/75)+b[2]*(2020-1850)-b[3]*((2020-1850)*np.minimum(1.5,(gy['Age Germination']/75))))-2
plt.plot(gy['Age Germination'],fG+gy['Net Growth C'],'g--')


#%%








# #%% Model
# def model(theta):
#     global ikp
#     t0=theta[0]
#     t1=theta[1]
#     m0=theta[2]
#     m1_y=theta[3]
#     m1_o=theta[4]
#     Age_y=theta[5]
#     Age_o=theta[6]
#     mSI=theta[7]
  
#     gym0=copy.deepcopy(gym)
  
#     m1=(m1_y+m1_o)-gu.Clamp(m0+(Age_y-gym0['Age'])/(Age_y-Age_o),m1_o,m1_y)
#     fG=mSI*(gu.Clamp(m0+(m1-m0)/(t1-t0)*(gym0['Time']-t0),m0,m1))
#     #fG=1.0
#     gym0['Net Growth C']=fG*gym0['Net Growth C']
#     gym0['Csw']=np.cumsum(gym0['Net Growth C'],axis=0)
#     #plt.plot(gym0['Csw'][:,0])
#     #plt.plot(gym0['Net Growth C'][:,0])
    
#     yG=np.zeros(ikp.size)
#     yB=np.zeros(ikp.size)
#     for i in range(ikp.size):
#         yG[i]=np.mean(gym0['Net Growth C'][ d['iAI'][i],i ])
#         yB[i]=gym0['Csw'][ d['iA0'][i],i ]
#     #print(np.mean(yG))
#     #print(np.mean(yB))
    
#     # Find an appropriate adjustment factor for biomass so that the likelihood
#     # is similar to that of net growth
#     sfBiomass=np.nanmean(yG)/np.nanmean(yB)
#     #np.nanmean(gpt['Ctot L t0'])*sfBiomass
    
#     # Error
#     E={}   
#     # Biomass    
#     yO=d['Csw L t0']
#     yM=yB
#     E['Biomass'],txt=gu.GetRegStats(yO,yM)    
#     y=sfBiomass*yO
#     yhat=sfBiomass*yM
    
#     # Net growth
#     yO=d['Csw Net']
#     yM=yG
#     E['Net Growth'],txt=gu.GetRegStats(yO,yM)    
#     y=np.append(y,yO)
#     yhat=np.append(yhat,yM)
    
#     return yhat,y,E

# #%% Run
# global chain,ll,ar,E
# thetaInit=np.array([1900,2010,1.0,1.0,1.0,10,100,1.0,0.00001]).reshape(1,-1)
# thetaMin=np.array([1850,1990,0.4,1.0,1.0, 0, 50,0.5,0.0000001])
# thetaMax=np.array([1945,2020,1.0,3.0,3.0,40,200,1.5,100.0])

# thetaScale=0.5*thetaInit[0,:]
# #thetaScale=np.array([10,10,0.15,1.0,1.0,10,10,0.2,2]) # Very fast to converge
# #thetaScale=np.array([6,6,0.1,0.75,0.75,6,6,0.15,0.25]) # 
# #thetaScale=np.array([1,1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]) # Very slow to converge
# SearchWidth=0.1*thetaScale
# #SearchWidth=1.5*thetaScale # NOt mixing any better
# N_Iter=10000
# chain=np.zeros((N_Iter,thetaInit.size))
# chain[0]=thetaInit
# ll=np.zeros(N_Iter)
# ar=np.zeros(N_Iter)
# E=[None]*N_Iter
# iStart=0
# chain,ll,ar,E=umc.MetropolisHastings2b(thetaInit,thetaScale,thetaMin,thetaMax,SearchWidth,[],model,N_Iter,chain,ll,ar,E,iStart)




