
"""
PSP - DESCRIPTIVE STATISTICS
"""

#%% Import modules

import numpy as np
import gc as garc
import time
import matplotlib.pyplot as plt
import geopandas as gpd
import pandas as pd
import fiona
import statsmodels.formula.api as smf
import fcgadgets.macgyver.utilities_general as gu
import fcgadgets.macgyver.utilities_gis as gis
import fcexplore.psp.Processing.psp_utilities as ugp
from fcgadgets.bc1ha import bc1ha_utilities as u1ha

#%% Import data

gp=gu.SetGraphics('Manuscript')

metaGP,gplt=ugp.ImportPSPs(type='Stand')

lut_1ha=u1ha.Import_BC1ha_LUTs()

zBGCZ=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\becz.tif')

zLCC1=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\LandUseLandCover\LandCoverClass1.tif')


#%%

plt.plot(gplt['Age VRI t0'],gplt['Age Max t0'],'b.')


#%% Age responses

ind=np.where( (gplt['PTF CN']==1) )[0]

#ind=np.where( (gplt['PTF YSM']==1) )[0]

# Coast
#ind=np.where( (gplt['PTF CN']==1) & (gplt['Ecozone BC L1']==metaGP['LUT']['Ecozone BC L1']['CWH']) | (gplt['PTF CN']==1) & (gplt['Ecozone BC L1']==metaGP['LUT']['Ecozone BC L1']['MH']) )[0]

#ind=np.where( (gplt['PTF CN']==1) & (gplt['Ecozone BC L1']==metaGP['LUT']['Ecozone BC L1']['ICH']) )[0]

#ind=np.where( (gplt['PTF CN']==1) & (gplt['Ecozone BC L1']==metaGP['LUT']['Ecozone BC L1']['ESSF']) )[0]

# Interior
#ind=np.where( (gplt['PTF CN']==1) & (gplt['Ecozone BC L1']!=metaGP['LUT']['Ecozone BC L1']['CWH']) & (gplt['Ecozone BC L1']!=metaGP['LUT']['Ecozone BC L1']['MH']) & (gplt['Ecozone BC L1']!=metaGP['LUT']['Ecozone BC L1']['ICH']) )[0]

#plt.plot(gplt['Age VRI t0'][ind],gplt['Ctot L t0'][ind],'b.')
#plt.close('all')
#plt.plot(gplt['Age VRI t0'][ind],gplt['Age Max t0'][ind],'b.')
#plt.plot(gplt['Age VRI t0'][ind],gplt['Age Min t0'][ind],'r.')

#plt.close('all')
#plt.plot(gplt['Age Med t0'][ind],gplt['Age Max t0'][ind],'b.')

x=gplt['Age Med t0'][ind]
bw=25; bin=np.arange(bw,300,bw)
xhat=np.arange(1,301,1)

lw=1; ms=5; mec='b'; mfc='w'; cl='b';
plt.close('all'); fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(15.5,11))
y=gplt['Ctot G Surv'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
ax[0,0].plot(bin,mu,'ko',ms=ms,lw=lw,color=cl,mfc=mfc,mec=mec)
yhat=np.interp(xhat,bin,mu)
ax[0,0].plot(xhat,yhat,'b-',lw=0.5)
ax[0,0].set(ylabel='Survivor growth (MgC ha$^-$$^1$ yr$^-$$^1$)',xlabel='Age, years',xlim=[0,300])
ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=gp['tickl'])

# #ind2=np.where( (gplt['PTF YSM']==1) & (gplt['Ecozone BC L1']!=metaGP['LUT']['Ecozone BC L1']['CWH']) & (gplt['Ecozone BC L1']!=metaGP['LUT']['Ecozone BC L1']['MH']) & (gplt['Ecozone BC L1']!=metaGP['LUT']['Ecozone BC L1']['ICH']) )[0]
# ind2=np.where( (gplt['PTF YSM']==1) & (gplt['Ecozone BC L1']==metaGP['LUT']['Ecozone BC L1']['CWH']) | (gplt['PTF YSM']==1) & \
#     (gplt['Ecozone BC L1']==metaGP['LUT']['Ecozone BC L1']['MH']) )[0]
# x=gplt['Age Med t0'][ind2]
# y=gplt['Ctot G Surv'][ind2]
# N,mu,med,sig,se=gu.discres(x,y,bw,bin)
# ax[0,0].plot(bin,mu,'ks',ms=ms,lw=lw,color='k',mfc='w',mec='r')

y=gplt['Ctot G Recr'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
ax[0,1].plot(bin,mu,'ko',ms=ms,lw=lw,color=cl,mfc=mfc,mec=mec)
ax[0,1].set(ylabel='Recruitment growth (MgC ha$^-$$^1$ yr$^-$$^1$)',xlabel='Age, years',xlim=[0,300])
ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=gp['tickl'])

y=gplt['Ctot Mort+Lost'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
ax[1,0].plot(bin,mu,'ko',ms=ms,lw=lw,color=cl,mfc=mfc,mec=mec,label='Total')
y=gplt['Ctot Mort+Lost'][ind]-gplt['Ctot Mort+Lost Fire'][ind]-gplt['Ctot Mort+Lost Insect'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
ax[1,0].plot(bin,mu,'ks',ms=ms,lw=lw,color='g',mfc=mfc,mec='g',label='W/O fire and insects')
ax[1,0].set(ylabel='Mortality (MgC ha$^-$$^1$ yr$^-$$^1$)',xlabel='Age, years',xlim=[0,300])
ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=gp['tickl'])
leg=ax[1,0].legend(loc='lower center',frameon=False,facecolor=None,edgecolor='w')

y=gplt['Ctot Net'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
ax[1,1].plot(bin,mu,'ks',ms=ms,lw=lw,color=cl,mfc=mfc,mec=mec,label='Total')
yhat_tot=np.interp(xhat,bin[np.isnan(mu)==False],mu[np.isnan(mu)==False])
ax[1,1].plot(xhat,yhat_tot,'b-',lw=0.5)

y=gplt['Ctot Net'][ind]+gplt['Ctot Mort+Lost Fire'][ind]+gplt['Ctot Mort+Lost Insect'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
ax[1,1].plot(bin,mu,'ko',ms=ms,lw=lw,color='g',mfc=mfc,mec='g',label='W/O fire and insects')
yhat_wofi=np.interp(xhat,bin[np.isnan(mu)==False],mu[np.isnan(mu)==False])
ax[1,1].plot(xhat,yhat_wofi,'g--',lw=0.5)

ax[1,1].set(ylabel='Net growth (MgC ha$^-$$^1$ yr$^-$$^1$)',xlabel='Age, years',xlim=[0,300])
ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=gp['tickl'])
leg=ax[1,1].legend(loc='upper right',frameon=False,facecolor=None,edgecolor='w')
gu.axletters(ax,plt,0.03,0.9,FontColor=gp['cla'],LetterStyle='Caps',FontWeight='Bold')

# Plot biomass

#ind=np.where( (gplt['PTF YSM']==1) )[0]

#ind=np.where( (gplt['PTF CN']==1) )[0]

#ind=np.where( (gplt['PTF CN']==1) & (gplt['Ecozone BC L1']==metaGP['LUT']['Ecozone BC L1']['ESSF']) )[0]

#ind=np.where( (gplt['PTF CNY']==1) & (gplt['Ecozone BC L1']!=metaGP['LUT']['Ecozone BC L1']['CWH']) & (gplt['Ecozone BC L1']!=metaGP['LUT']['Ecozone BC L1']['MH']) & \
#             (gplt['Ecozone BC L1']!=metaGP['LUT']['Ecozone BC L1']['ICH']) )[0]

x=gplt['Age Med t0'][ind]
y=gplt['Ctot L t0'][ind]
bw=10; bin=np.arange(bw,300,bw)

lw=1; ms=5; mec='b'; mfc='w'; cl='b';
fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(15.5,11))
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
ax[0,0].plot(bin,mu,'ko',ms=ms,lw=lw,color='k',mfc='w',mec='k')
ax[0,0].plot(xhat,np.cumsum(yhat_wofi),'g-',ms=ms,lw=lw,color='g',mfc=mfc,mec=mec)
ax[0,0].plot(xhat,np.cumsum(yhat_tot),'b--',ms=ms,lw=lw,color='b',mfc=mfc,mec=mec)
ax[0,0].set(ylabel='Biomass (MgC ha$^-$$^1$)',xlabel='Age, years',xlim=[0,300])
ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=gp['tickl'])

ind=np.where( (gplt['PTF YSM']==1) )[0]
x=gplt['Age Med t0'][ind]
y=gplt['Ctot L t0'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
ax[0,0].plot(bin,mu,'ks',ms=ms,lw=lw,color='k',mfc='w',mec='r')

#%% Biomass age response in areas with no history of major disturbance

zRef=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\BC_Land_Mask.tif')
zH=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\VEG_CONSOLIDATED_CUT_BLOCKS_SP_All.tif')
zF=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\PROT_HISTORICAL_FIRE_POLYS_SP_All.tif')
zI=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\IBM_Mask_ExcTrace.tif')

indG=gis.GetGridIndexToPoints(zRef,gplt['X'],gplt['Y'])

flgH=np.zeros(gplt['PTF CN'].size);
ind=np.where(zH['Data'][indG]>0)[0]
flgH[ind]=1

flgF=np.zeros(gplt['PTF CN'].size);
ind=np.where(zF['Data'][indG]>0)[0]
flgF[ind]=1

flgI=np.zeros(gplt['PTF CN'].size);
ind=np.where(zI['Data'][indG]>0)[0]
flgI[ind]=1

ind=np.where( (gplt['PTF CNV']==1) & (gplt['Ecozone BC L1']!=metaGP['LUT']['Ecozone BC L1']['CWH']) & (gplt['Ecozone BC L1']!=metaGP['LUT']['Ecozone BC L1']['MH']) & \
             (gplt['Ecozone BC L1']!=metaGP['LUT']['Ecozone BC L1']['ICH']) & \
            (flgH==0) & (flgF==0) & (flgI==0) )[0]

x1=gplt['Age VRI t0'][ind]
y1=gplt['Ctot L t0'][ind]
N,mu,med,sig,se=gu.discres(x1,y1,bw,bin)
ax[0,0].plot(bin,mu,'r^',ms=ms,lw=lw,color='k',mfc='w',mec='r')

ind=np.where( (gplt['PTF CN']==1) & (gplt['Ecozone BC L1']!=metaGP['LUT']['Ecozone BC L1']['CWH']) & (gplt['Ecozone BC L1']!=metaGP['LUT']['Ecozone BC L1']['MH']) & \
             (gplt['Ecozone BC L1']!=metaGP['LUT']['Ecozone BC L1']['ICH']) & (flgI==1) )[0]
x1=gplt['Age VRI t0'][ind]
y1=gplt['Ctot L t0'][ind]
N,mu,med,sig,se=gu.discres(x1,y1,bw,bin)
ax[0,0].plot(bin,mu,'rd',ms=ms,lw=lw,color='k',mfc='c',mec='c')









#%% Age responses (by climate class)

cc=2

ind=np.where( (gplt['PTF CN']==1) & (gplt['Climate Class']>=cc) )[0]

x=gplt['Age VRI t0'][ind]
bw=25; bin=np.arange(bw,300,bw)
xhat=np.arange(1,301,1)

lw=1; ms=5; mec='b'; mfc='w'; cl='b';
plt.close('all'); fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(15.5,11))
y=gplt['Ctot G Surv'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
ax[0,0].plot(bin,mu,'ko',ms=ms,lw=lw,color=cl,mfc=mfc,mec=mec)
yhat=np.interp(xhat,bin,mu)
ax[0,0].plot(xhat,yhat,'b-',lw=0.5)
ax[0,0].set(ylabel='Survivor growth (MgC ha$^-$$^1$ yr$^-$$^1$)',xlabel='Age, years',xlim=[0,300])
ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=gp['tickl'])

y=gplt['Ctot G Recr'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
ax[0,1].plot(bin,mu,'ko',ms=ms,lw=lw,color=cl,mfc=mfc,mec=mec)
ax[0,1].set(ylabel='Recruitment growth (MgC ha$^-$$^1$ yr$^-$$^1$)',xlabel='Age, years',xlim=[0,300])
ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=gp['tickl'])

y=gplt['Ctot Mort+Lost'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
ax[1,0].plot(bin,mu,'ko',ms=ms,lw=lw,color=cl,mfc=mfc,mec=mec,label='Total')
y=gplt['Ctot Mort+Lost'][ind]-gplt['Ctot Mort+Lost Fire'][ind]-gplt['Ctot Mort+Lost Insect'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
ax[1,0].plot(bin,mu,'ks',ms=ms,lw=lw,color='g',mfc=mfc,mec='g',label='W/O fire and insects')
ax[1,0].set(ylabel='Mortality (MgC ha$^-$$^1$ yr$^-$$^1$)',xlabel='Age, years',xlim=[0,300])
ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=gp['tickl'])
leg=ax[1,0].legend(loc='lower center',frameon=False,facecolor=None,edgecolor='w')

y=gplt['Ctot Net'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
ax[1,1].plot(bin,mu,'ks',ms=ms,lw=lw,color=cl,mfc=mfc,mec=mec,label='Total')
yhat_tot=np.interp(xhat,bin[np.isnan(mu)==False],mu[np.isnan(mu)==False])
ax[1,1].plot(xhat,yhat_tot,'b-',lw=0.5)

y=gplt['Ctot Net'][ind]+gplt['Ctot Mort+Lost Fire'][ind]+gplt['Ctot Mort+Lost Insect'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
ax[1,1].plot(bin,mu,'ko',ms=ms,lw=lw,color='g',mfc=mfc,mec='g',label='W/O fire and insects')
yhat_wofi=np.interp(xhat,bin[np.isnan(mu)==False],mu[np.isnan(mu)==False])
ax[1,1].plot(xhat,yhat_wofi,'g--',lw=0.5)

ax[1,1].set(ylabel='Net growth (MgC ha$^-$$^1$ yr$^-$$^1$)',xlabel='Age, years',xlim=[0,300])
ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=gp['tickl'])
leg=ax[1,1].legend(loc='upper right',frameon=False,facecolor=None,edgecolor='w')
gu.axletters(ax,plt,0.03,0.9,FontColor=gp['cla'],LetterStyle='Caps',FontWeight='Bold')

# Plot biomass

ind=np.where( (gplt['PTF CN']==1) & (gplt['Climate Class']>=cc) )[0]

x=gplt['Age VRI t0'][ind]
y=gplt['Ctot L t0'][ind]
bw=10; bin=np.arange(bw,300,bw)

lw=1; ms=5; mec='b'; mfc='w'; cl='b';
fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(15.5,11))
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
ax[0,0].plot(bin,mu,'ko',ms=ms,lw=lw,color='k',mfc='w',mec='k')
ax[0,0].plot(xhat,np.cumsum(yhat_wofi),'g-',ms=ms,lw=lw,color='g',mfc=mfc,mec=mec)
ax[0,0].plot(xhat,np.cumsum(yhat_tot),'b--',ms=ms,lw=lw,color='b',mfc=mfc,mec=mec)
ax[0,0].set(ylabel='Biomass (MgC ha$^-$$^1$)',xlabel='Age, years',xlim=[0,300])
ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=gp['tickl'])


#%%













#%% Modelling

#fAg=inline(['max(0,b(1).*(1+((b(2).*(x(:,1)./b(3)).^b(4)-1)./(exp(x(:,1)./b(3))))))'],'b','x'); % Weibull
#fAm=inline(['max(0,b(1).*(1+((b(2).*((x-b(5))./b(3)).^b(4)-1)./(exp((x-b(5))./b(3))))))'],'b','x');

def fun(x,a,b,c,d):
    yhat=d+a*((1-np.exp(-b*x))**(1/(1-c))) # From Dzierzon and Mason 2006
    return yhat
xhat=np.arange(0,100,1)

indM=np.where( (tl['DBH']>0) & (np.isnan(tl['H'])==True) | (tl['DBH']>0) & (tl['H']<=0) )[0]
indM.size/tl['DBH'].size

indG=np.where( (tl['DBH']>0) & (tl['H']>0) )[0]
indG.size/tl['DBH'].size

# Global model
ikp=np.where( (tl['DBH']>0) & (tl['H']>0) & (tl['Vital Status']==meta['LUT']['Vital Status']['Live']) & (tl['ID DA1']==meta['LUT']['Damage Agents']['None']) & (tl['Stature']==meta['LUT']['Stature']['Standing']) )[0]
x=tl['DBH'][ikp]
y=tl['H'][ikp]
popt_glob0=[26,0.1,0.66,2]
popt_glob,pcov=curve_fit(fun,x,y,popt_glob0)
#yhat=fun(xhat,popt_glob[0],popt_glob[1],popt_glob[2],popt_glob[3])
rs_glob,txt=gu.GetRegStats(x,y)



#%% Aridity

#ind=np.where( (gplt['PTF CN']==1) & (gplt['Ecozone BC L1']==metaGP['LUT']['Ecozone CA L1']['Pacific_Maritime']) )[0]
ind=np.where( (gplt['PTF CN']==1) )[0]
#ind=np.where( (gplt['PTF CN']==1) & (gplt['Ecozone BC L1']==metaGP['LUT']['Ecozone BC L1']['SBPS']) )[0]

x=gplt['Climate Class'][ind]
y=gplt['Ctot G Surv'][ind]

bw=1; bin=np.arange(1,6,bw)
N,mu,med,sig,se=gu.discres(x,y,bw,bin)

plt.close('all')
plt.plot(bin,mu,'bo')

