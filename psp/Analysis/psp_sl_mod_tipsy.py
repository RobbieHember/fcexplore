
"""
PSP - MODIFY TIPSY
"""

#%% Import modules

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
import fcgadgets.macgyver.utilities_general as gu
import fcexplore.psp.Processing.psp_utilities as ugp

# import pymc3 as pm
# import arviz as az
# import autograd.numpy as np
# from autograd import grad
# import matplotlib
# #import seaborn as sns
# import scipy.stats as stats
# import scipy.optimize as optimize

#%% Import data

gp=gu.SetGraphics('Manuscript')
metaGP,gplt=ugp.ImportPSPs(type='Stand')

# Import TIPSY
tipsy=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\TIPSY\tipsy_sw_carbon_plfd_si18_sph1600_oafdef.xlsx')
tipsy['Gsw']=gu.movingave(tipsy['Gsw'],10,'center')

# Add TIPSY to ground plot data
gplt['Gsw TIPSY']=np.nan*np.ones(gplt['Age VRI t0'].size)
for i in range(gplt['Age VRI t0'].size):
    if gplt['Csw G Surv'][i]>0:
        ind=np.where( (tipsy['Age']>=gplt['Age VRI t0'][i]) & (tipsy['Age']<gplt['Age VRI t1'][i]) )[0]
        gplt['Gsw TIPSY'][i]=np.mean(tipsy['Gsw'][ind])

gplt['Time']=(gplt['Year t0']+gplt['Year t1'])/2

#%%

# Filter
ikp=np.where( (gplt['PTF CN']==1) & \
    (gplt['Age VRI t0']>0) & \
    (gplt['Age VRI t0']<=100) & \
    (gplt['Csw Net']>0) )[0]
#    (gplt['Spc1 L ID t0']!=metaGP['LUT']['Species']['AS']) & \
# (gplt['Ecozone BC L1']!=metaGP['LUT']['Ecozone BC L1']['CWH']) & (gplt['Ecozone BC L1']!=metaGP['LUT']['Ecozone BC L1']['MH']) & (gplt['Ecozone BC L1']!=metaGP['LUT']['Ecozone BC L1']['ICH']) & \

#%%

AgeS=gplt['Age VRI t0']-50
TimeS=gplt['Time']-2010
ATS=AgeS*TimeS
y=gplt['Csw G Surv']
x=np.column_stack([gplt['Gsw TIPSY'],AgeS,TimeS,ATS])

ikp=np.where( (gplt['PTF CN']==1) &
    (np.isnan(np.sum(x,axis=1)+y)==False) &
    (gplt['Age VRI t0']>0) & (gplt['Age VRI t0']<=300) & \
    (gplt['Time']>0) & (gplt['Time']<=2040) & \
    (gplt['Csw Net']>-20) & (gplt['Csw Net']<20) & \
    (gplt['Gsw TIPSY']>0) & (gplt['Gsw TIPSY']<20) & \
    (ATS>-1000) & (ATS<1000) )[0]

x1=sm.tools.tools.add_constant(x)
md=sm.OLS(y[ikp],x1[ikp,:]).fit()
md.summary()
#xhat=np.linspace(np.min(x1[:,1]),np.max(x1[:,1]),10)
#yhat=md.predict(np.c_[np.ones(xhat.size),xhat])


#%% Model

def age_modifier(b0,b1,b2,A_o,G_tipsy):
    G_hat = G_tipsy * ( b0+(b1-b0)*1/(1+np.exp(-b2*(A_o-50))) )
    return G_hat

#%%

plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(20,9))
ax[0].plot(tipsy['Age'],0.5*tipsy['Gsw'],'b-')
ax[0].plot(tipsy['Age'],tipsy['Gsw'],'g--')

#%%

t=np.arange(1850,2051,1)
t0=1900
t1=2020
m0=0.5
m1_young=1.1
m1_old=0.7

plt.close('all')
Age=120
Age0=75
Age1=0
m1=gu.Clamp(m0+(Age1-Age0)/(Age1-Age0)*(Age-Age0),m1_old,m1_young)
ft=gu.Clamp(m0+(m1-m0)/(t1-t0)*(t-t0),m0,m1)
plt.plot(t,ft,'b-')
Age=20
Age0=75
Age1=0
m1=gu.Clamp(m0+(Age1-Age0)/(Age1-Age0)*(Age-Age0),m1_old,m1_young)
ft=gu.Clamp(m0+(m1-m0)/(t1-t0)*(t-t0),m0,m1)
plt.plot(t,ft,'g--')

#%%

def modifier(b0,b1,b2,b3,Time,Age,G_tipsy):
    #fY=np.maximum(1.0,1+(b0/100)*(t-1900))
    #fA=( b1+(b2-b1)*1/(1+np.exp(-b3*(Age-40))) )
    #fE=fY * fA
    TimeS=Time-2010
    AgeS=Age-60
    G_hat = G_tipsy + b0*TimeS + b1*AgeS + b2*(TimeS*AgeS)
    return G_hat

b0=1.0
b1=0.25
b2=1.0
b3=-0.1

t=np.arange(1900,2021,1)

plt.close('all')
Age=20*np.ones(t.size)
plt.plot(t,modifier(b0,b1,b2,b3,t,Age,1.1),'g-')
Age=120*np.ones(t.size)
plt.plot(t,modifier(b0,b1,b2,b3,t,Age,1.1),'b--')






#%%

b0=0.5
b1=2
b2=-0.1
xhat=np.arange(1,301,1)

plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(20,9))
bw=10; bin=np.arange(bw,200+bw,bw)
N,mu,med,sig,se=gu.discres(A_o,G_o,bw,bin)
ax[0].plot(bin,mu,'bo',ms=6)
ax[0].plot(tipsy['Age'],tipsy['Gsw'],'r-')

yhat=age_modifier(b0,b1,b2,tipsy['Age'],tipsy['Gsw'])
ax[0].plot(tipsy['Age'],yhat,'c--')

N,mu,med,sig,se=gu.discres(A_o,gplt['Gsw TIPSY'],bw,bin)
ax[0].plot(bin,mu,'gs',ms=6)

bw=20; bin=np.arange(bw,200+bw,bw)
N,mu,med,sig,se=gu.discres(A_o,C_o,bw,bin)
ax[1].plot(bin,mu,'bo',ms=6)
ax[1].plot(tipsy['Age'],tipsy['Csw'],'r-')




