"""
SAMPLING POWER
"""

#%% Import modules
import numpy as np
import gc as garc
import time
import matplotlib.pyplot as plt
import geopandas as gpd
import pandas as pd
import scipy.stats as sp
import statsmodels.formula.api as smf
from scipy.optimize import curve_fit
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcexplore.field_plots.Processing.psp_util as ugp
from fcgadgets.bc1ha import bc1ha_util as u1ha
gp=gu.SetGraphics('Manuscript')

#%% Import data
meta=u1ha.Init()
meta,gpt=ugp.ImportGroundPlotData(meta,type='Stand')
#d=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\GroundPlots\YSM and CMI Grid Coordinates\ground_plot_grid.xlsx')

# Function (tree height as a function of basal area)
def func(x,a,b):
    #y=a+(b-a)*np.exp(-c*x)
    #y=a*np.exp(-b*x)
    y=a*x**b
    return y

#%% Tree biomass

#wth=10000; x0=1380000; y0=800000
#wth=10000; x0=1440000; y0=700000

# Southern interior
wth=10000; x0=1585000; y0=550000
#wth=20000; x0=1585000; y0=550000
#wth=30000; x0=1585000; y0=550000

# N. Vancouver Island
#wth=11183; x0=983000; y0=526000

# Area of ROI
A=(wth*2)**2/10000/1e3 # Kha

# Index to plots within ROI
ind=np.where( (np.abs(gpt['X']-x0)<wth) & (np.abs(gpt['Y']-y0)<wth) )[0]
N_All=ind.size

# Gradient in sample size ranging from 5 plots to the full number of available plots
N_Samp=np.arange(5,N_All,5)

# Sampling density
rho_Samp=N_Samp/A
print(A) # Kha
#print(N_All) # Number plots
#print(N_All/A) # Sample frequency (plots/Kha)

plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(7.8,7))
ax.plot(gpt['X'],gpt['Y'],'bo',ms=5,mfc=None)
ax.plot(gpt['X'][ind],gpt['Y'][ind],'r.')

n=300
mu=np.zeros((N_Samp.size,n))
se=np.zeros((N_Samp.size,n))
cv=np.zeros((N_Samp.size,n))
for i in range(N_Samp.size):
    for j in range(n):
        r=np.random.random(ind.size)
        ir=np.argsort(r)[0:N_Samp[i]]
        mu[i,j]=np.mean(gpt['Ctot L t0'][ind[ir]])
        se[i,j]=2*np.std(gpt['Ctot L t0'][ind[ir]])/np.sqrt(N_Samp[i])
        cv[i,j]=np.std(gpt['Ctot L t0'][ind[ir]])/np.mean(gpt['Ctot L t0'][ind[ir]])*100

E=np.mean(se,axis=1)/np.mean(mu)*100

# Fitting
xhat=np.arange(1,1000,1)
p,pcov=curve_fit(func,N_Samp,E,[45,-0.4])
#print(p)
yhat=func(xhat,p[0],p[1])

# Find sample size that matches 5% error
ind2=np.where(np.abs(yhat-20)==np.min(np.abs(yhat-20)))[0]

plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(7.8,7))
ax.plot(xhat/A,20+0*xhat,'r--',label='Standard')
ax.plot(N_Samp/A,E,'ko',ms=5,mec='k',mfc='w',label='Observations (FAIB 2023)')
#ax.plot(xhat,func(xhat,5,25,0.03),'k-',label='OLS fit')
ax.plot(xhat/A,yhat,'k-',label='OLS fit')
ax.text(xhat[ind2]/A,20+2,str(np.round(xhat[ind2[0]]/A,decimals=1)),fontsize=12,ha='center')
ax.set(ylabel='Relative error (%)',xlabel='Sample frequency (plots Kha$^{-1}$)',xlim=[0,1.5],ylim=[0,50])
ax.legend(loc='upper center',facecolor=[1,1,1],frameon=False);
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=2)
plt.tight_layout()
#gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Biomass\Biomass_SamplingPower','png',900)

#%%

res=sp.normaltest(gpt['Ctot L t0'][ind])
#res=sp.normaltest(np.log(gpt['Ctot L t0'][ind]))

#s=0.54; mean,var,skew,kurt=sp.lognorm.stats(s,moments='mvsk')

x=np.linspace(sp.lognorm.ppf(0.01,s),sp.lognorm.ppf(0.99,s),100)
fig,ax=plt.subplots(1,1)
ax.plot(x,sp.lognorm.pdf(x,s),'r-', lw=5, alpha=0.6, label='lognorm pdf')


#shape,loc,scale=sp.lognorm.fit(gpt['Ctot L t0'][ind])
sp.lognorm.mean(shape,loc,scale)
sp.lognorm.std(shape,loc,scale)
#sp.lognorm.stats(shape,loc,scale,moments='mv')
#sp.lognorm.interval(0.05,shape,loc,scale)

x=np.arange(0,250,1)
plt.close('all'); fig,ax=plt.subplots(1,1)
ax.hist(gpt['Ctot L t0'][ind],bins=np.arange(0,250,10),density=True)
ax.plot(x,sp.lognorm.pdf(x,s,loc,scale),'r-',lw=5,alpha=0.6,label='lognorm pdf')


#%% Soil organic carbon

# Import Canadian Upland forest database
soc=gu.ipickle(r'C:\Users\rhember\Documents\Data\Soils\Shaw et al 2018 Database\SITES.pkl')

wth=150000
A=(wth*2)**2/10000/1e6
ind=np.where( (np.abs(soc['x']-1380000)<wth) & (np.abs(soc['y']-800000)<wth) )[0]
print(ind.size)
bin=np.arange(5,ind.size,5)
rho=bin/A

n=1000
mu=np.zeros((bin.size,n))
se=np.zeros((bin.size,n))
cv=np.zeros((bin.size,n))
for i in range(bin.size):
    for j in range(n):
        r=np.random.random(ind.size)
        ir=np.argsort(r)[0:bin[i]]
        mu[i,j]=np.mean(soc['TOT_C_THA'][ind[ir]])
        se[i,j]=np.std(soc['TOT_C_THA'][ind[ir]])/np.sqrt(bin[i])
        cv[i,j]=np.std(soc['TOT_C_THA'][ind[ir]])/np.mean(soc['TOT_C_THA'][ind[ir]])*100

E=np.mean(se,axis=1)/np.mean(mu)*100

# Fitting
xhat=np.arange(0,200,1)
p,pcov=curve_fit(func,bin,E,[45,-0.4])
print(p)
yhat=func(xhat,p[0],p[1])

# Find sample size that matches 5% error
ind=np.where(np.abs(yhat-5)==np.min(np.abs(yhat-5)))[0]

plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(7.8,7))
ax.plot(xhat/A,5+0*xhat,'r--',label='Standard')
ax.plot(bin/A,E,'ko',ms=5,mec='k',mfc='w',label='Observations (Shaw et al. 2018)')
#ax.plot(xhat,func(xhat,5,25,0.03),'k-',label='OLS fit')
ax.plot(xhat/A,yhat,'k-',label='OLS fit')
ax.text(xhat[ind]/A,5+2,str(np.round(xhat[ind[0]]/A,decimals=1)),fontsize=12,ha='center')
ax.set(ylabel='Relative error (%)',xlabel='Sample density (plots per Mha)',xlim=[0,20],ylim=[0,40])
ax.legend(loc='upper center',facecolor=[1,1,1],frameon=False);
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=2)
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Soils\SamplingPower','png',900)


#%% Map

# *** See PSP plot map script ***


