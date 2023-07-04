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
from scipy import stats
import statsmodels.formula.api as smf
from scipy.optimize import curve_fit
import fcgadgets.macgyver.utilities_general as gu
import fcgadgets.macgyver.utilities_gis as gis
import fcexplore.psp.Processing.psp_utilities as ugp
from fcgadgets.bc1ha import bc1ha_utilities as u1ha

#%% Set figure properties

gp=gu.SetGraphics('Manuscript')

#%% Import data

metaGP,gplt=ugp.ImportPSPs(type='Stand')

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
#wth=20000; x0=1585000; y0=550000
wth=11183; x0=1585000; y0=550000

# N. Vancouver Island
#wth=11183; x0=983000; y0=526000

A=(wth*2)**2/10000/1e3 # Kha
ind=np.where( (np.abs(gplt['X']-x0)<wth) & (np.abs(gplt['Y']-y0)<wth) )[0]
bin=np.arange(5,ind.size,5)
rho=bin/A
print(A) # Kha
print(ind.size) # Number plots
print(ind.size/A) # Sample frequency (plots/Kha)

plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(7.8,7))
ax.plot(gplt['X'],gplt['Y'],'bo',ms=5,mfc=None)
ax.plot(gplt['X'][ind],gplt['Y'][ind],'r.')

n=300
mu=np.zeros((bin.size,n))
se=np.zeros((bin.size,n))
cv=np.zeros((bin.size,n))
for i in range(bin.size):
    for j in range(n):
        r=np.random.random(ind.size)
        ir=np.argsort(r)[0:bin[i]]
        mu[i,j]=np.mean(gplt['Ctot L t0'][ind[ir]])
        se[i,j]=2*np.std(gplt['Ctot L t0'][ind[ir]])/np.sqrt(bin[i])
        cv[i,j]=np.std(gplt['Ctot L t0'][ind[ir]])/np.mean(gplt['Ctot L t0'][ind[ir]])*100

E=np.mean(se,axis=1)/np.mean(mu)*100

# Fitting
xhat=np.arange(1,1000,1)
p,pcov=curve_fit(func,bin,E,[45,-0.4])
#print(p)
yhat=func(xhat,p[0],p[1])

# Find sample size that matches 5% error
ind2=np.where(np.abs(yhat-20)==np.min(np.abs(yhat-20)))[0]

#plt.close('all');
fig,ax=plt.subplots(1,figsize=gu.cm2inch(7.8,7))
ax.plot(xhat/A,20+0*xhat,'r--',label='Standard')
ax.plot(bin/A,E,'ko',ms=5,mec='k',mfc='w',label='Observations (FAIB 2023)')
#ax.plot(xhat,func(xhat,5,25,0.03),'k-',label='OLS fit')
ax.plot(xhat/A,yhat,'k-',label='OLS fit')
ax.text(xhat[ind2]/A,20+2,str(np.round(xhat[ind2[0]]/A,decimals=1)),fontsize=12,ha='center')
ax.set(ylabel='Relative error (%)',xlabel='Sample frequency (plots Kha$^{-1}$)',xlim=[0,1.5],ylim=[0,50])
ax.legend(loc='upper center',facecolor=[1,1,1],frameon=False);
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=2)
plt.tight_layout()
#gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Biomass\Biomass_SamplingPower','png',900)

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

plt.close('all');
fig,ax=plt.subplots(1,figsize=gu.cm2inch(7.8,7))
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


