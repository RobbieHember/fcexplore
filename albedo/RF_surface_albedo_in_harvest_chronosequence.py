'''
Summary of emission intensity from harvesting
'''

#%% Import modules

from os import listdir
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from fcgadgets.macgyver import utilities_general as gu
import scipy.io
import copy
from scipy.optimize import curve_fit
import scipy.io as spio
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats, linalg
from matplotlib.patches import Rectangle

#%% Graphics parameters

fs=7
params={'font.sans-serif':'Arial','font.size':fs,'axes.edgecolor':'black','axes.labelsize':fs,'axes.labelcolor':'black','axes.titlesize':fs,'axes.linewidth':0.5,'lines.linewidth':0.5,
        'text.color':'black','xtick.color':'black','xtick.labelsize':fs,'xtick.major.width':0.5,'xtick.major.size':3,'xtick.direction':'in','ytick.color':'black','ytick.labelsize':fs,
        'ytick.major.width':0.5,'ytick.major.size':3,'ytick.direction':'in','legend.fontsize':fs,'savefig.dpi':300,'savefig.transparent':True}
plt.rcParams.update(params)


#%% Plot

# Age
Age=2006-np.array([1949,1988,2000])

# Shortwave albedo
swa=np.array([0.085,0.138,0.141])

# Net radiation
nr=np.array([2.318,2.09,1.85])

plt.close('all'); 
fig,ax=plt.subplots(1,figsize=gu.cm2inch(8.5,5.5)); lw=1; ms=3.5; cl1=[0.4,0.4,0.4]; cl2=[0.4,0.6,1]; cl3=[1,0.25,0];
ax.plot(Age,swa,'-bo',color=cl2,mec=cl2,mfc=cl2,lw=lw,ms=ms)
ax.set(position=[0.14,0.16,0.7,0.81],ylim=[0,0.25],yticks=np.arange(0,0.5,0.05), \
       xlim=[0,65],xticks=np.arange(0,75,10),xlabel='Stand age, years',ylabel='Surface albedo')
#ax.yaxis.set_ticks_position('both'); 
#ax.xaxis.set_ticks_position('both')
ax.spines['top'].set_visible(False); 
ax.spines['right'].set_visible(False); 
#ax.spines['left'].set_visible(False); 
ax.tick_params(axis='y',direction='out',colors=cl2,length=3)
ax.tick_params(axis='x',direction='out',colors=cl1,length=3)
ax.spines['left'].set_color(cl2)
ax.spines['bottom'].set_color(cl1)
ax.xaxis.label.set_color(cl1)
ax.yaxis.label.set_color(cl2)

ax2=ax.twinx()
ax2.plot(Age,nr,'-bs',color=cl3,mec=cl3,mfc=cl3,lw=lw,ms=ms)
ax2.set_ylabel('Net radiation (GJ m$^{-2}$ year$^{-1}$)',color=cl3)
ax2.set(position=[0.14,0.16,0.7,0.81],ylim=[1.6,2.4],yticks=np.arange(1.6,2.6,0.2), \
       xlim=[0,65],xlabel='')
ax2.spines['left'].set_visible(False); 
ax2.spines['right'].set_color(cl3)
ax2.tick_params(direction='out',colors=cl3,length=3)


gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Harvest\Albedo_DF_Chronosequence_Jassal2009','png',900)
fig.savefig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Harvest\Albedo_DF_Chronosequence_Jassal2009.svg',format='svg',dpi=1200)

