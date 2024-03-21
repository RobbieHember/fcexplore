"""
PSP - NET GROWTH OF NON-STEMWOOD BIOMASS FROM STEMWOOD BIOMASS
"""
#%% Import modules
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import copy
from scipy import stats
import statsmodels.formula.api as smf
from scipy.optimize import curve_fit
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.bc1ha.bc1ha_util as u1ha
import fcexplore.field_plots.Processing.psp_util as ugp
import fcexplore.field_plots.Processing.psp_plot as pgp
gp=gu.SetGraphics('Manuscript')

#%% Import data
meta=u1ha.Init()
meta,gpt,soc=ugp.ImportGroundPlotData(meta,type='Stand',include_soil='True')

gpt['Cbk Net']=gpt['Cbk G Surv']+gpt['Cbk G Recr']-gpt['Cbk Mort']
gpt['Cbr Net']=gpt['Cbr G Surv']+gpt['Cbr G Recr']-gpt['Cbr Mort']
gpt['Cf Net']=gpt['Cf G Surv']+gpt['Cf G Recr']-gpt['Cf Mort']

#%% Import Raster grids
# zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
# iGrd=gis.GetGridIndexToPoints(zRef,gpt['X'],gpt['Y'])
# vList=['tdc_vri15','lc_comp1_2019','bgcz']
# roi={'points':{'x':soc['x'],'y':soc['y']}}
# zS=u1ha.Import_Raster(meta,roi,vList)

#%% Allometric equation (Landsberg and Sands 2011)
def fun(x,b0,b1,b2):
	y=x[:,0]*b0+(b1-b0)*np.exp(-b2*x[:,1])
	return y

#%% Plot relationships

b0=[0.04,1.01,0.08]
Age=np.arange(1,250)
bw=20; bin=np.arange(5,250,bw)

plt.close('all'); fig,ax=plt.subplots(1,3,figsize=gu.cm2inch(16,5.8))
# Bark
ikp=np.where( (gpt['PTF CN']==1) & (gpt['Ecozone BC L1']!=meta['LUT']['GP']['Ecozone BC L1']['CWH']) & (gpt['Csw Net']>-2) & (gpt['Cbk Net']>-2) & (gpt['Age Mean t0']>0) )[0]
ax[0].plot([0,250],[0,0],'k-')
y=gpt['Cbk Net']
x=np.column_stack((gpt['Csw Net'],gpt['Age Mean t0']))
N,mu,med,sig,se=gu.discres(gpt['Age Mean t0'][ikp],y[ikp]/gpt['Csw Net'][ikp],bw,bin)
p,pcov=curve_fit(fun,x[ikp,:],y[ikp],b0)
yhat=p[0]+(p[1]-p[0])*np.exp(-p[2]*Age)
ax[0].plot(Age,yhat,'r-')
ax[0].plot(bin,med,'bo',mec='w',mfc=[0,0,0],ms=5)
ax[0].set(xticks=np.arange(0,550,25),xlabel='Age, years',ylabel='Bark-stemwood ratio',xlim=[0,250],ylim=[-0.5,1.2])
ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].tick_params(length=meta['Graphics']['gp']['tickl'])
# Branch
ikp=np.where( (gpt['PTF CN']==1) & (gpt['Ecozone BC L1']!=meta['LUT']['GP']['Ecozone BC L1']['CWH']) & (gpt['Csw Net']>-2) & (gpt['Cbr Net']>-2) & (gpt['Age Mean t0']>0) )[0]
ax[1].plot([0,250],[0,0],'k-')
y=gpt['Cbr Net']
x=np.column_stack((gpt['Csw Net'],gpt['Age Mean t0']))
N,mu,med,sig,se=gu.discres(gpt['Age Mean t0'][ikp],y[ikp]/gpt['Csw Net'][ikp],bw,bin)
p,pcov=curve_fit(fun,x[ikp,:],y[ikp],b0)
yhat=p[0]+(p[1]-p[0])*np.exp(-p[2]*Age)
ax[1].plot(Age,yhat,'r-')
ax[1].plot(bin,med,'bo',mec='w',mfc=[0,0,0],ms=5)
ax[1].set(xticks=np.arange(0,550,25),xlabel='Age, years',ylabel='Branch-stemwood ratio',xlim=[0,250],ylim=[-0.5,1.2])
ax[1].yaxis.set_ticks_position('both'); ax[1].xaxis.set_ticks_position('both'); ax[1].tick_params(length=meta['Graphics']['gp']['tickl'])
# Foliage
ikp=np.where( (gpt['PTF CN']==1) & (gpt['Ecozone BC L1']!=meta['LUT']['GP']['Ecozone BC L1']['CWH']) & (gpt['Csw Net']>-2) & (gpt['Cf Net']>-2) & (gpt['Age Mean t0']>0) )[0]
ax[2].plot([0,250],[0,0],'k-')
y=gpt['Cf Net']
x=np.column_stack((gpt['Csw Net'],gpt['Age Mean t0']))
N,mu,med,sig,se=gu.discres(gpt['Age Mean t0'][ikp],y[ikp]/gpt['Csw Net'][ikp],bw,bin)
p,pcov=curve_fit(fun,x[ikp,:],y[ikp],b0)
yhat=p[0]+(p[1]-p[0])*np.exp(-p[2]*Age)
ax[2].plot(Age,yhat,'r-')
ax[2].plot(bin,med,'bo',mec='w',mfc=[0,0,0],ms=5)
ax[2].set(xticks=np.arange(0,550,25),xlabel='Age, years',ylabel='Foliage-stemwood ratio',xlim=[0,250],ylim=[-0.5,1.2])
ax[2].yaxis.set_ticks_position('both'); ax[2].xaxis.set_ticks_position('both'); ax[2].tick_params(length=meta['Graphics']['gp']['tickl'])
plt.tight_layout()
gu.axletters(ax,plt,0.04,0.92,FontColor=meta['Graphics']['gp']['cla'],LetterStyle=meta['Graphics']['Modelling']['AxesLetterStyle'],FontWeight=meta['Graphics']['Modelling']['AxesFontWeight'])
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Allometric Equations\QA_Allometry_RelationshipsInterior','png',900)

#%% Branch
ikp=np.where( (gpt['PTF CN']==1) & (gpt['Ecozone BC L1']!=meta['LUT']['GP']['Ecozone BC L1']['CWH']) & (gpt['Csw Net']>-2) & (gpt['Cbr Net']>-2) & (gpt['Age Mean t0']>0) )[0]

y=gpt['Cbr Net']
x=np.column_stack((gpt['Csw Net'],gpt['Age Mean t0']))
b0=[0.04,1.01,0.08]
p,pcov=curve_fit(fun,x[ikp,:],y[ikp],b0)
Age=np.arange(1,250)

bw=20; bin=np.arange(5,250,bw)
N,mu,med,sig,se=gu.discres(gpt['Age Mean t0'][ikp],y[ikp]/gpt['Csw Net'][ikp],bw,bin)

plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(7.8,7))
ax.plot([0,250],[0,0],'k-')
#p=[0.04,1.01,0.08]
yhat=p[0]+(p[1]-p[0])*np.exp(-p[2]*Age)
ax.plot(Age,yhat,'r-')
plt.plot(bin,med,'bo',mec='w',mfc=[0,0,0],ms=5)
ax.set(xticks=np.arange(0,550,25),xlabel='Age, years',ylabel='Branch-stemwood ratio',xlim=[0,250],ylim=[-1,1.5])
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])

#%%

plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(7.8,7))
ax.plot(Age,yhat,'k-',linewidth=0.75,label='Default model')

ax.set(position=[0.11,0.11,0.88,0.88],xlim=[0,500],xticks=np.arange(0,550,50),xlabel='Age, years',ylabel='Annual probability of breakup')
ax.legend(loc='upper left',bbox_to_anchor=(0.06,0.92),frameon=False,facecolor='w')



#%% Foliage

ikp=np.where( (gpt['PTF CN']==1) & (gpt['Ecozone BC L1']!=meta['LUT']['GP']['Ecozone BC L1']['CWH']) & (gpt['Csw Net']>-2) & (gpt['Cf Net']>-2) & (gpt['Age Mean t0']>0) )[0]

bw=10; bin=np.arange(0,250,bw)
N,mu,med,sig,se=gu.discres(gpt['Age Mean t0'][ikp],gpt['Cf Net'][ikp],bw,bin)
plt.close('all');plt.plot(bin,med,'bo')

bw=0.25; bin=np.arange(0,5,bw)
N,mu,med,sig,se=gu.discres(gpt['Csw Net'][ikp],gpt['Cf Net'][ikp],bw,bin)
plt.close('all');plt.plot(bin,med,'bo')
