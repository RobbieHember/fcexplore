
"""
INITIAL TREE DENSITY

"""

#%% IMPORT MODULES

import numpy as np
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
import fcgadgets.macgyver.utilities_general as gu
import fcexplore.psp.Processing.psp_utilities as utl_gp

#%% Set figure properties

gp=gu.SetGraphics('Manuscript')

#%% Import data

meta={}
meta['Paths']={}
meta['Paths']['DB']=r'C:\Users\rhember\Documents\Data\GroundPlots\PSP-NADB2'
meta=utl_gp.ImportParameters(meta)
d=gu.ipickle(meta['Paths']['DB'] + '\\Processed\\L2\\L2_BC.pkl')
sl=d['sobs'].copy()
del d

#%% Filter

ind=np.where( (sl['Lat']>0) & (sl['Lon']!=0) )[0]
for k in sl.keys():
    sl[k]=sl[k][ind]


#%% Stem density age response

bw=10; bin=np.arange(10,260,bw)

plt.close('all')

ikp=np.where(sl['ClimateClass']==1)[0]
N,mu,med,sig,se=gu.discres(sl['Age t0'][ikp],sl['N L t0'][ikp],bw,bin)
plt.plot(bin,med,'-bo')

ikp=np.where(sl['ClimateClass']==2)[0]
N,mu,med,sig,se=gu.discres(sl['Age t0'][ikp],sl['N L t0'][ikp],bw,bin)
plt.plot(bin,med,'-go')

ikp=np.where(sl['ClimateClass']==3)[0]
N,mu,med,sig,se=gu.discres(sl['Age t0'][ikp],sl['N L t0'][ikp],bw,bin)
plt.plot(bin,med,'-ro')

#%% Plot by BGC zone

vL=['N L t0']

u=np.unique(sl['Ecozone BC L1'])
u=u[u>0]
lab=np.array(['' for _ in range(u.size)],dtype=object)

d={}
for v in vL:
    d[v]={}
    d[v]['N']=np.zeros(u.size)
    d[v]['mu']=np.zeros(u.size)
    d[v]['sd']=np.zeros(u.size)
    d[v]['se']=np.zeros(u.size)

for i in range(u.size):
    lab[i]=utl_gp.lut_id2cd(meta,'Ecozone BC L1',u[i])
    for v in vL:
        ind=np.where( (sl['Ecozone BC L1']==u[i]) & (sl['Age t0']>30) & (sl['Age t0']<60) )[0]
        d[v]['N'][i]=ind.size
        d[v]['mu'][i]=np.nanmean(sl[v][ind])
        d[v]['sd'][i]=np.nanstd(sl[v][ind])
        #d[v]['se'][i]=np.nanstd(sl[v][ind])/np.sqrt(ind[0].size)

# Remove nans
ind=np.where(np.isnan(d[v]['mu'])==False)[0]
u=u[ind]
lab=lab[ind]
for v in d:
    for k in d[v].keys():
        d[v][k]=d[v][k][ind]

# Put in order
ord=np.argsort(d['N L t0']['mu'])
uo=u[ord]
lab=np.flip(lab[ord])
for v in d:
    for k in d[v].keys():
        d[v][k]=np.flip(d[v][k][ord])

# Plot
cl=np.array([[0.27,0.49,0.77],[0.25,0.14,0.05],[1,0.5,0],[0.45,0.75,1],[0.6,1,0]])
plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,5.5))
ax.bar(np.arange(u.size),d['N L t0']['mu'],facecolor=cl[0,:],label='Stemwood')
ax.set(position=[0.09,0.15,0.9,0.8],xlim=[-0.5,u.size-0.5],ylim=[0,2600],yticks=np.arange(0,3200,500),xticks=np.arange(u.size),
       xticklabels=lab,ylabel='Tree density (stems ha$^{-1}$)',xlabel='Biogeoclimatic zone')
#plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Biomass\InitialTreeDensity_ByGBCZone','png',900)


