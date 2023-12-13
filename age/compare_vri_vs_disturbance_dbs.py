#%% Import modules
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn import linear_model
import statsmodels.formula.api as smf
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.bc1ha.bc1ha_util as u1ha
gp=gu.SetGraphics('Manuscript')

#%% Import data
meta=u1ha.Init()
z0=u1ha.Import_Raster(meta,[],['refg','lc_comp1_2019','age_vri','age_vri02','harv_yr_cc','fire_yr','ibm_yr'],'Extract Grid')

#%%

ikp=np.where( (z0['lc_comp1_2019']==1) & (z0['age_vri']>=0) & (z0['harv_yr_cc']>0) )
x=2022-z0['age_vri'][ikp].flatten()
y=z0['harv_yr_cc'][ikp].flatten()
d=x-y

plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(10,6))
ax.hist(d,np.arange(-200,25,5))
ax.set(xlabel='$\Delta$ time since stand-replacing harvest (VRI minus CCDB)',ylabel='Frequency')
txt='Mean difference = ' + str(np.round(np.mean(d))) + '\nMedian difference = ' + str(np.round(np.median(d)))
ax.text(-150,4000000,txt)

ind=np.where(np.abs(d)>5)[0]
print(ind.size/ikp[0].size*100)



#%% Subsample
for k in z0.keys():
    z0[k]=z0[k][0::10,0::10]

#%% Harvest

ikp=np.where( (z0['lc_comp1_2019']==1) & (z0['harv_yr_cc']>0) )
#plt.close('all');plt.plot(2023-z0['age_vri'][ikp],z0['harv_yr_cc'][ikp],'.')

x=2023-z0['age_vri'][ikp]
y=z0['harv_yr_cc'][ikp]
binX,binY,Z=gu.ScatterDensity(x,y,[1900,2023],[1900,2023],60,60)

plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(10,10))
ax.matshow(np.flip(Z.T,axis=0),clim=[0,200],extent=[np.min(binX),np.max(binX),np.min(binY),np.max(binY)],aspect='auto')
ax.set(position=[0.1,0.1,0.8,0.8],xlabel='VRI',ylabel='Wildfire DB')#,xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
#ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

#%% WIldfire

ikp=np.where( (z0['lc_comp1_2019']==1) & (z0['fire_yr']>0) )
#plt.close('all');plt.plot(2023-z0['age_vri'][ikp],z0['harv_yr_cc'][ikp],'.')

x=2023-z0['age_vri'][ikp]
y=z0['fire_yr'][ikp]
binX,binY,Z=gu.ScatterDensity(x,y,[1900,2023],[1900,2023],60,60)

plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(10,10))
ax.matshow(np.flip(Z.T,axis=0),clim=[0,200],extent=[np.min(binX),np.max(binX),np.min(binY),np.max(binY)],aspect='auto')
ax.set(position=[0.1,0.1,0.8,0.8],xlabel='VRI',ylabel='Wildfire DB')#,xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
#ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])


#%% IBM
ikp=np.where( (z0['lc_comp1_2019']==1) & (z0['ibm_yr']>0) )
#plt.close('all');plt.plot(2023-z0['age_vri'][ikp],z0['harv_yr_cc'][ikp],'.')

x=2023-z0['age_vri'][ikp]
y=z0['ibm_yr'][ikp]
binX,binY,Z=gu.ScatterDensity(x,y,[1900,2023],[1900,2023],60,60)

plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(10,10))
ax.matshow(np.flip(Z.T,axis=0),clim=[0,200],extent=[np.min(binX),np.max(binX),np.min(binY),np.max(binY)],aspect='auto')
ax.set(position=[0.1,0.1,0.8,0.8],xlabel='VRI',ylabel='AOS DB')#,xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
#ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])

