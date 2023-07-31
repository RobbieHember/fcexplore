
#%% Import modules

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import fcgadgets.macgyver.utilities_general as gu
import fcgadgets.macgyver.utilities_gis as gis
import fcgadgets.bc1ha.bc1ha_utilities as u1ha
import fcexplore.psp.Processing.psp_utilities as utl_gp

#%% Import data

gp=gu.SetGraphics('Manuscript')

lut_1ha=u1ha.Import_BC1ha_LUTs()

zLC2=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\lc2.tif')
zLC4=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\lc4.tif')
zLC5=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\lc5.tif')

zH=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\VEG_CONSOLIDATED_CUT_BLOCKS_SP_All.tif')
zWF=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\PROT_HISTORICAL_FIRE_POLYS_SP_All.tif')

zBGCZ=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\becz.tif')
lutBGC=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\becz_lut.xlsx')

zSPH=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\sphlive.tif')
zCC=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\crownc.tif')

# Ground plots
meta={}
meta['Paths']={}
meta['Paths']['DB']=r'C:\Users\rhember\Documents\Data\GroundPlots\PSP-NADB2'
meta['Paths']['Figs']=r'C:\Users\rhember\OneDrive - Government of BC\Figures\Biomass'
meta=utl_gp.ImportParameters(meta)
d=gu.ipickle(meta['Paths']['DB'] + '\\Processed\\L2\\L2_BC.pkl')
sl=d['sobs'].copy()
del d

#%% Plot tree density class by BGC zone

u=np.unique(zBGCZ['Data'])
u=u[(u>0) & (u<255)]
lab=np.array(['' for _ in range(u.size)],dtype=object)

vL=['Dense','Open','Sparse','All']

d={}
for v in vL:
    d[v]={}
    d[v]['N']=np.zeros(u.size)
    d[v]['mu']=np.zeros(u.size)
    d[v]['sd']=np.zeros(u.size)
    d[v]['se']=np.zeros(u.size)

for i in range(u.size):
    ind=np.where(lutBGC['VALUE']==u[i])[0]
    lab[i]=lutBGC['ZONE'][ind[0]]
    for v in vL:
        if v=='All':
            ind=np.where( (zBGCZ['Data']==u[i]) & (zLC2['Data']==lut_1ha['lc2']['Treed']) & (zH['Data']==0) & (zWF['Data']==0) )
            d[v]['N'][i]=ind[0].size
        else:
            ind=np.where( (zBGCZ['Data']==u[i]) & (zLC2['Data']==lut_1ha['lc2']['Treed']) & (zLC5['Data']==lut_1ha['lc5'][v]) & (zH['Data']==0) & (zWF['Data']==0) )
            d[v]['N'][i]=ind[0].size

# Calculate percent
for v in d:
    d[v]['P']=d[v]['N']/(d['Dense']['N']+d['Sparse']['N']+d['Open']['N'])*100

# Put in order
ord=np.argsort(d['Dense']['P']+d['Open']['P'])
uo=u[ord]
lab=np.flip(lab[ord])
for v in d:
    for k in d[v].keys():
        d[v][k]=np.flip(d[v][k][ord])

# Plot
plt.close('all'); cl=np.array([[0,0.3,0],[0.5,0.7,0.3],[0.8,0.8,0.5],[0.45,0.75,1],[0.6,1,0]])
fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,6.5));
ax.bar(np.arange(u.size),d['Dense']['P'],facecolor=cl[0,:],label='')
ax.bar(np.arange(u.size),d['Open']['P'],bottom=d['Dense']['P'],facecolor=cl[1,:],label='')
ax.bar(np.arange(u.size),d['Sparse']['P'],bottom=d['Dense']['P']+d['Open']['P'],facecolor=cl[2,:],label='')
ax.set(position=[0.08,0.075,0.9,0.91],yticks=np.arange(0,110,10),xticks=np.arange(u.size),xticklabels=lab,
       ylabel='Percent (%)',xlim=[-0.5,u.size-0.5],ylim=[0,100])
#plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Tree Density\DensityClassByGBCZone_WithoutHWF','png',900)

#%% Tree Density Class

vL1=['Frequency (%)','Crown cover from VRI (%)'] # ,'Tree density from VRI (stems/ha)'
vL2=['Dense','Open','Sparse']

d={}
for v1 in vL1:
    d[v1]={}
    for v2 in vL2:
        d[v1][v2]=0.0

# Frequency
vn='Frequency (%)'
n=0
for v2 in vL2:
    ind=np.where( (zLC2['Data']==lut_1ha['lc2']['Treed']) & (zLC5['Data']==lut_1ha['lc5'][v2]) ); d[vn][v2]=ind[0].size; n=n+ind[0].size

for k in d[vn].keys():
    d[vn][k]=np.round(d[vn][k]/n*100,decimals=1)

# # SPH by density class
# vn='Tree density from VRI (stems/ha)'
# for v2 in vL2:
#     ind=np.where( (zLC2['Data']==lut_1ha['lc2']['Treed']) & (zLC5['Data']==lut_1ha['lc5'][v2]) ); d[vn][v2]=np.round(np.nanmean(zSPH['Data'][ind]),decimals=0)

## CC by density class
vn='Crown cover from VRI (%)'
for v2 in vL2:
    ind=np.where( (zLC2['Data']==lut_1ha['lc2']['Treed']) & (zLC5['Data']==lut_1ha['lc5'][v2]) ); d[vn][v2]=np.round(np.nanmean(zCC['Data'][ind]),decimals=0)

# PSP stuff moved to PSP summary script

#%% OAF1 (systematic effect
vn='OAF1 (%)'
d[vn]={'Sparse':45,'Open':30,'Dense':15}

#%% Add labels and convert to dataframe

#lab=[' (61-100%)',' (26-60%)',' (10-25%)']
df=pd.DataFrame(d)
df.index.name='Density Class'

df.to_excel(r'C:\Users\rhember\Documents\Code_Python\fcgadgets\cbrunner\Parameters\Parameters_TreeDensityClassStatistics.xlsx')

#%% Plot

# lab=['Sparse\n(10-25%)','Open\n(26-60%)','Dense\n(61-100%)']

# plt.close('all'); cl=np.array([[0.8,0.8,0.8],[0.5,0.7,0.3],[0.8,0.8,0.5],[0.45,0.75,1],[0.6,1,0]])
# fig,ax=plt.subplots(1,figsize=gu.cm2inch(12,6.5));
# ax.bar(np.arange(mu.size),mu,facecolor=cl[0,:],label='')
# ax.set(position=[0.08,0.13,0.9,0.82],xticks=np.arange(mu.size),xticklabels=lab,ylabel='Tree biomass (tC ha$^{-1}$)',xlim=[-0.5,mu.size-0.5])
# #plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
# ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])
# gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Tree Density\TreeBiomass By DensityClass','png',900)

#%% Crown cover - comparison from various sources

plt.close('all')

ikp=np.where(sl['Flag Plot Types To Keep']==1)[0]
kd=gu.ksdensity(sl['Crown Closure (VRI)'][ikp])
plt.plot(kd,'b-')
print(np.mean(sl['Crown Closure (VRI)'][ikp]))

ikp=np.where(sl['Flag Keep CMI+NFI']==1)[0]
kd=gu.ksdensity(sl['Crown Closure (VRI)'][ikp])
plt.plot(kd,'g-')
print(np.mean(sl['Crown Closure (VRI)'][ikp]))

zLC2=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\lc2.tif')
zCC=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\crownc.tif')
ikp=np.where( (zLC2['Data'][0::10,0::10].flatten()==4) )[0]
kd=gu.ksdensity(zCC['Data'][0::10,0::10].flatten()[ikp])
plt.plot(kd,'c--')
print(np.mean(zCC['Data'][0::10,0::10].flatten()[ikp]))