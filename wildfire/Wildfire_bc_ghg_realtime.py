'''
BC WILDFIRE GHG EMISSIONS - PRELIMINARY (REAL-TIME) ESTIMATES

# Emission factors provided by CFS (from NIR average)

'''
#%% Import modules
import numpy as np
import time
import matplotlib.pyplot as plt
import pandas as pd
import fcgadgets.macgyver.utilities_general as gu
import fcgadgets.macgyver.utilities_gis as gis
import fcgadgets.bc1ha.bc1ha_utilities as u1ha
import fcgadgets.macgyver.utilities_fcs_graphs as ufcs

#%% Import paths and look-up-tables
meta=u1ha.Init()
gp=gu.SetGraphics('Manuscript')

#%% Define reference grid
z=u1ha.Import_Raster(meta,[],['refg','lc2','lcc1_c','ezcan'])

nir=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\Emission Reduction Projections\Summary of Reporting Initiatives.xlsx')

#%% Calculate time series of emissions
flg=0
if flg==1:
    uE=np.unique(z['ezcan']['Data'][z['ezcan']['Data']>0])
    tv=np.arange(1917,2024,1)
    A=np.zeros((tv.size,uE.size))
    E=np.zeros((tv.size,uE.size))
    for iT in range(tv.size):
        print(tv[iT])
        zF=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PROT_HISTORICAL_FIRE_POLYS_SP\\FIRE_YEAR_' + str(tv[iT]) + '.tif')['Data']
        #iFire=np.where( (zF>0) & (z['lc2']['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_2']['T']) )
        iFire=np.where( (zF>0) & (z['lcc1_c']['Data']==meta['LUT']['Derived']['lcc1']['Forest']) )
        zE=z['ezcan']['Data'][iFire]
        for iE in range(uE.size):
            iEF=np.where( (zE==uE[iE]) )
            A[iT,iE]=iEF[0].size
            iLUT=np.where(meta['LUT']['Raw']['ezcan']['ID']==uE[iE])[0]
            E[iT,iE]=iEF[0].size*meta['LUT']['Raw']['ezcan']['Emissions per hectare (tCO2e/ha) From Bhudson 23'][iLUT]

    d={'tv':tv,'A':A,'E':E}
    #gu.opickle(r'C:\Users\rhember\OneDrive - Government of BC\Requests\2023-07-21 Wildfire Emissions\Emissions_niref_lc2.pkl',d)
    gu.opickle(r'C:\Users\rhember\OneDrive - Government of BC\Requests\2023-07-21 Wildfire Emissions\Emissions_niref_lcc1.pkl',d)
else:
    #d=gu.ipickle(r'C:\Users\rhember\OneDrive - Government of BC\Requests\2023-07-21 Wildfire Emissions\Emissions_niref_lc2.pkl')
    d=gu.ipickle(r'C:\Users\rhember\OneDrive - Government of BC\Requests\2023-07-21 Wildfire Emissions\Emissions_niref_lcc1.pkl')

#%% Plot area affected
plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,6))
ax.bar(d['tv'],np.sum(d['A'],axis=1)/1e6,0.75,facecolor=[0.7,0,0])
ax.set(xticks=np.arange(1920,2040,10),ylabel='Affected area (Million hectares/year)',xlabel='Time, years',yticks=np.arange(0,1.8,0.2),xlim=[1919.25,2023.75])
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=1.5)
plt.tight_layout()
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Wildfire\RealTimeWildfireGHG_AreaAffected_lcc1c','png',900)

#%% Plot GHG emissions
iT=np.where( (d['tv']>=2017) & (d['tv']<=2023) )[0]
txt='Total over 2017-2023 = ' + str(int(np.sum(d['E'][iT,:])/1e6)) + '\n(tonnes carbon dioxide)'
plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,6))
ax.bar(d['tv'],np.sum(d['E'],axis=1)/1e6,0.75,facecolor=[0.7,0,0],label='BC estimate')
ax.plot(nir['Year'],nir['BC FLFL Wildfire (NIR22)'],'ko',mfc='w',label='NIR 2022')
ax.set(xticks=np.arange(1920,2040,10),ylabel='GHG emissions (MtCO2e/year)',xlabel='Time, years',yticks=np.arange(0,260,20),xlim=[1919.25,2023.75])
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=1.5)
ax.text(1987,100,txt,fontsize=7)
plt.legend(frameon=False,loc='upper left',facecolor=[1,1,1],labelspacing=0.25,ncol=1)
plt.tight_layout()
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Wildfire\RealTimeWildfireGHG_GHGEmissions_EF_lcc1c_WithNIR','png',900)
