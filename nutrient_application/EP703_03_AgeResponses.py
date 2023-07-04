'''============================================================================

Age responses at EP703

============================================================================'''

#%% Import modules

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyproj
import geopandas as gpd
from scipy.optimize import curve_fit
import scipy.io as spio
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats, linalg
from fcgadgets.macgyver import utilities_general as gu

#%% Set figure properties

gp=gu.SetGraphics('Manuscript')

#%% Import paths

PathProject=r'C:\Users\rhember\Documents\Data\EP703'
PathFigures=r'C:\Users\rhember\OneDrive - Government of BC\Figures\EP703'
PathManuscriptSI=r'G:\My Drive\Manuscripts\CCIPB_FertilizationGHGBenefit_EP703\SI'

#%% Import data

# Import site information
dS=gu.ReadExcel(PathProject + '\\EP703_SiteInfo.xlsx')

# Unique sites
uSite=np.unique(dS['ID_Instal'])

# Import stand-level data
sobs=gu.ipickle(PathProject + '\\Processed\\EP703_SL.pkl')

# Model
tipsy={}
tipsy['FD']=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\tipsy_fdc.xlsx')
tipsy['HW']=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\tipsy_hw31.xlsx')

#%% Derived variables and adjustments

# Create an adjusted value of TSF_t0 that will ensure time intervals are not double
# counting the year of measurement
sobs['TSF_t0_adj']=sobs['TSF_t0']+1

sobs['Bsw_Net_RGR']=(np.log(sobs['Bsw_t1'])-np.log(sobs['Bsw_t0']))/sobs['DT']

# Initial biomass
u=np.unique(np.column_stack((sobs['ID_Site'],sobs['ID_Plot'])).astype(float),axis=0)
sobs['Bsw_Init']=np.nan*np.ones(sobs['Bsw_t0'].size)
for i in range(u.shape[0]):
    ind=np.where( (sobs['ID_Site']==u[i,0]) & (sobs['ID_Plot']==u[i,1]) )[0]
    atsf=np.abs(sobs['TSF_t0'][ind])
    if ind.size!=0:
        ind2=np.where( atsf==np.min(atsf) )[0]
        sobs['Bsw_Init'][ind]=sobs['Bsw_t0'][ind[ind2[0]]]

# Relative change in stand density (%/yr)
sobs['dN_rel']=sobs['N_Net']/sobs['N_t0']*100

sobs['RGR']=(np.log(sobs['Bsw_t1'])-np.log(sobs['Bsw_t0']))/sobs['DT']

sobs['Bsw_G_ave']=sobs['Bsw_G']/sobs['N_t1']*1000

#%% Age response

tv=np.arange(1,100,1)
uSpc=['FD','HW']
uDose=[225,450]
uVar=['dN_rel','RGR','Age_t0','N_t0','N_Net','BA_t0','H_obs_t0','Bsw_t0','Bsw_Net','Bsw_G','Bsw_M']
rts={}
for iSpc in range(len(uSpc)):
    rts[uSpc[iSpc]]={}
    for iDose in range(len(uDose)):
        rts[uSpc[iSpc]][uDose[iDose]]={}
        rts[uSpc[iSpc]][uDose[iDose]]['y_c']={}
        rts[uSpc[iSpc]][uDose[iDose]]['y_f']={}
        for iVar in range(len(uVar)):
            rts[uSpc[iSpc]][uDose[iDose]]['y_c'][uVar[iVar]]=np.nan*np.ones((tv.size,sobs['Year_t0'].size))
            rts[uSpc[iSpc]][uDose[iDose]]['y_f'][uVar[iVar]]=np.nan*np.ones((tv.size,sobs['Year_t0'].size))
            cnt_c=0;
            cnt_f=0;
            for iYr in range(sobs['Age_t0'].size):
                #if (sobs['ID_Site'][iYr]==29) | (sobs['ID_Site'][iYr]==71) | (sobs['ID_Site'][iYr]==77):
                #    continue
                #if sobs['Age_t0'][iYr]>80:
                #    continue
                it=np.where( (tv>=sobs['Age_t0'][iYr]) & (tv<=sobs['Age_t1'][iYr]) )[0]
                if (sobs['Spc1_ID_t0'][iYr]==uSpc[iSpc]) & (sobs['N_Dose'][iYr]==0):
                    #it=np.where( (tv>=sobs['Year_t0'][iYr]) & (tv<=sobs['Year_t1'][iYr]) )[0]
                    rts[uSpc[iSpc]][uDose[iDose]]['y_c'][uVar[iVar]][it,cnt_c]=sobs[uVar[iVar]][iYr];
                    cnt_c=cnt_c+1
                if (sobs['Spc1_ID_t0'][iYr]==uSpc[iSpc]) & (sobs['N_Dose'][iYr]==uDose[iDose]):
                    #it=np.where( (tv>=sobs['Year_t0'][iYr]) & (tv<=sobs['Year_t1'][iYr]) )[0]
                    rts[uSpc[iSpc]][uDose[iDose]]['y_f'][uVar[iVar]][it,cnt_f]=sobs[uVar[iVar]][iYr];
                    cnt_f=cnt_f+1

#%% Age

Dose=225

ms=3
plt.close('all'); fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(15.5,10))

y=np.nanmean(rts['FD'][Dose]['y_c']['Bsw_Net'],axis=1)
ax[0,0].plot(tv,y,'-bo',ms=ms,label='Ground plots')
ax[0,0].plot(tipsy['FD']['Age'],gu.movingave(0.5*0.5*tipsy['FD']['G Total'],15,'center'),'r-',label='TIPSY')
ax[0,0].set(ylabel='Net growth (MgC ha$^-$$^1$ yr$^-$$^1$)',xlabel='Age, years',xlim=[0,100])
ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=gp['tickl'])


y=np.nanmean(rts['HW'][Dose]['y_c']['Bsw_Net'],axis=1)
ax[0,1].plot(tv,y,'-bo',ms=ms)
ax[0,1].plot(tipsy['HW']['Age'],gu.movingave(0.5*0.5*tipsy['HW']['G Total'],15,'center'),'r-')
ax[0,1].set(ylabel='Net growth (MgC ha$^-$$^1$ yr$^-$$^1$)',xlabel='Age, years',xlim=[0,100])
ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=gp['tickl'])

y=np.nanmean(rts['FD'][Dose]['y_c']['Bsw_t0'],axis=1)
ax[1,0].plot(tv,y,'-bo',ms=ms)
ax[1,0].plot(tipsy['FD']['Age'],0.5*0.5*tipsy['FD']['V Total'],'r-')
ax[1,0].set(ylabel='Biomass (MgC ha$^-$$^1$)',xlabel='Age, years',xlim=[0,100])
ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=gp['tickl'])

y=np.nanmean(rts['HW'][Dose]['y_c']['Bsw_t0'],axis=1)
ax[1,1].plot(tv,y,'-bo',ms=ms)
ax[1,1].plot(tipsy['HW']['Age'],0.5*0.5*tipsy['HW']['V Total'],'r-')
ax[1,1].set(ylabel='Biomass (MgC ha$^-$$^1$)',xlabel='Age, years',xlim=[0,100])
ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=gp['tickl'])

# Add modified
b=-0.1
ymn=0.45
ymx=1.7
A_inf=50
fG=ymn+(ymx-ymn)*1/(1+np.exp(-b*(tipsy['FD']['Age']-A_inf)))
ax[0,0].plot(tipsy['FD']['Age'],fG*gu.movingave(0.5*0.5*tipsy['FD']['G Total'],15,'center'),'r--',label='TIPSY (modified)')
ax[1,0].plot(tipsy['FD']['Age'],np.cumsum(fG*gu.movingave(0.5*0.5*tipsy['FD']['G Total'],15,'center')),'r--',label='TIPSY Modified')

fG=ymn+(ymx-ymn)*1/(1+np.exp(-b*(tipsy['HW']['Age']-A_inf)))
ax[0,1].plot(tipsy['HW']['Age'],fG*gu.movingave(0.5*0.5*tipsy['HW']['G Total'],15,'center'),'r--',label='TIPSY Modified')
ax[1,1].plot(tipsy['HW']['Age'],np.cumsum(fG*gu.movingave(0.5*0.5*tipsy['HW']['G Total'],15,'center')),'r--',label='TIPSY Modified')

leg=ax[0,0].legend(loc='upper right',frameon=False,facecolor=None,edgecolor='w')
gu.axletters(ax,plt,0.03,0.9,FontColor=gp['cla'],Labels=['Douglas-fir','Hemlock','Douglas-fir','Hemlock'],LetterStyle='Caps',FontWeight='Bold')
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\EP703\EP703_AgeResponseObsCompWithTIPSY','png',900)

#%% Time series analysis

tv=np.arange(1971,2012,1)
uSpc=['FD','HW']
uDose=[225,450]
uVar=['dN_rel','RGR','Age_t0','N_t0','N_Net','BA_t0','H_obs_t0','Bsw_t0','Bsw_Net','Bsw_G','Bsw_M']
rts={}
for iSpc in range(len(uSpc)):
    rts[uSpc[iSpc]]={}
    for iDose in range(len(uDose)):
        rts[uSpc[iSpc]][uDose[iDose]]={}
        rts[uSpc[iSpc]][uDose[iDose]]['y_c']={}
        rts[uSpc[iSpc]][uDose[iDose]]['y_f']={}
        for iVar in range(len(uVar)):
            rts[uSpc[iSpc]][uDose[iDose]]['y_c'][uVar[iVar]]=np.nan*np.ones((tv.size,sobs['Year_t0'].size))
            rts[uSpc[iSpc]][uDose[iDose]]['y_f'][uVar[iVar]]=np.nan*np.ones((tv.size,sobs['Year_t0'].size))
            cnt_c=0;
            cnt_f=0;
            for iYr in range(sobs['Year_t0'].size):
                #if (sobs['ID_Site'][iYr]==29) | (sobs['ID_Site'][iYr]==71) | (sobs['ID_Site'][iYr]==77):
                #    continue
                #if sobs['Age_t0'][iYr]>80:
                #    continue
                it=np.where( (tv>=sobs['Year_t0'][iYr]) & (tv<=sobs['Year_t1'][iYr]) )[0]
                if (sobs['Spc1_ID_t0'][iYr]==uSpc[iSpc]) & (sobs['N_Dose'][iYr]==0):
                    #it=np.where( (tv>=sobs['Year_t0'][iYr]) & (tv<=sobs['Year_t1'][iYr]) )[0]
                    rts[uSpc[iSpc]][uDose[iDose]]['y_c'][uVar[iVar]][it,cnt_c]=sobs[uVar[iVar]][iYr];
                    cnt_c=cnt_c+1
                if (sobs['Spc1_ID_t0'][iYr]==uSpc[iSpc]) & (sobs['N_Dose'][iYr]==uDose[iDose]):
                    #it=np.where( (tv>=sobs['Year_t0'][iYr]) & (tv<=sobs['Year_t1'][iYr]) )[0]
                    rts[uSpc[iSpc]][uDose[iDose]]['y_f'][uVar[iVar]][it,cnt_f]=sobs[uVar[iVar]][iYr];
                    cnt_f=cnt_f+1

#%% Time series

Dose=225

plt.close('all')
y=np.nanmean(rts['FD'][Dose]['y_c']['Bsw_Net'],axis=1)
plt.plot(tv,y,'-gs')
plt.plot(1950+tipsy['FD']['Age'],gu.movingave(0.5*0.5*tipsy['FD']['G Total'],15,'center'),'r-')

plt.close('all')
y=np.nanmean(rts['HW'][Dose]['y_c']['Bsw_Net'],axis=1)
plt.plot(tv,y,'-gs')
plt.plot(1942+tipsy['HW']['Age'],gu.movingave(0.5*0.5*tipsy['HW']['G Total'],15,'center'),'r-')


