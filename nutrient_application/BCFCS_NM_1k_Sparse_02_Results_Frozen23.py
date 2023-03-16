
#%% Import modules

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import gc as garc
import copy
import time
import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib import animation
import matplotlib as mpl
from fcgadgets.macgyver import utilities_general as gu
from fcgadgets.cbrunner import cbrun_utilities as cbu
from fcgadgets.macgyver import utilities_demo as udem
from fcgadgets.macgyver import utilities_fcs_graphs as utl

#%% Import data

metaC=gu.ipickle(r'D:\Data\FCI_Projects\BCFCS_NM_1k_Sparse\Inputs\Metadata.pkl')
metaC['Paths']['Figures']=r'C:\Users\rhember\OneDrive - Government of BC\Figures\BCFCS_NM_1k_Sparse'
tvC=np.arange(metaC['Project']['Year Start Saving'],metaC['Project']['Year End']+1,1)
geosC=gu.ipickle(metaC['Paths']['Project'] + '\\Geospatial\\geos.pkl')
ailC=gu.ipickle(metaC['Paths']['Project'] + '\\Inputs\\AIL.pkl')

# Area expansion factor
aefC=geosC['rgsf']**2

mosC=cbu.Import_MOS_ByScnAndStrata_GHGEcon_FromPoints(metaC)

# Define scenario comparisons
mosC['Delta']={}
mosC['Delta']['NM']={'iB':0,'iP':1}
mosC['Delta']['NMHH']={'iB':2,'iP':3}
mosC=cbu.Import_MOS_ByScnComparisonAndStrata_FromPoints(metaC,mosC)

# Import areas
mosC=cbu.Import_MOS_ByScnAndStrata_Area_FromPoints(metaC,mosC)

list(mosC['Scenarios'][0]['Sum'].keys())

# Indices
nPS='All'
nSS='All'
iPS=np.where(metaC['Project']['Strata']['Project']['Unique CD']==nPS)[0][0]
iSS=np.where(metaC['Project']['Strata']['Spatial']['Unique CD']==nSS)[0][0]
iT=np.where( (tvC>=1850) & (tvC<=2150) )[0]
cmp='NM'

#%% Import future simulations

metaF=gu.ipickle(r'D:\Data\FCI_Projects\BCFCS_NM_20k_Future\Inputs\Metadata.pkl')
#metaF['Paths']['Figures']=r'C:\Users\rhember\OneDrive - Government of BC\Figures\BCFCS_NM_1k_Sparse'
tvF=np.arange(metaF['Project']['Year Start Saving'],metaF['Project']['Year End']+1,1)
geosF=gu.ipickle(metaF['Paths']['Project'] + '\\Geospatial\\geos.pkl')

# Area expansion factor
aefF=55000000/metaF['Project']['N Stand']

mosF=cbu.Import_MOS_ByScnAndStrata_GHGEcon_FromPoints(metaF)

# Define scenario comparisons
mosF['Delta']={}
mosF['Delta']['NM']={'iB':0,'iP':1}
mosF['Delta']['NMHH']={'iB':2,'iP':3}
mosF=cbu.Import_MOS_ByScnComparisonAndStrata_FromPoints(metaF,mosF)

# Import areas
mosF=cbu.Import_MOS_ByScnAndStrata_Area_FromPoints(metaF,mosF)

#%% Combine completed and future

mosCF=copy.deepcopy(mosC)
iT=np.where( tvC>=tvF[0])[0]
for iScn in range(metaC['Project']['N Scenario']):
    for k1 in mosC['Scenarios'][iScn]['Sum'].keys():
        for k2 in mosC['Scenarios'][iScn]['Sum'][k1].keys():
            mosCF['Scenarios'][iScn]['Sum'][k1][k2]=mosCF['Scenarios'][iScn]['Sum'][k1][k2]*aefC
            mosCF['Scenarios'][iScn]['Mean'][k1][k2]=mosCF['Scenarios'][iScn]['Mean'][k1][k2]
            mosCF['Scenarios'][iScn]['Sum'][k1][k2][iT,iPS,iSS]=mosCF['Scenarios'][iScn]['Sum'][k1][k2][iT,iPS,iSS]+mosF['Scenarios'][iScn]['Sum'][k1][k2][:,iPS,iSS]*aefF
            mosCF['Scenarios'][iScn]['Mean'][k1][k2][iT,iPS,iSS]=mosCF['Scenarios'][iScn]['Mean'][k1][k2][iT,iPS,iSS]+mosF['Scenarios'][iScn]['Mean'][k1][k2][:,iPS,iSS]
for cmp in mosC['Delta'].keys():
    for k1 in mosC['Delta'][cmp]['ByStrata']['Sum'].keys():
        for k2 in mosC['Delta'][cmp]['ByStrata']['Sum'][k1].keys():
            mosCF['Delta'][cmp]['ByStrata']['Sum'][k1][k2]=mosCF['Delta'][cmp]['ByStrata']['Sum'][k1][k2]*aefC
            mosCF['Delta'][cmp]['ByStrata']['Mean'][k1][k2]=mosCF['Delta'][cmp]['ByStrata']['Mean'][k1][k2]
            mosCF['Delta'][cmp]['ByStrata']['Sum'][k1][k2][iT,iPS,iSS]=mosCF['Delta'][cmp]['ByStrata']['Sum'][k1][k2][iT,iPS,iSS]+mosF['Delta'][cmp]['ByStrata']['Sum'][k1][k2][:,iPS,iSS]*aefF
            mosCF['Delta'][cmp]['ByStrata']['Mean'][k1][k2][iT,iPS,iSS]=mosCF['Delta'][cmp]['ByStrata']['Mean'][k1][k2][iT,iPS,iSS]+mosF['Delta'][cmp]['ByStrata']['Mean'][k1][k2][:,iPS,iSS]

#%% Plot

plt.close('all')
iScn=1
plt.plot(tvC,mosC['Scenarios'][iScn]['Sum']['Area_Fertilization Aerial']['Ensemble Mean'][:,iPS,iSS]*aefC/1000,'-bo')
plt.plot(tvF,mosF['Scenarios'][iScn]['Sum']['Area_Fertilization Aerial']['Ensemble Mean'][:,iPS,iSS]*aefF/1000,'-gs')
plt.plot(ailC['Year'],ailC['Area Total']/1000,'rs')
print(np.mean(mosF['Scenarios'][iScn]['Sum']['Area_Fertilization Aerial']['Ensemble Mean'][:,iPS,iSS]*aefF/1000))

#%% Plot

plt.close('all')
iScn=1
plt.plot(tvF,mosF['Scenarios'][1]['Sum']['Area_Harvest']['Ensemble Mean'][:,iPS,iSS]*aefF/1e6,'-bo')
plt.plot(tvF,mosF['Scenarios'][3]['Sum']['Area_Harvest']['Ensemble Mean'][:,iPS,iSS]*aefF/1e6,'-gs')
print(np.mean(mosF['Scenarios'][1]['Sum']['Area_Harvest']['Ensemble Mean'][:,iPS,iSS]*aefF/1e6))
print(np.mean(mosF['Scenarios'][3]['Sum']['Area_Harvest']['Ensemble Mean'][:,iPS,iSS]*aefF/1e6))
#print(np.mean(mosF['Scenarios'][iScn]['Sum']['Area_Fertilization Aerial']['Ensemble Mean'][:,iPS,iSS]*aefF/1000))

#%%

plt.close('all')
plt.plot(tvC,mosC['Delta'][cmp]['ByStrata']['Sum']['E_CO2e_AGHGB_WSub']['Ensemble Mean'][:,iPS,iSS]*aefC/1e6,'-bo')
plt.plot(tvF,mosF['Delta'][cmp]['ByStrata']['Sum']['E_CO2e_AGHGB_WSub']['Ensemble Mean'][:,iPS,iSS]*aefF/1e6,'-gs')

#%%

plt.close('all')
plt.plot(tvC,mosCF['Delta'][cmp]['ByStrata']['Sum']['E_CO2e_AGHGB_WSub']['Ensemble Mean'][:,iPS,iSS]/1e6,'-bo')

#%%

plt.close('all')
iScn=1
plt.plot(tvC,gu.movingave(mosCF['Delta']['NM']['ByStrata']['Sum']['V_ToMillMerchTotal']['Ensemble Mean'][:,iPS,iSS]/1000,30,'historical'),'-bo')
plt.plot(tvC,gu.movingave(mosCF['Delta']['NMHH']['ByStrata']['Sum']['V_ToMillMerchTotal']['Ensemble Mean'][:,iPS,iSS]/1000,30,'historical'),'-gs')

#%%

plt.close('all')
iScn=1
plt.plot(tvC,gu.movingave(mosCF['Delta']['NM']['ByStrata']['Sum']['Cost Total_cumu']['Ensemble Mean'][:,iPS,iSS]/1e6,30,'historical'),'-bo')
plt.plot(tvC,gu.movingave(mosCF['Delta']['NMHH']['ByStrata']['Sum']['Cost Total_cumu']['Ensemble Mean'][:,iPS,iSS]/1e6,30,'historical'),'-gs')

#%% Plot mitigation value

plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(11,10)); fs2=7

# Raidial lines and CPT text
yrs=[1990,2020,2050,2070,2100]
ha=['right','right','left','left','left','left']
xa=[1.8,1.35,1.2,1.06,1.035,1.07]

ax.plot([0,0],[-2000,2000],'k-',lw=0.5,color='k')
ax.plot([-10000,10000],[0,0],'k-',lw=0.5,color='k')
rings=np.arange(100,1200,100)
for ring in rings:
    ax.plot(0,0,'o',ms=ring,mec=[0.8,0.8,0.8],mfc='none',mew=0.5)
cpt=np.array([160,80,40,20,10,5,1],dtype=float)
x_lab=np.array([-3000,-2800,-2400,-1625,-910,-470,-300])
Slope=np.zeros(cpt.size)
for iB in range(cpt.size):
    ax.plot([0,-10000],[0,-10000/cpt[iB]],'k-',color=[0.8,0.8,0.8])
#ax.text(xC,-y_new,int(cpt[iB]),fontsize=fs2)

cmp='NM'
dghgCF=mosCF['Delta'][cmp]['ByStrata']['Sum']['E_CO2e_AGHGB_WSub_cumu']['Ensemble Mean'][:,iPS,iSS]/1e6
dcostCF=mosCF['Delta'][cmp]['ByStrata']['Sum']['Cost Nutrient Management_cumu']['Ensemble Mean'][:,iPS,iSS]/1e6
drnCF=mosCF['Delta'][cmp]['ByStrata']['Sum']['Revenue Net_cumu']['Ensemble Mean'][:,iPS,iSS]/1e6

ax.plot(-1*dcostCF,dghgCF,'r-',color=[0.27,0.49,0.77],lw=1.5,label='Cost')
ax.plot(drnCF,dghgCF,'b-',color=[0.5,0.85,0],lw=1.5,label='Cost minus gross revenue')

# cmp='NMHH'
# dghgCF=mosCF['Delta'][cmp]['ByStrata']['Sum']['E_CO2e_AGHGB_WSub_cumu']['Ensemble Mean'][:,iPS,iSS]/1e6
# dcostCF=mosCF['Delta'][cmp]['ByStrata']['Sum']['Cost Total_cumu']['Ensemble Mean'][:,iPS,iSS]/1e6
# drnCF=mosCF['Delta'][cmp]['ByStrata']['Sum']['Revenue Net_cumu']['Ensemble Mean'][:,iPS,iSS]/1e6
# ax.plot(-1*dcostCF,dghgCF,'r--',color=[0.27,0.49,0.77],lw=1,label='Cost')
# ax.plot(drnCF,dghgCF,'b--',color=[0.5,0.85,0],lw=1,label='Cost minus gross revenue')

ax.set(xticks=np.arange(-12000,18000,500),xlabel='Cumulative $\Delta$ cost (CAD x 1 Million)',
       ylabel='Cumulative $\Delta$ GHG balance (MtCO$_2$e)',xlim=[-2500,200],ylim=[-110,10])
ax.legend(loc='center left',frameon=False,facecolor='w') # ,bbox_to_anchor=(0.66,1)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')
ax.get_xaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
plt.tight_layout()

flg=0
if flg==1:
    def onclick(event):
        global ix, iy
        ix,iy = event.xdata, event.ydata
        global xt,yt
        xt.append(ix)
        yt.append(iy)
        return xt,yt

    xt=[]; yt=[]
    for i in range(7):
        cid=fig.canvas.mpl_connect('button_press_event', onclick)

    for i in range(len(xt)):
        ax.text(xt[i],yt[i],int(cpt[i]),fontsize=fs2,color=[0,0,0])

gu.PrintFig(metaC['Paths']['Figures'] + '\\All_All_MitigationValue','png',900)

#%% Save tabular summary

flg=1
if flg==1:
    t_start=2021
    t_end=2100
    sum_mult=1.0/1e6
    df=udem.ExportSummariesByScenario(metaC,mosCF,t_start,t_end,sum_mult=sum_mult)

#%% Timber harvesting landbase

flg=1
if flg==1:
    iScn=1
    utl.Plot_THLB(meta,tv,aefC,thlb,iScn)

#%% Carbon pool change

flg=1
if flg==1:
    iScn=1
    utl.PlotCarbonPoolBarChart(meta,tv,mos,aefC,iScn,iPS,iSS)

#%% GHG flux change

flg=1
if flg==1:
    iScn=1
    utl.PlotFluxesBarChart(meta,tv,mos,aefC,iScn,iPS,iSS)

#%% Harvest removals

flg=1
if flg==1:
    iScn=1
    flg_ann=1
    flg_dead=0
    flg_hbs=1
    utl.PlotHarvestVolume(meta,mos,tv,aefC,iScn,iT,iPS,iSS,flg_ann,flg_dead,flg_hbs)

#%% GHG balance

flg=1
if flg==1:
    iScn=1
    utl.PlotGHGBalance(meta,mos,tv,aefC,iScn,iT,iPS,iSS)

#%% GHG balance (simple)

flg=1
if flg==1:
    iScn=1
    utl.PlotGHGBalanceSimple(meta,mos,tv,aefC,iScn,iT,iPS,iSS)

#%% Annual GHG balance and benefit

flg=1
if flg==1:
    utl.PlotGHGBenefit(meta,mos,tv,aefC,cmp,iT,iPS,iSS)

#%% GHG balance and benefit (annual and Cumulative)

flg=1
if flg==1:
    utl.PlotGHGBalanceAndBenefitWithCumulative(meta,mos,tv,aefC,cmp,iT,iPS,iSS)

#%% Forcing from scenario comparison

flg=1
if flg==1:
    utl.PlotForcingBarChart(meta,tv,mos,aefC,cmp,iPS,iSS)

#%% Comparison with PIR

flg=1
if flg==1:
    iScn=1
    utl.PlotComparisonWithPIR(meta,mos,tv,aefC,iScn,iT,iPS,iSS)

#%% Area disturbance and management

flg=1
if flg==1:
    iScn=1
    ivlT=1
    utl.PlotAreaDisturbed(meta,mos,tv,aefC,ivlT,iScn,iT,iPS,iSS)

#%% Mortality summary

flg=1
if flg==1:
    iScn=1
    ivlT=1
    utl.MortalitySummary(meta,mos,C_M_ByAgent,tv,iScn,iT,ivlT,iPS,iSS)

#%% Net biomass growth

flg=1
if flg==1:
    iScn=0
    utl.PlotNetGrowthTS(meta,mos,tv,C_M_ByAgent,iT,iScn,iPS,iSS)

#%% Harvest volume/ha

flg=1
if flg==1:
    iScn=1
    iT2=np.where( (tv>=1960) & (tv<=2025) )[0]
    utl.PlotVolumePerHectare(meta,mos,tv,iScn,iT2,iPS,iSS,aefC)

#%% Import map of mean for specified time period (for mapping)

flg=1
if flg==1:

    iScn=1

    vL=['A','C_Biomass_Tot','C_Litter_Tot','C_Soil_Tot','C_DeadWood_Tot','C_G_Gross_Tot',
        'C_G_Net_Tot','C_M_Reg_Tot','C_M_Dist','C_LF_Tot','C_ToMillMerch',
        'C_ToMillNonMerch', 'C_ToMillSnagStem', 'C_ToSlashpileBurnTot']

    flg_FromScratch=1

    if flg_FromScratch==1:
        mu_mod=cbu.Calc_MOS_MapMean(meta,iScn,[2010,2020],VariablesToKeep=vL)
        gu.opickle(meta['Paths']['Project'] + '\\Outputs\\MapData_iScn' + str(iScn) + '.pkl',mu_mod)
    else:
        mu_mod=gu.ipickle(meta['Paths']['Project'] + '\\Outputs\\MapData_iScn' + str(iScn) + '.pkl')

#%% QA Compare biomass by BGC Zone

flg=1
if flg==1:
    utl.QA_Biomass_ByBGCZone(meta,mu_mod)

#%% QA Compare soil organic carbon by BGC Zone

flg=1
if flg==1:
    soc=utl.QA_SOC_ByBGCZone(meta,mu_mod)

#%% QA Compare average biomass dynamics

flg=1
if flg==1:
    iScn=1
    utl.Plot_AverageBiomassDynamics(meta,mos,iScn,iPS,iSS,iT)

#%% Comparison of anthropgenic forcing

flg=1
if flg==1:
    utl.CompareAnthropogenicComponentWithCFS(meta,mos,iT,iPS,iSS,aefC)

#%% Summary of existing approaches

flg=1
if flg==1:
    utl.SummarizeExistingInitiatives(meta)

#%% Age class distribution

flg=1
if flg==1:
    iScn=1
    utl.Plot_AgeClassDist(meta,iScn,iPS,iSS)

#%% Plot all time series

flg=0
if flg==1:
    utl.PlotAllVariableTimeSeries(meta,mos,tv,iScn,iT,iPS,iSS)

#%% Radiative forcing

gp=gu.SetGraphics('Manuscript')

# Import
v0=cbu.LoadSingleOutputFile(meta,1,0,0)
tv=np.arange(meta['Project']['Year Start Saving'],meta['Project']['Year End']+1,1)

# Extract biophysical parameters
bB=meta['Param']['BE']['Biophysical']

# Import atmospheric GHG concentration data
dA=gu.ReadExcel(meta['Paths']['Model Code'] + '\\Parameters\\Parameters_RCP_GHG_Abundance_AR5.xlsx')
for k in dA.keys():
    if k=='Year':
        continue
    dA[k]=np.interp(tv,dA['Year'],dA[k])

#%% Radiative forcing

c_in=np.sum(v0['Atm_CO2_In'],axis=1)
c_out=np.sum(v0['Atm_CO2_Out'],axis=1)

t=np.arange(0,c_in.size)
r=0.001

c_out_decay=np.zeros( (c_in.size,c_in.size) )
for i in range(c_in.size):
    t=np.arange(0,c_in.size-i)
    c_out_decay[i:,i]=-np.append(0,np.diff(c_in[i]*np.exp(-r*t)))
c_out_decay=np.sum(c_out_decay,axis=1)

plt.close('all')
plt.plot(tv,c_in-c_out,'b-')
plt.plot(tv,c_in-c_out-c_out_decay,'r--')

#%% Radiative forcing

for i in range(c_in.size):
    #plt.plot(np.diff(c_out_decay))
    c_atm[i:,i]=c_atm[i:,i]+c_out_decay

# Remove outputs
#c_atm=np.mean(c_atm,axis=1)#-c_out

plt.close('all')
plt.plot(tv,c_atm,'b-')
plt.plot(tv,np.cumsum(c_in-c_out),'r--')
#plt.plot(tv,c_atm-c_out,'r--')

#plt.close('all')
#plt.plot(tv,c_in-c_out,'r-')
#plt.plot(tv,np.cumsum(c_out),'g-')

#%% Radiative forcing
#
plt.plot(c1[:,0]/c1[0,0],'b-')
#plt.plot(a[:,0],'b-')


#bB['Ratio_C_to_CO2']*
#/bB['Atmospheric Density CO2']


#%% Calibrate SI

# Not working well

# from scipy.optimize import curve_fit
# from fcexplore.psp.Processing import psp_utilities as utl_gp

# gp=gu.SetGraphics('Manuscript')

# # Import ground plot data
# metaGP={}
# metaGP['Paths']={}
# metaGP['Paths']['DB']=r'C:\Users\rhember\Documents\Data\GroundPlots\PSP-NADB2'
# metaGP=utl_gp.ImportParameters(metaGP)
# d=gu.ipickle(metaGP['Paths']['DB'] + '\\Processed\\L2\\L2_BC.pkl')
# sl=d['sobs']
# del d

# # Filter
# sl['pt_ind']=np.zeros(sl['ID Plot'].size)
# ind=np.where( (sl['Plot Type']==metaGP['LUT']['Plot Type BC']['CMI']) | (sl['Plot Type']==metaGP['LUT']['Plot Type BC']['NFI']) & (sl['Lat']>0) & (sl['Lon']!=0) )[0]
# sl['pt_ind'][ind]=1

# bBGC=gu.ReadExcel(r'C:\Users\rhember\Documents\Code_Python\fcgadgets\cbrunner\Parameters\Parameters_ClimateByBGC.xlsx')

# # Calculate stats by BGC zone
# vL=['Age t0','Cbk L t0','Cbr L t0','Cf L t0','Csw L t0','Cr L t0','Cag L t0','Ctot L t0']
# u=np.unique(sl['Ecozone BC L1'][sl['Ecozone BC L1']>0])
# lab=np.array(['' for _ in range(u.size)],dtype=object)
# d={}
# for v in vL:
#     d[v]={}
#     d[v]['N']=np.zeros(u.size)
#     d[v]['mu']=np.zeros(u.size)
#     d[v]['sd']=np.zeros(u.size)
#     d[v]['se']=np.zeros(u.size)

# for v in bBGC.keys():
#     d[v]=np.zeros(u.size)

# for i in range(u.size):
#     lab[i]=utl_gp.lut_id2cd(metaGP,'Ecozone BC L1',u[i])

#     # Observed values
#     for v in vL:
#         ind=np.where( (sl['Ecozone BC L1']==u[i]) &
#                      (sl['pt_ind']==1) &
#                      (sl['Cbk L t0']>=0) & (sl['Cbk L t0']<2000) &
#                      (sl['Cbr L t0']>=0) & (sl['Cbr L t0']<2000) &
#                      (sl['Cf L t0']>=0) & (sl['Cf L t0']<2000) &
#                      (sl['Cr L t0']>=0) & (sl['Cr L t0']<2000) &
#                      (sl['Csw L t0']>=0) & (sl['Csw L t0']<2000) &
#                      (sl['Ctot L t0']>=0) & (sl['Ctot L t0']<10000))[0]
#         d[v]['N'][i]=ind.size
#         d[v]['mu'][i]=np.nanmean(sl[v][ind])
#         d[v]['sd'][i]=np.nanstd(sl[v][ind])
#         #d[v]['se'][i]=np.nanstd(sl[v][ind])/np.sqrt(ind[0].size)

#     ind=np.where(bBGC['Name']==lab[i])[0]
#     for v in bBGC.keys():
#         if v=='Name':
#             continue
#         d[v][i]=bBGC[v][ind[0]]

# d['Ctot L t0']['mu']=d['Cbk L t0']['mu']+d['Cbr L t0']['mu']+d['Cf L t0']['mu']+d['Cr L t0']['mu']+d['Csw L t0']['mu']

# #%%

# #funT_txt='max(0,b(2).^((x(:,1)-10)./10)).*(1./(1+exp(-b(3).*(x(:,1)-16))))';
# #funN_txt='max(0,1-exp(-b(4).*(x(:,2)-b(5))))';
# #funC_txt='(1+b(6).*log(x(:,3)./280))';

# def funT(x,a,b):
#     yhat=a*b**((x-10)/10)
#     return yhat

# def funW(x,a,b,c):
#     yhat=a*(1/(1+np.exp(-b*(x-c))))
#     return yhat

# def func(x,b,c,d):
#     yhat=150.0*b**((x[:,0]-10)/10)*(1/(1+np.exp(-c*(x[:,1]-d))))
#     return yhat

# That=np.arange(-1,12,0.1)
# What=np.arange(0,201,1)

# plt.close('all')
# plt.plot(That,funT(That,150,3.3),'b-')

# plt.close('all')
# plt.plot(What,funW(What,150,0.023,-10.7),'b-')

# x=np.column_stack((d['MAT'],d['WS']))
# y=d['Ctot L t0']['mu']

# ikp=np.where(np.isin(lab,['CDF','BAFA','PP','BG','CMA'])==False)[0]

# popt,pcov=curve_fit(func,x[ikp,:],y[ikp],[3.3,0.04,110])
# yhat=func(x,popt[0],popt[1],popt[2])

# plt.close('all')
# plt.plot(y[ikp],yhat[ikp],'ko')

