
"""
PSP - DESCRIPTIVE STATISTICS
"""

#%% Import modules

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf
import fcgadgets.macgyver.utilities_general as gu
import fcgadgets.macgyver.utilities_gis as gis
import fcgadgets.bc1ha.bc1ha_utilities as u1ha
import fcexplore.psp.Processing.psp_utilities as ugp

gp=gu.SetGraphics('Manuscript')

#%% Import data

meta=u1ha.Init()
meta=u1ha.ImportLUTs(meta)
meta,gpt=ugp.ImportPlotData(meta,type='Stand')

#%% Import Raster grids

vList=['refg','lcc1_c','bgcz','harv_yr_con1','fire_yr','ibm_yr'] #'lcc1_c','gfcly','gfcly_filt',
z=u1ha.Import_Raster(meta,[],vList)
iGrd=gis.GetGridIndexToPoints(z['refg'],gpt['X'],gpt['Y'])

# zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])
# z['lcc1_c']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCoverClass1_Current.tif')
# z['bgcz']=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\BEC_BIOGEOCLIMATIC_POLY\\ZONE.tif')
# zAge=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\PROJ_AGE_1.tif')

# zH=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VEG_CONSOLIDATED_CUT_BLOCKS_SP\\HARVEST_YEAR_YearLast.tif')
# zF=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PROT_HISTORICAL_FIRE_POLYS_SP\\FIRE_YEAR_MaskAll.tif')
# zI=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\PEST_SEVERITY_CODE_IBM_MaskAll.tif')

#%% Export summary by Plot Type

d={}
ind=np.where( (gpt['Plot Type']==meta['LUT']['GP']['Plot Type BC']['VRI']) )[0]
for k in gpt.keys():
    d[k]=np.round(np.nanmean(gpt[k][ind]),decimals=2)
ind=np.where( (gpt['Plot Type']==meta['LUT']['GP']['Plot Type BC']['YSM']) )[0]
for k in gpt.keys():
    d[k]=np.append(d[k],np.round(np.nanmean(gpt[k][ind]),decimals=2))
ind=np.where( (gpt['PTF CN']==1) )[0]
for k in gpt.keys():
    d[k]=np.append(d[k],np.round(np.nanmean(gpt[k][ind]),decimals=2))
df=pd.DataFrame(d,index=[0,1,2])
df.to_excel(meta['Paths']['GP']['DB'] + '\\Processed\\SummarySL_ByPlotType.xlsx')

#%% QA - compare standing biomass with standing volume whole stem

plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(7.8,7.2))
ind=np.where( (gpt['PTF CN']==1) & (gpt['Vws L t0']>0) & (gpt['Ctot L t0']>0) )[0]
ax.plot(gpt['Vws L t0'][ind],gpt['Ctot L t0'][ind],'b.',ms=4,mec='w',mfc=[0.27,0.44,0.79],markeredgewidth=0.25)
rs,txt=gu.GetRegStats(gpt['Vws L t0'][ind],gpt['Ctot L t0'][ind])
ax.plot(rs['xhat'],rs['yhat'],'-k')
ax.set(ylabel='Tree biomass (MgC ha$^{-1}$)',xlabel='Whole-stem volume (m$^{3}$ ha$^{-1}$)',xlim=[0,1800],ylim=[0,600])
ax.text(900,100,txt,fontsize=8)
plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])
#gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_QA_BiomassVsVolume','png',900)

#%% Average biomass dynammics

ikp=np.where( (gpt['PTF CN']==1) & (gpt['Year t1']>0) )[0]
vL=['Year t0','Year t1','Ctot G Surv','Ctot G Recr','Ctot Mort Nat','Ctot Net','Ctot Mort Harv']
sts={}
for v in vL:
    sts['mu ' + v]=np.nanmean(gpt[v][ikp])
    sts['se ' + v]=np.nanstd(gpt[v][ikp])/np.sqrt(ikp.size)

print( str(np.nanpercentile(gpt['Year t0'],25)) + ' ' + str(np.nanpercentile(gpt['Year t1'],75)) )
print(sts['mu Ctot Net'])

cl=np.array([[0.75,0.75,0.75],[0.24,0.49,0.77],[0.6,1,0]])
cle=[0,0,0]#[0.05,0.2,0.45]
barw=0.6
lab=['Survivor\ngrowth','Recruitment\ngrowth','Natural\nmortality','Harvest\nmortality','Net\ngrowth'] #

plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(7.8,6))
ax.plot([0,6],[0,0],'-k',color=gp['cla'],lw=0.5)
ax.bar(1,sts['mu Ctot G Surv'],barw,facecolor=cl[0,:],label='Growth survivors')
ax.bar(2,sts['mu Ctot G Recr'],barw,facecolor=cl[0,:],label='Growth recruitment')
ax.bar(3,-sts['mu Ctot Mort Nat'],barw,facecolor=cl[0,:],label='Natural\nmortality')
ax.bar(4,-sts['mu Ctot Harv'],barw,facecolor=cl[0,:],label='Harvest')
ax.bar(5,sts['mu Ctot Net'],barw,facecolor=cl[0,:],label='Net')
ax.errorbar(1,sts['mu Ctot G Surv'],yerr=sts['se Ctot G Surv'],color=cle,fmt='none',capsize=2,lw=0.5)
ax.errorbar(2,sts['mu Ctot G Recr'],yerr=sts['se Ctot G Recr'],color=cle,fmt='none',capsize=2,lw=0.5)
ax.errorbar(3,-sts['mu Ctot Mort Nat'],yerr=sts['se Ctot Mort Nat'],color=cle,fmt='none',capsize=2,lw=0.5)
ax.errorbar(4,-sts['mu Ctot Harv'],yerr=sts['se Ctot Harv'],color=cle,fmt='none',capsize=2,lw=0.5)
ax.errorbar(5,sts['mu Ctot Net'],yerr=sts['se Ctot Net'],color=cle,fmt='none',capsize=2,lw=0.5)
ax.set(xticks=np.arange(1,len(lab)+1),xticklabels=lab,ylabel='Carbon balance of trees (MgC ha${^-1}$ yr$^{-1}$)',xlim=[0.5,5.5],ylim=[-1.5,2])
#plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])
plt.tight_layout()
gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_BiomassDynamics_Mean_CN','png',900)

#%% Average biomass dynammics (totals in CO2e)

# Area of forest
ind=np.where(z['lcc1_c']['Data']==meta['LUT']['Derived']['lcc1']['Forest'])
A=ind[0].size

ikp=np.where( (gpt['PTF CN']==1) & (gpt['Year t1']>0) )[0]
vL=['Year t0','Year t1','Ctot G Surv','Ctot G Recr','Ctot Mort Nat','Ctot Net','Ctot Mort Harv']
sts={}
for v in vL:
    sts['mu ' + v]=A*np.nanmean(gpt[v][ikp])/1e6*3.667
    sts['se ' + v]=A*np.nanstd(gpt[v][ikp])/np.sqrt(ikp.size)/1e6*3.667

print( str(np.nanpercentile(gpt['Year t0'],25)) + ' ' + str(np.nanpercentile(gpt['Year t1'],75)) )
print(sts['mu Ctot Net'])
print(sts['mu Ctot Harv'])

cl=np.array([[0.75,0.75,0.75],[0.24,0.49,0.77],[0.6,1,0]])
cle=[0,0,0]#[0.05,0.2,0.45]
barw=0.6
lab=['Survivor\ngrowth','Recruitment\ngrowth','Natural\nmortality','Harvest\nmortality','Net\ngrowth'] #

plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(7.8,6))
ax.plot([0,6],[0,0],'-k',color=gp['cla'],lw=0.5)
ax.bar(1,sts['mu Ctot G Surv'],barw,facecolor=cl[0,:],label='Growth survivors')
ax.bar(2,sts['mu Ctot G Recr'],barw,facecolor=cl[0,:],label='Growth recruitment')
ax.bar(3,-sts['mu Ctot Mort Nat'],barw,facecolor=cl[0,:],label='Mortality')
ax.bar(4,-sts['mu Ctot Harv'],barw,facecolor=cl[0,:],label='Mortality')
ax.bar(5,sts['mu Ctot Net'],barw,facecolor=cl[0,:],label='Mortality')
ax.errorbar(1,sts['mu Ctot G Surv'],yerr=sts['se Ctot G Surv'],color=cle,fmt='none',capsize=2,lw=0.5)
ax.errorbar(2,sts['mu Ctot G Recr'],yerr=sts['se Ctot G Recr'],color=cle,fmt='none',capsize=2,lw=0.5)
ax.errorbar(3,-sts['mu Ctot Mort Nat'],yerr=sts['se Ctot Mort Nat'],color=cle,fmt='none',capsize=2,lw=0.5)
ax.errorbar(4,-sts['mu Ctot Harv'],yerr=sts['se Ctot Harv'],color=cle,fmt='none',capsize=2,lw=0.5)
ax.errorbar(5,sts['mu Ctot Net'],yerr=sts['se Ctot Net'],color=cle,fmt='none',capsize=2,lw=0.5)
ax.set(xticks=np.arange(1,len(lab)+1),xticklabels=lab,ylabel='Carbon balance of trees (MtCO$_{2}$e yr$^{-1}$)',xlim=[0.5,5.5],ylim=[-350,450])
#plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])
plt.tight_layout()
gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_BiomassDynamics_Sum_CN','png',900)

#%% Average volume dynammics

ikp=np.where( (gpt['PTF CN']==1) & (gpt['Year t1']>0) )[0]

vL=['Year t0','Year t1','Vws G Surv','Vws G Recr','Vws Mort Nat','Vws Harv','Vws Mort Net'] #,
sts={}
for v in vL:
    sts['mu ' + v]=np.nanmean(gpt[v][ikp])
    sts['se ' + v]=np.nanstd(gpt[v][ikp])/np.sqrt(ikp.size)

print( str(np.nanpercentile(gpt['Year t0'],25)) + ' ' + str(np.nanpercentile(gpt['Year t1'],75)) )
print(sts['mu Vws Net'])
#print(sts['mu Vws Mort Harv'])

cl=np.array([[0.8,0.8,0.8],[0.6,1,0]])
cle=[0,0,0]
barw=0.5
lab=['Survivor\ngrowth','Recruitment\ngrowth','Natural\nmortality','Harvest\nmortality','Net\ngrowth'] #

plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(7.8,6))
ax.plot([0,6],[0,0],'-k',color=gp['cla'],lw=0.5)
ax.bar(1,sts['mu Vws G Surv'],barw,facecolor=cl[0,:],label='Growth survivors')
ax.bar(2,sts['mu Vws G Recr'],barw,facecolor=cl[0,:],label='Growth recruitment')
ax.bar(3,-sts['mu Vws Mort Nat'],barw,facecolor=cl[0,:],label='Mortality')
ax.bar(4,-sts['mu Vws Harv'],barw,facecolor=cl[0,:],label='Mortality')
ax.bar(5,sts['mu Vws Net'],barw,facecolor=cl[0,:],label='Mortality')
ax.errorbar(1,sts['mu Vws G Surv'],yerr=sts['se Vws G Surv'],color=cle,fmt='none',capsize=2,lw=0.5)
ax.errorbar(2,sts['mu Vws G Recr'],yerr=sts['se Vws G Recr'],color=cle,fmt='none',capsize=2,lw=0.5)
ax.errorbar(3,-sts['mu Vws Mort Nat'],yerr=sts['se Vws Mort Nat'],color=cle,fmt='none',capsize=2,lw=0.5)
ax.errorbar(4,-sts['mu Vws Harv'],yerr=sts['se Vws Harv'],color=cle,fmt='none',capsize=2,lw=0.5)
ax.errorbar(5,sts['mu Vws Net'],yerr=sts['se Vws Net'],color=cle,fmt='none',capsize=2,lw=0.5)
ax.set(xticks=np.arange(1,len(lab)+1),xticklabels=lab,ylabel='Whole stem volume periodic increment (m$^{3}$ ha$^{-1}$ yr$^{-1}$)',xlim=[0.5,5.5],ylim=[-4,5])
#plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])
plt.tight_layout()
gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_VolumeDynamics_Mean_CN','png',900)

#%% Fire profiles

tvFire=np.arange(1990,2022,1)
fy=z['fire_yr']['Data'][iGrd]
tsf=gpt['Year t0']-fy

bw=5; bin=np.arange(-12.5,20,bw)
prfl={}
for k in gpt.keys():
    prfl[k]={}
    prfl[k]['N'],prfl[k]['mu'],prfl[k]['med'],prfl[k]['sig'],prfl[k]['se']=gu.discres(tsf,gpt[k],bw,bin)

# Plot
v='Age Med t0'
#v='Ctot L t0'
#v='Ctot G Recr'
#v='Ctot Mort'
#v='N Net (%)'
it0=np.where(bin<0)[0]
it1=np.where( (bin>0) & (bin<=15) )[0]
plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(7.8,6));
ax.plot([0,0],[0,1.07*np.nanmax(prfl[v]['mu'])],'k-',lw=2,color=[0.75,0.75,0.75])
ax.plot(bin,prfl[v]['mu'],'ko',ms=3,mfc='k',mec=cl[0,:],label='Live')
ax.errorbar(bin,prfl[v]['mu'],yerr=prfl[v]['se'],color='k',fmt='none',capsize=1.5,lw=0.5)
mu0=np.nanmean(prfl[v]['mu'][it0])
mu1=np.nanmean(prfl[v]['mu'][it1])
ax.plot(bin[it0],mu0*np.ones(it0.size),'k--',color='k',lw=0.5)
ax.plot(bin[it1],mu1*np.ones(it1.size),'k-.',color='k',lw=0.5)
ax.text(np.mean(bin[it1]),1.15*mu1,'Change = ' + str(np.round(mu1-mu0,decimals=1)) + '\n(' + str(np.round( (mu1-mu0)/mu0*100 ,decimals=1)) + '%)',ha='center',fontsize=7,style='normal',color='k')
ax.set(ylabel='Tree carbon (tC ha$^1$)',xlabel='Time since detection, years',ylim=[0,1.07*np.nanmax(prfl[v]['mu'])])
#ax[0].legend(loc='upper right',frameon=False,facecolor='w',edgecolor='w')
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])
plt.tight_layout()
#gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\QA_Wildfire_Profiles','png',900)

#%% IBM profiles

tso=-100*np.ones(gpt['Year t0'].size)
for iY in range(10):
    zY=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\PEST_SEVERITY_CODE_IBM_' + str(iY+1) + '_Year.tif')['Data'][iGrd]
    zS=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\PEST_SEVERITY_CODE_IBM_' + str(iY+1) + '_Severity.tif')['Data'][iGrd]
    ind=np.where( (zS==meta['LUT']['PEST_INFESTATION_POLY']['PEST_SEVERITY_CODE']['S']) & (gpt['Year t0']-zY>tso) )[0]
    tso[ind]=gpt['Year t0'][ind]-zY[ind]

#%%

bw=5; bin=np.arange(-15,25,bw)
prfl={}
for k in gpt.keys():
    prfl[k]={}
    prfl[k]['N'],prfl[k]['mu'],prfl[k]['med'],prfl[k]['sig'],prfl[k]['se']=gu.discres(tso,gpt[k],bw,bin)

# Plot
#v='Ctot L t0'
#v='Ctot D t0'
#v='Ctot G Surv'
#v='Ctot G Recr'
v='Ctot Mort Nat'
#v='N Recr (%)'
it0=np.where(bin<0)[0]
it1=np.where( (bin>0) & (bin<=15) )[0]
plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(7.8,6));
ax.plot([0,0],[0,1.07*np.nanmax(prfl[v]['mu'])],'k-',lw=2,color=[0.75,0.75,0.75])
ax.plot(bin,prfl[v]['mu'],'ko',ms=3,mfc='k',mec=cl[0,:],label='Live')
ax.errorbar(bin,prfl[v]['mu'],yerr=prfl[v]['se'],color='k',fmt='none',capsize=1.5,lw=0.5)
mu0=np.nanmean(prfl[v]['mu'][it0])
mu1=np.nanmean(prfl[v]['mu'][it1])
ax.plot(bin[it0],mu0*np.ones(it0.size),'k--',color='k',lw=0.5)
ax.plot(bin[it1],mu1*np.ones(it1.size),'k-.',color='k',lw=0.5)
ax.text(np.mean(bin[it1]),1.15*mu1,'Change = ' + str(np.round(mu1-mu0,decimals=1)) + '\n(' + str(np.round( (mu1-mu0)/mu0*100 ,decimals=1)) + '%)',ha='center',fontsize=7,style='normal',color='k')
ax.set(ylabel='Tree carbon (tC ha$^1$)',xlabel='Time since detection, years',ylim=[0,1.07*np.nanmax(prfl[v]['mu'])])
#ax[0].legend(loc='upper right',frameon=False,facecolor='w',edgecolor='w')
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])
plt.tight_layout()
#gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\QA_Wildfire_Profiles','png',900)

#%% Mortality regression analysis

# Modelled harvest mortality is too high because salvage harvest is being counted as mortality
# even though something else killed most of the trees.

gplt2=copy.deepcopy(gplt)
for k in gplt2.keys():
    try:
        gplt2[k]=gplt2[k][(gplt['PTF CNY']==1)]
    except:
        pass

# Look at number of occurrences
print(np.where(gplt2['IBB Occurrence']==1)[0].size)
print(np.where(gplt2['IBD Occurrence']==1)[0].size)
print(np.where(gplt2['IBM Occurrence']==1)[0].size)
print(np.where(gplt2['IBS Occurrence']==1)[0].size)
print(np.where(gplt2['Wildfire Occurrence']==1)[0].size)
print(np.where(gplt2['Harv Occurrence']==1)[0].size)

df=pd.DataFrame.from_dict(gplt2)
df.columns=df.columns.str.replace(' ','_')
#x=df[['debt_ratio','industry']]
#y=df['cash_flow']

#frm='Ctot_Mort ~ Age_t0 + C(IBB_Occurrence) + C(IBD_Occurrence) + C(IBM_Occurrence) + C(IBS_Occurrence) + C(IDW_Occurrence) + C(Wildfire_Occurrence) + C(Harv_Occurrence)'
frm1='Ctot_Mort ~ Ctot_L_t0 + C(IBB_Occurrence) + C(IBM_Occurrence) + C(IBS_Occurrence) + C(IDW_Occurrence) + C(Wildfire_Occurrence) + C(Harv_Occurrence)'
md1=smf.ols(formula=frm1,data=df)
mr1=md1.fit()
print(mr1.summary())

frm2='Mod_C_M_Tot ~ Mod_C_Biomass_Tot_t0 + C(IBB_Occurrence) + C(IBM_Occurrence) + C(IBS_Occurrence) + C(IDW_Occurrence) + C(Wildfire_Occurrence) + C(Harv_Occurrence)'
md2=smf.ols(formula=frm2,data=df)
mr2=md2.fit()
print(mr2.summary())

b1=mr1.params.values
b2=mr2.params.values
se1=mr1.bse.values
se2=mr2.bse.values
lab=['IBB','IBM','IBS','IDW','Wildfire','Harvest'] #'Int', ,'Initial\nbiomass'

# Correct harvest to account for salvage
b2[-2]=0.65*b2[-2]

barw=0.32
cl=np.array([[0.9,0.9,1.0],[0.8,0.9,0.8]])

plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(11,9))
ax.plot([-1,b1.size+1],[0,0],'k-',lw=0.5)
ax.bar(np.arange(b1.size)-barw/2-0.01,b1,barw,facecolor=cl[0,:],label='Observations')
ax.bar(np.arange(b1.size)+barw/2+0.01,b2,barw,facecolor=cl[1,:],label='Predictions (FCS)')
ax.errorbar(np.arange(b1.size)-barw/2-0.01,b1,yerr=se1,color=gp['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)
ax.errorbar(np.arange(b1.size)+barw/2+0.01,b2,yerr=se2,color=gp['cla'],fmt='none',capsize=1.5,lw=0.25,markeredgewidth=0.5)

ax.set(position=[0.08,0.12,0.9,0.86],xticks=np.arange(1,b1.size-1,1),xticklabels=lab,
       ylabel='Regression coefficient (MgC ha$^{-1}$ yr$^{-1}$)',
       xlim=[-0.5+1,b1.size-1-0.5],ylim=[-4,10.5])
plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])
#gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\QA_Mortality_Regression','png',900)

#%% Mountain Pine Beetle

# Import AOS data
ind=gis.GetGridIndexToPoints(zRef,gplt['X'],gplt['Y'])
tvMPB=np.arange(1999,2022,1)
ios=np.zeros( (tvMPB.size,ind[0].size) ,dtype='int8')
for iT in range(tvMPB.size):
    z=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\PEST_INFESTATION_POLY_IBM_SeverityClass_' + str(tvMPB[iT]) + '.tif')
    ios[iT,:]=z['Data'][ind].copy()

# Calculate average outbreak profiles

ios_mx=np.max(ios,axis=0)
ind=np.where(ios_mx>=2)[0]
sts={}
sts['tso']=np.zeros(ind.size)
for v in gplt.keys():
    sts[v]=np.zeros(ind.size)
    for iStand in range(ind.size):
        indOutbreak=np.where(ios[:,ind[iStand]]>0)[0]
        sts['tso'][iStand]=gplt['Year t0'][ind[iStand]]-tvMPB[indOutbreak[0]]
        sts[v][iStand]=gplt[v][ind[iStand]]

bw=4; bin=np.arange(-10,20,bw)
prfl={}
for k in sts.keys():
    prfl[k]={}
    prfl[k]['N'],prfl[k]['mu'],prfl[k]['med'],prfl[k]['sig'],prfl[k]['se']=gu.discres(sts['tso'],sts[k],bw,bin)

# Plot

it0=np.where(bin<0)[0]
it1=np.where(bin>=10)[0]
it2=np.where(bin>0)[0]

# Plot observed profile

plt.close('all')
fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(15.5,7)); ms=4; cl=np.array([[0.8,0.8,0.8],[0.27,0.49,0.77],[0.75,0,0],[0.67,0.89,0.95],[1,0.6,0.6]])
ax[0].plot([0,0],[0,130],'k-',lw=2,color=[0.75,0.75,0.75])
ax[0].plot(bin,prfl['Ctot L t0']['mu']+prfl['Cdw t0']['mu'],'k^',ms=ms,mfc=cl[0,:],mec=cl[0,:],label='Live + dead')

mu0=np.mean(prfl['Ctot L t0']['mu'][it0])
mu1=np.mean(prfl['Ctot L t0']['mu'][it1])
ax[0].plot(bin[it0],mu0*np.ones(it0.size),'k-',color=cl[3,:],lw=2)
ax[0].plot(bin[it1],mu1*np.ones(it1.size),'k-',color=cl[3,:],lw=2)
ax[0].text(np.mean(bin[it1]),mu1+5,str(np.round(mu1-mu0,decimals=1)) + ' (' + str(np.round( (mu1-mu0)/mu0*100 ,decimals=1)) + '%)',ha='center',fontsize=10,style='normal',weight='bold',color=cl[1,:])

mu0=np.mean(prfl['Cdw t0']['mu'][it0])
ax[0].plot(bin[it0],mu0*np.ones(it0.size),'k-',color=cl[4,:],lw=2)
it1b=np.where(prfl['Cdw t0']['mu']==np.max(prfl['Cdw t0']['mu']) )[0]
mu1=prfl['Cdw t0']['mu'][it1b[0]]
ax[0].plot(bin[it2],mu1*np.ones(it2.size),'k-',color=cl[4,:],lw=2)
ax[0].text(np.mean(bin[it2]),mu1+5,str(np.round(mu1-mu0,decimals=1)) + ' (' + str(np.round( (mu1-mu0)/mu0*100 ,decimals=1)) + '%)',ha='center',fontsize=10,style='normal',weight='bold',color=cl[2,:])

ax[0].plot(bin,prfl['Ctot L t0']['mu'],'go',ms=ms,mfc=cl[1,:],mec=cl[1,:],label='Live')
ax[0].plot(bin,prfl['Cdw t0']['mu'],'rs',ms=ms,mfc='w',mec=cl[2,:],label='Dead')

ax[0].set(ylim=[0,140],ylabel='Tree carbon (MtCO$_2$e yr$^-$$^1$)',xlabel='Time since detection, years')
ax[0].legend(loc='upper right',frameon=False,facecolor='w',edgecolor='w')
ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].tick_params(length=gp['tickl'])

# Plot modelled profile
#fig,ax=plt.subplots(1,figsize=gu.cm2inch(8.5,7.5)); ms=4; cl=np.array([[0.8,0.8,0.8],[0.27,0.49,0.77],[0.75,0,0],[0.67,0.89,0.95],[1,0.6,0.6]])
ax[1].plot([0,0],[0,250],'k-',lw=2,color=[0.75,0.75,0.75])
ax[1].plot(bin,prfl['Mod C_Biomass_Tot t0']['mu']+prfl['Mod C_DeadWood_Tot t0']['mu'],'k^',ms=ms,mfc=cl[0,:],mec=cl[0,:],label='Live + dead')

mu0=np.mean(prfl['Mod C_Biomass_Tot t0']['mu'][it0])
mu1=np.mean(prfl['Mod C_Biomass_Tot t0']['mu'][it1])
ax[1].plot(bin[it0],mu0*np.ones(it0.size),'k-',color=cl[3,:],lw=2)
ax[1].plot(bin[it1],mu1*np.ones(it1.size),'k-',color=cl[3,:],lw=2)
ax[1].text(np.mean(bin[it1]),mu1+5,str(np.round(mu1-mu0,decimals=1)) + ' (' + str(np.round( (mu1-mu0)/mu0*100 ,decimals=1)) + '%)',ha='center',fontsize=10,style='normal',weight='bold',color=cl[1,:])

mu0=np.mean(prfl['Mod C_DeadWood_Tot t0']['mu'][it0])
ax[1].plot(bin[it0],mu0*np.ones(it0.size),'k-',color=cl[4,:],lw=2)
it1b=np.where(prfl['Mod C_DeadWood_Tot t0']['mu']==np.max(prfl['Mod C_DeadWood_Tot t0']['mu']) )[0]
mu1=prfl['Mod C_DeadWood_Tot t0']['mu'][it1b[0]]
ax[1].plot(bin[it2],mu1*np.ones(it2.size),'k-',color=cl[4,:],lw=2)
ax[1].text(np.mean(bin[it2]),mu1+5,str(np.round(mu1-mu0,decimals=1)) + ' (' + str(np.round( (mu1-mu0)/mu0*100 ,decimals=1)) + '%)',ha='center',fontsize=10,style='normal',weight='bold',color=cl[2,:])

ax[1].plot(bin,prfl['Mod C_Biomass_Tot t0']['mu'],'go',ms=ms,mfc=cl[1,:],mec=cl[1,:],label='Live')
ax[1].plot(bin,prfl['Mod C_DeadWood_Tot t0']['mu'],'rs',ms=ms,mfc='w',mec=cl[2,:],label='Dead')

ax[1].set(ylim=[0,140],ylabel='Tree carbon (MtCO$_2$e yr$^-$$^1$)',xlabel='Time since detection, years')
ax[1].legend(loc='upper right',frameon=False,facecolor='w',edgecolor='w')
ax[1].yaxis.set_ticks_position('both'); ax[1].xaxis.set_ticks_position('both'); ax[1].tick_params(length=gp['tickl'])

gu.axletters(ax,plt,0.045,0.92,FontColor=gp['cla'],LetterStyle='Caps',FontWeight='Bold')
plt.tight_layout()
gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\QA_IBM_Profiles','png',900)


#%% Plot Stocks by BGC zone

vL=['Age VRI t0','Cbk L t0','Cbr L t0','Cf L t0','Csw L t0','Cr L t0','Cag L t0','Ctot L t0','Ctot G Surv','Ctot G Recr','Ctot Mort+Lost','Ctot Net']

u=np.unique(gpt['Ecozone BC L1'])
u=u[u>0]
lab=np.array(['' for _ in range(u.size)],dtype=object)

d={}
for v in vL:
    d[v]={}
    d[v]['N']=np.zeros(u.size)
    d[v]['mu']=np.zeros(u.size)
    d[v]['sd']=np.zeros(u.size)
    d[v]['se']=np.zeros(u.size)
data=[None]*u.size
for i in range(u.size):
    lab[i]=ugp.lut_id2cd(meta['LUT']['GP']['Ecozone BC L1'],u[i])
    for v in vL:
        ind=np.where( (gpt['Ecozone BC L1']==u[i]) &
                     (gpt['PTF CNV']==1) &
                     (gpt['Cbk L t0']>=0) & (gpt['Cbk L t0']<2000) &
                     (gpt['Cbr L t0']>=0) & (gpt['Cbr L t0']<2000) &
                     (gpt['Cf L t0']>=0) & (gpt['Cf L t0']<2000) &
                     (gpt['Cr L t0']>=0) & (gpt['Cr L t0']<2000) &
                     (gpt['Csw L t0']>=0) & (gpt['Csw L t0']<2000) &
                     (gpt['Ctot L t0']>=0) & (gpt['Ctot L t0']<10000))[0]
        d[v]['N'][i]=ind.size
        d[v]['mu'][i]=np.nanmean(gpt[v][ind])
        d[v]['sd'][i]=np.nanstd(gpt[v][ind])
        #d[v]['se'][i]=np.nanstd(gpt[v][ind])/np.sqrt(ind[0].size)
    ind=np.where( (gpt['Ecozone BC L1']==u[i]) & (gpt['PTF CNV']==1) & (gpt['Ctot L t0']>=0) & (gpt['Ctot L t0']<10000))[0]
    data[i]=gpt['Ctot L t0'][ind]

# Put in order
d['Ctot L t0']['mu']=d['Cbk L t0']['mu']+d['Cbr L t0']['mu']+d['Cf L t0']['mu']+d['Cr L t0']['mu']+d['Csw L t0']['mu']
ord=np.argsort(d['Ctot L t0']['mu'])
uo=u[ord]
lab=np.flip(lab[ord])
for v in d:
    for k in d[v].keys():
        d[v][k]=np.flip(d[v][k][ord])
data2=[None]*len(data)
for i in range(len(data)):
    data2[i]=data[np.flip(ord)[i]]

# Plot
cl=np.array([[0.65,0.85,0.05],[0.75,0.5,0.95],[0.6,0.85,1],[0,0.6,0],[0.75,0.5,0]])
plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,8))
ax.bar(np.arange(u.size),d['Csw L t0']['mu'],facecolor=cl[0,:],label='Stemwood')
ax.bar(np.arange(u.size),d['Cbk L t0']['mu'],bottom=d['Csw L t0']['mu'],facecolor=cl[1,:],label='Bark')
ax.bar(np.arange(u.size),d['Cbr L t0']['mu'],bottom=d['Csw L t0']['mu']+d['Cbk L t0']['mu'],facecolor=cl[2,:],label='Branches')
ax.bar(np.arange(u.size),d['Cf L t0']['mu'],bottom=d['Csw L t0']['mu']+d['Cbk L t0']['mu']+d['Cbr L t0']['mu'],facecolor=cl[3,:],label='Foliage')
ax.bar(np.arange(u.size),d['Cr L t0']['mu'],bottom=d['Csw L t0']['mu']++d['Cbk L t0']['mu']+d['Cbr L t0']['mu']+d['Cf L t0']['mu'],facecolor=cl[4,:],label='Roots')
ax.set(position=[0.08,0.065,0.9,0.92],yticks=np.arange(0,550,20),xticks=np.arange(u.size),
       xticklabels=lab,ylabel='Biomass (MgC ha$^{-1}$)',xlim=[-0.5,u.size-0.5],ylim=[0,180])
plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])

flg=0
if flg==1:
    vio=ax.violinplot(data2,np.arange(u.size),widths=0.7,showmeans=False,showextrema=True,showmedians=True)
    for pc in vio['bodies']:
        pc.set_facecolor('none')
        pc.set_edgecolor([0.25,0.25,0.25])
        pc.set_linewidth(0.5)
        #pc.set_alpha(1)
    for partname in ('cbars','cmins','cmaxes','cmedians'):
        vp=vio[partname]
        vp.set_edgecolor([0.75,0.75,0.75])
        vp.set_linewidth(0.5)
        vp.set_alpha(0.75)

    vp=vio['cbars']
    vp.set_alpha(0)

    vp=vio['cmedians']
    vp.set_edgecolor([0.75,0.75,0.75])
    vp.set_linewidth(2)
    #vp.set_alpha(0.75)

for i in range(u.size):
    ax.text(i,6,str(d['Csw L t0']['N'][i].astype(int)),color='k',ha='center',fontsize=7,fontweight='normal')

gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_Biomass_ByGBCZone_CNV','png',900)

#%% Plot net biomass production by BGC zone (CN)

vL=['Ctot G Surv','Ctot G Recr','Ctot Mort Nat','Ctot Mort Harv','Ctot Net']

u=np.unique(gpt['Ecozone BC L1'])
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
    lab[i]=ugp.lut_id2cd(meta['LUT']['GP']['Ecozone BC L1'],u[i])
    for v in vL:
        ind=np.where( (gpt['Ecozone BC L1']==u[i]) &
                     (gpt['PTF CN']==1) &
                     (gpt['Ctot Net']>=-1000) & (gpt['Ctot Net']<1000) )[0]

                     #(gpt['Ctot Mort Harv']==0) &
        d[v]['N'][i]=ind.size
        d[v]['mu'][i]=np.nanmean(gpt[v][ind])
        d[v]['sd'][i]=np.nanstd(gpt[v][ind])
        d[v]['se'][i]=np.nanstd(gpt[v][ind])/np.sqrt(ind.size)

# Remove classes with inadequate data
ind=np.where(d['Ctot Net']['N']>=3)[0]
for v in vL:
    for k in d[v].keys():
        d[v][k]=d[v][k][ind]
u=u[ind]
lab=lab[ind]

# Put in order
ord=np.argsort(d['Ctot Net']['mu'])
uo=u[ord]
lab=np.flip(lab[ord])
for v in d:
    for k in d[v].keys():
        d[v][k]=np.flip(d[v][k][ord])

# Area weighting
d['Area']=np.zeros(lab.size)
for i in range(lab.size):
    ind2=np.where( (z['bgcz']['Data']==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'][lab[i]]) & (z['lcc1_c']['Data']==meta['LUT']['Derived']['lcc1']['Forest']) )
    d['Area'][i]=ind2[0].size/1e6

wa=np.sum(d['Ctot Net']['mu']*d['Area'])/np.sum(d['Area'])

lab2=np.append(lab,'Area\nweighted')

plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(7.8,6)); cl=np.array([[0.75,0.75,0.75],[0.5,0.5,0.5],[1,0.5,0],[0.45,0.75,1],[0.6,1,0]])
ax.bar(np.arange(u.size),d['Ctot Net']['mu'],facecolor=cl[0,:],label='')
for i in range(u.size):
    ax.errorbar(i,d['Ctot Net']['mu'][i],yerr=d['Ctot Net']['se'][i],color=gp['cla'],fmt='none',capsize=2,lw=0.5)
ax.bar(u.size,wa,facecolor=cl[1,:],label='Weighted average')
ax.errorbar(u.size,wa,yerr=np.mean(d['Ctot Net']['se']),color=gp['cla'],fmt='none',capsize=2,lw=0.5)
ax.plot([-1,20],[0,0],'-k',color=gp['cla'],lw=0.5)
ax.set(yticks=np.arange(-5,3.5,0.5),xticks=np.arange(u.size+1),xticklabels=lab2,
       ylabel='Net biomass production (MgC ha$^{-1}$ yr$^{-1}$)',
       xlim=[-0.5,u.size+1-0.5],ylim=[-1.5,3])
#plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])
for i in range(u.size):
    if d['Ctot Net']['mu'][i]>0:
        adj=0.1
    else:
        adj=-0.1
    ax.text(i+0.24,d['Ctot Net']['mu'][i]+adj,str(d['Ctot Net']['N'][i].astype(int)),color='k',ha='center',va='center',fontsize=5,fontweight='normal')
plt.tight_layout()
gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_NetGrowth_ByGBCZone_CN','png',900)

#%% Plot net biomass production by BGC zone (YSM)

vL=['Ctot G Surv','Ctot G Recr','Ctot Mort+Lost','Ctot Mort+Lost Harv','Ctot Net']

u=np.unique(gpt['Ecozone BC L1'])
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
    lab[i]=ugp.lut_id2cd(meta['LUT']['GP']['Ecozone BC L1'],u[i])
    for v in vL:
        ind=np.where( (gpt['Ecozone BC L1']==u[i]) &
                     (gpt['PTF YSM']==1) &
                     (gpt['Ctot Net']>=-1000) & (gpt['Ctot Net']<1000) )[0]

                     #(gpt['Ctot Mort Harv']==0) &
        d[v]['N'][i]=ind.size
        d[v]['mu'][i]=np.nanmean(gpt[v][ind])
        d[v]['sd'][i]=np.nanstd(gpt[v][ind])
        d[v]['se'][i]=np.nanstd(gpt[v][ind])/np.sqrt(ind.size)

# Remove classes with inadequate data
ind=np.where(d['Ctot Net']['N']>=3)[0]
for v in vL:
    for k in d[v].keys():
        d[v][k]=d[v][k][ind]
u=u[ind]
lab=lab[ind]

# Put in order
ord=np.argsort(d['Ctot Net']['mu'])
uo=u[ord]
lab=np.flip(lab[ord])
for v in d:
    for k in d[v].keys():
        d[v][k]=np.flip(d[v][k][ord])

# Area weighting
d['Area']=np.zeros(lab.size)
for i in range(lab.size):
    ind2=np.where( (z['bgcz']['Data']==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'][lab[i]]) & (z['lcc1_c']['Data']==meta['LUT']['Derived']['lcc1']['Forest Land']) )
    d['Area'][i]=ind2[0].size/1e6

wa=np.sum(d['Ctot Net']['mu']*d['Area'])/np.sum(d['Area'])

lab2=np.append(lab,'Area\nweighted')

# Plot

plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,9)); cl=np.array([[0.75,0.75,0.75],[0.5,0.5,0.5],[1,0.5,0],[0.45,0.75,1],[0.6,1,0]])
ax.bar(np.arange(u.size),d['Ctot Net']['mu'],facecolor=cl[0,:],label='')
for i in range(u.size):
    ax.errorbar(i,d['Ctot Net']['mu'][i],yerr=d['Ctot Net']['se'][i],color=gp['cla'],fmt='none',capsize=2,lw=0.5)
ax.bar(u.size,wa,facecolor=cl[1,:],label='Weighted average')
ax.errorbar(u.size,wa,yerr=np.mean(d['Ctot Net']['se']),color=gp['cla'],fmt='none',capsize=2,lw=0.5)
ax.plot([-1,20],[0,0],'-k',color=gp['cla'],lw=0.5)
ax.set(position=[0.08,0.075,0.9,0.91],yticks=np.arange(-5,10,0.5),xticks=np.arange(u.size+1),xticklabels=lab2,
       ylabel='Net biomass production (MgC ha$^{-1}$ yr$^{-1}$)',
       xlim=[-0.5,u.size+1-0.5],ylim=[-0.5,5.25])
#plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])
for i in range(u.size):
    if d['Ctot Net']['mu'][i]>0:
        adj=0.08
    else:
        adj=-0.08
    ax.text(i+0.2,d['Ctot Net']['mu'][i]+adj,str(d['Ctot Net']['N'][i].astype(int)),color='k',ha='center',va='center',fontsize=7,fontweight='normal')
gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_NetGrowth_ByGBCZone_YSM','png',900)



#%% Age distributions

x=np.arange(0,401,1)
y=zAge['Data'][0::5,0::5].flatten()
y=y[y>=0]
kde=stats.gaussian_kde(y)
p1=kde(x)

plt.close('all'); lw=1.25
fig,ax=plt.subplots(1,figsize=gu.cm2inch(14.5,7))
#plt.plot(x,p1/np.sum(p1)*100,'b-',lw=lw,color=[0.27,0.49,0.77],label='VRI R1 Comp Polygons')

ind=np.where( (gpt['Plot Type']==meta['LUT']['GP']['Plot Type BC']['VRI']) & (gpt['Age Mean t0']>=0) )[0]
kde=stats.gaussian_kde(gpt['Age Mean t0'][ind])
p=kde(x)
plt.plot(x,p/np.sum(p)*100,'r-',lw=lw,color=[0.27,0.49,0.77],label='VRI ground plots')

ind=np.where( (gpt['PTF CN']==1) & (gpt['Age Mean t0']>=0) )[0]
kde=stats.gaussian_kde(gpt['Age Mean t0'][ind])
p=kde(x)
plt.plot(x,p/np.sum(p)*100,'r--',lw=lw,color=[1,0.5,0],label='CMI + NFI cores (mean)')

kde=stats.gaussian_kde(gpt['Age Min t0'][ind])
p=kde(x)
plt.plot(x,p/np.sum(p)*100,'r-.',lw=lw,color=[0.5,1,0],label='CMI + NFI cores (min)')

kde=stats.gaussian_kde(gpt['Age Max t0'][ind])
p=kde(x)
plt.plot(x,p/np.sum(p)*100,'r:',lw=lw,color=[0,0.4,0],label='CMI + NFI cores (max)')

# ind=np.where( (gpt['PTF YSM']==1) & (gpt['Age Mean t0']>=0) )[0]
# kde=stats.gaussian_kde(gpt['Age Mean t0'][ind])
# p=kde(x)
# plt.plot(x,p,'g--')
ax.set(ylabel='Frequency (%)',xlabel='Stand age, years',xlim=[0,400],ylim=[0,0.85])
plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])
plt.tight_layout()
gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_AgeDistribution','png',900)

#%% Age responses

# All
#ind=np.where( (gpt['PTF CN']==1) )[0]
#ind=np.where( (gpt['PTF YSM']==1) )[0]

# Coast
#ind=np.where( (gpt['PTF CN']==1) & (gpt['Ecozone BC L1']==meta['LUT']['GP']['Ecozone BC L1']['CWH']) | (gpt['PTF CN']==1) & (gpt['Ecozone BC L1']==meta['LUT']['GP']['Ecozone BC L1']['MH']) )[0]
ind=np.where( (gpt['PTF CN']==1) & (gpt['Ecozone BC L1']==meta['LUT']['GP']['Ecozone BC L1']['ICH']) )[0]
#ind=np.where( (gpt['PTF CN']==1) & (gpt['Ecozone BC L1']==meta['LUT']['GP']['Ecozone BC L1']['ESSF']) )[0]

# Interior
#ind=np.where( (gpt['PTF CN']==1) & (gpt['Ecozone BC L1']!=meta['LUT']['GP']['Ecozone BC L1']['CWH']) & (gpt['Ecozone BC L1']!=meta['LUT']['GP']['Ecozone BC L1']['MH']) & (gpt['Ecozone BC L1']!=meta['LUT']['GP']['Ecozone BC L1']['ICH']) )[0]

x=gpt['Age Med t0'][ind]
bw=25; bin=np.arange(bw,300,bw)
xhat=np.arange(1,301,1)

lw=0.75; ms=4; mec=[0.27,0.49,0.77]; mfc='w'; cl=[0.27,0.49,0.77];
plt.close('all'); fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(15.5,10))
y=gpt['Ctot G Surv'][ind]+gpt['Ctot G Recr'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
yhat=np.interp(xhat,bin,mu)
#ax[0,0].plot(xhat,yhat,'b-',lw=0.5)
ax[0,0].plot(bin,mu,'ko',ms=ms,lw=lw,color=cl,mfc=mfc,mec=mec)
ax[0,0].set(ylabel='Survivor growth (MgC ha$^-$$^1$ yr$^-$$^1$)',xlabel='Age, years',xlim=[0,300],ylim=[0,5])
ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=gp['tickl'])

y=gpt['Ctot G Recr'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
yhat=np.interp(xhat,bin,mu)
#ax[0,1].plot(xhat,yhat,'b-',lw=0.5)
ax[0,1].plot(bin,mu,'ko',ms=ms,lw=lw,color=cl,mfc=mfc,mec=mec)
ax[0,1].set(ylabel='Ingrowth (MgC ha$^-$$^1$ yr$^-$$^1$)',xlabel='Age, years',xlim=[0,300],ylim=[0,5])
ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=gp['tickl'])

y=gpt['Ctot Mort+Lost'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
yhat=np.interp(xhat,bin,mu)
#ax[1,0].plot(xhat,yhat,'b-',lw=0.5)
ax[1,0].plot(bin,mu,'ko',ms=ms,lw=lw,color=cl,mfc=mfc,mec=mec,label='Total')

y=gpt['Ctot Mort+Lost'][ind]-gpt['Ctot Mort+Lost Fire'][ind]-gpt['Ctot Mort+Lost Insect'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
yhat=np.interp(xhat,bin,mu)
#ax[1,0].plot(xhat,yhat,'g--',lw=0.5)
ax[1,0].plot(bin,mu,'ks',ms=ms,lw=lw,color='g',mfc=mfc,mec='g',label='W/O fire and insects')
ax[1,0].set(ylabel='Mortality (MgC ha$^-$$^1$ yr$^-$$^1$)',xlabel='Age, years',xlim=[0,300],ylim=[0,5])
ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=gp['tickl'])
leg=ax[1,0].legend(loc='upper left',frameon=False,facecolor=None,edgecolor='w')

ax[1,1].plot(xhat,0*xhat,'k-',lw=1.5,color=[0.8,0.8,0.8])
y=gpt['Ctot Net'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
yhat_tot=np.interp(xhat,bin[np.isnan(mu)==False],mu[np.isnan(mu)==False])
#ax[1,1].plot(xhat,yhat_tot,'b-',lw=0.5)
ax[1,1].plot(bin,mu,'ks',ms=ms,lw=lw,color=cl,mfc=mfc,mec=mec,label='Total')

y=gpt['Ctot Net'][ind]+gpt['Ctot Mort+Lost Fire'][ind]+gpt['Ctot Mort+Lost Insect'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
yhat_wofi=np.interp(xhat,bin[np.isnan(mu)==False],mu[np.isnan(mu)==False])
#ax[1,1].plot(xhat,yhat_wofi,'g--',lw=0.5)
ax[1,1].plot(bin,mu,'ko',ms=ms,lw=lw,color='g',mfc=mfc,mec='g',label='W/O fire and insects')

ax[1,1].set(ylabel='Net growth (MgC ha$^-$$^1$ yr$^-$$^1$)',xlabel='Age, years',xlim=[0,300],ylim=[-2,5])
ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=gp['tickl'])
leg=ax[1,1].legend(loc='upper right',frameon=False,facecolor=None,edgecolor='w')
gu.axletters(ax,plt,0.03,0.9,FontColor=gp['cla'],LetterStyle='Caps',FontWeight='Bold')
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Ground Plots\AgeResponse_ICH','png',900)

# Plot biomass

x=gpt['Age Med t0'][ind]
y=gpt['Ctot L t0'][ind]
bw=10; bin=np.arange(bw,300,bw)

lw=1; ms=5; mec='b'; mfc='w'; cl='b';
fig,ax=plt.subplots(1,figsize=gu.cm2inch(15.5,11))
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
ax.plot(bin,mu,'ko',ms=ms,lw=lw,color='k',mfc='w',mec='k')
ax.plot(xhat,np.cumsum(yhat_wofi),'g-',ms=ms,lw=lw,color='g',mfc=mfc,mec=mec)
ax.plot(xhat,np.cumsum(yhat_tot),'b--',ms=ms,lw=lw,color='b',mfc=mfc,mec=mec)
ax.set(ylabel='Biomass (MgC ha$^-$$^1$)',xlabel='Age, years',xlim=[0,300])
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])

#ind2=np.where( (gpt['PTF YSM']==1) )[0]
#x=gpt['Age Med t0'][ind2]
#y=gpt['Ctot L t0'][ind2]
#N,mu,med,sig,se=gu.discres(x,y,bw,bin)
#ax.plot(bin,mu,'ks',ms=ms,lw=lw,color='k',mfc='w',mec='r')


#%% Biomass age response in areas with no history of major disturbance

indG=gis.GetGridIndexToPoints(zRef,gpt['X'],gpt['Y'])

flgH=np.zeros(gpt['PTF CN'].size);
ind=np.where(zH['Data'][indG]>0)[0]
flgH[ind]=1

flgF=np.zeros(gpt['PTF CN'].size);
ind=np.where(zF['Data'][indG]>0)[0]
flgF[ind]=1

flgI=np.zeros(gpt['PTF CN'].size);
ind=np.where(zI['Data'][indG]>0)[0]
flgI[ind]=1

ind=np.where( (gpt['PTF CNV']==1) & (gpt['Ecozone BC L1']!=meta['LUT']['GP']['Ecozone BC L1']['CWH']) & (gpt['Ecozone BC L1']!=meta['LUT']['GP']['Ecozone BC L1']['MH']) & \
             (gpt['Ecozone BC L1']!=meta['LUT']['GP']['Ecozone BC L1']['ICH']) & \
            (flgH==0) & (flgF==0) & (flgI==0) )[0]

x1=gpt['Age VRI t0'][ind]
y1=gpt['Ctot L t0'][ind]
N,mu,med,sig,se=gu.discres(x1,y1,bw,bin)
ax[0,0].plot(bin,mu,'r^',ms=ms,lw=lw,color='k',mfc='w',mec='r')

ind=np.where( (gpt['PTF CN']==1) & (gpt['Ecozone BC L1']!=meta['LUT']['GP']['Ecozone BC L1']['CWH']) & (gpt['Ecozone BC L1']!=meta['LUT']['GP']['Ecozone BC L1']['MH']) & \
             (gpt['Ecozone BC L1']!=meta['LUT']['GP']['Ecozone BC L1']['ICH']) & (flgI==1) )[0]
x1=gpt['Age VRI t0'][ind]
y1=gpt['Ctot L t0'][ind]
N,mu,med,sig,se=gu.discres(x1,y1,bw,bin)
ax[0,0].plot(bin,mu,'rd',ms=ms,lw=lw,color='k',mfc='c',mec='c')

#%% Age responses

zone='IDF'
bw=25; bin=np.arange(bw,350+bw,bw)
xhat=np.arange(1,351,1)

lw=1; ms=5; mec='b'; mfc='w'; cl='b';
plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(15.5,10))

ind=np.where( (gpt['Ecozone BC L1']==meta['LUT']['GP']['Ecozone BC L1'][zone]) )[0]
x=gpt['Age Med t0'][ind]
y=gpt['Ctot L t0'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
yhat=np.interp(xhat,bin,mu)
#ax.plot(xhat,yhat,'b-',lw=0.5)
ax.plot(bin,mu,'bo',ms=ms,lw=lw,color='b',mfc=mfc,mec='b')

ind=np.where( (gpt['PTF CNY']==1) & (gpt['Ecozone BC L1']==meta['LUT']['GP']['Ecozone BC L1'][zone]) )[0]
x=gpt['Age Med t0'][ind]
y=gpt['Ctot L t0'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
yhat=np.interp(xhat,bin,mu)
#ax.plot(xhat,yhat,'r-',lw=0.5)
ax.plot(bin,mu,'rs',ms=ms,lw=lw,color='r',mfc=mfc,mec='r')

ax.set(ylabel='Survivor growth (MgC ha$^-$$^1$ yr$^-$$^1$)',xlabel='Age, years',xlim=[0,360])
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])

#%% Age responses (by climate class)

cc=2

ind=np.where( (gpt['PTF CN']==1) & (gpt['Climate Class']>=cc) )[0]

x=gpt['Age VRI t0'][ind]
bw=25; bin=np.arange(bw,300,bw)
xhat=np.arange(1,301,1)

lw=1; ms=5; mec='b'; mfc='w'; cl='b';
plt.close('all'); fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(15.5,11))
y=gpt['Ctot G Surv'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
ax[0,0].plot(bin,mu,'ko',ms=ms,lw=lw,color=cl,mfc=mfc,mec=mec)
yhat=np.interp(xhat,bin,mu)
ax[0,0].plot(xhat,yhat,'b-',lw=0.5)
ax[0,0].set(ylabel='Survivor growth (MgC ha$^-$$^1$ yr$^-$$^1$)',xlabel='Age, years',xlim=[0,300])
ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=gp['tickl'])

y=gpt['Ctot G Recr'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
ax[0,1].plot(bin,mu,'ko',ms=ms,lw=lw,color=cl,mfc=mfc,mec=mec)
ax[0,1].set(ylabel='Recruitment growth (MgC ha$^-$$^1$ yr$^-$$^1$)',xlabel='Age, years',xlim=[0,300])
ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=gp['tickl'])

y=gpt['Ctot Mort+Lost'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
ax[1,0].plot(bin,mu,'ko',ms=ms,lw=lw,color=cl,mfc=mfc,mec=mec,label='Total')
y=gpt['Ctot Mort+Lost'][ind]-gpt['Ctot Mort+Lost Fire'][ind]-gpt['Ctot Mort+Lost Insect'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
ax[1,0].plot(bin,mu,'ks',ms=ms,lw=lw,color='g',mfc=mfc,mec='g',label='W/O fire and insects')
ax[1,0].set(ylabel='Mortality (MgC ha$^-$$^1$ yr$^-$$^1$)',xlabel='Age, years',xlim=[0,300])
ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=gp['tickl'])
leg=ax[1,0].legend(loc='lower center',frameon=False,facecolor=None,edgecolor='w')

y=gpt['Ctot Net'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
ax[1,1].plot(bin,mu,'ks',ms=ms,lw=lw,color=cl,mfc=mfc,mec=mec,label='Total')
yhat_tot=np.interp(xhat,bin[np.isnan(mu)==False],mu[np.isnan(mu)==False])
ax[1,1].plot(xhat,yhat_tot,'b-',lw=0.5)

y=gpt['Ctot Net'][ind]+gpt['Ctot Mort+Lost Fire'][ind]+gpt['Ctot Mort+Lost Insect'][ind]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
ax[1,1].plot(bin,mu,'ko',ms=ms,lw=lw,color='g',mfc=mfc,mec='g',label='W/O fire and insects')
yhat_wofi=np.interp(xhat,bin[np.isnan(mu)==False],mu[np.isnan(mu)==False])
ax[1,1].plot(xhat,yhat_wofi,'g--',lw=0.5)

ax[1,1].set(ylabel='Net growth (MgC ha$^-$$^1$ yr$^-$$^1$)',xlabel='Age, years',xlim=[0,300])
ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=gp['tickl'])
leg=ax[1,1].legend(loc='upper right',frameon=False,facecolor=None,edgecolor='w')
gu.axletters(ax,plt,0.03,0.9,FontColor=gp['cla'],LetterStyle='Caps',FontWeight='Bold')

# Plot biomass

ind=np.where( (gpt['PTF CN']==1) & (gpt['Climate Class']>=cc) )[0]

x=gpt['Age VRI t0'][ind]
y=gpt['Ctot L t0'][ind]
bw=10; bin=np.arange(bw,300,bw)

lw=1; ms=5; mec='b'; mfc='w'; cl='b';
fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(15.5,11))
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
ax[0,0].plot(bin,mu,'ko',ms=ms,lw=lw,color='k',mfc='w',mec='k')
ax[0,0].plot(xhat,np.cumsum(yhat_wofi),'g-',ms=ms,lw=lw,color='g',mfc=mfc,mec=mec)
ax[0,0].plot(xhat,np.cumsum(yhat_tot),'b--',ms=ms,lw=lw,color='b',mfc=mfc,mec=mec)
ax[0,0].set(ylabel='Biomass (MgC ha$^-$$^1$)',xlabel='Age, years',xlim=[0,300])
ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=gp['tickl'])


#%%


#%% Modelling

#fAg=inline(['max(0,b(1).*(1+((b(2).*(x(:,1)./b(3)).^b(4)-1)./(exp(x(:,1)./b(3))))))'],'b','x'); % Weibull
#fAm=inline(['max(0,b(1).*(1+((b(2).*((x-b(5))./b(3)).^b(4)-1)./(exp((x-b(5))./b(3))))))'],'b','x');

def fun(x,a,b,c,d):
    yhat=d+a*((1-np.exp(-b*x))**(1/(1-c))) # From Dzierzon and Mason 2006
    return yhat
xhat=np.arange(0,100,1)

indM=np.where( (tl['DBH']>0) & (np.isnan(tl['H'])==True) | (tl['DBH']>0) & (tl['H']<=0) )[0]
indM.size/tl['DBH'].size

indG=np.where( (tl['DBH']>0) & (tl['H']>0) )[0]
indG.size/tl['DBH'].size

# Global model
ikp=np.where( (tl['DBH']>0) & (tl['H']>0) & (tl['Vital Status']==meta['LUT']['Vital Status']['Live']) & (tl['ID DA1']==meta['LUT']['Damage Agents']['None']) & (tl['Stature']==meta['LUT']['Stature']['Standing']) )[0]
x=tl['DBH'][ikp]
y=tl['H'][ikp]
popt_glob0=[26,0.1,0.66,2]
popt_glob,pcov=curve_fit(fun,x,y,popt_glob0)
#yhat=fun(xhat,popt_glob[0],popt_glob[1],popt_glob[2],popt_glob[3])
rs_glob,txt=gu.GetRegStats(x,y)



#%% Aridity

#ind=np.where( (gpt['PTF CN']==1) & (gpt['Ecozone BC L1']==meta['LUT']['GP']['Ecozone CA L1']['Pacific_Maritime']) )[0]
ind=np.where( (gpt['PTF CN']==1) )[0]
#ind=np.where( (gpt['PTF CN']==1) & (gpt['Ecozone BC L1']==meta['LUT']['GP']['Ecozone BC L1']['SBPS']) )[0]

x=gpt['Climate Class'][ind]
y=gpt['Ctot G Surv'][ind]

bw=1; bin=np.arange(1,6,bw)
N,mu,med,sig,se=gu.discres(x,y,bw,bin)

plt.close('all')
plt.plot(bin,mu,'bo')

