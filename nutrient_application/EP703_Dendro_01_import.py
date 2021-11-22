'''
EP703 RETROSPECTIVE MONITORING STUDY - IMPORT DATA
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

#%% Graphics parameters

fs=8
params={'font.sans-serif':'Arial','font.size':fs,'axes.edgecolor':'black','axes.labelsize':fs,'axes.labelcolor':'black','axes.titlesize':fs,'axes.linewidth':0.5,'lines.linewidth':0.5,
        'text.color':'black','xtick.color':'black','xtick.labelsize':fs,'xtick.major.width':0.5,'xtick.major.size':3,'xtick.direction':'in','ytick.color':'black','ytick.labelsize':fs,
        'ytick.major.width':0.5,'ytick.major.size':3,'ytick.direction':'in','legend.fontsize':fs,'savefig.dpi':300,'savefig.transparent':True}
plt.rcParams.update(params)

#%% Preapre project

meta={}
meta['Paths']={}
meta['Paths']['Inputs']=r'C:\Users\rhember\Documents\Data\TreeRings\EP703\Given\cross_dated_cores received 20211105'
meta['Paths']['Figures']=r'C:\Users\rhember\OneDrive - Government of BC\Figures\EP703\TreeRings'

#%%

# Import tree-level observations
dTL=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\TreeRings\EP703\Given\cross_dated_cores received 20211105\EP703_RM_Master_stacked_pith_est.xlsx')

# Import tree summary file 
dTS=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\TreeRings\EP703\Given\Field Data received 20210125\RM Tree Summary.xlsx')

# Import plot summary file
dPS=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\TreeRings\EP703\Given\Field Data received 20210125\RM Plot Summary.xlsx')

# Climate data
clm=gu.ImportMat(r'C:\Users\rhember\Documents\Data\TreeRings\EP703\Given\Climate\EP703_TR_Climate.mat','clm')

#%% Add derived variables

# Add unique tree ID
dTL['ID_Tree Unique']=np.zeros(dTL['TRW'].size,dtype=int)
cnt=0
u=np.unique(np.column_stack((dTL['ID_Inst'],dTL['ID_Plot'],dTL['ID_Tree'])),axis=0)
for iU in range(u.shape[0]):
    ind=np.where( (dTL['ID_Inst']==u[iU,0]) & (dTL['ID_Plot']==u[iU,1]) & (dTL['ID_Tree']==u[iU,2]) )[0]
    dTL['ID_Tree Unique'][ind]=cnt
    cnt=cnt+1

# Adjust units of BAI (cm2/yr)
dTL['BAI']=dTL['BAI']/100

# Add info from site and tree summary files and add some derived variables
dTL['BGC SS']=np.zeros(dTL['TRW'].size)
dTL['Dob 2020']=np.zeros(dTL['TRW'].size)
dTL['H 2020']=np.zeros(dTL['TRW'].size)
dTL['Year First']=np.zeros(dTL['TRW'].size)
dTL['Year Last']=np.zeros(dTL['TRW'].size)
dTL['Time since first ring']=np.zeros(dTL['TRW'].size)
dTL['Dib']=np.zeros(dTL['TRW'].size)
dTL['Dib Last']=np.zeros(dTL['TRW'].size)
dTL['Bsw']=np.zeros(dTL['TRW'].size)
dTL['Gsw']=np.zeros(dTL['TRW'].size)
dTL['RGR']=np.zeros(dTL['TRW'].size)
dTL['TRW Standardized RR']=np.zeros(dTL['TRW'].size)
dTL['BAI Standardized RR']=np.zeros(dTL['TRW'].size)
dTL['Bsw Standardized RR']=np.zeros(dTL['TRW'].size)
dTL['Gsw Standardized AD']=np.zeros(dTL['TRW'].size)
dTL['Gsw Standardized RR']=np.zeros(dTL['TRW'].size)
dTL['RGR Standardized RR']=np.zeros(dTL['TRW'].size)
dTL['Lat']=np.zeros(dTL['TRW'].size)
dTL['Lon']=np.zeros(dTL['TRW'].size)

uTID=np.unique(dTL['ID_Tree Unique'])
for iU in range(uTID.size):
    
    # Index to unique tree ID
    indUT=np.where(dTL['ID_Tree Unique']==uTID[iU])[0]
    
    # Define first and last year of measurement for each tree
    iWithData=np.where(dTL['Dib'][indUT]>-1.0)[0]
    dTL['Year First'][indUT]=dTL['Year'][indUT[iWithData[0]]]
    dTL['Year Last'][indUT]=dTL['Year'][indUT[iWithData[-1]]]
    dTL['Time since first ring'][indUT[iWithData]]=np.arange(1,iWithData.size+1,1)
    
    # Index to site summary info
    indPS=np.where( (dPS['ID_Inst']==dTL['ID_Inst'][indUT[0]]) & (dPS['ID_Plot']==dTL['ID_Plot'][indUT[0]]) )[0]
    
    # Site series
    dTL['BGC SS'][indUT]=dPS['BGC site series'][indPS]
    dTL['Lat'][indUT]=dPS['Lat'][indPS]
    dTL['Lon'][indUT]=dPS['Lon'][indPS]
    
    # Index to tree summary file
    indTS=np.where( (dTS['ID_Inst']==dTL['ID_Inst'][indUT[0]]) & (dTS['ID_Plot']==dTL['ID_Plot'][indUT[0]]) & (dTS['ID_Tree']==dTL['ID_Tree'][indUT[0]]) )[0]
        
    if indTS.size>0:
        try:
            dTL['Dob 2020'][indUT]=dTS['DBH_2020 (cm)'][indTS].astype(np.float)   
        except:
            print('Something odd about DBH 2020')
        dTL['H 2020'][indUT]=dTS['Ht_2020 (m)'][indTS][0]
    else:
        print('No crosswalk found with tree summary')
        continue

    # Inside-bark diameter (cm)
    dTL['Dib'][indUT]=0.1*np.cumsum(2*np.nan_to_num(dTL['TRW'][indUT]))
    
    # Inside-bark diameter for last year (cm) 
    # *** Only populated if last year >= 2019 (for comparison with 2020 tape measurement)
    if dTL['Year'][indUT[iWithData[-1]]]>=2019:
        dTL['Dib Last'][indUT[iWithData[-1]]]=dTL['Dib'][indUT[iWithData[-1]]]
    
    # Stemwood biomass from allometric function (kgDM) - DBH only (Ung et al. 2008)
    dTL['Bsw'][indUT]=0.0204*dTL['Dib'][indUT]**2.6974
    
    # Stemwood biomass growth (kgDM/yr)
    dTL['Gsw'][indUT[1:]]=np.diff(dTL['Bsw'][indUT])
    
    # Relative growth rate
    dTL['RGR'][indUT[1:]]=np.log(dTL['Bsw'][indUT[1:]])-np.log(dTL['Bsw'][indUT[0:-1]])

    # Standardized growth relative to mean growth during a period leading up to N application
    ind0=np.where( (dTL['ID_Tree Unique']==uTID[iU]) & (dTL['Year']>=1971-5) & (dTL['Year']<=1970) )[0]
    dTL['TRW Standardized RR'][indUT]=dTL['TRW'][indUT]/np.mean(dTL['TRW'][ind0])
    dTL['BAI Standardized RR'][indUT]=dTL['BAI'][indUT]/np.mean(dTL['BAI'][ind0])
    dTL['Bsw Standardized RR'][indUT]=dTL['Bsw'][indUT]/np.mean(dTL['Bsw'][ind0])
    dTL['Gsw Standardized AD'][indUT]=dTL['Gsw'][indUT]-np.mean(dTL['Gsw'][ind0])
    dTL['Gsw Standardized RR'][indUT]=dTL['Gsw'][indUT]/np.maximum(0.1,np.mean(dTL['Gsw'][ind0]))
    dTL['RGR Standardized RR'][indUT]=dTL['RGR'][indUT]/np.mean(dTL['RGR'][ind0])
    
# Change name of treatment
dTL['Treatment']=dTL['Treatment.EP703']

# Delete fields
del dTL['Unnamed: 0'],dTL['Treatment.EP703']

# Summarize QA flags
dTL['QA Summary']=np.zeros(dTL['TRW'].size)
ind=np.where( (dTL['QC_xdate_flags']=='Good_correlation') | (dTL['QC_xdate_flags']=='Weak_correlation') )[0]
dTL['QA Summary'][ind]=1
print('% retained = ' + str(100*np.where(dTL['QA Summary']==1)[0].size/dTL['QA Summary'].size))

# Define unique installations
uInst=np.unique(dTL['ID_Inst'])

# Define unique site series
uSS=np.unique(dTL['BGC SS'])

#%% Add climate data

# Define unique coordinates
uLL=np.unique(np.column_stack((dTL['ID_Inst'],dTL['ID_Plot'],dTL['Lat'],dTL['Lon'])),axis=0)

# Export to process climate data
flg=0
if flg==1:
    df=pd.DataFrame(uLL)
    df.to_excel(r'G:\My Drive\EP703 TR coordinates.xlsx')

# Define base period for normals
iTn_ann=np.where( (clm['tv'][:,0]>=1971) & (clm['tv'][:,0]<=2000) )[0]
iTn_ws=np.where( (clm['tv'][:,0]>=1971) & (clm['tv'][:,0]<=2000) & (clm['tv'][:,1]>=5) & (clm['tv'][:,0]<=9) )[0]

yrs=np.arange(1950,2022,1)

dTL['tmean_ann_n']=np.zeros(dTL['TRW'].size)
dTL['prcp_ann_n']=np.zeros(dTL['TRW'].size)
dTL['tmean_gs_r']=np.zeros(dTL['TRW'].size)
dTL['prcp_gs_r']=np.zeros(dTL['TRW'].size)
dTL['cwd_gs_r']=np.zeros(dTL['TRW'].size)
dTL['ws_gs_r']=np.zeros(dTL['TRW'].size)
for iLL in range(uLL.shape[0]):
    ind=np.where( (dTL['Lat']==uLL[iLL,2]) & (dTL['Lon']==uLL[iLL,3]) )[0]
    dTL['tmean_ann_n'][ind]=np.mean(clm['tmean'][iTn_ann,iLL])
    dTL['prcp_ann_n'][ind]=np.sum(clm['prcp'][iTn_ann,iLL])/(iTn_ann.size/12)
    for iY in range(yrs.size):
        iT=np.where( (clm['tv'][:,0]==yrs[iY]) & (clm['tv'][:,1]>=4) & (clm['tv'][:,1]<=9) )[0]
        iT2=np.where(dTL['Year'][ind]==yrs[iY])[0]
        dTL['tmean_gs_r'][ind[iT2]]=np.mean(clm['tmean'][iT,iLL])
        dTL['prcp_gs_r'][ind[iT2]]=np.mean(clm['prcp'][iT,iLL])
        dTL['cwd_gs_r'][ind[iT2]]=np.mean(clm['etp_tmw'][iT,iLL]-clm['eta_tmw'][iT,iLL])
        dTL['ws_gs_r'][ind[iT2]]=np.mean(clm['ws_tmw'][iT,iLL])

uC=np.unique(np.column_stack((dTL['tmean_ann_n'],dTL['prcp_ann_n'])),axis=0)

#plt.plot(uC[:,0],uC[:,1],'.')

#%% Get stand-level information from EP703 dataset

sobs=gu.ipickle(r'C:\Users\rhember\Documents\Data\EP703\Processed\EP703_SL.pkl')
sobs['RGR']=(np.log(sobs['Bsw_t1'])-np.log(sobs['Bsw_t0']))/sobs['DT']


dByInst={}
dByInst['DA_BA_G']=np.zeros(uInst.size)
dByInst['DR_BA_G']=np.zeros(uInst.size)
dByInst['DA_Hobs_G']=np.zeros(uInst.size)
dByInst['DR_Hobs_G']=np.zeros(uInst.size)
dByInst['DA_Bsw_G']=np.zeros(uInst.size)
dByInst['DR_Bsw_G']=np.zeros(uInst.size)
dByInst['DA_RGR']=np.zeros(uInst.size)
dByInst['DR_RGR']=np.zeros(uInst.size)
for iInst in range(uInst.size):
    iC=np.where( (sobs['ID_Site']==uInst[iInst]) & (sobs['N_Dose']==0) & (sobs['TSF_t0']>=0) & (sobs['TSF_t1']<=9) )[0]
    iF=np.where( (sobs['ID_Site']==uInst[iInst]) & (sobs['N_Dose']==225) & (sobs['TSF_t0']>=0) & (sobs['TSF_t1']<=9) )[0]
    dByInst['DA_BA_G'][iInst]=np.mean(sobs['BA_G'][iF]-sobs['BA_G'][iC])
    dByInst['DR_BA_G'][iInst]=np.mean( (sobs['BA_G'][iF]-sobs['BA_G'][iC])/sobs['BA_G'][iC]*100 )
    dByInst['DA_Hobs_G'][iInst]=np.nanmean(sobs['H_obs_G'][iF]-sobs['H_obs_G'][iC])
    dByInst['DR_Hobs_G'][iInst]=np.nanmean( (sobs['H_obs_G'][iF]-sobs['H_obs_G'][iC])/sobs['H_obs_G'][iC]*100 )
    dByInst['DA_Bsw_G'][iInst]=np.mean(sobs['Bsw_G'][iF]-sobs['Bsw_G'][iC])
    dByInst['DR_Bsw_G'][iInst]=np.mean( (sobs['Bsw_G'][iF]-sobs['Bsw_G'][iC])/sobs['Bsw_G'][iC]*100 )
    dByInst['DA_RGR'][iInst]=np.median(sobs['RGR'][iF]-sobs['RGR'][iC])
    dByInst['DR_RGR'][iInst]=np.median( (sobs['RGR'][iF]-sobs['RGR'][iC])/sobs['RGR'][iC]*100 )
    
    #plt.close('all')
    #plt.plot(sobs['TSF_t0'][iC],sobs['N_t0'][iC],'bo')
    #plt.plot(sobs['TSF_t0'][iF],sobs['N_t0'][iF],'rs')
    
    #plt.plot(sobs['TSF_t0'][iC],sobs['Bsw_G'][iC],'bo')
    #plt.plot(sobs['TSF_t0'][iF],sobs['Bsw_G'][iF],'rs')


#%% Comparison between Dob and Dib from cores

# Isolate data
ikp=np.where( (dTL['Dob 2020']>0) & (dTL['Dib Last']>0) )[0]
#ikp=np.where( (dTL['Dob 2020']>0) & (dTL['Dib Last']>0) & (dTL['Pith_status']=='Yes') )[0]

# Fit linear best fit relationship
y=dTL['Dib Last'][ikp]
x=dTL['Dob 2020'][ikp]
x1=sm.tools.tools.add_constant(x)
md=sm.OLS(y,x1).fit()
md.summary()
xhat=np.linspace(np.min(x1[:,1]),np.max(x1[:,1]),10)
yhat=md.predict(np.c_[np.ones(xhat.size),xhat])


# Plot
plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(8.5,8.5))
ax.plot([0,90],[0,90],'k-',lw=2,color=[0.8,0.8,0.8])
ax.plot(x,y,'o',ms=4,mec='w',mfc='k')
ax.plot(xhat,yhat,'r-',lw=1.25)
ax.set(xlim=[0,80],ylim=[0,80],xlabel='Tape-based Dob in 2020 (cm)',ylabel='Core-based Dib in 2019 (cm)')
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')
plt.tight_layout()
#gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\EP703\TreeRings\DiameterComparison','png',150)

#%% Look at individual trees

ind=np.where( (dTL['ID_Tree Unique']==13) )[0]

print(dTL['Est/mm'][ind[0]])

plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(8.5,8.5))
#ax.plot([0,90],[0,90],'k-',lw=2,color=[0.8,0.8,0.8])
#ax.plot(dTL['Time since first ring'][ind],dTL['TRW'][ind],'-ko',ms=3,mec='k',mfc='k')
ax.plot(dTL['Year'][ind],dTL['TRW'][ind],'-ko',ms=3,mec='k',mfc='k')

#plt.close('all')
#fig,ax=plt.subplots(1,figsize=gu.cm2inch(8.5,8.5))
#ax.plot([0,90],[0,90],'k-',lw=2,color=[0.8,0.8,0.8])
#ax.plot(dTL['TRW'][ind],dTL['Gsw'][ind],'-ko',ms=3,mec='k',mfc='k')