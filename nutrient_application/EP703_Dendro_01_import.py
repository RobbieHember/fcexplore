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
from matplotlib.patches import Rectangle

#%% Graphics parameters

fs=8
params={'font.sans-serif':'Arial','font.size':fs,'axes.edgecolor':'black','axes.labelsize':fs,'axes.labelcolor':'black','axes.titlesize':fs,'axes.linewidth':0.5,'lines.linewidth':0.5,
        'text.color':'black','xtick.color':'black','xtick.labelsize':fs,'xtick.major.width':0.5,'xtick.major.size':3,'xtick.direction':'in','ytick.color':'black','ytick.labelsize':fs,
        'ytick.major.width':0.5,'ytick.major.size':3,'ytick.direction':'in','legend.fontsize':fs,'savefig.dpi':300,'savefig.transparent':True}
plt.rcParams.update(params)

#%% Preapre project

meta={}
meta['Paths']={}
meta['Paths']['Figures']=r'C:\Users\rhember\OneDrive - Government of BC\Figures\EP703\TreeRings'

#%%

# Import tree-level observations
dTL=gu.ReadExcel(r'C:\Users\rhember\OneDrive - Government of BC\Shared\EP703\Retrospective Monitoring Study\Data\Given\Core Data\Received 20211105 from JBerg\EP703_RM_Master_stacked_pith_est.xlsx')

# Import tree summary file 
dTS=gu.ReadExcel(r'C:\Users\rhember\OneDrive - Government of BC\Shared\EP703\Retrospective Monitoring Study\Data\Given\Field Data\Received 20210125 from JBapty\RM Tree Summary.xlsx')

# Import plot summary file
dPS=gu.ReadExcel(r'C:\Users\rhember\OneDrive - Government of BC\Shared\EP703\Retrospective Monitoring Study\Data\Given\Field Data\Received 20210125 from JBapty\RM Plot Summary.xlsx')

# Climate data
clm=gu.ImportMat(r'C:\Users\rhember\OneDrive - Government of BC\Shared\EP703\Retrospective Monitoring Study\Data\Given\Climate\Received 20211024 from RHember\EP703_TR_Climate.mat','clm')

# Import EP703 tree-level data
dEPT=gu.ipickle(r'C:\Users\rhember\Documents\Data\EP703\Processed\EP703_TL.pkl')

#%% Custom functions

# Function (tree height as a function of basal area)
def funcH(x,a,b):
    return a*x**b

#%% Add derived variables

# Add unique tree ID
dTL['ID_Tree_Unique']=np.zeros(dTL['TRW'].size,dtype=int)
cnt=0
u=np.unique(np.column_stack((dTL['ID_Inst'],dTL['ID_Plot'],dTL['ID_Tree'])),axis=0)
for iU in range(u.shape[0]):
    ind=np.where( (dTL['ID_Inst']==u[iU,0]) & (dTL['ID_Plot']==u[iU,1]) & (dTL['ID_Tree']==u[iU,2]) )[0]
    dTL['ID_Tree_Unique'][ind]=cnt
    cnt=cnt+1

# Add unique plot ID
dTL['ID_Plot_Unique']=np.zeros(dTL['TRW'].size,dtype=int)
cnt=0
u=np.unique(np.column_stack( (dTL['ID_Inst'],dTL['ID_Plot']) ),axis=0)
for iU in range(u.shape[0]):
    ind=np.where( (dTL['ID_Inst']==u[iU,0]) & (dTL['ID_Plot']==u[iU,1]) )[0]
    dTL['ID_Plot_Unique'][ind]=cnt
    cnt=cnt+1

# Adjust units of BAI (cm2/yr)
dTL['BAI']=dTL['BAI']/100

# Add info from site and tree summary files and add some derived variables
dTL['BGC SS']=np.zeros(dTL['TRW'].size)
dTL['Dob 2020']=np.zeros(dTL['TRW'].size)
dTL['H 2020']=np.zeros(dTL['TRW'].size)
dTL['Year First']=np.zeros(dTL['TRW'].size)
dTL['Year Last']=np.zeros(dTL['TRW'].size)
dTL['Age']=np.zeros(dTL['TRW'].size)
dTL['Dib']=np.zeros(dTL['TRW'].size)
dTL['Dib Last']=np.zeros(dTL['TRW'].size)
dTL['H obs EP']=np.zeros(dTL['TRW'].size)
dTL['Dob obs EP']=np.zeros(dTL['TRW'].size)
dTL['Bsw']=np.zeros(dTL['TRW'].size)
dTL['Gsw']=np.zeros(dTL['TRW'].size)
dTL['RGR']=np.zeros(dTL['TRW'].size)
dTL['TRW RR']=np.zeros(dTL['TRW'].size)
dTL['BAI RR']=np.zeros(dTL['TRW'].size)
dTL['Bsw RR']=np.zeros(dTL['TRW'].size)
dTL['Gsw AD']=np.zeros(dTL['TRW'].size)
dTL['Gsw RR']=np.zeros(dTL['TRW'].size)
dTL['RGR RR']=np.zeros(dTL['TRW'].size)
dTL['Lat']=np.zeros(dTL['TRW'].size)
dTL['Lon']=np.zeros(dTL['TRW'].size)

uTID=np.unique(dTL['ID_Tree_Unique'])
for iU in range(uTID.size):
    
    # Index to unique tree ID
    indUT=np.where(dTL['ID_Tree_Unique']==uTID[iU])[0]
    
    # Define first and last year of measurement for each tree
    iWithData=np.where(dTL['Dib'][indUT]>-1.0)[0]
    dTL['Year First'][indUT]=dTL['Year'][indUT[iWithData[0]]]
    dTL['Year Last'][indUT]=dTL['Year'][indUT[iWithData[-1]]]
    dTL['Age'][indUT[iWithData]]=np.arange(1,iWithData.size+1,1)
    
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
    ind0=np.where( (dTL['ID_Tree_Unique']==uTID[iU]) & (dTL['Year']>=1971-5) & (dTL['Year']<=1970) )[0]
    dTL['TRW RR'][indUT]=dTL['TRW'][indUT]/np.mean(dTL['TRW'][ind0])
    dTL['BAI RR'][indUT]=dTL['BAI'][indUT]/np.mean(dTL['BAI'][ind0])
    dTL['Bsw RR'][indUT]=dTL['Bsw'][indUT]/np.mean(dTL['Bsw'][ind0])
    dTL['Gsw AD'][indUT]=dTL['Gsw'][indUT]-np.mean(dTL['Gsw'][ind0])
    dTL['Gsw RR'][indUT]=dTL['Gsw'][indUT]/np.maximum(0.1,np.mean(dTL['Gsw'][ind0]))
    dTL['RGR RR'][indUT]=dTL['RGR'][indUT]/np.mean(dTL['RGR'][ind0])
    
# Change name of treatment
dTL['Treatment']=dTL['Treatment.EP703']

# Delete fields
del dTL['Unnamed: 0'],dTL['Treatment.EP703']

# Summarize QA flags
dTL['QA Summary']=np.zeros(dTL['TRW'].size)

# Cross dating
ind=np.where( (dTL['QC_xdate_flags']=='Good_correlation') | (dTL['QC_xdate_flags']=='Weak_correlation') )[0]
dTL['QA Summary'][ind]=1

# Unrealistic size removed
ind=np.where( (dTL['Dib Last']/dTL['Dob 2020']>1.0) )[0]
dTL['QA Summary'][ind]=0

print('% retained = ' + str(100*np.where(dTL['QA Summary']==1)[0].size/dTL['QA Summary'].size))

# Define unique installations
uInst=np.unique(dTL['ID_Inst'])

# Define unique site series
uSS=np.unique(dTL['BGC SS'])

# Combine SS 5 and 6 (may not be a good idea for all analysis!)
dTL['BGC SS Comb']=dTL['BGC SS'].copy()
ind=np.where( (dTL['BGC SS']==5) | (dTL['BGC SS']==6) )[0]
dTL['BGC SS Comb'][ind]=99


#%% Add detrended ring width (from J Axelson in dplR)

# Path to data
pthin=r'C:\Users\rhember\OneDrive - Government of BC\Shared\EP703\Retrospective Monitoring Study\Data\Given\Core Data\Received 20211222 from JAxelson'

# Initialize
dTL['TRW S NE DPLR']=np.zeros(dTL['Year'].size)

# Loop through installations
for iI in range(uInst.size):
    
    if uInst[iI]<10:
        nI='0' + str(uInst[iI])
    else:
        nI=str(uInst[iI])
    
    # Determine the letter associated with this installation
    letL=['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p']    
    for let in letL:
        try:
            d0=gu.ReadExcel(pthin + '\\703_' + nI + '\\Series\\703_' + nI + let + '_rwi_NE.xlsx')
            break
        except:
            pass
    
    # Pull year out
    Year=d0['Unnamed: 0']; del d0['Unnamed: 0']
    
    # Loop through trees
    for k in d0.keys():

        # Get tree and plot IDs
        ind=0
        for i in range(len(k)):
            if k[i]=='_':
                ind=i
                break
        ID_Plot=int(k[0:ind])
        ID_Tree=int(k[ind+1:])
        
        # Index to tree-level structure
        ind0=np.where( (dTL['ID_Inst']==uInst[iI]) & (dTL['ID_Plot']==ID_Plot) & (dTL['ID_Tree']==ID_Tree) & (dTL['Year']>=Year[0]) & (dTL['Year']<=Year[-1]) )[0]
        if ind0.size==0:
            print('Tree not found in DB')
            continue
        #ind0=np.where( (dTL['ID_Inst']==uInst[iI]) & (dTL['ID_Plot']==ID_Plot) )[0]
        #np.unique(dTL['ID_Tree'][ind0])
        #ind0.size
        
        # Populate
        dTL['TRW S NE DPLR'][ind0]=d0[k]

# Standardize the detrended time series
dTL['TRW S NE DPLR RR']=np.zeros(dTL['TRW S NE DPLR'].size)
for iU in range(uTID.size):
    
    # Index to unique tree ID
    indUT=np.where(dTL['ID_Tree_Unique']==uTID[iU])[0]
    
    # Standardized growth relative to mean growth during a period leading up to N application
    ind0=np.where( (dTL['ID_Tree_Unique']==uTID[iU]) & (dTL['Year']>=1971-5) & (dTL['Year']<=1970) )[0]
    dTL['TRW S NE DPLR RR'][indUT]=dTL['TRW S NE DPLR'][indUT]/np.mean(dTL['TRW S NE DPLR'][ind0])

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

#%% Get tree-level EP703 data
 
dTL['H obs EP']=np.zeros(dTL['TRW'].size)
dTL['Dob obs EP']=np.zeros(dTL['TRW'].size)
for iInst in range(uInst.size):    
    ind0=np.where( (dTL['ID_Inst']==uInst[iInst]) & (dTL['Treatment']=='C') | (dTL['ID_Inst']==uInst[iInst]) & (dTL['Treatment']=='F1') )[0]    
    uPlot=np.unique(dTL['ID_Plot'][ind0])    
    for iPlot in range(uPlot.size):        
        ind1=np.where( (dTL['ID_Inst']==uInst[iInst]) & (dTL['ID_Plot']==uPlot[iPlot]) & (dTL['QA Summary']==1) )[0]        
        uTree=np.unique(dTL['ID_Tree'][ind1])        
        for iTree in range(uTree.size):            
            ind2=np.where( (dEPT['ID_Site']==uInst[iInst]) & (dEPT['ID_Plot']==uPlot[iPlot]) & (dEPT['ID_Tree']==uTree[iTree]) & (dEPT['H_obs_t0']>0) )[0]            
            for iY in range(ind2.size):            
                ind3=np.where( (dTL['ID_Inst']==uInst[iInst]) & (dTL['ID_Plot']==uPlot[iPlot]) & (dTL['ID_Tree']==uTree[iTree]) & (dTL['Year']==dEPT['Year_t0'][ind2[iY]]) )[0]                
                dTL['H obs EP'][ind3]=dEPT['H_obs_t0'][ind2[iY]]
                dTL['Dob obs EP'][ind3]=dEPT['D_t0'][ind2[iY]]

##--------------------------------------------------------------------------
## Gap fill heights
## *** This is a mess ***
##--------------------------------------------------------------------------
#
## Add 2020 measurements to "obs EP" variables
#ind=np.where( (dTL['H 2020']>0) & (dTL['Dob 2020']>0) & (dTL['Year']==2020) )[0]
#dTL['H obs EP'][ind]=dTL['H 2020'][ind]
#dTL['Dob obs EP'][ind]=dTL['Dob 2020'][ind]
#
#def funcH2(x,a,b,c,d):
#    return a*( 1-np.exp(-b*(x-c)) )**d
#
#plt.close('all')
##plt.plot(dTL['Dob 2020'],dTL['H 2020'],'gs')
#plt.plot(dTL['Dob obs EP'],dTL['H obs EP'],'k.')
#
#iFit1=np.where( (np.isnan(dTL['Dob obs EP']+dTL['H obs EP'])==False) & (dTL['H obs EP']>0) & (dTL['Dob obs EP']>0) )[0]
#x=dTL['Dob obs EP'][iFit1]
#y=dTL['H obs EP'][iFit1]
#poptG,pcovG=curve_fit(funcH2,x,y,[55,0.045,10,1.2])
#yhat=funcH2(xhat,poptG[0],poptG[1],poptG[2],poptG[3])
#yhat=funcH2(xhat,55,0.045,10,1.2)
#plt.plot(xhat,yhat,'r-')        
#
#
#xhat=np.arange(0,60,1)
#
#dTL['H gf']=np.zeros(dTL['TRW'].size)
#for iInst in range(uInst.size):    
#    ind0=np.where( (dTL['ID_Inst']==uInst[iInst]) & (dTL['Treatment']=='C') | (dTL['ID_Inst']==uInst[iInst]) & (dTL['Treatment']=='F1') )[0]    
#    uPlot=np.unique(dTL['ID_Plot'][ind0])    
#    for iPlot in range(uPlot.size):        
#        ind1=np.where( (dTL['ID_Inst']==uInst[iInst]) & (dTL['ID_Plot']==uPlot[iPlot]) & (dTL['QA Summary']==1) )[0]        
#        uTree=np.unique(dTL['ID_Tree'][ind1])        
#        
#        # Global fit
#        iFit1=np.where((np.isnan(dTL['Dob obs EP'][ind1]+dTL['H obs EP'][ind1])==False))[0]
#        x=dTL['Dob obs EP'][ind1[iFit1]]
#        y=dTL['H obs EP'][ind1[iFit1]]
#        poptG,pcovG=curve_fit(funcH,x,y,[5,0.25])
#        yhat=funcH(xhat,poptG[0],poptG[1])
#        
#        plt.close('all')
#        plt.plot(x,y,'ko')
#        plt.plot(xhat,yhat,'r-')
#        
#        for iTree in range(uTree.size):            
#            ind2=np.where( (dEPT['ID_Site']==uInst[iInst]) & (dEPT['ID_Plot']==uPlot[iPlot]) & (dEPT['ID_Tree']==uTree[iTree]) & (dEPT['H_obs_t0']>0) )[0]            
#            for iY in range(ind2.size):            
#                ind3=np.where( (dTL['ID_Inst']==uInst[iInst]) & (dTL['ID_Plot']==uPlot[iPlot]) & (dTL['ID_Tree']==uTree[iTree]) & (dTL['Year']==dEPT['Year_t0'][ind2[iY]]) )[0]                
#                dTL['H gf'][ind3]=dEPT['H_obs_t0'][ind2[iY]]
#                dTL['Dob obs EP'][ind3]=dEPT['D_t0'][ind2[iY]]
#
#uPlotS=np.unique(dInst.plot)
#for iPlotS in range(len(uPlotS)):
#    
#    # Index to plot
#    indP=np.where((dInst.plot==uPlotS[iPlotS]))[0]        
#    idP_all=dInst.tree[indP]
#    hP_all=dInst.ht[indP]
#    baP_all=dInst.ba[indP]
#    
#    # Global fit
#    indP_Fit=np.where((np.isnan(hP_all+baP_all)==False))[0]
#    x=baP_all[indP_Fit].astype(float); 
#    y=hP_all[indP_Fit].astype(float)        
#    poptG,pcovG=curve_fit(funcH,x,y,[5,0.25])
#    
#    flg=0
#    if flg==1:
#        yhat=funcH(x,poptG[0],poptG[1])
#        N_Fit=indP_Fit.size
#        b_Perf=np.polyfit(y,yhat,1)
#        r_Perf=np.corrcoef(y,yhat)[0,1]        
#        plt.close('all')
#        fig,ax=plt.subplots(1,2)
#        ax[0].plot(x,y,'.')
#        ax[0].plot(x,yhat,'.')
#        ax[1].plot(y,yhat,'.')        
#    
#        #x2=np.c_[np.ones(x.size),np.log(x)]
#        #md=sm.OLS(np.log(y),x2).fit()
#        #md.summary()        
#    
#    indGF=np.where( (dInst.plot==uPlotS[iPlotS]) & (np.isnan(dInst.ba)==False) )[0]
#    dInst.ht_mod[indGF]=funcH(dInst.ba[indGF],poptG[0],poptG[1])
#    
#    # Individual tree fit
#    uTree=np.unique(idP_all[indP_Fit])
#    for k in range(uTree.size):
#        indCal=np.where( (idP_all==uTree[k]) & (np.isnan(hP_all+baP_all)==False) )[0]
#        if indCal.size>2:
#            x=baP_all[indCal].astype(float)
#            y=hP_all[indCal].astype(float)                      
#            try:
#                poptI,pcovI=curve_fit(funcH,x,y,[5,0.25])
#                yhat=funcH(x,poptI[0],poptI[1])
#                ax[0].plot(x,yhat,'-')
#                indGF=np.where( (dInst.plot==uPlotS[j]) & (dInst.tree==uTree[k]) & (np.isnan(dInst.ba)==False) )[0]
#                dInst.ht_mod[indGF]=funcH(dInst.ba[indGF],poptI[0],poptI[1])
#                dInst.H_Fit_ITM[indGF]=1
#            except:
#                continue        
#    
## Create gap-filled variable        
#ind=np.where( (np.isnan(dInst.ht)==True) & (np.isnan(dInst.ht_mod)==False) )[0]
#dInst.ht_gf[ind]=dInst.ht_mod[ind]

#%% Comparison between Dob and Dib from cores

def CompareDiameters(dTL,meta):
    
    # Isolate data
    ikp=np.where( (dTL['Dob 2020']>0) & (dTL['Dib Last']>0) & (dTL['QA Summary']==1) )[0]
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
    gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\EP703\TreeRings\DiameterComparison','png',150)

    return

#%% Look at individual trees

def LookAtIndividualTrees():

    ind=np.where( (dTL['ID_Tree_Unique']==25) )[0]

    print(dTL['Est/mm'][ind[0]])
    
    plt.close('all')
    fig,ax=plt.subplots(1,figsize=gu.cm2inch(8.5,8.5))
    #ax.plot([0,90],[0,90],'k-',lw=2,color=[0.8,0.8,0.8])
    ax.plot(dTL['Age'][ind],dTL['TRW'][ind],'-ko',ms=3,mec='k',mfc='k')
    #ax.plot(dTL['Year'][ind],dTL['TRW'][ind],'-ko',ms=3,mec='k',mfc='k')
    
    
    fig,ax=plt.subplots(1,figsize=gu.cm2inch(8.5,8.5))
    #ax.plot([0,90],[0,90],'k-',lw=2,color=[0.8,0.8,0.8])
    ax.plot(dTL['Age'][ind],dTL['BAI'][ind],'-ko',ms=3,mec='k',mfc='k')
    #ax.plot(dTL['Year'][ind],dTL['TRW'][ind],'-ko',ms=3,mec='k',mfc='k')
    
    #plt.close('all')
    #fig,ax=plt.subplots(1,figsize=gu.cm2inch(8.5,8.5))
    #ax.plot([0,90],[0,90],'k-',lw=2,color=[0.8,0.8,0.8])
    #ax.plot(dTL['TRW'][ind],dTL['Gsw'][ind],'-ko',ms=3,mec='k',mfc='k')
    
    return

#%% Gap-fill early missing rings using Weibul function

# Weibul distribution function
def func(x,a,b,c,d):
    y=a*(1+((b*(x/c)**d-1)/(np.exp(x/c))))
    return y

def GapFill():
    
    uTID=np.unique(dTL['ID_Tree_Unique'])
    
    # Initialize new structure
    dGF={}
    ExclusionList=['TRW S NE DPLR', 'TRW S NE DPLR RR', 'tmean_ann_n', 'prcp_ann_n', 'tmean_gs_r', 'prcp_gs_r', 'cwd_gs_r', 'ws_gs_r']
    for k in dTL.keys():
        if np.isin(k,ExclusionList)==True:
            continue
        if (dTL[k].dtype=='float') | (dTL[k].dtype=='int32'):
            dGF[k]=np.zeros(200000)
        else:
            dGF[k]=np.array(['' for _ in range(200000)],dtype=dTL[k].dtype)
    
    cnt=0
    
    for iU in range(uTID.size):
        
        # Index to unique tree ID
        iTree=np.where( (dTL['ID_Tree_Unique']==uTID[iU]) & (dTL['QA Summary']==1) )[0]
        
        TRW=dTL['TRW'][iTree]
        Year=dTL['Year'][iTree]
        Dob20=dTL['Dob 2020'][iTree]
        
        iBad=np.where( (np.isnan(TRW)==True) )[0]
        
        # Initialize full age series
        Age=np.arange(0,iTree.size,1)
        
        Dcore=0.1*np.cumsum(2*np.nan_to_num(TRW))
        
        
        D=Dob20[-1]-Dcore[-1]
        Dcore=Dcore+D
        Dcore[iBad]=0
        
        rw=np.append(0,np.diff(Dcore))
        rw[0:iBad[-1]+2]=np.nan
        
        plt.close('all')
        plt.plot(Age,rw,'.b-')
        
        rw_gf=rw.copy()
        iFill=np.where(rw==0)[0]
        p0=[0.5,2.5,5,3]
        rw_gf=func(Age,p0[0],p0[1],p0[2],p0[3])  
        #plt.plot(Age,rw_gf,'g--')
        
        iFit=np.where( (np.isnan(rw)==False) )[0]
        x=np.append(np.arange(0,5,1),Age[iFit])
        y=np.append(np.zeros(5),rw[iFit])
        p,pcov=curve_fit(func,x,y,p0)
        rw_gf=func(Age,p[0],p[1],p[2],p[3])  
        plt.plot(Age,rw_gf,'g--',lw=1.5)
        
        
        
      
        
        
        
        
        
        
        # Prediction
        x0=np.append(np.zeros(1000),dTL['Age'][iGood])
        y0=np.append(np.zeros(1000),dTL['TRW'][iGood])    
        binFit=np.arange(0.05,1.6,0.2)
        dD=np.zeros(binFit.size)
        par=[None]*binFit.size
        for iFit in range(binFit.size):        
            x=np.append(1*np.ones(1000),x0); 
            y=np.append(binFit[iFit]*np.ones(1000),y0)
            try:
                p,pcov=curve_fit(func,x,y,[1,1,25,0.5])
            except:
                continue        
            TRW_full=-999*np.ones(A_full.size)
            TRW_full[ia]=dTL['TRW'][iGood]
            iFill=np.where(TRW_full==-999)[0]
            TRW_full[iFill]=func(A_full[iFill],p[0],p[1],p[2],p[3])        
            Dib_full=0.1*np.cumsum(2*TRW_full)
            dD[iFit]=Dib_full[-1]/dTL['Dob 2020'][iGood[-1]]
            par[iFit]=p
        
        # Gap-fill with the best fit to Dob at time of measurement
        iMinD=np.where( (np.abs(1-dD)==np.min(np.abs(1-dD))) )[0]
        
        #iMinD=np.array([0])
        
        if (iMinD.size>0) & (np.sum(dD)>0) & (dTL['Dob 2020'][iGood[-1]]>0):
            TRW_full=-999*np.ones(A_full.size)
            TRW_full[ia]=dTL['TRW'][iGood]
            iFill=np.where(TRW_full==-999)[0]
            p=par[iMinD[0]]
            TRW_full[iFill]=func(A_full[iFill],p[0],p[1],p[2],p[3])  
            Dib_full=0.1*np.cumsum(2*TRW_full)
        else:
            TRW_full=np.zeros(A_full.size)
            TRW_full[ia]=dTL['TRW'][iGood]
            iFill=np.where(TRW_full==0)[0]
            Dib_full=0.1*np.cumsum(2*TRW_full)
        
        flg=0
        if flg==1:
            p=par[0]
            plt.close('all'); 
            plt.plot(dTL['Age'][iGood],dTL['TRW'][iGood],'k-',color=np.random.random(3),lw=1)
            plt.plot(A_full,func(A_full,p[0],p[1],p[2],p[3]),'k--')
    
        #TRW_full[iFill]=func(A_full[iFill],p[0],p[1],p[2],p[3])
        #plt.plot(A_full,TRW_full,'ko')
        flg=0
        if flg==1:
            Dib_given=0.1*np.cumsum(2*dTL['TRW'][iGood])
            plt.close('all'); 
            plt.plot(dTL['Age'][iGood],Dib_given,'k-',color='b',lw=1)
            plt.plot(A_full,Dib_full,'k-',color='g',lw=1)
            plt.plot(dTL['Age'][iAll[-1]],dTL['Dob 2020'][iGood[-1]],'ks')
            #gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\EP703\TreeRings\GapFill_' + str(uTID[iU]),'png',200)
        
        # Index to final data dictionary
        ind=np.arange(cnt,cnt+Year_full.size,1)
        
        # Populate scalar variables
        for k in dGF.keys():
            dGF[k][ind]=dTL[k][iFirst]
            
        #Filler=np.ones(Year_full.size)
        
        dGF['Year'][ind]=Year_full
        dGF['Age'][ind]=A_full
        dGF['TRW'][ind]=TRW_full
        dGF['Dib'][ind]=Dib_full
        dGF['Dib Last'][ind]=np.max(Dib_full)
        
        BA=np.pi*(dGF['Dib'][ind]/2)**2
        dGF['BAI'][ind]=np.append(0,np.diff(BA))
        
        # Stemwood biomass from allometric function (kgDM) - DBH only (Ung et al. 2008)
        dGF['Bsw'][ind]=0.0204*dGF['Dib'][ind]**2.6974
        
        # Stemwood biomass growth (kgDM/yr)
        dGF['Gsw'][ind]=np.append(0,np.diff(dGF['Bsw'][ind]))
        
        # Relative growth rate
        dGF['RGR'][ind]=np.append(0,np.log(dGF['Bsw'][ind[1:]])-np.log(dGF['Bsw'][ind[0:-1]]))
    
        # Standardized growth relative to mean growth during a period leading up to N application
        ind0=np.where( (dGF['ID_Tree_Unique']==dGF['ID_Tree_Unique'][ind[0]]) & (dGF['Year']>=1971-5) & (dGF['Year']<=1970) )[0]
        dGF['TRW RR'][ind]=dGF['TRW'][ind]/np.mean(dGF['TRW'][ind0])
        dGF['BAI RR'][ind]=dGF['BAI'][ind]/np.mean(dGF['BAI'][ind0])
        dGF['Bsw RR'][ind]=dGF['Bsw'][ind]/np.mean(dGF['Bsw'][ind0])
        dGF['Gsw AD'][ind]=dGF['Gsw'][ind]-np.mean(dGF['Gsw'][ind0])
        dGF['Gsw RR'][ind]=dGF['Gsw'][ind]/np.maximum(0.1,np.mean(dGF['Gsw'][ind0]))
        dGF['RGR RR'][ind]=dGF['RGR'][ind]/np.mean(dGF['RGR'][ind0])
        
        # Update counter
        cnt=cnt+Year_full.size
    
    # Remove excess data
    for k in dGF.keys():
        dGF[k]=dGF[k][0:cnt]
    


#%% Comparison between Dob and Dib from cores (after gap-filling)

def CompareDiametersAfterGF():
    # Isolate data
    ikp=np.where( (dGF['Dob 2020']>0) & (dGF['Dib Last']>0) & (dGF['QA Summary']==1) )[0]
    #ikp=np.where( (dGF['Dob 2020']>0) & (dGF['Dib Last']>0) & (dGF['Pith_status']=='Yes') )[0]
    
    # Fit linear best fit relationship
    y=dGF['Dib Last'][ikp]
    x=dGF['Dob 2020'][ikp]
    x1=sm.tools.tools.add_constant(x)
    md=sm.OLS(y,x1).fit()
    md.summary()
    xhat=np.linspace(np.min(x1[:,1]),np.max(x1[:,1]),10)
    yhat=md.predict(np.c_[np.ones(xhat.size),xhat])
    
    # Plot
    #plt.close('all')
    fig,ax=plt.subplots(1,figsize=gu.cm2inch(8.5,8.5))
    ax.plot([0,90],[0,90],'k-',lw=2,color=[0.8,0.8,0.8])
    ax.plot(x,y,'o',ms=4,mec='w',mfc='k')
    ax.plot(xhat,yhat,'r-',lw=1.25)
    ax.set(xlim=[0,80],ylim=[0,80],xlabel='Tape-based Dob in 2020 (cm)',ylabel='Core-based Dib in 2019 (cm)')
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')
    plt.tight_layout()
    #gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\EP703\TreeRings\DiameterComparison_GF','png',150)
