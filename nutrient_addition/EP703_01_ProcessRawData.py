''' 

PROCESS EP703 DATA

Notes:
    1) Thinned plots are excluded
    2) Missing tree heights are gap-filled
    3) Climate data are added at the end

'''

#%% Import modules

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyproj
from scipy.optimize import curve_fit
import scipy.io as spio
import fcgadgets.macgyver.utilities_general as gu
import warnings

warnings.filterwarnings('ignore')

#%% Configure metadata

meta={}
meta['Paths']={}
meta['Paths']['Project']=r'I:\My Drive\Data\EP703'
meta['Paths']['Inputs']=meta['Paths']['Project'] + '\\Received20190515'
meta['Paths']['Outputs']=meta['Paths']['Project'] + '\\Processed'
meta['Paths']['NACID']=r'E:\Data\Climate\NACID'

#%% Define functions

# Function (tree height as a function of basal area)
def funcH(x,a,b):
    return a*x**b

#%% Import site data
# This spreadsheet has extracted info for the unthinned plots.

pthin=meta['Paths']['Project'] + '\\EP703_SiteInfo.xlsx'
dSite=pd.read_excel(pthin,sheet_name='Sheet1')

# Convert to data structure
dSite=gu.DataFrameToDataStruct(dSite)

# Unique sites
uSite=np.unique(dSite.ID_Instal)

#%% Define parameters

# Allometry - DBH only (Ung et al. 2008)
pA=list()
d2b={'Species':'FD','Stemwood':[0.0204,2.6974],'Foliage':[0.1233,1.6636],'Branch':[0.0404,2.1388],'Bark':[0.0069,2.5462]}
pA.append(d2b.copy())
d2b={'Species':'HW','Stemwood':[0.0141,2.8668],'Foliage':[0.1676,1.4339],'Branch':[0.0703,1.9547],'Bark':[0.0025,2.8062]}
pA.append(d2b.copy())
d2b={'Species':'CW','Stemwood':[0.0111,2.8027],'Foliage':[0.1233,1.5152],'Branch':[0.1158,1.7196],'Bark':[0.0003,3.2721]}
pA.append(d2b.copy())
d2b={'Species':'BA','Stemwood':[0.0424,2.4289],'Foliage':[0.0645,1.9400],'Branch':[0.0322,2.1313],'Bark':[0.0057,2.4786]}
pA.append(d2b.copy())

# Allometry - DBH plus height (Ung et al. 2008)
pA2=list()
d2b={'Species':'FD','Stemwood':[0.0191,1.5365,1.3634],'Foliage':[0.0718,2.2935,-0.4744],'Branch':[0.0351,2.24811,0],'Bark':[0.0083,2.4811,0]}
pA2.append(d2b.copy())
d2b={'Species':'HW','Stemwood':[0.0113,1.9332,1.1125],'Foliage':[0.2656,2.0107,-0.7963],'Branch':[0.0609,2.0021,0],'Bark':[0.0019,2.3356,0.6371]}
pA2.append(d2b.copy())
d2b={'Species':'CW','Stemwood':[0.0188,1.3376,1.5293],'Foliage':[0.1097,1.553,0],'Branch':[0.0611,1.9208,0],'Bark':[0.0002,2.4369,1.1315]}
pA2.append(d2b.copy())
d2b={'Species':'BA','Stemwood':[0.0315,1.8297,0.8056],'Foliage':[0.0453,2.4867,-0.4982],'Branch':[0.0420,2.0313,0],'Bark':[0.0067,2.6970,-0.3105]}
pA2.append(d2b.copy())

# Convert kg to Mg
kg2Mg=0.001

# Convert dry mass to carbon
dm2c=0.5

#%% List of damage agents

da_list=['A', 'AB', 'AD', 'AE', 'AX', 'BT', 'CK', 'D', 'DDF', 'DDP', 'DDS',
       'DM', 'DN', 'DR', 'DRA', 'DS', 'DT', 'FCK', 'FK', 'MT', 'N', 'NGC',
       'NW', 'NWS', 'NX', 'NY', 'PC', 'SC', 'SCB', 'SCM', 'SCT']

#%% Initialize stand and tree level data

# Initialize stand-level dictionary
varS=['ID_Site','ID_Plot','Lat','Lon','BEC','SI','EstabType','AreaPlot', \
    'N_Init','pH_Min_Init','pH_Humus_Init','P_Min_Init','P_Humus_Init','Ca_Min_Init', \
    'N_Dose','N_Type','N_AppMeth','NumTreesMeas','Num_H_G_NegOrMis', \
    'Spc1_ID_t0','Spc1_Pct_t0','Spc2_ID_t0','Spc2_Pct_t0', \
    'Year_t0','Year_t1','DT','Age_t0','Age_t1','TSF_t0','TSF_t1', \
    'N_t0','N_t1','N_Net','Dam_t0','Dam_t1', \
    'BA_t0','BA_t1','BA_G','H_obs_t0','H_obs_t1','H_obs_G', \
    'H_gf_t0','H_gf_t1','H_gf_G', \
    'Bsw_t0','Bsw_t1','Bsw_G','Bsw_M','Bsw_R','Bsw_Net', \
    'Bf_t0','Bf_t1','Bf_G','Bf_M','Bf_R','Bf_Net', \
    'Bbr_t0','Bbr_t1','Bbr_G','Bbr_M','Bbr_R','Bbr_Net', \
    'Bbk_t0','Bbk_t1','Bbk_G','Bbk_M','Bbk_R','Bbk_Net', \
    'Bsw_G_cc1','Bsw_G_cc2','Bsw_G_cc3','Bsw_G_cc4','Bsw_G_cc5','Bsw_G_cc6', \
    'Bsw_M_cc1','Bsw_M_cc2','Bsw_M_cc3','Bsw_M_cc4','Bsw_M_cc5','Bsw_M_cc6', \
    'Bsw_Net_cc1','Bsw_Net_cc2','Bsw_Net_cc3','Bsw_Net_cc4','Bsw_Net_cc5','Bsw_Net_cc6', \
    'DA_A_t0','DA_AB_t0','DA_AD_t0','DA_AE_t0','DA_AX_t0','DA_BT_t0', \
    'DA_CK_t0','DA_D_t0','DA_DDF_t0','DA_DDP_t0','DA_DDS_t0', \
    'DA_DM_t0','DA_DN_t0','DA_DR_t0','DA_DRA_t0','DA_DS_t0','DA_DT_t0', \
    'DA_FCK_t0','DA_FK_t0','DA_MT_t0','DA_N_t0','DA_NGC_t0','DA_NW_t0','DA_NWS_t0', \
    'DA_NX_t0','DA_NY_t0','DA_PC_t0','DA_SC_t0','DA_SCB_t0','DA_SCM_t0','DA_SCT_t0', \
    'DA_A_t1','DA_AB_t1','DA_AD_t1','DA_AE_t1','DA_AX_t1','DA_BT_t1', \
    'DA_CK_t1','DA_D_t1','DA_DDF_t1','DA_DDP_t1','DA_DDS_t1', \
    'DA_DM_t1','DA_DN_t1','DA_DR_t1','DA_DRA_t1','DA_DS_t1','DA_DT_t1', \
    'DA_FCK_t1','DA_FK_t1','DA_MT_t1','DA_N_t1','DA_NGC_t1','DA_NW_t1','DA_NWS_t1', \
    'DA_NX_t1','DA_NY_t1','DA_PC_t1','DA_SC_t1','DA_SCB_t1','DA_SCM_t1','DA_SCT_t1']

sobs={}
for k in varS:
    sobs[k]=np.zeros(3000)

# Initialize a counter for stand-level dictionary
cntS=0

sobs['BEC']=sobs['BEC'].astype('U')
sobs['EstabType']=sobs['EstabType'].astype('U')
sobs['N_Type']=sobs['N_Type'].astype('U')
sobs['N_AppMeth']=sobs['N_AppMeth'].astype('U')
sobs['Spc1_ID_t0']=sobs['Spc1_ID_t0'].astype('U')
sobs['Spc2_ID_t0']=sobs['Spc2_ID_t0'].astype('U')

# Initialize tree-level dictionary
varT=['ID_Site','ID_Plot','ID_Tree','Lat','Lon','BEC','SI','AEF','N_Init','pH_Min_Init', \
    'N_Dose','Year_t0','Year_t1','DT', \
    'TSF_t0','TSF_t1','SA_t0','SA_t1','SN_t0','SN_t1', \
    'SBsw_t0','SBsw_t1','SBswLT_t0','SBswLT_t1', \
    'ID_Spc_t0','Mort','Rec','TreeClass_t0','CrownClass_t0','TopDead_t0','TopDead_t1', \
    'TopBroke_t0','TopBroke_t1','DA_t0','DA_t1','BA_t0','BA_t1','BA_G', \
    'H_obs_t0','H_obs_t1','H_obs_G','H_mod_t0','H_mod_t1','H_mod_G','H_gf_t0','H_gf_t1','H_gf_G', \
    'Bsw_t0','Bsw_t1','Bsw_G','Bsw_mod_t0','Bsw_mod_t1','H_Fit_ITM_t0','H_Fit_ITM_t1']
tobs={}
for k in varT:
    tobs[k]=np.zeros(250000)

# Initialize a counter for tree-level dictionary
cntT=0

tobs['BEC']=tobs['BEC'].astype('U')
tobs['DA_t0']=tobs['DA_t0'].astype('U')
tobs['DA_t1']=tobs['DA_t1'].astype('U')

#%% Loop through installations (40 sec)

N=0
for i in range(uSite.size):
    
    # Exclude sites that are missing
    if (uSite[i]==13) | (uSite[i]==62) | (uSite[i]==82): 
        continue
    
    #--------------------------------------------------------------------------
    # Import tree-level data for site
    #--------------------------------------------------------------------------
    
    if uSite[i]<10: 
        fnam='0' + str(uSite[i])
    else: 
        fnam=str(uSite[i])
    
    pthin=meta['Paths']['Project']['Inputs'] + '\\EP703_inst' + fnam + '.csv'
        
    dfInst=pd.read_csv(pthin)
    
    # Convert to data structure
    dInst=gu.DataFrameToDataStruct(dfInst)
    
    #N=N+dInst.tree'].unique().size
    
    # Calculate basal area (cm2)
    dInst.ba=np.pi*(dInst.dbh/2)**2
    
    # Vital status
    dInst.VS=np.ones(dInst.tree.size)
    dInst.VS[dInst.mort!='nan']=0
    
    # Initialize modelled and gap-filled tree height
    dInst.ht_mod=np.nan*np.ones(dInst.tree.size)
    dInst.ht_gf=dInst.ht.copy()
    
    # Indicator
    dInst.H_Fit_ITM=np.zeros(dInst.tree.size)
    
    #--------------------------------------------------------------------------
    # Site information
    #--------------------------------------------------------------------------
    
    # Index to site table
    iSite=np.where(dSite.ID_Instal==uSite[i])[0]
    
    # Unique plots within this installation
    uPlot=dSite.ID_Plot[iSite]
    
    # Unique years
    uYear=np.sort(np.unique(dInst.year))
    
    # Area expansion factor
    Area=dSite.Plot_area[iSite[0]]
    AEF=1/Area
    
    # Age 
    SA=dSite.StandAge[iSite[0]]
    
    # Year of establishment
    YrEst=dSite.Year_Est[iSite[0]]
    
    # Geographic coordinates
    Lat=dSite.Lat[iSite[0]]+dSite.Lat_1[iSite[0]]/60+dSite.Lat_2[iSite[0]]/(60*60)
    Lon=-dSite.Lon[iSite[0]]-dSite.Lon_1[iSite[0]]/60-dSite.Lon_2[iSite[0]]/(60*60)
    
    #--------------------------------------------------------------------------
    # Gap fill heights
    #--------------------------------------------------------------------------
    
    uPlotS=np.unique(dInst.plot)
    for iPlotS in range(len(uPlotS)):
        
        # Index to plot
        indP=np.where((dInst.plot==uPlotS[iPlotS]))[0]        
        idP_all=dInst.tree[indP]
        hP_all=dInst.ht[indP]
        baP_all=dInst.ba[indP]
        
        # Global fit
        indP_Fit=np.where((np.isnan(hP_all+baP_all)==False))[0]
        x=baP_all[indP_Fit].astype(float); 
        y=hP_all[indP_Fit].astype(float)        
        poptG,pcovG=curve_fit(funcH,x,y,[5,0.25])
        
        flg=0
        if flg==1:
            yhat=funcH(x,poptG[0],poptG[1])
            N_Fit=indP_Fit.size
            b_Perf=np.polyfit(y,yhat,1)
            r_Perf=np.corrcoef(y,yhat)[0,1]        
            plt.close('all')
            fig,ax=plt.subplots(1,2)
            ax[0].plot(x,y,'.')
            ax[0].plot(x,yhat,'.')
            ax[1].plot(y,yhat,'.')        
        
            #x2=np.c_[np.ones(x.size),np.log(x)]
            #md=sm.OLS(np.log(y),x2).fit()
            #md.summary()        
        
        indGF=np.where( (dInst.plot==uPlotS[iPlotS]) & (np.isnan(dInst.ba)==False) )[0]
        dInst.ht_mod[indGF]=funcH(dInst.ba[indGF],poptG[0],poptG[1])
        
        # Individual tree fit
        uTree=np.unique(idP_all[indP_Fit])
        for k in range(uTree.size):
            indCal=np.where( (idP_all==uTree[k]) & (np.isnan(hP_all+baP_all)==False) )[0]
            if indCal.size>2:
                x=baP_all[indCal].astype(float)
                y=hP_all[indCal].astype(float)                      
                try:
                    poptI,pcovI=curve_fit(funcH,x,y,[5,0.25])
                    yhat=funcH(x,poptI[0],poptI[1])
                    ax[0].plot(x,yhat,'-')
                    indGF=np.where( (dInst.plot==uPlotS[j]) & (dInst.tree==uTree[k]) & (np.isnan(dInst.ba)==False) )[0]
                    dInst.ht_mod[indGF]=funcH(dInst.ba[indGF],poptI[0],poptI[1])
                    dInst.H_Fit_ITM[indGF]=1
                except:
                    continue        
        
    # Create gap-filled variable        
    ind=np.where( (np.isnan(dInst.ht)==True) & (np.isnan(dInst.ht_mod)==False) )[0]
    dInst.ht_gf[ind]=dInst.ht_mod[ind]
    
    #--------------------------------------------------------------------------
    # Loop through plots and years
    #--------------------------------------------------------------------------
    
    for j in range(len(uPlot)):
        for k in range(len(uYear)-1):
            
            print(str(i) + ' ' + str(j) + ' ' + str(k))
            
            yr0=uYear[k]
            yr1=uYear[k+1]
            dt=yr1-yr0
            
            ind0=np.where((dInst.plot==uPlot[j]) & (dInst.year==uYear[k]))[0]
            
            if ind0.size==0:
                continue
            
            ind1=np.where((dInst.plot==uPlot[j]) & (dInst.year==uYear[k+1]))[0]
            
            #------------------------------------------------------------------    
            # Extract variables for all trees in plot j and year k
            #------------------------------------------------------------------
            
            # First measurement
            id0_all=dInst.tree[ind0]
            sp0_all=dInst.spp[ind0]
            vs0_all=dInst.VS[ind0]
            dbh0_all=dInst.dbh[ind0]
            h_obs0_all=dInst.ht[ind0]
            h_mod0_all=dInst.ht_mod[ind0]
            h_gf0_all=dInst.ht_gf[ind0]            
            ba0_all=np.pi*(dbh0_all/2)**2
            TreeClass0_all=dInst.TrClass[ind0]
            CrownClass0_all=dInst.CrClass[ind0]
            btop0_all=dInst.bTop[ind0]
            dtop0_all=dInst.dTop[ind0]
            da0_all=dInst.DmAg1[ind0]
            H_Fit_ITM0_all=dInst.H_Fit_ITM[ind0] 
            
            Bsw0_all=pA2[0]['Stemwood'][0]*dbh0_all**pA2[0]['Stemwood'][1]*h_gf0_all**pA2[0]['Stemwood'][2]
            Bsw0_all[sp0_all=='FD']=pA2[0]['Stemwood'][0]*dbh0_all[sp0_all=='FD']**pA2[0]['Stemwood'][1]*h_gf0_all[sp0_all=='FD']**pA2[0]['Stemwood'][2]
            Bsw0_all[sp0_all=='HW']=pA2[1]['Stemwood'][0]*dbh0_all[sp0_all=='HW']**pA2[1]['Stemwood'][1]*h_gf0_all[sp0_all=='HW']**pA2[1]['Stemwood'][2]
            Bsw0_all[sp0_all=='CW']=pA2[2]['Stemwood'][0]*dbh0_all[sp0_all=='CW']**pA2[2]['Stemwood'][1]*h_gf0_all[sp0_all=='CW']**pA2[2]['Stemwood'][2]
            Bsw0_all[sp0_all=='BA']=pA2[3]['Stemwood'][0]*dbh0_all[sp0_all=='BA']**pA2[3]['Stemwood'][1]*h_gf0_all[sp0_all=='BA']**pA2[3]['Stemwood'][2]
            
            Bsw_mod0_all=pA2[0]['Stemwood'][0]*dbh0_all**pA2[0]['Stemwood'][1]*h_mod0_all**pA2[0]['Stemwood'][2]
            Bsw_mod0_all[sp0_all=='FD']=pA2[0]['Stemwood'][0]*dbh0_all[sp0_all=='FD']**pA2[0]['Stemwood'][1]*h_mod0_all[sp0_all=='FD']**pA2[0]['Stemwood'][2]
            Bsw_mod0_all[sp0_all=='HW']=pA2[1]['Stemwood'][0]*dbh0_all[sp0_all=='HW']**pA2[1]['Stemwood'][1]*h_mod0_all[sp0_all=='HW']**pA2[1]['Stemwood'][2]
            Bsw_mod0_all[sp0_all=='CW']=pA2[2]['Stemwood'][0]*dbh0_all[sp0_all=='CW']**pA2[2]['Stemwood'][1]*h_mod0_all[sp0_all=='CW']**pA2[2]['Stemwood'][2]
            Bsw_mod0_all[sp0_all=='BA']=pA2[3]['Stemwood'][0]*dbh0_all[sp0_all=='BA']**pA2[3]['Stemwood'][1]*h_mod0_all[sp0_all=='BA']**pA2[3]['Stemwood'][2]
            
            Bf0_all=pA2[0]['Foliage'][0]*dbh0_all**pA2[0]['Foliage'][1]*h_gf0_all**pA2[0]['Foliage'][2]
            Bf0_all[sp0_all=='FD']=pA2[0]['Foliage'][0]*dbh0_all[sp0_all=='FD']**pA2[0]['Foliage'][1]*h_gf0_all[sp0_all=='FD']**pA2[0]['Foliage'][2]
            Bf0_all[sp0_all=='HW']=pA2[1]['Foliage'][0]*dbh0_all[sp0_all=='HW']**pA2[1]['Foliage'][1]*h_gf0_all[sp0_all=='HW']**pA2[1]['Foliage'][2]
            Bf0_all[sp0_all=='CW']=pA2[2]['Foliage'][0]*dbh0_all[sp0_all=='CW']**pA2[2]['Foliage'][1]*h_gf0_all[sp0_all=='CW']**pA2[2]['Foliage'][2]
            Bf0_all[sp0_all=='BA']=pA2[3]['Foliage'][0]*dbh0_all[sp0_all=='BA']**pA2[3]['Foliage'][1]*h_gf0_all[sp0_all=='BA']**pA2[3]['Foliage'][2]
            
            Bbr0_all=pA2[0]['Branch'][0]*dbh0_all**pA2[0]['Branch'][1]*h_gf0_all**pA2[0]['Branch'][2]
            Bbr0_all[sp0_all=='FD']=pA2[0]['Branch'][0]*dbh0_all[sp0_all=='FD']**pA2[0]['Branch'][1]*h_gf0_all[sp0_all=='FD']**pA2[0]['Branch'][2]
            Bbr0_all[sp0_all=='HW']=pA2[1]['Branch'][0]*dbh0_all[sp0_all=='HW']**pA2[1]['Branch'][1]*h_gf0_all[sp0_all=='HW']**pA2[1]['Branch'][2]
            Bbr0_all[sp0_all=='CW']=pA2[2]['Branch'][0]*dbh0_all[sp0_all=='CW']**pA2[2]['Branch'][1]*h_gf0_all[sp0_all=='CW']**pA2[2]['Branch'][2]
            Bbr0_all[sp0_all=='BA']=pA2[3]['Branch'][0]*dbh0_all[sp0_all=='BA']**pA2[3]['Branch'][1]*h_gf0_all[sp0_all=='BA']**pA2[3]['Branch'][2]
            
            Bbk0_all=pA2[0]['Bark'][0]*dbh0_all**pA2[0]['Bark'][1]*h_gf0_all**pA2[0]['Bark'][2]
            Bbk0_all[sp0_all=='FD']=pA2[0]['Bark'][0]*dbh0_all[sp0_all=='FD']**pA2[0]['Bark'][1]*h_gf0_all[sp0_all=='FD']**pA2[0]['Bark'][2]
            Bbk0_all[sp0_all=='HW']=pA2[1]['Bark'][0]*dbh0_all[sp0_all=='HW']**pA2[1]['Bark'][1]*h_gf0_all[sp0_all=='HW']**pA2[1]['Bark'][2]
            Bbk0_all[sp0_all=='CW']=pA2[2]['Bark'][0]*dbh0_all[sp0_all=='CW']**pA2[2]['Bark'][1]*h_gf0_all[sp0_all=='CW']**pA2[2]['Bark'][2]
            Bbk0_all[sp0_all=='BA']=pA2[3]['Bark'][0]*dbh0_all[sp0_all=='BA']**pA2[3]['Bark'][1]*h_gf0_all[sp0_all=='BA']**pA2[3]['Bark'][2]
            
            # Second measurement
            id1_all=dInst.tree[ind1]
            sp1_all=dInst.spp[ind1]
            vs1_all=dInst.VS[ind1]
            dbh1_all=dInst.dbh[ind1]
            h_obs1_all=dInst.ht[ind1]
            h_mod1_all=dInst.ht_mod[ind1]
            h_gf1_all=dInst.ht_gf[ind1]            
            H_Fit_ITM1_all=dInst.H_Fit_ITM[ind1] 
            ba1_all=np.pi*(dbh1_all/2)**2
            TreeClass1_all=dInst.TrClass[ind1]
            CrownClass1_all=dInst.CrClass[ind1]
            btop1_all=dInst.bTop[ind1]
            dtop1_all=dInst.dTop[ind1]
            da1_all=dInst.DmAg1[ind1]
            
            Bsw1_all=pA2[0]['Stemwood'][0]*dbh1_all**pA2[0]['Stemwood'][1]*h_gf1_all**pA2[0]['Stemwood'][2]
            Bsw1_all[sp1_all=='FD']=pA2[0]['Stemwood'][0]*dbh1_all[sp1_all=='FD']**pA2[0]['Stemwood'][1]*h_gf1_all[sp1_all=='FD']**pA2[0]['Stemwood'][2]
            Bsw1_all[sp1_all=='HW']=pA2[1]['Stemwood'][0]*dbh1_all[sp1_all=='HW']**pA2[1]['Stemwood'][1]*h_gf1_all[sp1_all=='HW']**pA2[1]['Stemwood'][2]
            Bsw1_all[sp1_all=='CW']=pA2[2]['Stemwood'][0]*dbh1_all[sp1_all=='CW']**pA2[2]['Stemwood'][1]*h_gf1_all[sp1_all=='CW']**pA2[2]['Stemwood'][2]
            Bsw1_all[sp1_all=='BA']=pA2[3]['Stemwood'][0]*dbh1_all[sp1_all=='BA']**pA2[3]['Stemwood'][1]*h_gf1_all[sp1_all=='BA']**pA2[3]['Stemwood'][2]
            
            Bsw_mod1_all=pA2[0]['Stemwood'][0]*dbh1_all**pA2[0]['Stemwood'][1]*h_mod1_all**pA2[0]['Stemwood'][2]
            Bsw_mod1_all[sp1_all=='FD']=pA2[0]['Stemwood'][0]*dbh1_all[sp1_all=='FD']**pA2[0]['Stemwood'][1]*h_mod1_all[sp1_all=='FD']**pA2[0]['Stemwood'][2]
            Bsw_mod1_all[sp1_all=='HW']=pA2[1]['Stemwood'][0]*dbh1_all[sp1_all=='HW']**pA2[1]['Stemwood'][1]*h_mod1_all[sp1_all=='HW']**pA2[1]['Stemwood'][2]
            Bsw_mod1_all[sp1_all=='CW']=pA2[2]['Stemwood'][0]*dbh1_all[sp1_all=='CW']**pA2[2]['Stemwood'][1]*h_mod1_all[sp1_all=='CW']**pA2[2]['Stemwood'][2]
            Bsw_mod1_all[sp1_all=='BA']=pA2[3]['Stemwood'][0]*dbh1_all[sp1_all=='BA']**pA2[3]['Stemwood'][1]*h_mod1_all[sp1_all=='BA']**pA2[3]['Stemwood'][2]
                        
            Bf1_all=pA2[0]['Foliage'][0]*dbh1_all**pA2[0]['Foliage'][1]*h_gf1_all**pA2[0]['Foliage'][2]
            Bf1_all[sp1_all=='FD']=pA2[0]['Foliage'][0]*dbh1_all[sp1_all=='FD']**pA2[0]['Foliage'][1]*h_gf1_all[sp1_all=='FD']**pA2[0]['Foliage'][2]
            Bf1_all[sp1_all=='HW']=pA2[1]['Foliage'][0]*dbh1_all[sp1_all=='HW']**pA2[1]['Foliage'][1]*h_gf1_all[sp1_all=='HW']**pA2[1]['Foliage'][2]
            Bf1_all[sp1_all=='CW']=pA2[2]['Foliage'][0]*dbh1_all[sp1_all=='CW']**pA2[2]['Foliage'][1]*h_gf1_all[sp1_all=='CW']**pA2[2]['Foliage'][2]
            Bf1_all[sp1_all=='BA']=pA2[3]['Foliage'][0]*dbh1_all[sp1_all=='BA']**pA2[3]['Foliage'][1]*h_gf1_all[sp1_all=='BA']**pA2[3]['Foliage'][2]
            
            Bbr1_all=pA2[0]['Branch'][0]*dbh1_all**pA2[0]['Branch'][1]*h_gf1_all**pA2[0]['Branch'][2]
            Bbr1_all[sp1_all=='FD']=pA2[0]['Branch'][0]*dbh1_all[sp1_all=='FD']**pA2[0]['Branch'][1]*h_gf1_all[sp1_all=='FD']**pA2[0]['Branch'][2]
            Bbr1_all[sp1_all=='HW']=pA2[1]['Branch'][0]*dbh1_all[sp1_all=='HW']**pA2[1]['Branch'][1]*h_gf1_all[sp1_all=='HW']**pA2[1]['Branch'][2]
            Bbr1_all[sp1_all=='CW']=pA2[2]['Branch'][0]*dbh1_all[sp1_all=='CW']**pA2[2]['Branch'][1]*h_gf1_all[sp1_all=='CW']**pA2[2]['Branch'][2]
            Bbr1_all[sp1_all=='BA']=pA2[3]['Branch'][0]*dbh1_all[sp1_all=='BA']**pA2[3]['Branch'][1]*h_gf1_all[sp1_all=='BA']**pA2[3]['Branch'][2]
            
            Bbk1_all=pA2[0]['Bark'][0]*dbh1_all**pA2[0]['Bark'][1]*h_gf1_all**pA2[0]['Bark'][2]
            Bbk1_all[sp1_all=='FD']=pA2[0]['Bark'][0]*dbh1_all[sp1_all=='FD']**pA2[0]['Bark'][1]*h_gf1_all[sp1_all=='FD']**pA2[0]['Bark'][2]
            Bbk1_all[sp1_all=='HW']=pA2[1]['Bark'][0]*dbh1_all[sp1_all=='HW']**pA2[1]['Bark'][1]*h_gf1_all[sp1_all=='HW']**pA2[1]['Bark'][2]
            Bbk1_all[sp1_all=='CW']=pA2[2]['Bark'][0]*dbh1_all[sp1_all=='CW']**pA2[2]['Bark'][1]*h_gf1_all[sp1_all=='CW']**pA2[2]['Bark'][2]
            Bbk1_all[sp1_all=='BA']=pA2[3]['Bark'][0]*dbh1_all[sp1_all=='BA']**pA2[3]['Bark'][1]*h_gf1_all[sp1_all=='BA']**pA2[3]['Bark'][2]
            
            #------------------------------------------------------------------ 
            # Indices
            #------------------------------------------------------------------ 
            
            # Find indices to trees found at each measurement
            c,ia,ib=gu.intersect(id0_all,id1_all)
            
            # Indices to trees that were lost (ia2) and those that appeared (ib2)
            #c2,ia2,ib2=setxor(id0_all,id1_all)            
            ib2=np.empty(0,dtype=int)
            for s in range(len(id1_all)):
                ind=np.where(id0_all==id1_all[s])[0]
                if ind.size==0:
                    ib2=np.append(ib2,s)
            
            # Index to all live trees at t0
            ind0_live=np.where(vs0_all==1)[0]
            
            # Index to all live trees at t1
            ind1_live=np.where(vs1_all==1)[0]
            
            # Index to trees that survived
            ind_surv=np.where((vs0_all[ia]==1) & (vs1_all[ib]==1))[0]
            
            ind_surv_cc1=np.where((vs0_all[ia]==1) & (vs1_all[ib]==1) & (CrownClass0_all[ia]==1))[0]
            ind_surv_cc2=np.where((vs0_all[ia]==1) & (vs1_all[ib]==1) & (CrownClass0_all[ia]==2))[0]
            ind_surv_cc3=np.where((vs0_all[ia]==1) & (vs1_all[ib]==1) & (CrownClass0_all[ia]==3))[0]
            ind_surv_cc4=np.where((vs0_all[ia]==1) & (vs1_all[ib]==1) & (CrownClass0_all[ia]==4))[0]
            ind_surv_cc5=np.where((vs0_all[ia]==1) & (vs1_all[ib]==1) & (CrownClass0_all[ia]==5))[0]
            ind_surv_cc6=np.where((vs0_all[ia]==1) & (vs1_all[ib]==1) & (CrownClass0_all[ia]==6))[0]
            
            # Index to trees that died
            ind_mort=np.where((vs0_all[ia]==1) & (vs1_all[ib]==0))[0]
            
            ind_mort_cc1=np.where((vs0_all[ia]==1) & (vs1_all[ib]==0) & (CrownClass0_all[ia]==1))[0]
            ind_mort_cc2=np.where((vs0_all[ia]==1) & (vs1_all[ib]==0) & (CrownClass0_all[ia]==2))[0]
            ind_mort_cc3=np.where((vs0_all[ia]==1) & (vs1_all[ib]==0) & (CrownClass0_all[ia]==3))[0]
            ind_mort_cc4=np.where((vs0_all[ia]==1) & (vs1_all[ib]==0) & (CrownClass0_all[ia]==4))[0]
            ind_mort_cc5=np.where((vs0_all[ia]==1) & (vs1_all[ib]==0) & (CrownClass0_all[ia]==5))[0]
            ind_mort_cc6=np.where((vs0_all[ia]==1) & (vs1_all[ib]==0) & (CrownClass0_all[ia]==6))[0]

            # Index to recruited trees that died during interval
            if ib2.size!=0:
                ind_rec=np.where((vs1_all[ib2]==1))[0]

            #------------------------------------------------------------------ 
            # Calculate stand-level variables
            #------------------------------------------------------------------ 

            # Assign each species with a numerical ID
            ID_Spc_all_t0=np.zeros((id0_all.shape))
            ind=np.where(sp0_all=='FD')[0]; ID_Spc_all_t0[ind]=1
            ind=np.where(sp0_all=='HW')[0]; ID_Spc_all_t0[ind]=2
            ind=np.where(sp0_all=='CW')[0]; ID_Spc_all_t0[ind]=3
            ind=np.where(sp0_all=='BA')[0]; ID_Spc_all_t0[ind]=4            
            ID_Spc_all_t1=np.zeros((id1_all.shape))
            ind=np.where(sp1_all=='FD')[0]; ID_Spc_all_t1[ind]=1
            ind=np.where(sp1_all=='HW')[0]; ID_Spc_all_t1[ind]=2
            ind=np.where(sp1_all=='CW')[0]; ID_Spc_all_t1[ind]=3
            ind=np.where(sp1_all=='BA')[0]; ID_Spc_all_t1[ind]=4

            # Species composition
            uSpc=np.unique(sp0_all)
            nSpc=np.zeros((uSpc.shape[0],2))
            for s in range(len(uSpc)):
                ind=np.where(sp0_all==uSpc[s])[0]
                nSpc[s,0]=s
                nSpc[s,1]=ind.size
            nSpc=np.flip(nSpc[nSpc[:,1].argsort(),:],axis=0)
            Spc1_ID_t0=uSpc[int(nSpc[0,0])]
            Spc1_Pct_t0=nSpc[0,1]/ind0.size*100
            if uSpc.size>1:
                Spc2_ID_t0=uSpc[int(nSpc[1,0])]
                Spc2_Pct_t0=nSpc[1,1]/ind0.size*100
            else:
                Spc2_ID_t0=' '
                Spc2_Pct_t0=0
            
            # Number of height measurements
            Num_H_Obs_t0=np.where(h_gf0_all[ind0_live]>0)[0].size
            
            # Number of missing or negative height growth
            H_G=(h_gf1_all[ib[ind_surv]]-h_gf0_all[ia[ind_surv]])/dt
            ind=np.where((H_G<0) | (np.isnan(H_G==True)))[0]
            Num_H_G_NegOrMis=ind.size
            
            # Stand age
            SA_t0=SA+uYear[k]-YrEst
            SA_t1=SA+uYear[k+1]-YrEst
            
            # Stand density (stems ha-1)
            SN_t0=AEF*ind0_live.size
            SN_t1=AEF*ind1_live.size
            SN_Net=(SN_t1-SN_t0)/dt
            
            # Mean diameter (cm)
            SDam_t0=np.nanmean(dbh0_all[ind0_live])
            SDam_t1=np.nanmean(dbh1_all[ind1_live])

            # Basal area (m2/ha)
            SBA_t0=AEF*np.nansum(ba0_all[ind0_live])/10000
            SBA_t1=AEF*np.nansum(ba1_all[ind1_live])/10000
            SBA_G=(SBA_t1-SBA_t0)/dt
            
            # Mean tree height (m)
            SH_obs_t0=np.nanmean(h_obs0_all[ind0_live])
            SH_obs_t1=np.nanmean(h_obs1_all[ind1_live])
            SH_obs_G=(SH_obs_t1-SH_obs_t0)/dt
            
            # Mean tree height (m)
            SH_gf_t0=np.nanmean(h_gf0_all[ind0_live])
            SH_gf_t1=np.nanmean(h_gf1_all[ind1_live])
            SH_gf_G=(SH_gf_t1-SH_gf_t0)/dt

            # Biomass, live (MgC/ha)
            SBsw_t0=dm2c*kg2Mg*AEF*np.nansum(Bsw0_all[ind0_live])
            SBf_t0=dm2c*kg2Mg*AEF*np.nansum(Bf0_all[ind0_live])
            SBbr_t0=dm2c*kg2Mg*AEF*np.nansum(Bbr0_all[ind0_live])
            SBbk_t0=dm2c*kg2Mg*AEF*np.nansum(Bbk0_all[ind0_live])
            SBsw_t1=dm2c*kg2Mg*AEF*np.nansum(Bsw1_all[ind1_live])
            SBf_t1=dm2c*kg2Mg*AEF*np.nansum(Bf1_all[ind1_live])
            SBbr_t1=dm2c*kg2Mg*AEF*np.nansum(Bbr1_all[ind1_live])
            SBbk_t1=dm2c*kg2Mg*AEF*np.nansum(Bbk1_all[ind1_live])
            
            # Biomass growth (MgC/ha/yr) 
            tmp0=np.nansum(Bsw0_all[ia[ind_surv]].flatten())
            tmp1=np.nansum(Bsw1_all[ib[ind_surv]].flatten())
            SBsw_G=dm2c*kg2Mg*AEF*(tmp1-tmp0)/dt
            
            tmp0=np.nansum(Bsw0_all[ia[ind_surv_cc1]].flatten())
            tmp1=np.nansum(Bsw1_all[ib[ind_surv_cc1]].flatten())
            SBsw_G_cc1=dm2c*kg2Mg*AEF*(tmp1-tmp0)/dt
            
            tmp0=np.nansum(Bsw0_all[ia[ind_surv_cc2]].flatten())
            tmp1=np.nansum(Bsw1_all[ib[ind_surv_cc2]].flatten())
            SBsw_G_cc2=dm2c*kg2Mg*AEF*(tmp1-tmp0)/dt
            
            tmp0=np.nansum(Bsw0_all[ia[ind_surv_cc3]].flatten())
            tmp1=np.nansum(Bsw1_all[ib[ind_surv_cc3]].flatten())
            SBsw_G_cc3=dm2c*kg2Mg*AEF*(tmp1-tmp0)/dt
            
            tmp0=np.nansum(Bsw0_all[ia[ind_surv_cc4]].flatten())
            tmp1=np.nansum(Bsw1_all[ib[ind_surv_cc4]].flatten())
            SBsw_G_cc4=dm2c*kg2Mg*AEF*(tmp1-tmp0)/dt
            
            tmp0=np.nansum(Bsw0_all[ia[ind_surv_cc5]].flatten())
            tmp1=np.nansum(Bsw1_all[ib[ind_surv_cc5]].flatten())
            SBsw_G_cc5=dm2c*kg2Mg*AEF*(tmp1-tmp0)/dt
            
            tmp0=np.nansum(Bsw0_all[ia[ind_surv_cc6]].flatten())
            tmp1=np.nansum(Bsw1_all[ib[ind_surv_cc6]].flatten())
            SBsw_G_cc6=dm2c*kg2Mg*AEF*(tmp1-tmp0)/dt
            
            tmp0=np.nansum(Bf0_all[ia[ind_surv]].flatten())
            tmp1=np.nansum(Bf1_all[ib[ind_surv]].flatten())
            SBf_G=dm2c*kg2Mg*AEF*(tmp1-tmp0)/dt
            
            tmp0=np.nansum(Bbk0_all[ia[ind_surv]].flatten())
            tmp1=np.nansum(Bbk1_all[ib[ind_surv]].flatten())
            SBbk_G=dm2c*kg2Mg*AEF*(tmp1-tmp0)/dt
            
            tmp0=np.nansum(Bbr0_all[ia[ind_surv]].flatten())
            tmp1=np.nansum(Bbr1_all[ib[ind_surv]].flatten())
            SBbr_G=dm2c*kg2Mg*AEF*(tmp1-tmp0)/dt
            
            # Biomass mortality (MgC/ha/yr)
            SBsw_M=dm2c*kg2Mg*AEF*np.nansum(Bsw1_all[ib[ind_mort]].flatten())/dt
            
            SBsw_M_cc1=dm2c*kg2Mg*AEF*np.nansum(Bsw1_all[ib[ind_mort_cc1]].flatten())/dt
            SBsw_M_cc2=dm2c*kg2Mg*AEF*np.nansum(Bsw1_all[ib[ind_mort_cc2]].flatten())/dt
            SBsw_M_cc3=dm2c*kg2Mg*AEF*np.nansum(Bsw1_all[ib[ind_mort_cc3]].flatten())/dt
            SBsw_M_cc4=dm2c*kg2Mg*AEF*np.nansum(Bsw1_all[ib[ind_mort_cc4]].flatten())/dt
            SBsw_M_cc5=dm2c*kg2Mg*AEF*np.nansum(Bsw1_all[ib[ind_mort_cc5]].flatten())/dt
            SBsw_M_cc6=dm2c*kg2Mg*AEF*np.nansum(Bsw1_all[ib[ind_mort_cc6]].flatten())/dt
            
            SBsw_Net_cc1=SBsw_G_cc1-SBsw_M_cc1
            SBsw_Net_cc2=SBsw_G_cc2-SBsw_M_cc2
            SBsw_Net_cc3=SBsw_G_cc3-SBsw_M_cc3
            SBsw_Net_cc4=SBsw_G_cc4-SBsw_M_cc4
            SBsw_Net_cc5=SBsw_G_cc5-SBsw_M_cc5
            SBsw_Net_cc6=SBsw_G_cc6-SBsw_M_cc6
            
            SBf_M=dm2c*kg2Mg*AEF*np.nansum(Bf1_all[ib[ind_mort]].flatten())/dt
            SBbr_M=dm2c*kg2Mg*AEF*np.nansum(Bbr1_all[ib[ind_mort]].flatten())/dt
            SBbk_M=dm2c*kg2Mg*AEF*np.nansum(Bbk1_all[ib[ind_mort]].flatten())/dt

            # Biomass recruitment (MgC/ha/yr)
            if ib2.size!=0:
                SBsw_R=dm2c*kg2Mg*AEF*np.nansum(Bsw1_all[ib2[ind_rec]].flatten())/dt
                SBf_R=dm2c*kg2Mg*AEF*np.nansum(Bf1_all[ib2[ind_rec]].flatten())/dt
                SBbr_R=dm2c*kg2Mg*AEF*np.nansum(Bbr1_all[ib2[ind_rec]].flatten())/dt
                SBbk_R=dm2c*kg2Mg*AEF*np.nansum(Bbk1_all[ib2[ind_rec]].flatten())/dt
            else:
                SBsw_R=float(0.0)
                SBf_R=float(0.0)
                SBbr_R=float(0.0)
                SBbk_R=float(0.0)
            
            # Damage agents
            da0_t={}
            da1_t={}
            for da in da_list:
                ind=np.where(da0_all==da)[0]
                da0_t['DA_' + str(da) + '_t0']=0
                if ind.size>0:
                    da0_t['DA_' + str(da) + '_t0']=ind.size 
                ind=np.where(da1_all==da)[0]
                da1_t['DA_' + str(da) + '_t1']=0
                if ind.size>0:
                    da1_t['DA_' + str(da) + '_t1']=ind.size     
            
            #------------------------------------------------------------------
            # Add to stand-level dataframe
            #------------------------------------------------------------------
            
            DataToAdd={'ID_Site':uSite[i],
               'ID_Plot':uPlot[j],
               'Lat':Lat,
               'Lon':Lon,
               'BEC':dSite.BGC[iSite[0]],
               'SI':dSite.SI[iSite[0]],
               'EstabType':dSite.EstTyp_Dist[iSite[0]],
               'AreaPlot':Area,
               'N_Init':dSite.StandDensity[iSite[0]].astype(float),
               'pH_Min_Init':dSite.pH_Min[iSite[0]],
               'pH_Humus_Init':dSite.pH_Humus[iSite[0]],
               'P_Min_Init':dSite.P_Min[iSite[0]],
               'P_Humus_Init':dSite.P_Humus[iSite[0]],
               'Ca_Min_Init':dSite.Ca_Min[iSite[0]],
               'N_Dose':dSite.N_Dose[iSite[j]].astype(float),
               'N_Type':dSite.N_Type[iSite[j]],
               'N_AppMeth':dSite.N_AppMeth[iSite[j]],
               'NumTreesMeas':np.array(ia.size+ib.size).astype(float),
               'Num_H_G_NegOrMis':Num_H_G_NegOrMis,
               'Spc1_ID_t0':Spc1_ID_t0,
               'Spc1_Pct_t0':Spc1_Pct_t0,
               'Spc2_ID_t0':Spc2_ID_t0,
               'Spc2_Pct_t0':Spc2_Pct_t0,
               'Year_t0':yr0,
               'Year_t1':yr1,
               'DT':dt,
               'Age_t0':SA_t0.astype(float),
               'Age_t1':SA_t1.astype(float),
               'TSF_t0':yr0-YrEst,
               'TSF_t1':yr1-YrEst,
               'N_t0':SN_t0.astype(float),
               'N_t1':SN_t1.astype(float),
               'N_Net':SN_Net.astype(float),
               'Dam_t0':SDam_t0,
               'Dam_t1':SDam_t1,
               'BA_t0':SBA_t0,
               'BA_t1':SBA_t1,
               'BA_G':SBA_G,
               'H_obs_t0':SH_obs_t0,
               'H_obs_t1':SH_obs_t1,
               'H_obs_G':SH_obs_G,
               'H_gf_t0':SH_gf_t0,
               'H_gf_t1':SH_gf_t1,
               'H_gf_G':SH_gf_G,
               'Bsw_t0':SBsw_t0,
               'Bsw_t1':SBsw_t1,
               'Bsw_G':SBsw_G,
               'Bsw_M':SBsw_M,
               'Bsw_R':SBsw_R,
               'Bsw_Net':SBsw_R+SBsw_G-SBsw_M,
               'Bf_t0':SBf_t0,
               'Bf_t1':SBf_t1,
               'Bf_G':SBf_G,
               'Bf_M':SBf_M,
               'Bf_R':SBf_R,
               'Bf_Net':SBf_R+SBf_G-SBf_M,
               'Bbr_t0':SBbr_t0,
               'Bbr_t1':SBbr_t1,
               'Bbr_G':SBbr_G,
               'Bbr_M':SBbr_M,
               'Bbr_R':SBbr_R,
               'Bbr_Net':SBbr_R+SBbr_G-SBbr_M,
               'Bbk_t0':SBbk_t0,
               'Bbk_t1':SBbk_t1,
               'Bbk_G':SBbk_G,
               'Bbk_M':SBbk_M,
               'Bbk_R':SBbk_R,
               'Bbk_Net':SBbk_R+SBbk_G-SBbk_M,
               'Bsw_G_cc1':SBsw_G_cc1,
               'Bsw_G_cc2':SBsw_G_cc2,
               'Bsw_G_cc3':SBsw_G_cc3,
               'Bsw_G_cc4':SBsw_G_cc4,
               'Bsw_G_cc5':SBsw_G_cc5,
               'Bsw_G_cc6':SBsw_G_cc6,
               'Bsw_M_cc1':SBsw_M_cc1,
               'Bsw_M_cc2':SBsw_M_cc2,
               'Bsw_M_cc3':SBsw_M_cc3,
               'Bsw_M_cc4':SBsw_M_cc4,
               'Bsw_M_cc5':SBsw_M_cc5,
               'Bsw_M_cc6':SBsw_M_cc6,
               'Bsw_Net_cc1':SBsw_Net_cc1,
               'Bsw_Net_cc2':SBsw_Net_cc2,
               'Bsw_Net_cc3':SBsw_Net_cc3,
               'Bsw_Net_cc4':SBsw_Net_cc4,
               'Bsw_Net_cc5':SBsw_Net_cc5,
               'Bsw_Net_cc6':SBsw_Net_cc6,
               'DA_A_t0':da0_t['DA_A_t0'],
               'DA_AB_t0':da0_t['DA_AB_t0'],
               'DA_AD_t0':da0_t['DA_AD_t0'],
               'DA_AE_t0':da0_t['DA_AE_t0'],
               'DA_AX_t0':da0_t['DA_AX_t0'],
               'DA_BT_t0':da0_t['DA_BT_t0'],
               'DA_CK_t0':da0_t['DA_CK_t0'],
               'DA_D_t0':da0_t['DA_D_t0'],
               'DA_DDF_t0':da0_t['DA_DDF_t0'],
               'DA_DDP_t0':da0_t['DA_DDP_t0'],
               'DA_DDS_t0':da0_t['DA_DDS_t0'],
               'DA_DM_t0':da0_t['DA_DM_t0'],
               'DA_DN_t0':da0_t['DA_DN_t0'],
               'DA_DR_t0':da0_t['DA_DR_t0'],
               'DA_DRA_t0':da0_t['DA_DRA_t0'],
               'DA_DS_t0':da0_t['DA_DS_t0'],
               'DA_DT_t0':da0_t['DA_DT_t0'],
               'DA_FCK_t0':da0_t['DA_FCK_t0'],
               'DA_FK_t0':da0_t['DA_FK_t0'],
               'DA_MT_t0':da0_t['DA_MT_t0'],
               'DA_N_t0':da0_t['DA_N_t0'],
               'DA_NGC_t0':da0_t['DA_NGC_t0'],
               'DA_NW_t0':da0_t['DA_NW_t0'],
               'DA_NWS_t0':da0_t['DA_NWS_t0'],
               'DA_NX_t0':da0_t['DA_NX_t0'],
               'DA_NY_t0':da0_t['DA_NY_t0'],
               'DA_PC_t0':da0_t['DA_PC_t0'],
               'DA_SC_t0':da0_t['DA_SC_t0'],
               'DA_SCB_t0':da0_t['DA_SCB_t0'],
               'DA_SCM_t0':da0_t['DA_SCM_t0'],
               'DA_SCT_t0':da0_t['DA_SCT_t0'],
               'DA_A_t1':da1_t['DA_A_t1'],
               'DA_AB_t1':da1_t['DA_AB_t1'],
               'DA_AD_t1':da1_t['DA_AD_t1'],
               'DA_AE_t1':da1_t['DA_AE_t1'],
               'DA_AX_t1':da1_t['DA_AX_t1'],
               'DA_BT_t1':da1_t['DA_BT_t1'],
               'DA_CK_t1':da1_t['DA_CK_t1'],
               'DA_D_t1':da1_t['DA_D_t1'],
               'DA_DDF_t1':da1_t['DA_DDF_t1'],
               'DA_DDP_t1':da1_t['DA_DDP_t1'],
               'DA_DDS_t1':da1_t['DA_DDS_t1'],
               'DA_DM_t1':da1_t['DA_DM_t1'],
               'DA_DN_t1':da1_t['DA_DN_t1'],
               'DA_DR_t1':da1_t['DA_DR_t1'],
               'DA_DRA_t1':da1_t['DA_DRA_t1'],
               'DA_DS_t1':da1_t['DA_DS_t1'],
               'DA_DT_t1':da1_t['DA_DT_t1'],
               'DA_FCK_t1':da1_t['DA_FCK_t1'],
               'DA_FK_t1':da1_t['DA_FK_t1'],
               'DA_MT_t1':da1_t['DA_MT_t1'],
               'DA_N_t1':da1_t['DA_N_t1'],
               'DA_NGC_t1':da1_t['DA_NGC_t1'],
               'DA_NW_t1':da1_t['DA_NW_t1'],
               'DA_NWS_t1':da1_t['DA_NWS_t1'],
               'DA_NX_t1':da1_t['DA_NX_t1'],
               'DA_NY_t1':da1_t['DA_NY_t1'],
               'DA_PC_t1':da1_t['DA_PC_t1'],
               'DA_SC_t1':da1_t['DA_SC_t1'],
               'DA_SCB_t1':da1_t['DA_SCB_t1'],
               'DA_SCM_t1':da1_t['DA_SCM_t1'],
               'DA_SCT_t1':da1_t['DA_SCT_t1']}   
            
            for k in DataToAdd:
                #if (type(DataToAdd[k])==np.str_):
                #    sobs[k]=sobs[k].astype('U')
                sobs[k][cntS]=DataToAdd[k]
            cntS=cntS+1
            #sobs=sobs.append(DataToAdd.copy(),ignore_index=True)
            
            #------------------------------------------------------------------
            # Add survivor data to tree-level dataframe
            #------------------------------------------------------------------
            
            ID_Spc_t0=ID_Spc_all_t0[ia[ind_surv]]
            
            ID_Tree=id0_all[ia[ind_surv]]
            
            VS_t0=vs0_all[ia[ind_surv]]
                        
            BA_t0=ba0_all[ia[ind_surv]]
            BA_t1=ba1_all[ib[ind_surv]]
            BA_G=(BA_t1-BA_t0)/dt
            
            H_obs_t0=h_obs0_all[ia[ind_surv]]
            H_obs_t1=h_obs1_all[ib[ind_surv]]
            H_obs_G=(H_obs_t1-H_obs_t0)/dt
            
            H_mod_t0=h_mod0_all[ia[ind_surv]]
            H_mod_t1=h_mod1_all[ib[ind_surv]]
            H_mod_G=(H_mod_t1-H_mod_t0)/dt
            
            H_gf_t0=h_gf0_all[ia[ind_surv]]
            H_gf_t1=h_gf1_all[ib[ind_surv]]
            H_gf_G=(H_gf_t1-H_gf_t0)/dt
            
            H_Fit_ITM_t0=H_Fit_ITM0_all[ia[ind_surv]]
            H_Fit_ITM_t1=H_Fit_ITM1_all[ib[ind_surv]]
            
            Bsw_t0=dm2c*Bsw0_all[ia[ind_surv]]
            Bsw_t1=dm2c*Bsw1_all[ib[ind_surv]]
            Bsw_G=dm2c*(Bsw1_all[ib[ind_surv]]-Bsw0_all[ia[ind_surv]])/dt
            
            Bsw_mod_t0=dm2c*Bsw_mod0_all[ia[ind_surv]]
            Bsw_mod_t1=dm2c*Bsw_mod1_all[ib[ind_surv]]
              
            TreeClass_t0=TreeClass0_all[ia[ind_surv]]
            TreeClass_t1=TreeClass1_all[ib[ind_surv]]
            
            CrownClass_t0=CrownClass0_all[ia[ind_surv]]
            CrownClass_t1=CrownClass1_all[ib[ind_surv]]
            
            TopDead_t0=dtop0_all[ia[ind_surv]]
            TopDead_t1=dtop1_all[ib[ind_surv]]
            
            TopBroke_t0=btop0_all[ia[ind_surv]]
            TopBroke_t1=btop1_all[ib[ind_surv]]
            
            DA_t0=da0_all[ia[ind_surv]]
            DA_t1=da1_all[ib[ind_surv]]
            
            # Stand biomass of larger trees
            tlist=np.zeros((ind_surv.size,3))
            tlist[:,0]=np.arange(0,ind_surv.size,1)
            tlist[:,1]=Bsw_t0
            tlist=np.nan_to_num(tlist)
            tlist=np.flip(tlist[tlist[:,1].argsort(),:],axis=0)
            tlist[:,2]=np.cumsum(tlist[:,1]*kg2Mg*AEF)            
            SBswLT_t0=np.zeros(ind_surv.shape)
            SBswLT_t0[tlist[:,0].astype(int)]=tlist[:,2]
            #plt.plot(Bsw_t0,SBswLT,'.')
            
            for v in range(ind_surv.size):
                DataToAdd={'ID_Site':uSite[i],
                   'ID_Plot':uPlot[j],
                   'ID_Tree':ID_Tree[v],
                   'Lat':Lat,
                   'Lon':Lon,'BEC':dSite.BGC[iSite[0]],'SI':dSite.SI[iSite[0]],
                   'AEF':AEF,'N_Init':dSite.StandDensity[iSite[0]],
                   'pH_Min_Init':dSite.pH_Min[iSite[0]],
                   'N_Dose':dSite.N_Dose[iSite[j]].astype(float),'Year_t0':yr0,'Year_t1':yr1,'DT':dt,
                   'TSF_t0':yr0-YrEst,'TSF_t1':yr1-YrEst,
                   'SA_t0':SA_t0.astype(float),'SA_t1':SA_t1.astype(float),
                   'SN_t0':SN_t0,'SN_t1':SN_t1,
                   'SBsw_t0':SBsw_t0,'SBsw_t1':SBsw_t1,'SBswLT_t0':SBswLT_t0[v],'SBswLT_t1':0,
                   'ID_Spc_t0':ID_Spc_t0[v],'Mort':0,'Rec':0,'TreeClass_t0':TreeClass_t0[v],
                   'CrownClass_t0':CrownClass_t0[v],'TopDead_t0':TopDead_t0[v],'TopDead_t1':TopDead_t1[v],
                   'TopBroke_t0':TopBroke_t0[v],'TopBroke_t1':TopBroke_t1[v],'DA_t0':DA_t0[v],'DA_t1':DA_t1[v],
                   'BA_t0':BA_t0[v],'BA_t1':BA_t1[v],'BA_G':BA_G[v],
                   'H_obs_t0':H_obs_t0[v],'H_obs_t1':H_obs_t1[v],'H_obs_G':H_obs_G[v],
                   'H_mod_t0':H_mod_t0[v],'H_mod_t1':H_mod_t1[v],'H_mod_G':H_mod_G[v],
                   'H_gf_t0':H_gf_t0[v],'H_gf_t1':H_gf_t1[v],'H_gf_G':H_gf_G[v],
                   'Bsw_t0':Bsw_t0[v],'Bsw_t1':Bsw_t1[v],'Bsw_G':Bsw_G[v],
                   'Bsw_mod_t0':Bsw_mod_t0[v],'Bsw_mod_t1':Bsw_mod_t1[v],
                   'H_Fit_ITM_t0':H_Fit_ITM_t0[v],'H_Fit_ITM_t1':H_Fit_ITM_t1[v]}
                for k in DataToAdd:
                    tobs[k][cntT]=DataToAdd[k]
                cntT=cntT+1

            #------------------------------------------------------------------
            # Add mortality data to tree-level dataframe
            #------------------------------------------------------------------
            
            if ind_mort.size>0:
                
                ID_Spc_t0=ID_Spc_all_t0[ia[ind_mort]]
            
                ID_Tree=id0_all[ia[ind_mort]]
            
                VS_t0=vs0_all[ia[ind_mort]]
                        
                BA_t0=ba0_all[ia[ind_mort]]
                BA_t1=ba1_all[ib[ind_mort]]
                BA_G=(BA_t1-BA_t0)/dt
            
                H_obs_t0=h_obs0_all[ia[ind_mort]]
                H_obs_t1=h_obs1_all[ib[ind_mort]]
                H_obs_G=np.nan*np.ones(ind_mort.size,dtype=float)
                
                H_mod_t0=h_mod0_all[ia[ind_mort]]
                H_mod_t1=h_mod1_all[ib[ind_mort]]
                H_mod_G=np.nan*np.ones(ind_mort.size,dtype=float)
                
                H_gf_t0=h_gf0_all[ia[ind_mort]]
                H_gf_t1=h_gf1_all[ib[ind_mort]]
                H_gf_G=np.nan*np.ones(ind_mort.size,dtype=float)                
            
                H_Fit_ITM_t0=H_Fit_ITM0_all[ia[ind_mort]]
                H_Fit_ITM_t1=H_Fit_ITM1_all[ib[ind_mort]]
            
                Bsw_t0=dm2c*Bsw0_all[ia[ind_mort]]
                Bsw_t1=dm2c*Bsw1_all[ib[ind_mort]]
                Bsw_G=dm2c*(Bsw1_all[ib[ind_mort]]-Bsw0_all[ia[ind_mort]])/dt
                
                TreeClass_t0=TreeClass0_all[ia[ind_mort]]
                TreeClass_t1=TreeClass1_all[ib[ind_mort]]
            
                CrownClass_t0=CrownClass0_all[ia[ind_mort]]
                CrownClass_t1=CrownClass1_all[ib[ind_mort]]
            
                TopDead_t0=dtop0_all[ia[ind_mort]]
                TopDead_t1=dtop1_all[ib[ind_mort]]
            
                TopBroke_t0=btop0_all[ia[ind_mort]]
                TopBroke_t1=btop1_all[ib[ind_mort]]
            
                DA_t0=da0_all[ia[ind_mort]]
                DA_t1=da1_all[ib[ind_mort]]
                
                # Stand biomass of larger trees
                tlist=np.zeros((ind_mort.size,3))
                tlist[:,0]=np.arange(0,ind_mort.size,1)
                tlist[:,1]=Bsw_t0
                tlist=np.nan_to_num(tlist)
                tlist=np.flip(tlist[tlist[:,1].argsort(),:],axis=0)
                tlist[:,2]=np.cumsum(tlist[:,1]*kg2Mg*AEF)            
                SBswLT_t0=np.zeros(ind_mort.shape)
                SBswLT_t0[tlist[:,0].astype(int)]=tlist[:,2]
                #plt.plot(Bsw_t0,SBswLT,'.')                
                
                for v in range(ind_mort.size):
                    DataToAdd={'ID_Site':uSite[i],
                               'ID_Plot':uPlot[j],
                               'ID_Tree':ID_Tree[v],
                               'Lat':Lat,
                               'Lon':Lon,'BEC':dSite.BGC[iSite[0]],'SI':dSite.SI[iSite[0]],
                               'AEF':AEF,'N_Init':dSite.StandDensity[iSite[0]],
                               'pH_Min_Init':dSite.pH_Min[iSite[0]],
                               'N_Dose':dSite.N_Dose[iSite[j]].astype(float),'Year_t0':yr0,'Year_t1':yr1,'DT':dt,
                               'TSF_t0':yr0-YrEst,'TSF_t1':yr1-YrEst,
                               'SA_t0':SA_t0.astype(float),'SA_t1':SA_t1.astype(float),
                               'SN_t0':SN_t0,'SN_t1':SN_t1,
                               'SBsw_t0':SBsw_t0,'SBsw_t1':SBsw_t1,'SBswLT_t0':SBswLT_t0[v],'SBswLT_t1':0,
                               'ID_Spc_t0':ID_Spc_t0[v],'Mort':1,'Rec':0,'TreeClass_t0':TreeClass_t0[v],
                               'CrownClass_t0':CrownClass_t0[v],'TopDead_t0':TopDead_t0[v],'TopDead_t1':TopDead_t1[v],
                               'TopBroke_t0':TopBroke_t0[v],'TopBroke_t1':TopBroke_t1[v],'DA_t0':DA_t0[v],'DA_t1':DA_t1[v],
                               'BA_t0':BA_t0[v],'BA_t1':BA_t1[v],'BA_G':BA_G[v],
                               'H_obs_t0':H_obs_t0[v],'H_obs_t1':H_obs_t1[v],'H_obs_G':H_obs_G[v],
                               'H_mod_t0':H_mod_t0[v],'H_mod_t1':H_mod_t1[v],'H_mod_G':H_mod_G[v],
                               'H_gf_t0':H_gf_t0[v],'H_gf_t1':H_gf_t1[v],'H_gf_G':H_gf_G[v],
                               'Bsw_t0':Bsw_t0[v],'Bsw_t1':Bsw_t1[v],'Bsw_G':Bsw_G[v],
                               'Bsw_mod_t0':np.nan,'Bsw_mod_t1':np.nan,
                               'H_Fit_ITM_t0':H_Fit_ITM_t0[v],'H_Fit_ITM_t1':H_Fit_ITM_t1[v]}
                    for k in DataToAdd:
                        tobs[k][cntT]=DataToAdd[k]
                    cntT=cntT+1 

            #------------------------------------------------------------------
            # Add recruitment data to tree-level dataframe
            #------------------------------------------------------------------
            
            if ib2.size>0:
                
                n=ib2[ind_rec].size
                
                ID_Spc_t1=ID_Spc_all_t1[ib2[ind_rec]]
            
                ID_Tree=id1_all[ib2[ind_rec]]
            
                VS_t1=vs1_all[ib2[ind_rec]]
                        
                BA_t0=np.nan*np.ones(n,dtype=float)
                BA_t1=ba1_all[ib2[ind_rec]]
                BA_G=np.nan*np.ones(n,dtype=float)
            
                H_obs_t0=np.nan*np.ones(n,dtype=float)
                H_obs_t1=h_obs1_all[ib2[ind_rec]]
                H_obs_G=np.nan*np.ones(n,dtype=float)
            
                H_mod_t0=np.nan*np.ones(n,dtype=float)
                H_mod_t1=h_mod1_all[ib2[ind_rec]]
                H_mod_G=np.nan*np.ones(n,dtype=float)
                
                H_gf_t0=np.nan*np.ones(n,dtype=float)
                H_gf_t1=h_gf1_all[ib2[ind_rec]]
                H_gf_G=np.nan*np.ones(n,dtype=float)
            
                H_Fit_ITM_t0=np.nan*np.ones(n,dtype=float)
                H_Fit_ITM_t1=H_Fit_ITM1_all[ib2[ind_rec]]
            
                Bsw_t0=np.nan*np.ones(n,dtype=float)
                Bsw_t1=dm2c*Bsw1_all[ib2[ind_rec]]
                Bsw_G=np.nan*np.ones(n,dtype=float)
                
                TreeClass_t0=np.nan*np.ones(n,dtype=float)
                TreeClass_t1=TreeClass1_all[ib2[ind_rec]]
            
                CrownClass_t0=np.nan*np.ones(n,dtype=float)
                CrownClass_t1=CrownClass1_all[ib2[ind_rec]]
            
                TopDead_t0=np.nan*np.ones(n,dtype=float)
                TopDead_t1=dtop1_all[ib2[ind_rec]]
            
                TopBroke_t0=np.nan*np.ones(n,dtype=float)
                TopBroke_t1=btop1_all[ib2[ind_rec]]
            
                DA_t0=np.nan*np.ones(n,dtype=float)
                DA_t1=da1_all[ib2[ind_rec]]
                
                # Stand biomass of larger trees
                tlist=np.zeros((ind_rec.size,3))
                tlist[:,0]=np.arange(0,ind_rec.size,1)
                tlist[:,1]=Bsw_t1
                tlist=np.nan_to_num(tlist)
                tlist=np.flip(tlist[tlist[:,1].argsort(),:],axis=0)
                tlist[:,2]=np.cumsum(tlist[:,1]*kg2Mg*AEF)            
                SBswLT_t1=np.zeros(ind_rec.shape)
                SBswLT_t1[tlist[:,0].astype(int)]=tlist[:,2]
                #plt.plot(Bsw_t0,SBswLT,'.')
            
                for v in range(ind_rec.size):
                    DataToAdd={'ID_Site':uSite[i],
                               'ID_Plot':uPlot[j],
                               'ID_Tree':ID_Tree[v],
                               'Lat':Lat,
                               'Lon':Lon,'BEC':dSite.BGC[iSite[0]],'SI':dSite.SI[iSite[0]],
                               'AEF':AEF,'N_Init':dSite.StandDensity[iSite[0]],
                               'pH_Min_Init':dSite.pH_Min[iSite[0]],
                               'N_Dose':dSite.N_Dose[iSite[j]].astype(float),'Year_t0':yr0,'Year_t1':yr1,'DT':dt,
                               'TSF_t0':yr0-YrEst,'TSF_t1':yr1-YrEst,
                               'SA_t0':SA_t0.astype(float),'SA_t1':SA_t1.astype(float),
                               'SN_t0':SN_t0,'SN_t1':SN_t1,
                               'SBsw_t0':SBsw_t0,'SBsw_t1':SBsw_t1,'SBswLT_t0':0,'SBswLT_t1':SBswLT_t1[v],
                               'ID_Spc_t0':ID_Spc_t1[v],'Mort':0,'Rec':1,'TreeClass_t0':TreeClass_t0[v],
                               'CrownClass_t0':CrownClass_t0[v],'TopDead_t0':TopDead_t0[v],'TopDead_t1':TopDead_t1[v],
                               'TopBroke_t0':TopBroke_t0[v],'TopBroke_t1':TopBroke_t1[v],'DA_t0':DA_t0[v],'DA_t1':DA_t1[v],
                               'BA_t0':BA_t0[v],'BA_t1':BA_t1[v],'BA_G':BA_G[v],
                               'H_obs_t0':H_obs_t0[v],'H_obs_t1':H_obs_t1[v],'H_obs_G':H_obs_G[v],
                               'H_mod_t0':H_mod_t0[v],'H_mod_t1':H_mod_t1[v],'H_mod_G':H_mod_G[v],
                               'H_gf_t0':H_gf_t0[v],'H_gf_t1':H_gf_t1[v],'H_gf_G':H_gf_G[v],
                               'Bsw_t0':Bsw_t0[v],'Bsw_t1':Bsw_t1[v],'Bsw_G':Bsw_G[v],
                               'Bsw_mod_t0':np.nan,'Bsw_mod_t1':np.nan,
                               'H_Fit_ITM_t0':H_Fit_ITM_t0[v],'H_Fit_ITM_t1':H_Fit_ITM_t1[v]}
                    for k in DataToAdd:
                        tobs[k][cntT]=DataToAdd[k]
                    cntT=cntT+1                

#------------------------------------------------------------------------------
# Remove excess zeros
#------------------------------------------------------------------------------

for k in sobs:
    sobs[k]=sobs[k][0:cntS]
for k in tobs:
    tobs[k]=tobs[k][0:cntT]

#%% Convert to float
              
tobs['SA_t0']=tobs['SA_t0'].astype(float)
tobs['SBswLT_t0']=tobs['SBswLT_t0'].astype(float)
tobs['TSF_t0']=tobs['TSF_t0'].astype(float)
tobs['TSF_t1']=tobs['TSF_t1'].astype(float)
tobs['N_Dose']=tobs['N_Dose'].astype(float)

#%% Save to file

gu.opickle(meta['Paths']['Project']['Outputs'] + '\\EP703_SL.pkl',sobs)
gu.opickle(meta['Paths']['Project']['Outputs'] + '\\EP703_TL.pkl',tobs)


#%% Add environmental data

# Import data
sobs=gu.ipickle(meta['Paths']['Project']['Outputs'] + '\\EP703_SL.pkl')
tobs=gu.ipickle(meta['Paths']['Project']['Outputs'] + '\\EP703_TL.pkl')

# Unique location of installations
uLL=np.unique(np.column_stack([sobs['Lat'],sobs['Lon']]),axis=0)

# Spatial reference system of NACID
SRS_NACID=pyproj.Proj('+proj=lcc +lat_1=49 +lat_2=77 +lat_0=40 +lon_0=-100 +x_0=0 +y_0=0 +no_defs +a=6378137 +rf=298.257222101 +to_meter=1')

# Map of plots
plt.close('all')
fig,ax=plt.subplots(figsize=(9,9))
#gdf1.plot(ax=ax,color=[0.5,0.5,0.5],edgecolor='k',label='Land')
plt.plot(uLL[:,1],uLL[:,0],'go')
ax.grid(color='k', linestyle='-', linewidth=0.25)

# Calculate projected coordinates
x,y=SRS_NACID(uLL[:,1],uLL[:,0])

# Import NACID X and Y grids
mat=spio.loadmat(meta['Paths']['NACID'] + '\\grid.mat',squeeze_me=True)
ix=np.where(np.asarray(mat['grd'].dtype.names)=='xL')[0][0]
iy=np.where(np.asarray(mat['grd'].dtype.names)=='yL')[0][0]
X=mat['grd'][()][ix].astype(float)[0,:]
Y=mat['grd'][()][iy].astype(float)[:,0]

# Index to NACID grid
ixy=np.zeros((x.shape[0],2),dtype=int)
for i in range(len(ixy)):
    adx=np.abs(x[i]-X)
    ady=np.abs(y[i]-Y)
    ixy[i,0]=int(np.where(adx==np.min(adx))[0][0])
    ixy[i,1]=int(np.where(ady==np.min(ady))[0][0])

# Check that we are indexing the right location
#plt.close('all')
#plt.imshow(tmean0,vmin=0,vmax=20)
#for i in range(len(x)):
#    plt.plot(ixy[i,0],ixy[i,1],'yo')

# Define a time period for the annual data
tv=np.arange(1970,2016,1)

# Import annual total nitrogen deposition
ndep=np.zeros((len(tv),len(x)))
for i in range(len(tv)):    
    if tv[i]>2013:
        continue # No data past 2013
    pthin=meta['Paths']['NACID'] + '\\Grids\\NACID_ndep_tot_ann_abso_comp_hist_v1\\NACID_ndep_tot_ann_abso_comp_hist_v1_' + str(tv[i]) + '.mat'
    z=spio.loadmat(pthin,squeeze_me=True)
    ndep0=z['z'][()][1].astype(float)*z['z'][()][0]
    for j in range(len(x)):
        print(str(i) + ' ' + str(j))
        ndep[i,j]=ndep0[ixy[j,1],ixy[j,0]]
#plt.plot(tv,np.mean(ndep,axis=1))

# Import mean annual temperature normal
tmean_ann_n=np.zeros((1,len(x)))
for i in range(12):
    pthin=meta['Paths']['NACID'] + '\\Grids\\NACID_tmean_mon_norm_1971to2000_si_hist_v1\\NACID_tmean_mon_norm_1971to2000_si_hist_v1_' + str(i+1) + '.mat'
    z=spio.loadmat(pthin,squeeze_me=True)
    idat=np.where(np.asarray(z['z'].dtype.names)=='Data')[0][0]
    iSF=np.where(np.asarray(z['z'].dtype.names)=='ScaleFactor')[0][0]
    tmean0=z['z'][()][idat].astype(float)*z['z'][()][iSF]    
    for j in range(len(x)):
        tmean_ann_n[0,j]=tmean_ann_n[0,j]+tmean0[ixy[j,1],ixy[j,0]]
tmean_ann_n=np.tile(tmean_ann_n/12,(tv.size,1))
 
# Import mean annual precipitation normal
prcp_ann_n=np.zeros((1,len(x)))
for i in range(12):
    pthin=meta['Paths']['NACID'] + '\\Grids\\NACID_prcp_mon_norm_1971to2000_si_hist_v1\\NACID_prcp_mon_norm_1971to2000_si_hist_v1_' + str(i+1) + '.mat'
    z=spio.loadmat(pthin,squeeze_me=True)
    idat=np.where(np.asarray(z['z'].dtype.names)=='Data')[0][0]
    iSF=np.where(np.asarray(z['z'].dtype.names)=='ScaleFactor')[0][0]
    prcp0=z['z'][()][idat].astype(float)*z['z'][()][iSF]    
    for j in range(len(x)):
        prcp_ann_n[0,j]=prcp_ann_n[0,j]+prcp0[ixy[j,1],ixy[j,0]]
prcp_ann_n=np.tile(prcp_ann_n,(tv.size,1))       

# Import mean annual temperature actual
tmean_ann_r=np.zeros((len(tv),len(x)))
for i in range(len(tv)):    
    pthin=meta['Paths']['NACID'] + '\\Grids\\NACID_tmean_ann_abso_si_hist_v1\\NACID_tmean_ann_abso_si_hist_v1_' + str(tv[i]) + '.mat'
    z=spio.loadmat(pthin,squeeze_me=True)
    idat=np.where(np.asarray(z['z'].dtype.names)=='Data')[0][0]
    iSF=np.where(np.asarray(z['z'].dtype.names)=='ScaleFactor')[0][0]
    tmean0=z['z'][()][idat].astype(float)*z['z'][()][iSF]
    for j in range(len(x)):
        print(str(i) + ' ' + str(j))
        tmean_ann_r[i,j]=tmean0[ixy[j,1],ixy[j,0]]
        
# Import annual precipitation actual
prcp_ann_r=np.zeros((len(tv),len(x)))
for i in range(len(tv)):    
    pthin=meta['Paths']['NACID'] + '\\Grids\\NACID_prcp_ann_abso_si_hist_v1\\NACID_prcp_ann_abso_si_hist_v1_' + str(tv[i]) + '.mat'
    z=spio.loadmat(pthin,squeeze_me=True)
    idat=np.where(np.asarray(z['z'].dtype.names)=='Data')[0][0]
    iSF=np.where(np.asarray(z['z'].dtype.names)=='ScaleFactor')[0][0]
    prcp0=z['z'][()][idat].astype(float)*z['z'][()][iSF]
    for j in range(len(x)):
        print(str(i) + ' ' + str(j))
        prcp_ann_r[i,j]=prcp0[ixy[j,1],ixy[j,0]]

# Import warm-season mean temperature
tmean_ws_r=np.zeros((len(tv),len(x)))
for i in range(len(tv)):    
    pthin=meta['Paths']['NACID'] + '\\Grids\\NACID_tmean_gs_abso_si_hist_v1\\NACID_tmean_gs_abso_si_hist_v1_' + str(tv[i]) + '.mat'
    z=spio.loadmat(pthin,squeeze_me=True)
    idat=np.where(np.asarray(z['z'].dtype.names)=='Data')[0][0]
    iSF=np.where(np.asarray(z['z'].dtype.names)=='ScaleFactor')[0][0]
    tmean0=z['z'][()][idat].astype(float)*z['z'][()][iSF]
    for j in range(len(x)):
        print(str(i) + ' ' + str(j))
        tmean_ws_r[i,j]=tmean0[ixy[j,1],ixy[j,0]]

# Import warm-season mean soil water content
ws_ws_r=np.zeros((len(tv),len(x)))
for i in range(len(tv)):    
    pthin=meta['Paths']['NACID'] + '\\Grids\\NACID_ws_tmw_gs_abso_comp_hist_v1\\NACID_ws_tmw_gs_abso_comp_hist_v1_' + str(tv[i]) + '.mat'
    z=spio.loadmat(pthin,squeeze_me=True)
    idat=np.where(np.asarray(z['z'].dtype.names)=='Data')[0][0]
    iSF=np.where(np.asarray(z['z'].dtype.names)=='ScaleFactor')[0][0]
    ws0=z['z'][()][idat].astype(float)*z['z'][()][iSF]
    for j in range(len(x)):
        print(str(i) + ' ' + str(j))
        ws_ws_r[i,j]=ws0[ixy[j,1],ixy[j,0]]
        
# Import warm-season mean climatic water deficit
cwd_ws_r=np.zeros((len(tv),len(x)))
for i in range(len(tv)):    
    pthin=meta['Paths']['NACID'] + '\\Grids\\NACID_cwd_tmw_gs_abso_comp_hist_v1\\NACID_cwd_tmw_gs_abso_comp_hist_v1_' + str(tv[i]) + '.mat'
    z=spio.loadmat(pthin,squeeze_me=True)
    idat=np.where(np.asarray(z['z'].dtype.names)=='Data')[0][0]
    iSF=np.where(np.asarray(z['z'].dtype.names)=='ScaleFactor')[0][0]
    cwd0=z['z'][()][idat].astype(float)*z['z'][()][iSF]
    for j in range(len(x)):
        print(str(i) + ' ' + str(j))
        cwd_ws_r[i,j]=cwd0[ixy[j,1],ixy[j,0]] 

# Import warm-season mean climatic water deficit
etp_ws_r=np.zeros((len(tv),len(x)))
for i in range(len(tv)):    
    pthin=meta['Paths']['NACID'] + '\\Grids\\NACID_etp_tmw_gs_abso_comp_hist_v1\\NACID_etp_tmw_gs_abso_comp_hist_v1_' + str(tv[i]) + '.mat'
    z=spio.loadmat(pthin,squeeze_me=True)
    idat=np.where(np.asarray(z['z'].dtype.names)=='Data')[0][0]
    iSF=np.where(np.asarray(z['z'].dtype.names)=='ScaleFactor')[0][0]
    etp0=z['z'][()][idat].astype(float)*z['z'][()][iSF]
    for j in range(len(x)):
        print(str(i) + ' ' + str(j))
        etp_ws_r[i,j]=etp0[ixy[j,1],ixy[j,0]] 

# Add to SOBS structure
sobs['cwd_mjjas_r']=np.zeros(sobs['ID_Site'].size)
sobs['etp_mjjas_r']=np.zeros(sobs['ID_Site'].size)
sobs['ndep_ann_r']=np.zeros(sobs['ID_Site'].size)
sobs['prcp_ann_r']=np.zeros(sobs['ID_Site'].size)
sobs['prcp_ann_n']=np.zeros(sobs['ID_Site'].size)
sobs['tmean_ann_r']=np.zeros(sobs['ID_Site'].size)
sobs['tmean_ann_n']=np.zeros(sobs['ID_Site'].size)
sobs['tmean_mjjas_r']=np.zeros(sobs['ID_Site'].size)
sobs['ws_mjjas_r']=np.zeros(sobs['ID_Site'].size)
for i in range(len(x)):
    isite=np.where( (sobs['Lat']==uLL[i,0]) & (sobs['Lon']==uLL[i,1]) )[0]
    for j in range(len(isite)):
        it=np.where( (tv>=sobs['Year_t0'][isite[j]]) & (tv<=sobs['Year_t1'][isite[j]]) )[0]
        sobs['ndep_ann_r'][isite[j]]=np.mean(ndep[it,i])
        sobs['cwd_mjjas_r'][isite[j]]=np.mean(cwd_ws_r[it,i])
        sobs['etp_mjjas_r'][isite[j]]=np.mean(etp_ws_r[it,i])
        sobs['prcp_ann_n'][isite[j]]=np.mean(prcp_ann_n[it,i])
        sobs['prcp_ann_r'][isite[j]]=np.mean(prcp_ann_r[it,i])
        sobs['tmean_ann_n'][isite[j]]=np.mean(tmean_ann_n[it,i])
        sobs['tmean_ann_r'][isite[j]]=np.mean(tmean_ann_r[it,i])
        sobs['tmean_mjjas_r'][isite[j]]=np.mean(tmean_ws_r[it,i])        
        sobs['ws_mjjas_r'][isite[j]]=np.mean(ws_ws_r[it,i])


#%% Save to file 

gu.opickle(meta['Paths']['Project']['Outputs'] + '\\EP703_SL.pkl',sobs)
gu.opickle(meta['Paths']['Project']['Outputs'] + '\\EP703_TL.pkl',tobs)

