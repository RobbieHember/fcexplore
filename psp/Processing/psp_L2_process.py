
'''
PSP-NADB Process Level 2 Data

Import raw (Level 1) data for Canadian data sources and the national US
FIA database and produce Level 2 stand- and tree-level structures SOBS
and tobs['

Notes
This will process both PSPs and TSPs

PNW-IDB2 and Alaksa source data are organized so differently that they
needs to be run using a seperate script.

The order of fields in the SOBS and TOBS structures must match the
equivalent script for the PNW-IDB2

If there are no parameters from the Canadian National Biomass equations,
then the allometric equations default to the DBH-based parameters from
Jenkins (2004).

Discrepency between dN1 and dN2 in some jurisdictions (BC, QC). It is
because of live trees that remain alive and pass a utilization level.
The AEF decreases, affecting (N1-N0), but it has no effect oN Recr
and mortality.

% dN1=(sobs['N L t1'][cnt]-sobs['N L t0'][cnt])/sobs['Delta t'][cnt]
% dN2=sobs['N_R'][cnt]-sobs['N Mort'][cnt]-sobs['N_H'][cnt]
% if abs(dN1-dN2)>10
%   a=1
%   a0=[id0_all aef0_all vital_status0_all dbh0_all]
%   a1=[id1_all aef1_all vital_status1_all dbh1_all]
%   tmp=np.nan*np.ones(500,8)
%   tmp(1:length(a0),1:4)=a0
%   tmp(1:length(a1),5:end)=a1
%   xlswrite('E:\Data\ForestInventory\PSP-NADB\Data\prob',tmp)
% end

Litterfall was derived from constant functions of biomass
(see Kurz et al 2009 Li et al 2003 Kurz et al 1996).

Growth of trees that die is recorded, but it is not reliable. Sometimes
they don't measure the DBH and H of the trees that die or are harvested.
They might have DBH <=0 or NaN. In Quebec, they are -999.These were set
to the DBH and H of the previous measurement.

'''

#%% Import modules

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import gc as garc
import time
import warnings
from fcgadgets.macgyver import utilities_general as gu
from fcexplore.psp.Processing import psp_utilities as utl
warnings.filterwarnings("ignore")

#%% Import project info

meta={}
meta['Paths']={}
meta['Paths']['DB']=r'C:\Users\rhember\Documents\Data\GroundPlots\PSP-NADB2'

meta=utl.ImportParameters(meta)

#%% Parmaters and constants

meta['prm']={}
meta['prm']['Carbon Content']=0.5

#%% Select source databases to run

jurL=['BC']

#%% Import parameters

#[tmp,tmp,sdb]=xlsread([pthin 'Parameters\PSP-NADB_AllometricParameters.xlsx'],'LUT','A3:B17')
#ind=np.where(strcmp(db,sdb[:,1])==1 | strncmp(db,sdb[:,1],3)==1)
#SourceDB=cell2mat(sdb(ind,1))
#if isempty(SourceDB)==1 SourceDB=15 end

# State, prov, territory code
#[tmp,tmp,sdb]=xlsread([pthin 'Parameters\PSP-NADB_AllometricParameters.xlsx'],'LUT','C3:D63')
#if strncmp(db,'FIA',3)==1
#  ind=np.where(strcmp(db(end-1:end),sdb[:,1])==1)
#else
#  ind=np.where(strcmp(db,sdb[:,1])==1)
#end
#State=cell2mat(sdb(ind,1))

# Species codes from allometric equation table
#[spc_id,tmp,tmp]=xlsread([pthin 'Parameters\PSP-NADB_AllometricParameters.xlsx'],'Allometry','A5:A85')
#[tmp,spc_code,tmp]=xlsread([pthin 'Parameters\PSP-NADB_AllometricParameters.xlsx'],'Allometry','B5:B85') clear tmp
#[tmp,tmp,spc_code_replace]=xlsread([pthin 'Parameters\PSP-NADB_AllometricParameters.xlsx'],'Allometry','Q5:Q85') clear tmp

# Plant function type
#[pft_flg,tmp,tmp]=xlsread([pthin 'Parameters\PSP-NADB_AllometricParameters.xlsx'],'Allometry','U5:U85')

# Parameters (biomass):
#pt_amd2b_sw=xlsread([pthin 'Parameters\PSP-NADB_AllometricParameters.xlsx'],'Allometry','D5:F85')
#pt_amd2b_bk=xlsread([pthin 'Parameters\PSP-NADB_AllometricParameters.xlsx'],'Allometry','G5:I85')
#pt_amd2b_f=xlsread([pthin 'Parameters\PSP-NADB_AllometricParameters.xlsx'],'Allometry','J5:L85')
#pt_amd2b_br=xlsread([pthin 'Parameters\PSP-NADB_AllometricParameters.xlsx'],'Allometry','M5:O85')

#%% Gap-fill missing allometric parameters

# for i=1:length(pt_amd2b_sw)
#   if isnan(pt_amd2b_sw(i,1))==1
#     ind=np.where(spc_id==spc_code_replace{i})
#     pt_amd2b_sw(i,:)=pt_amd2b_sw(ind,:)
#     pt_amd2b_bk(i,:)=pt_amd2b_bk(ind,:)
#     pt_amd2b_br(i,:)=pt_amd2b_br(ind,:)
#     pt_amd2b_f(i,:)=pt_amd2b_f(ind,:)
#   end
# end

#%% QA flags

QA_Flag1=[]

#%% Loop through each source database

for iJur in range(len(jurL)):

    # Import Level 1 tree and plot data structures, TL and PLT
    d=gu.ipickle(meta['Paths']['DB'] + '\\Processed\\L1\\L1_' + jurL[iJur] + '.pkl')
    tl=d['tl']
    pl=d['pl']
    del d

    # Unique species
    uSpc=np.unique(tl['ID Species'])

    # Unique plots
    uP=np.unique(tl['ID Plot'])

    # Parameters based on DBH+H (Ung et al. 2008 Lambert et al. 2005)
    beta={}
    tissueL=['bk','br','f','sw']
    for tissue in tissueL:
        beta[tissue]=np.nan*np.ones( (tl['ID Species'].size,3) )

    tl['ID PFT']=np.zeros(tl['ID Species'].size,dtype=int)
    for iU in range(uSpc.size):
        ind_all0=np.where(tl['ID Species']==uSpc[iU])[0]
        ind_all1=np.where(meta['Allo B']['ID']==uSpc[iU])[0]
        tl['ID PFT'][ind_all0]=meta['Allo B']['PFT'][ind_all1]
        for tissue in tissueL:
            for iP in range(3):
                beta[tissue][ind_all0,iP]=meta['Allo B']['sw' + str(iP+1)][ind_all1]

    # Carbon (kgC tree-1)
    for tissue in tissueL:
        tl['C' + tissue]=meta['prm']['Carbon Content']*(beta[tissue][:,0]*tl['DBH']**beta[tissue][:,1]*tl['H']**beta[tissue][:,2])

    # Merchcantable stemwood
    tl['Csw125']=tl['Csw']
    ind=np.where(tl['DBH']<12.5)[0]
    tl['Csw125'][ind]=0

    # Aboveground biomass carbon  (kg C tree-1)
    tl['Cag']=tl['Csw']+tl['Cbk']+tl['Cbr']+tl['Cf']

    # Roots from Li et al (2003), coarse root litterfall of 2 % yr-1 from
    # Kurz et al. (1996). ***
    tl['Cr']=np.zeros(tl['Csw'].size)

    ind=np.where(tl['ID PFT']==meta['LUT']['PFT']['Deciduous'])[0]
    tl['Cr'][ind]=1.576*tl['Cag'][ind]**0.615

    ind=np.where(tl['ID PFT']==meta['LUT']['PFT']['Coniferous'])[0]
    tl['Cr'][ind_all1]=0.222*tl['Cag'][ind_all1]

    # Total biomass carbon  (kg C tree-1)
    tl['Ctot']=tl['Csw']+tl['Cbk']+tl['Cbr']+tl['Cf']+tl['Cr']

    # Biomass Litterfall due to litterfall
    # *** Aboveground components taken from Table 3 of Kurz et al (2009)

    FracLitterfallF=np.zeros(tl['Csw'].size)
    FracLitterfallBk=np.zeros(tl['Csw'].size)
    FracLitterfallBr=np.zeros(tl['Csw'].size)
    FracLitterfallR=np.zeros(tl['Csw'].size)

    ind=np.where(tl['ID PFT']==meta['LUT']['PFT']['Deciduous'])[0]
    FracLitterfallBk[ind]=0.035
    FracLitterfallBr[ind]=0.035
    FracLitterfallF[ind]=0.95

    ind=np.where(tl['ID PFT']==meta['LUT']['PFT']['Coniferous'])[0]
    FracLitterfallBk[ind_all1]=0.035
    FracLitterfallBr[ind_all1]=0.035
    FracLitterfallF[ind_all1]=0.10

    tl['Cbk Litterfall']=FracLitterfallBk*tl['Cbk']
    tl['Cbr Litterfall']=FracLitterfallBr*tl['Cbr']
    tl['Cf Litterfall']=FracLitterfallF*tl['Cf']

    pFine=0.072+0.354*np.exp(-0.06*(2*tl['Cr']))
    pCoarse=1.0-pFine
    tl['Cr Litterfall']=0.641*pFine*tl['Cr']+0.02*pCoarse*tl['Cr']

    # Initialize stand observation structure
    slvL=['ID Source','ID Plot','ID Visit','Lat','Lon','Elev','Aspect','Slope','Area',
         'Jourisiction','Ecozone CA L1','Ecozone BC L1','Ecozone BC L2','Plot Type',
         'Utilization Level','Management','SI','Year t0','Year t1','Month t0','Month t1',
         'Delta t','Age t0','Age t1','Spc1 ID','Spc1 %BA','Spc1 %N','Spc2 ID','Spc2 %BA',
         'Spc2 %N','Spc3 ID','Spc3 %BA','Spc3 %N','Spc4 ID','Spc4 %BA','Spc4 %N',
         'N L t0','N L t1','N Recr','N Mort','N Harv','N Net',
         'Dam L t0','Dam L t1','Dam G Surv',
         'Dqm L t0','Dqm L t1',
         'BA L t0','BA L t1','BA G Surv','BA Mort',
         'H L t0','H L t1','H G Surv',
         'H L max t0','H L max t1',
         'Cag L t0','Cag L t1','Cag G Surv','Cag G Recr','Cag Mort','Cag Harv','Cag D t0','Cag D t1',
         'Cbk L t0','Cbk L t1','Cbk G Surv','Cbk G Recr','Cbk Mort',
         'Cbr L t0','Cbr L t1','Cbr G Surv','Cbr G Recr','Cbr Mort',
         'Cf L t0','Cf L t1','Cf G Surv','Cf G Recr','Cf Mort',
         'Cr L t0','Cr L t1','Cr G Surv','Cr G Recr','Cr Mort',
         'Csw L t0','Csw L t1','Csw G Surv','Csw G Recr','Csw Mort',
         'Csw Mort Insect','Csw Mort Disease','Csw Mort Plant Competition','Csw Mort Animal','Csw Mort Fire','Csw Mort Frost','Csw Mort Drought','Csw Mort SnowAndIce','Csw Mort Wind','Csw Mort Silviculture','Csw Mort Other','Csw Mort Unknown',
         'Csw Harv','Csw D t0','Csw D t1','Csw125 L t0','Csw125 L t1',
         'Ctot L t0','Ctot L t1','Ctot G Surv','Ctot G Recr','Ctot Mort','Ctot Harv','Ctot Litf','Ctot D t0','Ctot D t1',
         'N Insect t0','N Instect t1',
         'N Disease t0','N Disease t1',
         'N Plant Competition t0','N Plant Competition t1',
         'N Animal t0','N Animal t1',
         'N Fire t0','N Fire t1',
         'N Frost t0','N Frost t1',
         'N Drought t0','N Drought t1',
         'N SnowAndIce t0','N SnowAndIce t1',
         'N Wind t0','N Wind t1',
         'N Silviculture t0','N Silviculture t1',
         'N Other t0','N Other t1',
         'N Unknown t0','N Unknown t1']

    sobs={}
    for v in slvL:
        sobs[v]=np.nan*np.ones(pl['ID Plot'].size)

    # Initialize fields in the Level 2 tree-level structure TOBS
    tlvL=['ID Source','ID Plot','ID Visit','ID Tree','ID Species','Lat','Lon','Elev','Aspect','Slope',
         'Ecozone BC L1','Ecozone BC L2','Ecozone CA L1','Plot Type','Stand Management',
         'Year t0','Year t1','Delta t','Stand Spc1 ID','Stand Spc1 %BA','Stand Spc1 %N','Stand Age t0','Stand Age t1',
         'Stand N L t0','Stand N L t1','Stand Cag L t0','Stand Cag L t1','Stand Cag L Larger t0','Stand BA L Larger t0',
         'Stand DA Insect t0','Stand DA Instect t1','Stand DA Disease t0','Stand DA Disease t1',
         'Stand DA Fire t0','Stand DA Fire t1','Stand DA Animal t0','Stand DA Animal t1','Stand DA Weather t0','Stand DA WEather t1',
         'Stand DA Silviculture t0','Stand DA Silviculture t1','Stand DA Mechanical t0','Stand DA Mechanical t1',
         'Stand DA SnowAndIce t0','Stand DA SnowAndIce t1','Stand DA Unknown t0','Stand DA Unknown t1',
         'AEF','DA t0','DA t1','DAI t0','DAI t1','DAP t0','DAP t1','DA IDW t0','DA IDW t1',
         'Mortality','Recruitment','DBH t0','DBH t1','DBH G','BA t0','BA t1','BA G',
         'H t0','H t1','H G','H Obs t0','H Obs t1',
         'Cag t0','Cag t1','Cag G','Csw t0','Csw t1','Csw G']

    tobs={}
    for v in tlvL:
        tobs[v]=np.nan*np.ones(tl['ID Plot'].size)

    # Initialize counter
    cnt=1

    # Loop through plots
    for iP in range(uP.size):

        # Index to plot from tree structure
        iPlot=np.where(tl['ID Plot']==uP[iP])[0]

        # Data in Level 1 PLT structure may be missing from Level 1 TL
        # structure, skip and record in QAQC accounting
        if iPlot.size==0:
            print('Warning, plot missing from tl structure!')
            sobs['ID Plot'][cnt]=np.nan #pl['ID Plot'][iPlot_t0]
            #sobs['QAQC_Error1'][cnt]=1
            continue

        # Index to inventory years for plot, I
        uVisit=np.unique(tl['ID Visit'][iPlot])

        # Loop through measurement intervals
        for iVisit in range(uVisit.size):

            # Index to plot structure
            iPlot_t0=np.where( (pl['ID Plot']==uP[iP]) & (pl['ID Visit']==uVisit[iVisit]) )[0]
            iPlot_t1=np.where( (pl['ID Plot']==uP[iP]) & (pl['ID Visit']==uVisit[iVisit]+1) )[0]

            # Empty field plot, record and continue
            if iPlot_t0.size==0:
                #print('Plot ID in TL file, but missing from PL file: ' + str(uP[iP]))
                continue

            # Flag indicating PSP (1) or TSP (0)
            if iPlot_t1.size==0:
                flg_PSP=0
            else:
                flg_PSP=1

            # Isolate trees from plot I and inventory year iVisit
            ind_all0=np.where( (tl['ID Plot']==uP[iP]) & (tl['ID Visit']==uVisit[iVisit]) )[0]

            if flg_PSP==1:

                # Isolate trees from plot I and inventory year iVisit+1
                ind_all1=np.where( (tl['ID Plot']==uP[iP]) & (tl['ID Visit']==uVisit[iVisit]+1) )[0]

                # Find indices to current and previous trees
                c,ia1,ib1=np.intersect1d(tl['ID Tree'][ind_all0],tl['ID Tree'][ind_all1],return_indices=True)

                # Find index to lost trees - trees that were counted in t0 and
                # disappear (ia2).
                ia2=np.where(np.isin(tl['ID Tree'][ind_all0],tl['ID Tree'][ind_all1])==False)[0]

                # Find trees that were Recruited - trees that were
                # not counted at t0, but were at t1 [ib2].
                ib2=np.where(np.isin(tl['ID Tree'][ind_all1],tl['ID Tree'][ind_all0])==False)[0]

                # All trees disappear, record and continue.
                if ia1.size==0:
                    print('Warning, PL file thinks it is a PSP, but there are no matching trees at t1.')
                    QA_Flag1.append(iPlot_t0)
                    #sobs['ID Plot'][cnt]=pl['ID Plot'][iPlot_t0]
                    #sobs['ID Plot'][cnt]=pl['ID Plot'][iPlot_t1]
                    continue

            # Populate stand-level site variables
            sobs['ID Plot'][cnt]=pl['ID Plot'][iPlot_t0]
            sobs['ID Visit'][cnt]=uVisit[iVisit]
            sobs['Lat'][cnt]=pl['Lat'][iPlot_t0]
            sobs['Lon'][cnt]=pl['Lon'][iPlot_t0]
            #sobs['Area'][cnt]=pl['Area'][iPlot_t0]
            #sobs['Utilization Level'][cnt]=np.nan
            sobs['Year t0'][cnt]=pl['Year'][iPlot_t0]
            sobs['Month t0'][cnt]=pl['Month'][iPlot_t0]
            sobs['Age t0'][cnt]=pl['Age'][iPlot_t0]
            #sobs['Management'][cnt]=pl['Management'][iPlot_t0]

            #sobs['StandOrigin'][cnt]=pl['StandOrigin'][iPlot_t0]
            sobs['Plot Type'][cnt]=pl['Plot Type'][iPlot_t0]
            sobs['Ecozone BC L1'][cnt]=pl['Ecozone BC L1'][iPlot_t0]
            sobs['Ecozone BC L2'][cnt]=pl['Ecozone BC L2'][iPlot_t0]
            #sobs['SI'][cnt]=pl['SI_SourceDB'][iPlot_t0]

            if flg_PSP==1:
                sobs['Year t1'][cnt]=pl['Year'][iPlot_t1]
                sobs['Month t1'][cnt]=pl['Month'][iPlot_t1]
                dt=sobs['Year t1'][cnt]-sobs['Year t0'][cnt]
                sobs['Delta t'][cnt]=dt
                sobs['Age t1'][cnt]=pl['Age'][iPlot_t1]

            # Extract trees

            # All trees at t0
            id0_all=tl['ID Tree'][ind_all0]
            dbh0_all=tl['DBH'][ind_all0]
            ba0_all=np.pi*(tl['DBH'][ind_all0]/2)**2
            h0_all=tl['H'][ind_all0]
            h0_obs_all=tl['H Obs'][ind_all0]
            vital_status0_all=tl['Vital Status'][ind_all0]
            aef0_all=tl['AEF'][ind_all0]
            spc0_all=tl['ID Species'][ind_all0]
            da0_all=tl['ID DA1'][ind_all0]

            Cag0_all=tl['Cag'][ind_all0]
            Cbk0_all=tl['Cbk'][ind_all0]
            Cbr0_all=tl['Cbr'][ind_all0]
            Cf0_all=tl['Cf'][ind_all0]
            Cr0_all=tl['Cr'][ind_all0]
            Csw0_all=tl['Csw'][ind_all0]
            Csw125_0_all=tl['Csw125'][ind_all0]
            Ctot0_all=tl['Ctot'][ind_all0]
            Cf_Litterfall0_all=tl['Cf Litterfall'][ind_all0]
            Cbk_Litterfall0_all=tl['Cbk Litterfall'][ind_all0]
            Cbr_Litterfall0_all=tl['Cbr Litterfall'][ind_all0]
            Cr_Litterfall0_all=tl['Cr Litterfall'][ind_all0]

            # All trees at t1
            if flg_PSP==1:
                id1_all=tl['ID Tree'][ind_all1]
                dbh1_all=tl['DBH'][ind_all1]
                ba1_all=np.pi*(tl['DBH'][ind_all1]/2)**2
                h1_all=tl['H'][ind_all1]
                h1_obs_all=tl['H Obs'][ind_all1]
                vital_status1_all=tl['Vital Status'][ind_all1]
                aef1_all=tl['AEF'][ind_all1]
                spc1_all=tl['ID Species'][ind_all1]
                da1_all=tl['ID DA1'][ind_all1]

                Cag1_all=tl['Cag'][ind_all1]
                Cbk1_all=tl['Cbk'][ind_all1]
                Cbr1_all=tl['Cbr'][ind_all1]
                Cf1_all=tl['Cf'][ind_all1]
                Csw1_all=tl['Csw'][ind_all1]
                Csw125_1_all=tl['Csw'][ind_all1]
                Ctot1_all=tl['Ctot'][ind_all1]
                Cf_Litterfall1_all=tl['Cf Litterfall'][ind_all1]
                Cbk_Litterfall1_all=tl['Cbk Litterfall'][ind_all1]
                Cbr_Litterfall1_all=tl['Cbr Litterfall'][ind_all1]
                Cr_Litterfall1_all=tl['Cr Litterfall'][ind_all1]

            # Define indices to trees with a specific status and fate

            # Index to trees that were alive in t0 (including those that died)
            ind_live0=np.where(vital_status0_all==1)[0]

            # Index to trees that were dead at t0
            ind_dead0=np.where(vital_status0_all==0)[0]

            if flg_PSP==1:

                # Index to trees that were alive in t1 (including those that were
                # Recrited)
                ind_live1=np.where(vital_status1_all==1)[0]

                # Index to trees that were dead at t1
                ind_dead1=np.where(vital_status1_all==0)[0]

                # Index to survivors
                ind_surv=np.where( (vital_status0_all[ia1]==1) & (vital_status1_all[ib1]==1) )[0]

                # Index to trees that were counted as dead and then counted as alive
                ind_reborn=np.where( (vital_status0_all[ia1]==0) & (vital_status1_all[ib1]==1) )[0]

                # Index to trees that died during interval
                ind_mort_v=np.where( (vital_status0_all[ia1]==1) & (vital_status1_all[ib1]==0) )[0]

                # Index to trees that died during interval
                #ind_mort_ins=np.where( (vital_status0_all[ia1]==1) & (vital_status1_all[ib1]==0) & (da0_all[ia1]==1) | (vital_status0_all[ia1]==1) & (vital_status1_all[ib1]==0) & (da1_all[ib1]==1) )[0]
                #ind_mort_dis=np.where( (vital_status0_all[ia1]==1) & (vital_status1_all[ib1]==0) & (da0_all[ia1]==2) | (vital_status0_all[ia1]==1) & (vital_status1_all[ib1]==0) & (da1_all[ib1]==2) )[0]
                #ind_mort_fir=np.where( (vital_status0_all[ia1]==1) & (vital_status1_all[ib1]==0) & (da0_all[ia1]==3) | (vital_status0_all[ia1]==1) & (vital_status1_all[ib1]==0) & (da1_all[ib1]==3) )[0]
                #ind_mort_ani=np.where( (vital_status0_all[ia1]==1) & (vital_status1_all[ib1]==0) & (da0_all[ia1]==4) | (vital_status0_all[ia1]==1) & (vital_status1_all[ib1]==0) & (da1_all[ib1]==4) )[0]
                #ind_mort_wea=np.where( (vital_status0_all[ia1]==1) & (vital_status1_all[ib1]==0) & (da0_all[ia1]==5) | (vital_status0_all[ia1]==1) & (vital_status1_all[ib1]==0) & (da1_all[ib1]==5) )[0]
                #ind_mort_sno=np.where( (vital_status0_all[ia1]==1) & (vital_status1_all[ib1]==0) & (da0_all[ia1]==10) | (vital_status0_all[ia1]==1) & (vital_status1_all[ib1]==0) & (da1_all[ib1]==10) )[0]

                # Index to Recrited trees that died during interval
                ind_mort_r=np.where(vital_status1_all[ib2]==0)[0]

                # Index to Recrited trees that remain alive
                ind_rec=np.where(vital_status1_all[ib2]==1)[0]

                # Index to harvested trees
                ind_harv=np.where( (vital_status0_all[ia1]==1) & (vital_status1_all[ib1]==2) )[0]

            # Cleaning of trees that die
            # *** Sometimes they don't measure the DBH and H of the trees that die
            # or are harvested. They might have DBH <=0 or NaN. In Quebec, they are
            # -999.These were set to the DBH and H of the previous measurement. ***
            if flg_PSP==1:

                # Mortality
                iBadD=np.where( (dbh0_all[ia1[ind_mort_v]]>0) & (dbh1_all[ib1[ind_mort_v]]<=0) | (dbh0_all[ia1[ind_mort_v]]>0) & (np.isnan(dbh1_all[ib1[ind_mort_v]])==1) )[0]
                if iBadD.size!=0:
                    aef1_all[ib1[ind_mort_v[iBadD]]]=aef0_all[ia1[ind_mort_v[iBadD]]]
                    dbh1_all[ib1[ind_mort_v[iBadD]]]=dbh0_all[ia1[ind_mort_v[iBadD]]]
                    ba1_all[ib1[ind_mort_v[iBadD]]]=ba0_all[ia1[ind_mort_v[iBadD]]]
                    h1_all[ib1[ind_mort_v[iBadD]]]=h0_all[ia1[ind_mort_v[iBadD]]]
                    h1_obs_all[ib1[ind_mort_v[iBadD]]]=h0_obs_all[ia1[ind_mort_v[iBadD]]]
                    Csw1_all[ib1[ind_mort_v[iBadD]]]=Csw0_all[ia1[ind_mort_v[iBadD]]]
                    Csw125_1_all[ib1[ind_mort_v[iBadD]]]=Csw125_0_all[ia1[ind_mort_v[iBadD]]]
                    Cag1_all[ib1[ind_mort_v[iBadD]]]=Cag0_all[ia1[ind_mort_v[iBadD]]]
                    Ctot1_all[ib1[ind_mort_v[iBadD]]]=Ctot0_all[ia1[ind_mort_v[iBadD]]]

                # Harvesting
                iBadD=np.where( (dbh0_all[ia1[ind_harv]]>0) & (dbh1_all[ib1[ind_harv]]<=0) | (dbh0_all[ia1[ind_harv]]>0) & (np.isnan(dbh1_all[ib1[ind_harv]])==1) )[0]
                if iBadD.size!=0:
                    aef1_all[ib1[ind_harv[iBadD]]]=aef0_all[ia1[ind_harv[iBadD]]]
                    spc1_all[ib1[ind_harv[iBadD]]]=spc0_all[ia1[ind_harv[iBadD]]]
                    dbh1_all[ib1[ind_harv[iBadD]]]=dbh0_all[ia1[ind_harv[iBadD]]]
                    ba1_all[ib1[ind_harv[iBadD]]]=ba0_all[ia1[ind_harv[iBadD]]]
                    h1_all[ib1[ind_harv[iBadD]]]=h0_all[ia1[ind_harv[iBadD]]]
                    h1_obs_all[ib1[ind_harv[iBadD]]]=h0_obs_all[ia1[ind_harv[iBadD]]]
                    Csw1_all[ib1[ind_harv[iBadD]]]=Csw0_all[ia1[ind_harv[iBadD]]]
                    Csw125_1_all[ib1[ind_harv[iBadD]]]=Csw125_0_all[ia1[ind_harv[iBadD]]]
                    Cag1_all[ib1[ind_harv[iBadD]]]=Cag0_all[ia1[ind_harv[iBadD]]]
                    Ctot1_all[ib1[ind_harv[iBadD]]]=Ctot0_all[ia1[ind_harv[iBadD]]]

            # Species composition at t0

            # Get a list of all unique species in the stand
            uSpc0=np.unique(spc0_all)

            SpcComp=np.zeros((uSpc0.size,3))
            SpcComp[:,0]=uSpc0
            for iSpc in range(uSpc0.size):
                indS=np.where(spc0_all==uSpc0[iSpc])[0]

                # Composition by basal area
                BA_tot=np.nansum((np.pi*(dbh0_all/200)**2)*aef0_all)
                BA_spc=np.nansum((np.pi*(dbh0_all[indS]/200)**2)*aef0_all[indS])
                SpcComp[iSpc,1]=BA_spc/BA_tot

                # Composition by tree density
                N_spc=np.sum(aef0_all[indS])
                N_tot=np.sum(aef0_all)
                SpcComp[iSpc,2]=N_spc/N_tot

            ord=np.argsort(SpcComp[:,1])
            SpcComp=np.flip(SpcComp[ord,:],axis=0)

            for iSpc in range( np.minimum(4,uSpc0.size) ):
                sobs['Spc' + str(iSpc+1) + ' ID'][cnt]=SpcComp[iSpc,0]
                sobs['Spc' + str(iSpc+1) + ' %BA'][cnt]=np.round(SpcComp[iSpc,1],decimals=3)
                sobs['Spc' + str(iSpc+1) + ' %N'][cnt]=np.round(SpcComp[iSpc,2],decimals=3)

            # Populate stand-level state variables

            # Stand density (stems ha-1)
            sobs['N L t0'][cnt]=np.round(np.sum(aef0_all[ind_live0]),decimals=2)

            #sobs['N_R'][cnt]=np.round(sum(aef1_all[ib2[ind_rec]])/dt*100)/100

            # Diameter, arithmetic mean (cm)
            sobs['Dam L t0'][cnt]=np.round(np.nanmean(dbh0_all[ind_live0]),decimals=1)

            # Diameter, quadratic mean (cm)
            sobs['Dqm L t0'][cnt]=np.round(np.sqrt(np.mean(dbh0_all[ind_live0]**2)),decimals=1)

            # Basal area (m2 ha-1)
            sobs['BA L t0'][cnt]=np.round(np.nansum((np.pi*(dbh0_all[ind_live0]/200)**2)*aef0_all[ind_live0]),decimals=2)

            # Height (m)
            sobs['H L t0'][cnt]=np.round(np.nanmean(h0_all[ind_live0]),decimals=2)

            # Aboveground carbon, live (Mg C ha-1)
            sobs['Cag L t0'][cnt]=np.nansum(aef0_all[ind_live0]*0.001*Cag0_all[ind_live0])

            # Bark carbon, live (Mg C ha-1)
            sobs['Cbk L t0'][cnt]=np.nansum(aef0_all[ind_live0]*0.001*Cbk0_all[ind_live0])

            # Branch carbon, live (Mg C ha-1)
            sobs['Cbr L t0'][cnt]=np.nansum(aef0_all[ind_live0]*0.001*Cbr0_all[ind_live0])

            # Foliage carbon, live (Mg C ha-1)
            sobs['Cf L t0'][cnt]=np.nansum(aef0_all[ind_live0]*0.001*Cf0_all[ind_live0])

            # Stemwood carbon, live (Mg C ha-1)
            sobs['Csw L t0'][cnt]=np.nansum(aef0_all[ind_live0]*0.001*Csw0_all[ind_live0])

            # Stemwood carbon, live (Mg C ha-1)
            sobs['Csw125 L t0'][cnt]=np.nansum(aef0_all[ind_live0]*0.001*Csw125_0_all[ind_live0])

            # Total carbon, live (Mg C ha-1)
            sobs['Ctot L t0'][cnt]=np.nansum(aef0_all[ind_live0]*0.001*Ctot0_all[ind_live0])

            # Stemwood carbon, dead (Mg C ha-1)
            sobs['Csw D t0'][cnt]=np.nansum(aef0_all[ind_dead0]*0.001*Csw0_all[ind_dead0])

            # Aboveground carbon, dead (Mg C ha-1)
            sobs['Cag D t0'][cnt]=np.nansum(aef0_all[ind_dead0]*0.001*Cag0_all[ind_dead0])

            # Total carbon, dead (Mg C ha-1)
            sobs['Ctot D t0'][cnt]=np.nansum(aef0_all[ind_dead0]*0.001*Ctot0_all[ind_dead0])

            # 'N Insect t0'
            # 'N Disease t0'
            # 'N Plant Competition t0'
            # 'N Animal t0'
            # 'N Fire t0'
            # 'N Frost t0'
            # 'N Drought t0'
            # 'N SnowAndIce t0'
            # 'N Wind t0'
            # 'N Silviculture t0'
            # 'N Other t0'
            # 'N Unknown t0'

            # Populate stand-level growth and litterfall of survivors

            if flg_PSP==1:

                # Stand density (stems ha-1)
                sobs['N L t1'][cnt]=np.sum(aef1_all[ind_live1])

                # Diameter, arithmetic mean (cm)
                sobs['Dam L t1'][cnt]=np.nanmean(dbh1_all[ind_live1])

                # Diameter, quadratic mean (cm)
                sobs['Dqm L t1'][cnt]=np.sqrt(np.mean(dbh1_all[ind_live1]**2))

                # Basal area (m2 ha-1)
                sobs['BA L t1'][cnt]=np.nansum((np.pi*(dbh1_all[ind_live1]/200)**2)*aef1_all[ind_live1])

                # Height (m)
                sobs['H L t1'][cnt]=np.nanmean(h1_all[ind_live1])

                # Stemwood carbon, live (Mg C ha-1)
                sobs['Csw L t1'][cnt]=np.nansum(aef1_all[ind_live1]*0.001*Csw1_all[ind_live1])

                # Stemwood carbon, live (Mg C ha-1)
                sobs['Csw125 L t1'][cnt]=np.nansum(aef1_all[ind_live1]*0.001*Csw125_1_all[ind_live1])

                # Foliage carbon, live (Mg C ha-1)
                sobs['Cf L t1'][cnt]=np.nansum(aef1_all[ind_live1]*0.001*Cf1_all[ind_live1])

                # Branch carbon, live (Mg C ha-1)
                sobs['Cbr L t1'][cnt]=np.nansum(aef1_all[ind_live1]*0.001*Cbr1_all[ind_live1])

                # Bark carbon, live (Mg C ha-1)
                sobs['Cbk L t1'][cnt]=np.nansum(aef1_all[ind_live1]*0.001*Cbk1_all[ind_live1])

                # Aboveground carbon, live (Mg C ha-1)
                sobs['Cag L t1'][cnt]=np.nansum(aef1_all[ind_live1]*0.001*Cag1_all[ind_live1])

                # total carbon, live (Mg C ha-1)
                sobs['Ctot L t1'][cnt]=np.nansum(aef1_all[ind_live1]*0.001*Ctot1_all[ind_live1])

                # Stemwood carbon, dead (Mg C ha-1)
                sobs['Csw D t1'][cnt]=np.nansum(aef1_all[ind_dead1]*0.001*Csw1_all[ind_dead1])

                # Aboveground carbon, dead (Mg C ha-1)
                sobs['Cag D t1'][cnt]=np.nansum(aef1_all[ind_dead1]*0.001*Cag1_all[ind_dead1])

                # Total carbon, dead (Mg C ha-1)
                sobs['Ctot D t1'][cnt]=np.nansum(aef1_all[ind_dead1]*0.001*Ctot1_all[ind_dead1])

                # Diameter, arithmetic mean, growth (cm yr-1)
                Dam_t0_surv=dbh0_all[ia1[ind_surv]]
                Dam_t1_surv=dbh1_all[ib1[ind_surv]]
                dDam_surv=(Dam_t1_surv-Dam_t0_surv)/dt
                sobs['Dam G Surv'][cnt]=np.nanmean(dDam_surv)

                # Basal area growth (m2 ha-1)
                BA_t0_surv=(np.pi*(dbh0_all[ia1[ind_surv]]/200)**2)*aef0_all[ia1[ind_surv]]
                BA_t1_surv=(np.pi*(dbh1_all[ib1[ind_surv]]/200)**2)*aef1_all[ib1[ind_surv]]
                dBA_surv=(BA_t1_surv-BA_t0_surv)/dt
                dBA_surv_sum=np.nansum(dBA_surv)
                sobs['BA G Surv'][cnt]=np.nansum(dBA_surv)

                # Height growth (m yr-1)
                H_t0_surv=h0_all[ia1[ind_surv]]
                H_t1_surv=h1_all[ib1[ind_surv]]
                dH_surv=(H_t1_surv-H_t0_surv)/dt
                sobs['H G Surv'][cnt]=np.nanmean(dH_surv)

                # Aboveground carbon growth (Mg C ha-1 yr-1)
                Cag_t0_surv=aef0_all[ia1[ind_surv]]*0.001*Cag0_all[ia1[ind_surv]]
                Cag_t1_surv=aef1_all[ib1[ind_surv]]*0.001*Cag1_all[ib1[ind_surv]]
                dCag_surv=(Cag_t1_surv-Cag_t0_surv)/dt
                sobs['Cag G Surv'][cnt]=np.nansum(dCag_surv)

                # Bark carbon growth (Mg C ha-1 yr-1)
                Cbk_t0_surv=aef0_all[ia1[ind_surv]]*0.001*Cbk0_all[ia1[ind_surv]]
                Cbk_t1_surv=aef1_all[ib1[ind_surv]]*0.001*Cbk1_all[ib1[ind_surv]]
                dCbk_surv=(Cbk_t1_surv-Cbk_t0_surv)/dt
                sobs['Cbk G Surv'][cnt]=np.nansum(dCbk_surv)

                # Branch carbon growth (Mg C ha-1 yr-1)
                Cbr_t0_surv=aef0_all[ia1[ind_surv]]*0.001*Cbr0_all[ia1[ind_surv]]
                Cbr_t1_surv=aef1_all[ib1[ind_surv]]*0.001*Cbr1_all[ib1[ind_surv]]
                dCbr_surv=(Cbr_t1_surv-Cbr_t0_surv)/dt
                sobs['Cbr G Surv'][cnt]=np.nansum(dCbr_surv)

                # Foliage carbon growth (Mg C ha-1 yr-1)
                Cf_t0_surv=aef0_all[ia1[ind_surv]]*0.001*Cf0_all[ia1[ind_surv]]
                Cf_t1_surv=aef1_all[ib1[ind_surv]]*0.001*Cf1_all[ib1[ind_surv]]
                dCf_surv=(Cf_t1_surv-Cf_t0_surv)/dt
                sobs['Cf G Surv'][cnt]=np.nansum(dCf_surv)

                # Stemwood carbon growth (Mg C ha-1 yr-1)
                Csw_t0_surv=aef0_all[ia1[ind_surv]]*0.001*Csw0_all[ia1[ind_surv]]
                Csw_t1_surv=aef1_all[ib1[ind_surv]]*0.001*Csw1_all[ib1[ind_surv]]
                dCsw_surv=(Csw_t1_surv-Csw_t0_surv)/dt
                sobs['Csw G Surv'][cnt]=np.nansum(dCsw_surv)

                # Total carbon growth (Mg C ha-1 yr-1)
                Ctot_t0_surv=aef0_all[ia1[ind_surv]]*0.001*Ctot0_all[ia1[ind_surv]]
                Ctot_t1_surv=aef1_all[ib1[ind_surv]]*0.001*Ctot1_all[ib1[ind_surv]]
                dCtot_surv=(Ctot_t1_surv-Ctot_t0_surv)/dt
                sobs['Ctot G Surv'][cnt]=np.nansum(dCtot_surv)

                # Litterfall Litterfall from foliage biomass (Mg C ha-1 yr-1)
                Cf_Litterfall_t0=aef0_all[ia1[ind_surv]]*0.001*Cf_Litterfall0_all[ia1[ind_surv]]
                Cf_Litterfall_t1=aef1_all[ib1[ind_surv]]*0.001*Cf_Litterfall1_all[ib1[ind_surv]]
                Cf_Litterfall_Ave=(Cf_Litterfall_t0+Cf_Litterfall_t1)/2
                Cf_Litterfall=np.nansum(Cf_Litterfall_Ave)

                # Litterfall Litterfall from bark biomass (Mg C ha-1 yr-1)
                Cbk_Litterfall_t0=aef0_all[ia1[ind_surv]]*0.001*Cbk_Litterfall0_all[ia1[ind_surv]]
                Cbk_Litterfall_t1=aef1_all[ib1[ind_surv]]*0.001*Cbk_Litterfall1_all[ib1[ind_surv]]
                Cbk_Litterfall_Ave=(Cbk_Litterfall_t0+Cbk_Litterfall_t1)/2
                Cbk_Litterfall=np.nansum(Cbk_Litterfall_Ave)

                # Litterfall Litterfall from branch biomass (Mg C ha-1 yr-1)
                Cbr_Litterfall_t0=aef0_all[ia1[ind_surv]]*0.001*Cbr_Litterfall0_all[ia1[ind_surv]]
                Cbr_Litterfall_t1=aef1_all[ib1[ind_surv]]*0.001*Cbr_Litterfall1_all[ib1[ind_surv]]
                Cbr_Litterfall_Ave=(Cbr_Litterfall_t0+Cbr_Litterfall_t1)/2
                Cbr_Litterfall=np.nansum(Cbr_Litterfall_Ave)

                # Litterfall Litterfall from root biomass (Mg C ha-1 yr-1)
                Cr_Litterfall_t0=aef0_all[ia1[ind_surv]]*0.001*Cr_Litterfall0_all[ia1[ind_surv]]
                Cr_Litterfall_t1=aef1_all[ib1[ind_surv]]*0.001*Cr_Litterfall1_all[ib1[ind_surv]]
                Cr_Litterfall_Ave=(Cr_Litterfall_t0+Cr_Litterfall_t1)/2
                Cr_Litterfall=np.nansum(Cr_Litterfall_Ave)

                # Litterfall total (Mg C ha-1 yr-1)
                sobs['Ctot Litf'][cnt]=Cf_Litterfall+Cbk_Litterfall+Cbr_Litterfall+Cr_Litterfall

                # Populate tree-level growth of survivors

                # Find where to add tree level data
                iStart=np.where(np.isnan(tobs['ID Tree'])==True)[0]
                iFill=np.arange(iStart[0],iStart[0]+ind_surv.size,1).astype(int)

                tobs['ID Source'][iFill]=sobs['ID Source'][cnt]
                tobs['ID Plot'][iFill]=pl['ID Plot'][iPlot_t0]
                tobs['ID Visit'][iFill]=sobs['ID Visit'][cnt]
                tobs['ID Tree'][iFill]=id0_all[ia1[ind_surv]]
                tobs['ID Species'][iFill]=spc0_all[ia1[ind_surv]]
                tobs['Lat'][iFill]=sobs['Lat'][cnt]
                tobs['Lon'][iFill]=sobs['Lon'][cnt]
                tobs['Year t0'][iFill]=sobs['Year t0'][cnt]
                tobs['Year t1'][iFill]=sobs['Year t1'][cnt]
                tobs['Delta t'][iFill]=sobs['Delta t'][cnt]
                tobs['Plot Type'][iFill]=sobs['Plot Type'][cnt]
                tobs['Ecozone BC L1'][iFill]=sobs['Ecozone BC L1'][cnt]
                tobs['Ecozone BC L2'][iFill]=sobs['Ecozone BC L2'][cnt]

                tobs['Stand Age t0'][iFill]=sobs['Age t0'][cnt]
                tobs['Stand Age t1'][iFill]=sobs['Age t1'][cnt]
                tobs['Stand Spc1 ID'][iFill]=sobs['Spc1 ID'][cnt]
                tobs['Stand Spc1 %BA'][iFill]=sobs['Spc1 %BA'][cnt]
                tobs['Stand Spc1 %N'][iFill]=sobs['Spc1 %N'][cnt]
                tobs['Stand N L t0'][iFill]=sobs['N L t0'][cnt]
                tobs['Stand N L t1'][iFill]=sobs['N L t1'][cnt]
                #tobs['Stand Dqm L t0'][iFill]=sobs['Dqm L t0'][cnt]
                tobs['Stand Cag L t0'][iFill]=sobs['Cag L t0'][cnt]
                tobs['Stand Cag L t1'][iFill]=sobs['Cag L t1'][cnt]
                #tobs['Management'][iFill]=sobs['Management'][cnt]

                # tobs['Stand DA Insect t0'][iFill]=sobs['DA Insect t0'][cnt]
                # tobs['Stand DA Disease t0'][iFill]=sobs['DA Disease t0'][cnt]
                # tobs['Stand DA Fire t0'][iFill]=sobs['DA Fire t0'][cnt]
                # tobs['Stand DA Animal t0'][iFill]=sobs['DA Animal t0'][cnt]
                # tobs['Stand DA Weather t0'][iFill]=sobs['DA Weather t0'][cnt]
                # tobs['Stand DA Silviculture t0'][iFill]=sobs['DA Silviculture t0'][cnt]
                # tobs['Stand DA Mechanical t0'][iFill]=sobs['DA Mechanical t0'][cnt]
                # tobs['Stand DA SnowAndIce t0'][iFill]=sobs['DA SnowAndIce t0'][cnt]
                # tobs['Stand DA Unknown t0'][iFill]=sobs['DA Unknown t0'][cnt]

                # tobs['Stand DA Insect_t1'][iFill]=sobs['DA Insect_t1'][cnt]
                # tobs['Stand DA Disease_t1'][iFill]=sobs['DA Disease_t1'][cnt]
                # tobs['Stand DA Fire_t1'][iFill]=sobs['DA Fire_t1'][cnt]
                # tobs['Stand DA AnimaL t1'][iFill]=sobs['DA AnimaL t1'][cnt]
                # tobs['Stand DA Weather_t1'][iFill]=sobs['DA Weather_t1'][cnt]
                # tobs['Stand DA Silviculture_t1'][iFill]=sobs['DA Silviculture_t1'][cnt]
                # tobs['Stand DA MechanicaL t1'][iFill]=sobs['DA MechanicaL t1'][cnt]
                # tobs['Stand DA SnowAndIce_t1'][iFill]=sobs['DA SnowAndIce_t1'][cnt]
                # tobs['Stand DA Unknown_t1'][iFill]=sobs['DA Unknown_t1'][cnt]

                tobs['AEF'][iFill]=aef0_all[ia1[ind_surv]]
                tobs['Mortality'][iFill]=0*np.ones(ind_surv.size)
                tobs['Recruitment'][iFill]=0*np.ones(ind_surv.size)
                tobs['DA t0'][iFill]=da0_all[ia1[ind_surv]]
                tobs['DA t1'][iFill]=da1_all[ib1[ind_surv]]

                tobs['DBH t0'][iFill]=dbh0_all[ia1[ind_surv]]
                tobs['DBH t1'][iFill]=dbh1_all[ib1[ind_surv]]
                Gd=(dbh1_all[ib1[ind_surv]]-dbh0_all[ia1[ind_surv]])/dt
                tobs['DBH G'][iFill]=Gd

                tobs['BA t0'][iFill]=ba0_all[ia1[ind_surv]]
                tobs['BA t1'][iFill]=ba1_all[ib1[ind_surv]]
                Gba=(ba1_all[ib1[ind_surv]]-ba0_all[ia1[ind_surv]])/dt
                tobs['BA G'][iFill]=Gba

                tobs['H t0'][iFill]=h0_all[ia1[ind_surv]]
                tobs['H t1'][iFill]=h1_all[ib1[ind_surv]]
                Gh=(h1_all[ib1[ind_surv]]-h0_all[ia1[ind_surv]])/dt
                tobs['H G'][iFill]=dH_surv

                tobs['H Obs t0'][iFill]=h0_obs_all[ia1[ind_surv]]
                tobs['H Obs t1'][iFill]=h1_obs_all[ib1[ind_surv]]

                tobs['Cag t0'][iFill]=Cag0_all[ia1[ind_surv]]
                tobs['Cag t1'][iFill]=Cag1_all[ib1[ind_surv]]
                Gag=(Cag1_all[ib1[ind_surv]]-Cag0_all[ia1[ind_surv]])/dt
                tobs['Cag G'][iFill]=Gag

                #tobs['Cbk t0'][iFill]=Cbk0_all[ia1[ind_surv]]
                #tobs['Cbk t1'][iFill]=Cbk1_all[ib1[ind_surv]]

                #tobs['Cbr t0'][iFill]=Cbr0_all[ia1[ind_surv]]
                #tobs['Cbr t1'][iFill]=Cbr1_all[ib1[ind_surv]]

                #tobs['Cf t0'][iFill]=Cf0_all[ia1[ind_surv]]
                #tobs['Cf t1'][iFill]=Cf1_all[ib1[ind_surv]]

                tobs['Csw t0'][iFill]=Csw0_all[ia1[ind_surv]]
                tobs['Csw t1'][iFill]=Csw1_all[ib1[ind_surv]]
                Gsw=(Csw1_all[ib1[ind_surv]]-Csw0_all[ia1[ind_surv]])/dt
                tobs['Csw G'][iFill]=Gsw

                # Competition indices
                for i in range(iFill.size):
                    iLarger=np.where( tobs['Cag t0'][iFill]>tobs['Cag t0'][iFill[i]] )[0]
                    tobs['Stand Cag L Larger t0'][iFill[i]]=np.sum(tobs['AEF'][iFill[iLarger]]*tobs['Cag t0'][iFill[iLarger]])*0.001
                    tobs['Stand BA L Larger t0'][iFill[i]]=np.sum(tobs['AEF'][iFill[iLarger]]*(np.pi*(tobs['DBH t0'][iFill[iLarger]]/200)**2))

                # Rank
                # DBH_ForRank=dbh0_all[ia]
                # Lx=length(DBH_ForRank)
                # rnk']=flipdim(sortrows([DBH_ForRank[1:Lx]']),1)
                # rnk[:,2]=linspace(1,0,Lx)
                # rnk2=zeros(Lx,1) rnk2(rnk[:,1])=rnk[:,2]
                # tobs['Rank'][iFill]=rnk2(ind_surv)

                # Populate stand-level biomass Litterfall due to natural mortality

                if ind_mort_v.size==0:

                    # Ensure that an absence of dead trees translates into a rate of 0 rather
                    # than NaN (which implies an absence of monitoring)
                    sobs['N Mort'][cnt]=0
                    sobs['Cag Mort'][cnt]=0
                    sobs['Csw Mort'][cnt]=0
                    sobs['Cbk Mort'][cnt]=0
                    sobs['Cbr Mort'][cnt]=0
                    sobs['Cf Mort'][cnt]=0
                    sobs['Ctot Mort'][cnt]=0
                    #sobs['Csw_M_Ins'][cnt]=0
                    #sobs['Csw_M_Dis'][cnt]=0
                    #sobs['Csw_M_Fir'][cnt]=0
                    #sobs['Csw_M_Ani'][cnt]=0
                    #sobs['Csw_M_Wea'][cnt]=0
                    #sobs['Csw_M_Sno'][cnt]=0

                else:

                    # Demographic mortality (trees yr-1)
                    sobs['N Mort'][cnt]=np.sum(aef0_all[ia1[ind_mort_v]])/dt
                    # sobs['N Mort'][cnt]=np.round(sum(aef1_all[ib1[ind_mort_v]])/dt*100)/100

                    # Demographic mortality of Recrited trees (that did not make it) (trees yr-1)
                    #sobs['N Mr'][cnt]=np.sum(aef1_all[ib2[ind_mort_r]])/dt

                    #Csw0_dying_NP=np.nansum(aef0_all[ia[ind_mort_v]]*0.001*Csw0_all[ia[ind_mort_v]])
                    #ind_NotPine=spc1_all[ib1[ind_mort_v]]~=48
                    #Csw1_dying_NotPine=np.nansum(aef1_all[ib1[ind_mort_v[ind_NotPine]]*0.001*Csw1_all[ib1[ind_mort_v[ind_NotPine]])
                    #sobs['Csw_M_NotPine'][cnt]=Csw1_dying_NotPine/dt

                    # Stemwood growth of trees that died
                    # *** Not reliable ***
                    #sobs['Csw_Gd'][cnt]=(Csw1_dying-Csw0_dying)/dt

                    # Aboveground carbon mortality (Mg C ha-1 yr-1)
                    Cag0_dying=np.nansum(aef0_all[ia1[ind_mort_v]]*0.001*Cag0_all[ia1[ind_mort_v]])
                    Cag1_dying=np.nansum(aef1_all[ib1[ind_mort_v]]*0.001*Cag1_all[ib1[ind_mort_v]])
                    sobs['Cag Mort'][cnt]=Cag1_dying/dt

                    #sobs['Cag_Gd'][cnt]=(Cag1_dying-Cag0_dying)/dt

                    # Branch carbon mortality (Mg C ha-1 yr-1)
                    Cbr0_dying=np.nansum(aef0_all[ia1[ind_mort_v]]*0.001*Cbr0_all[ia1[ind_mort_v]])
                    Cbr1_dying=np.nansum(aef1_all[ib1[ind_mort_v]]*0.001*Cbr1_all[ib1[ind_mort_v]])
                    sobs['Cbr Mort'][cnt]=Cbr1_dying/dt

                    # Bark carbon mortality (Mg C ha-1 yr-1)
                    Cbk0_dying=np.nansum(aef0_all[ia1[ind_mort_v]]*0.001*Cbk0_all[ia1[ind_mort_v]])
                    Cbk1_dying=np.nansum(aef1_all[ib1[ind_mort_v]]*0.001*Cbk1_all[ib1[ind_mort_v]])
                    sobs['Cbk Mort'][cnt]=Cbk1_dying/dt

                    # Foliage carbon mortality (Mg C ha-1 yr-1)
                    Cf0_dying=np.nansum(aef0_all[ia1[ind_mort_v]]*0.001*Cf0_all[ia1[ind_mort_v]])
                    Cf1_dying=np.nansum(aef1_all[ib1[ind_mort_v]]*0.001*Cf1_all[ib1[ind_mort_v]])
                    sobs['Cf Mort'][cnt]=Cf1_dying/dt

                    # Stemwood carbon mortality (Mg C ha-1 yr-1)
                    Csw0_dying=np.nansum(aef0_all[ia1[ind_mort_v]]*0.001*Csw0_all[ia1[ind_mort_v]])
                    Csw1_dying=np.nansum(aef1_all[ib1[ind_mort_v]]*0.001*Csw1_all[ib1[ind_mort_v]])
                    sobs['Csw Mort'][cnt]=Csw1_dying/dt

                    # Aboveground carbon mortality (Mg C ha-1 yr-1)
                    Ctot0_dying=np.nansum(aef0_all[ia1[ind_mort_v]]*0.001*Ctot0_all[ia1[ind_mort_v]])
                    Ctot1_dying=np.nansum(aef1_all[ib1[ind_mort_v]]*0.001*Ctot1_all[ib1[ind_mort_v]])
                    sobs['Ctot Mort'][cnt]=Ctot1_dying/dt

                    # Csw1_dying=np.nansum(aef1_all[ib1(ind_mort_ins))*0.001*Csw1_all[ib1(ind_mort_ins)))
                    # sobs['Csw_M_Ins'][cnt]=Csw1_dying/dt
                    # Csw1_dying=np.nansum(aef1_all[ib1(ind_mort_dis))*0.001*Csw1_all[ib1(ind_mort_dis)))
                    # sobs['Csw_M_Dis'][cnt]=Csw1_dying/dt
                    # Csw1_dying=np.nansum(aef1_all[ib1(ind_mort_fir))*0.001*Csw1_all[ib1(ind_mort_fir)))
                    # sobs['Csw_M_Fir'][cnt]=Csw1_dying/dt
                    # Csw1_dying=np.nansum(aef1_all[ib1(ind_mort_ani))*0.001*Csw1_all[ib1(ind_mort_ani)))
                    # sobs['Csw_M_Ani'][cnt]=Csw1_dying/dt
                    # Csw1_dying=np.nansum(aef1_all[ib1(ind_mort_wea))*0.001*Csw1_all[ib1(ind_mort_wea)))
                    # sobs['Csw_M_Wea'][cnt]=Csw1_dying/dt
                    # Csw1_dying=np.nansum(aef1_all[ib1(ind_mort_sno))*0.001*Csw1_all[ib1(ind_mort_sno)))
                    # sobs['Csw_M_Sno'][cnt]=Csw1_dying/dt

                # Populate tree-level natural mortality

                # Find where to add to tree-level data
                iStart=np.where(np.isnan(tobs['ID Tree'])==True)[0]
                iFillD=np.arange(iStart[0],iStart[0]+ind_mort_v.size,1).astype(int)

                tobs['ID Source'][iFillD]=sobs['ID Source'][cnt]
                tobs['ID Plot'][iFillD]=pl['ID Plot'][iPlot_t0]
                tobs['ID Visit'][iFillD]=sobs['ID Visit'][cnt]
                tobs['ID Tree'][iFillD]=id0_all[ia1[ind_mort_v]]
                tobs['ID Species'][iFillD]=spc0_all[ia1[ind_mort_v]]
                tobs['Lat'][iFillD]=sobs['Lat'][cnt]
                tobs['Lon'][iFillD]=sobs['Lon'][cnt]
                #tobs['Elevation_Given'][iFillD]=sobs['Elevation_Given'][cnt]
                #tobs['Aspect_Given'][iFillD]=sobs['Aspect_Given'][cnt]
                #tobs['Slope_Given'][iFillD]=sobs['Slope_Given'][cnt]

                tobs['Year t0'][iFillD]=sobs['Year t0'][cnt]
                tobs['Year t1'][iFillD]=sobs['Year t1'][cnt]
                tobs['Delta t'][iFillD]=sobs['Delta t'][cnt]
                tobs['Plot Type'][iFillD]=sobs['Plot Type'][cnt]
                tobs['Ecozone BC L1'][iFillD]=sobs['Ecozone BC L1'][cnt]
                tobs['Ecozone BC L2'][iFillD]=sobs['Ecozone BC L2'][cnt]

                tobs['Stand Age t0'][iFillD]=sobs['Age t0'][cnt]
                tobs['Stand Age t1'][iFillD]=sobs['Age t1'][cnt]
                tobs['Stand Spc1 ID'][iFillD]=sobs['Spc1 ID'][cnt]
                tobs['Stand Spc1 %BA'][iFillD]=sobs['Spc1 %BA'][cnt]
                tobs['Stand Spc1 %N'][iFillD]=sobs['Spc1 %N'][cnt]
                tobs['Stand N L t0'][iFillD]=sobs['N L t0'][cnt]
                tobs['Stand N L t1'][iFillD]=sobs['N L t1'][cnt]
                #tobs['Stand Dqm L t0'][iFillD]=sobs['Dqm L t0'][cnt]
                #tobs['Stand Csw L t0'][iFillD]=sobs['Csw L t0'][cnt]
                #tobs['Stand Csw L t1'][iFillD]=sobs['Csw L t1'][cnt]
                tobs['Stand Cag L t0'][iFillD]=sobs['Cag L t0'][cnt]
                tobs['Stand Cag L t1'][iFillD]=sobs['Cag L t1'][cnt]
                tobs['Stand Management'][iFillD]=sobs['Management'][cnt]

                # tobs['DA_Insect t0'][iFillD]=sobs['DA_Insect t0'][cnt]
                # tobs['DA_Disease t0'][iFillD]=sobs['DA_Disease t0'][cnt]
                # tobs['DA_Fire t0'][iFillD]=sobs['DA_Fire t0'][cnt]
                # tobs['DA_Animal t0'][iFillD]=sobs['DA_Animal t0'][cnt]
                # tobs['DA_Weather t0'][iFillD]=sobs['DA_Weather t0'][cnt]
                # tobs['DA_Silviculture t0'][iFillD]=sobs['DA_Silviculture t0'][cnt]
                # tobs['DA_Mechanical t0'][iFillD]=sobs['DA_Mechanical t0'][cnt]
                # tobs['DA_SnowAndIce t0'][iFillD]=sobs['DA_SnowAndIce t0'][cnt]
                # tobs['DA_Unknown t0'][iFillD]=sobs['DA_Unknown t0'][cnt]

                # tobs['DA_Insect_t1'][iFillD]=sobs['DA_Insect_t1'][cnt]
                # tobs['DA_Disease_t1'][iFillD]=sobs['DA_Disease_t1'][cnt]
                # tobs['DA_Fire_t1'][iFillD]=sobs['DA_Fire_t1'][cnt]
                # tobs['DA_AnimaL t1'][iFillD]=sobs['DA_AnimaL t1'][cnt]
                # tobs['DA_Weather_t1'][iFillD]=sobs['DA_Weather_t1'][cnt]
                # tobs['DA_Silviculture_t1'][iFillD]=sobs['DA_Silviculture_t1'][cnt]
                # tobs['DA_MechanicaL t1'][iFillD]=sobs['DA_MechanicaL t1'][cnt]
                # tobs['DA_SnowAndIce_t1'][iFillD]=sobs['DA_SnowAndIce_t1'][cnt]
                # tobs['DA_Unknown_t1'][iFillD]=sobs['DA_Unknown_t1'][cnt]

                tobs['AEF'][iFillD]=aef0_all[ia1[ind_mort_v]]
                tobs['Mortality'][iFillD]=np.ones(ind_mort_v.size)
                #tobs['MortSampYN'][iFillD]=np.ones(ind_mort_v.size)
                tobs['Recruitment'][iFillD]=0*np.ones(ind_mort_v.size)
                tobs['DA t0'][iFillD]=da0_all[ia1[ind_mort_v]]
                tobs['DA t1'][iFillD]=da1_all[ib1[ind_mort_v]]
                #tobs['DAI t0'][iFillD]=dai0_all[ia1[ind_mort_v]]
                #tobs['DAI_t1'][iFillD]=dai1_all[ib1[ind_mort_v]]
                #tobs['DAP t0'][iFillD]=dap0_all[ia1[ind_mort_v]]
                #tobs['DAP_t1'][iFillD]=dap1_all[ib1[ind_mort_v]]
                #tobs['DA_IDW t0'][iFillD]=da_idw0_all[ia1[ind_mort_v]]
                #tobs['DA_IDW_t1'][iFillD]=da_idw1_all[ib1[ind_mort_v]]

                tobs['DBH t0'][iFillD]=dbh0_all[ia1[ind_mort_v]]
                tobs['DBH t1'][iFillD]=dbh1_all[ib1[ind_mort_v]]
                tobs['BA t0'][iFillD]=ba0_all[ia1[ind_mort_v]]
                tobs['BA t1'][iFillD]=ba1_all[ib1[ind_mort_v]]
                tobs['H t0'][iFillD]=h0_all[ia1[ind_mort_v]]
                tobs['H t1'][iFillD]=h1_all[ib1[ind_mort_v]]
                tobs['H Obs t0'][iFillD]=h0_all[ia1[ind_mort_v]]
                tobs['H Obs t1'][iFillD]=h1_all[ib1[ind_mort_v]]
                tobs['Csw t0'][iFillD]=Csw0_all[ia1[ind_mort_v]]
                tobs['Csw t1'][iFillD]=Csw1_all[ib1[ind_mort_v]]
                tobs['Cag t0'][iFillD]=Cag0_all[ia1[ind_mort_v]]
                tobs['Cag t1'][iFillD]=Cag1_all[ib1[ind_mort_v]]

                # Rank (rnk variable created above)
                #tobs['Rank'][iFillD]=np.round(rnk2[ind_mort_v],decimals=3)

                # Competition indices
                for i in range(iFillD.size):
                    iLarger=np.where(tobs['Cag t0'][iFill]>tobs['Cag t0'][iFillD[i]])[0]
                    tobs['Stand Cag L Larger t0'][iFillD[i]]=np.sum(tobs['AEF'][iFill[iLarger]]*tobs['Cag t0'][iFill[iLarger]])*0.001
                    tobs['Stand BA L Larger t0'][iFillD[i]]=np.sum(tobs['AEF'][iFill[iLarger]]*(np.pi*(tobs['DBH t0'][iFill[iLarger]]/200)**2))

                # Populate stand-level Recruitment (ingrowth)

                if (ib2.size==0) | (ind_rec.size==0):

                    # Ensure that an absence of Recruitment translates into a rate of 0
                    # rather than NaN (which implies an absence of monitoring)
                    sobs['Csw G Recr'][cnt]=0
                    sobs['Cag G Recr'][cnt]=0
                    sobs['Ctot G Recr'][cnt]=0
                    #sobs['Ctot T Recr'][cnt]=0
                    sobs['N Recr'][cnt]=0

                else:

                    # Demographic Recruitment (trees yr-1)
                    sobs['N Recr'][cnt]=np.sum(aef1_all[ib2[ind_rec]])/dt

                    # Aboveground biomass growth of Recruitment
                    Cag_r=np.nansum(aef1_all[ib2[ind_rec]]*0.001*Cag1_all[ib2[ind_rec]])
                    sobs['Cag G Recr'][cnt]=Cag_r/dt

                    # Bark biomass growth of Recruitment
                    Cbk_r=np.nansum(aef1_all[ib2[ind_rec]]*0.001*Cbk1_all[ib2[ind_rec]])
                    sobs['Cbk G Recr'][cnt]=Cbk_r/dt

                    # Branch biomass growth of Recruitment
                    Cbr_r=np.nansum(aef1_all[ib2[ind_rec]]*0.001*Cbr1_all[ib2[ind_rec]])
                    sobs['Cbr G Recr'][cnt]=Cbr_r/dt

                    # Foliage biomass growth of Recruitment
                    Cf_r=np.nansum(aef1_all[ib2[ind_rec]]*0.001*Cf1_all[ib2[ind_rec]])
                    sobs['Cf G Recr'][cnt]=Cf_r/dt

                    # Stemwood biomass growth of Recruitment
                    Csw_r=np.nansum(aef1_all[ib2[ind_rec]]*0.001*Csw1_all[ib2[ind_rec]])
                    sobs['Csw G Recr'][cnt]=Csw_r/dt

                    # Total biomass growth of Recruitment
                    Ctot_r=np.nansum(aef1_all[ib2[ind_rec]]*0.001*Ctot1_all[ib2[ind_rec]])
                    sobs['Ctot G Recr'][cnt]=Ctot_r/dt

                    # Calculate litterfall Litterfall of Recrits
                    # *** Averaging cannot be done because biomass components at t0 are
                    # unknown. Apply a correction factor that assumes litterfall was 25
                    # percent lower at t0 relative to t1. ***

                    Cf_Litterfall_t1=aef1_all[ib2[ind_rec]]*0.001*Cf_Litterfall1_all[ib2[ind_rec]]
                    Cf_Litterfall_t0=0.75*Cf_Litterfall_t1
                    Cf_Litterfall_Ave=(Cf_Litterfall_t0+Cf_Litterfall_t1)/2
                    Cf_Litterfall=np.nansum(Cf_Litterfall_Ave)

                    Cbk_Litterfall_t1=aef1_all[ib2[ind_rec]]*0.001*Cbk_Litterfall1_all[ib2[ind_rec]]
                    Cbk_Litterfall_t0=0.75*Cbk_Litterfall_t1
                    Cbk_Litterfall_Ave=(Cbk_Litterfall_t0+Cbk_Litterfall_t1)/2
                    Cbk_Litterfall=np.nansum(Cbk_Litterfall_Ave)

                    Cbr_Litterfall_t1=aef1_all[ib2[ind_rec]]*0.001*Cbr_Litterfall1_all[ib2[ind_rec]]
                    Cbr_Litterfall_t0=0.75*Cbr_Litterfall_t1
                    Cbr_Litterfall_Ave=(Cbr_Litterfall_t0+Cbr_Litterfall_t1)/2
                    Cbr_Litterfall=np.nansum(Cbr_Litterfall_Ave)

                    Cr_Litterfall_t1=aef1_all[ib2[ind_rec]]*0.001*Cr_Litterfall1_all[ib2[ind_rec]]
                    Cr_Litterfall_t0=0.75*Cr_Litterfall_t1
                    Cr_Litterfall_Ave=(Cr_Litterfall_t0+Cr_Litterfall_t1)/2
                    Cr_Litterfall=np.nansum(Cr_Litterfall_Ave)

                    # Total Litterfall of Recrits
                    Ctot_Litterfall=Cf_Litterfall+Cbk_Litterfall+Cbr_Litterfall+Cr_Litterfall
                    sobs['Ctot Litf'][cnt]=sobs['Ctot Litf'][cnt]+Ctot_Litterfall

                # Populate tree-level Recruitment

                # Find where to add to tree-level data
                iStart=np.where(np.isnan(tobs['ID Tree'])==True)[0]
                iFill=np.arange(iStart[0],iStart[0]+ind_rec.size,1).astype(int)

                tobs['ID Source'][iFill]=sobs['ID Source'][cnt]
                tobs['ID Plot'][iFill]=pl['ID Plot'][iPlot_t0]
                tobs['ID Visit'][iFill]=sobs['ID Visit'][cnt]
                tobs['ID Tree'][iFill]=id1_all[ib2[ind_rec]]
                tobs['ID Species'][iFill]=spc1_all[ib2[ind_rec]]
                tobs['Lat'][iFill]=sobs['Lat'][cnt]
                tobs['Lon'][iFill]=sobs['Lon'][cnt]
                #tobs['Elevation_Given'][iFill]=sobs['Elevation_Given'][cnt]
                #tobs['Aspect_Given'][iFill]=sobs['Aspect_Given'][cnt]
                #tobs['Slope_Given'][iFill]=sobs['Slope_Given'][cnt]
                tobs['Year t0'][iFill]=sobs['Year t0'][cnt]
                tobs['Year t1'][iFill]=sobs['Year t1'][cnt]
                tobs['Delta t'][iFill]=sobs['Delta t'][cnt]
                tobs['Plot Type'][iFill]=sobs['Plot Type'][cnt]
                tobs['Ecozone BC L1'][iFill]=sobs['Ecozone BC L1'][cnt]
                tobs['Ecozone BC L2'][iFill]=sobs['Ecozone BC L2'][cnt]

                tobs['Stand Age t0'][iFill]=sobs['Age t0'][cnt]
                tobs['Stand Age t1'][iFill]=sobs['Age t1'][cnt]
                tobs['Stand Spc1 ID'][iFill]=sobs['Spc1 ID'][cnt]
                tobs['Stand Spc1 %BA'][iFill]=sobs['Spc1 %BA'][cnt]
                tobs['Stand Spc1 %N'][iFill]=sobs['Spc1 %N'][cnt]
                tobs['Stand N L t0'][iFill]=sobs['N L t0'][cnt]
                tobs['Stand N L t1'][iFill]=sobs['N L t1'][cnt]
                #tobs['Stand Dqm L t0'][iFill]=sobs['Dqm L t0'][cnt]
                #tobs['Stand Csw L t0'][iFill]=sobs['Csw L t0'][cnt]
                #tobs['Stand Csw L t1'][iFill]=sobs['Csw L t1'][cnt]
                tobs['Stand Cag L t0'][iFill]=sobs['Cag L t0'][cnt]
                tobs['Stand Cag L t1'][iFill]=sobs['Cag L t1'][cnt]
                tobs['Stand Management'][iFill]=sobs['Management'][cnt]

                # tobs['Stand DA Insect t0'][iFill]=sobs['DA Insect t0'][cnt]
                # tobs['Stand DA Disease t0'][iFill]=sobs['DA Disease t0'][cnt]
                # tobs['Stand DA Fire t0'][iFill]=sobs['DA Fire t0'][cnt]
                # tobs['Stand DA Animal t0'][iFill]=sobs['DA Animal t0'][cnt]
                # tobs['Stand DA Weather t0'][iFill]=sobs['DA Weather t0'][cnt]
                # tobs['Stand DA Silviculture t0'][iFill]=sobs['DA Silviculture t0'][cnt]
                # tobs['Stand DA Mechanical t0'][iFill]=sobs['DA Mechanical t0'][cnt]
                # tobs['Stand DA SnowAndIce t0'][iFill]=sobs['DA SnowAndIce t0'][cnt]
                # tobs['Stand DA Unknown t0'][iFill]=sobs['DA Unknown t0'][cnt]

                # tobs['Stand DA Insect_t1'][iFill]=sobs['DA Insect t1'][cnt]
                # tobs['Stand DA Disease t1'][iFill]=sobs['DA Disease t1'][cnt]
                # tobs['Stand DA Fire t1'][iFill]=sobs['DA Fire t1'][cnt]
                # tobs['Stand DA AnimaL t1'][iFill]=sobs['DA AnimaL t1'][cnt]
                # tobs['Stand DA Weather t1'][iFill]=sobs['DA Weather t1'][cnt]
                # tobs['Stand DA Silviculture t1'][iFill]=sobs['DA Silviculture t1'][cnt]
                # tobs['Stand DA MechanicaL t1'][iFill]=sobs['DA MechanicaL t1'][cnt]
                # tobs['Stand DA SnowAndIce t1'][iFill]=sobs['DA SnowAndIce t1'][cnt]
                # tobs['Stand DA Unknown t1'][iFill]=sobs['DA Unknown t1'][cnt]

                tobs['AEF'][iFill]=aef1_all[ib2[ind_rec]]
                tobs['Mortality'][iFill]=np.zeros(ind_rec.size)
                #tobs['MortSampYN'][iFill]=np.ones(ind_rec.size]
                tobs['Recruitment'][iFill]=np.ones(ind_rec.size)
                tobs['DA t0'][iFill]=np.nan
                tobs['DA t1'][iFill]=da1_all[ib2[ind_rec]]
                tobs['DBH t0'][iFill]=np.nan
                tobs['DBH t1'][iFill]=dbh1_all[ib2[ind_rec]]
                #tobs['DAI t0'][iFill]=np.nan
                #tobs['DAI t1'][iFill]=dai1_all[ib2[ind_rec]]
                #tobs['DAP t0'][iFill]=np.nan
                #tobs['DAP_t1'][iFill]=dap1_all[ib2[ind_rec]]
                #tobs['DA_IDW t0'][iFill]=np.nan
                #tobs['DA_IDW_t1'][iFill]=da_idw1_all[ib2[ind_rec]]
                tobs['BA t0'][iFill]=np.nan
                tobs['BA t1'][iFill]=ba1_all[ib2[ind_rec]]
                tobs['H t0'][iFill]=np.nan
                tobs['H t1'][iFill]=h1_all[ib2[ind_rec]]
                tobs['H Obs t0'][iFill]=np.nan
                tobs['H Obs t1'][iFill]=h1_all[ib2[ind_rec]]
                tobs['Csw t0'][iFill]=np.nan
                tobs['Csw t1'][iFill]=Csw1_all[ib2[ind_rec]]
                tobs['Cag t0'][iFill]=np.nan
                tobs['Cag t1'][iFill]=Cag1_all[ib2[ind_rec]]

                # Populate stand-level mortality due to harvesting (Removals)

                if ind_harv.size==0:

                    # Ensure that an absence of harvested trees translates into a rate of 0 rather
                    # than NaN (which implies an absence of monitoring)
                    sobs['N Harv'][cnt]=0
                    sobs['Csw Harv'][cnt]=0
                    sobs['Cag Harv'][cnt]=0
                    sobs['Ctot Harv'][cnt]=0

                else:

                    # Demographic harvesting (trees yr-1).
                    sobs['N Harv'][cnt]=np.nansum(aef0_all[ia1[ind_harv]])/dt

                    # Aboveground carbon harvesting (Mg C ha-1 yr-1)
                    Cag_h=np.nansum(aef0_all[ia1[ind_harv]]*0.001*Cag0_all[ia1[ind_harv]])
                    sobs['Cag Harv'][cnt]=Cag_h/dt

                    # Stemwood carbon harvesting (Mg C ha-1 yr-1)
                    Csw_h=np.nansum(aef0_all[ia1[ind_harv]]*0.001*Csw0_all[ia1[ind_harv]])
                    sobs['Csw Harv'][cnt]=Csw_h/dt

                    # Total carbon harvesting (Mg C ha-1 yr-1)
                    Ctot_h=np.nansum(aef0_all[ia1[ind_harv]]*0.001*Ctot0_all[ia1[ind_harv]])
                    sobs['Ctot Harv'][cnt]=Ctot_h/dt

            # Update stand level counter
            cnt=cnt+1

    # Save
    d1={'tobs':tobs,'sobs':sobs}
    gu.opickle(meta['Paths']['DB'] + '\\Processed\\L2\\L2_' + jurL[iJur] + '.pkl',d1)

