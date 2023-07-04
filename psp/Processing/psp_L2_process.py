
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
import geopandas as gpd
import gc as garc
import time
import warnings
from shapely.geometry import Point, Polygon
from shapely import geometry
from fcgadgets.macgyver import utilities_general as gu
from fcgadgets.macgyver import utilities_gis as gis
import fcgadgets.macgyver.utilities_query_gdb as qgdb
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

#%% QA flags

QA_Flag1=[]

#%% Import raster data

zRef=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Admin\BC_Land_Mask.tif')

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

    # Plant functional type
    tl['ID PFT']=np.zeros(tl['ID Species'].size,dtype=int)
    for iU in range(uSpc.size):
        ind_all0=np.where(tl['ID Species']==uSpc[iU])[0]
        ind_all1=np.where(meta['Allo B']['ID']==uSpc[iU])[0]
        tl['ID PFT'][ind_all0]=meta['Allo B']['PFT'][ind_all1]

    # Volume based on DBH+H (Nigh et al. 2016)
    # *** Using volume provided with dataset ***

    # Biomass based on DBH+H (Ung et al. 2008 Lambert et al. 2005)
    beta={}
    tissueL=['bk','br','f','sw']
    for tissue in tissueL:
        beta[tissue]=np.nan*np.ones( (tl['ID Species'].size,3) )

    for iU in range(uSpc.size):
        ind_all0=np.where(tl['ID Species']==uSpc[iU])[0]
        ind_all1=np.where(meta['Allo B']['ID']==uSpc[iU])[0]

        if np.isnan(meta['Allo B']['GapFillCode'][ind_all1])==False:
            # Missing equaiton, using an alternative species
            ind_all1=np.where(meta['Allo B']['ID']==meta['Allo B']['GapFillCode'][ind_all1])[0]

        for tissue in tissueL:
            for iP in range(3):
                beta[tissue][ind_all0,iP]=meta['Allo B'][tissue + str(iP+1)][ind_all1]

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
    tl['Cr'][ind]=0.222*tl['Cag'][ind]

    # Total biomass carbon  (kg C tree-1)
    tl['Ctot']=tl['Cag']+tl['Cr']

    # Biomass Litterfall due to litterfall
    # *** Aboveground components taken from Table 3 of Kurz et al (2009)

    FracLitterfallBk=np.zeros(tl['Csw'].size)
    FracLitterfallBr=np.zeros(tl['Csw'].size)
    FracLitterfallF=np.zeros(tl['Csw'].size)
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
    slvL=['ID Source','ID Plot','ID Visit','Lat','Lon','X','Y','Elev','Aspect','Slope','Area',
         'Jourisiction','Ecozone CA L1','Ecozone BC L1','Ecozone BC L2','Plot Type','Num Plots',
         'Utilization Level','Management','SI','Year t0','Year t1','Month t0','Month t1',
         'Delta t','Age VRI t0','Age VRI t1','Age Mean t0','Age Med t0','Age Min t0','Age Max t0','Spc1 L ID t0','Spc1 L %BA t0','Spc1 L %N t0','Spc2 L ID t0','Spc2 L %BA t0',
         'Spc2 L %N t0','Spc3 L ID t0','Spc3 L %BA t0','Spc3 L %N t0','Spc4 L ID t0','Spc4 L %BA t0','Spc4 L %N t0',
         'Conifer L %BA t0','Conifer L %N t0','Deciduous L %BA t0','Deciduous L %N t0',
         'N L t0','N L t1','N D t0','N D t1','N Recr','N Mort','N L Lost','N Net','N L PL t0',
         'N Recr (%)','N Mort (%)','N L Lost (%)','N Mort+Lost (%)','N Net (%)',
         'N D Lost','N D PL t0','N L PL t1','N D PL t1',
         'Dam L t0','Dam L t1','Dam G Surv',
         'Dqm L t0','Dqm L t1',
         'BA L t0','BA L t1','BA G Surv','BA Mort',
         'H L t0','H L t1','H G Surv',
         'H L max t0','H L max t1',
         'Vws L t0','Vws L t1','Vws G Surv','Vws G Recr','Vws Mort','Vws L Lost','Vws Mort+Lost','Vws Net','Vws Harv',
         'Vntwb L t0','Vntwb D t0',
         'Ctot L t0','Ctot L t1',
         'Ctot L Resid t0',
         'Ctot G Surv','Ctot G Recr','Ctot Mort','Ctot L Lost','Ctot Mort+Lost','Ctot Net',
         'Ctot Mort+Lost None','Ctot Mort+Lost Insect','Ctot Mort+Lost Disease','Ctot Mort+Lost Plant Competition','Ctot Mort+Lost Animal Browsing','Ctot Mort+Lost Fire','Ctot Mort+Lost Frost',
         'Ctot Mort+Lost Drought','Ctot Mort+Lost Snow and Ice','Ctot Mort+Lost Wind','Ctot Mort+Lost Silviculture','Ctot Mort+Lost Other','Ctot Mort+Lost Unknown','Ctot Mort+Lost Harv',
         'Ctot D t0','Ctot D t1','Ctot D Fallen t0','Ctot D Lost','Ctot Litf',
         'Ctot L PL t0','Ctot D PL t0','Ctot L PL t1','Ctot D PL t1',
         'Ctot L FD t0','Ctot D FD t0',
         'Ctot L BL t0','Ctot D BL t0',
         'Ctot L SW t0','Ctot D SW t0',
         'Ctot L SE t0','Ctot D SE t0',
         'Ctot L Decid t0','Ctot D Decid t0',
         'Ctot G Surv HGT10',
         'Cag L t0','Cag L t1','Cag D t0','Cag D t1','Cag L Fallen t0','Cag L Fallen t1','Cag D Fallen t0','Cag D Fallen t1','Cag G Surv','Cag G Recr','Cag Mort','Cag Harv',
         'Cbk L t0','Cbk L t1','Cbk G Surv','Cbk G Recr','Cbk Mort',
         'Cbr L t0','Cbr L t1','Cbr G Surv','Cbr G Recr','Cbr Mort',
         'Cf L t0','Cf L t1','Cf G Surv','Cf G Recr','Cf Mort',
         'Cr L t0','Cr L t1','Cr G Surv','Cr G Recr','Cr Mort',
         'Csw L t0','Csw L t1','Csw G Surv','Csw G Recr','Csw Mort','Csw L Lost','Csw Mort+Lost','Csw Net',
         'Csw Harv','Csw D t0','Csw D t1','Csw125 L t0','Csw125 L t1',
         'Csw Indiv Med t0','Csw Indiv Mean t0','Csw Indiv Max t0',
         'Cdw t0']

    sobs={}
    for v in slvL:
        sobs[v]=np.nan*np.ones(pl['ID Plot'].size)

    # Initialize fields in the Level 2 tree-level structure TOBS
    tlvL=['ID Source','ID Plot','ID Visit','ID Tree','ID Species','ID PFT','Lat','Lon','X','Y','Elev','Aspect','Slope',
         'Ecozone BC L1','Ecozone BC L2','Ecozone CA L1','Plot Type','Stand Management',
         'Year t0','Year t1','Delta t','Stand Spc1 ID','Stand Spc1 %BA','Stand Spc1 %N',
         'Stand Age VRI t0','Stand Age VRI t1','Stand Age Med t0',
         'Stand N L t0','Stand N L t1','Stand Cag L t0','Stand Cag L t1','Stand Cag L Larger t0','Stand BA L Larger t0',
         'Stand DA Insect t0','Stand DA Instect t1','Stand DA Disease t0','Stand DA Disease t1',
         'Stand DA Fire t0','Stand DA Fire t1','Stand DA Animal t0','Stand DA Animal t1','Stand DA Weather t0','Stand DA WEather t1',
         'Stand DA Silviculture t0','Stand DA Silviculture t1','Stand DA Mechanical t0','Stand DA Mechanical t1',
         'Stand DA SnowAndIce t0','Stand DA SnowAndIce t1','Stand DA Unknown t0','Stand DA Unknown t1',
         'Vital Status t0','Vital Status t1','Stature t0','Stature t1','AEF','Age t0','Resid','DA t0','DA t1','DAI t0','DAI t1','DAP t0','DAP t1','DA IDW t0','DA IDW t1',
         'Mortality','Recruitment','DBH t0','DBH t1','DBH G','BA t0','BA t1','BA G',
         'H t0','H t1','H G','H Obs t0','H Obs t1','Vws t0','Vws t1',
         'Cag t0','Cag t1','Cag G',
         'Cbk t0','Cbk t1','Cbk G',
         'Cbr t0','Cbr t1','Cbr G',
         'Cf t0','Cf t1','Cf G',
         'Cr t0','Cr t1','Cr G',
         'Csw t0','Csw t1','Csw G']

    tobs={}
    for v in tlvL:
        tobs[v]=np.nan*np.ones(tl['ID Plot'].size)

    # Initialize counter
    cnt=0

    # Loop through plots
    for iP in range(uP.size):

        #if uP[iP]==1100981:
        #    break

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
            if (iPlot_t1.size==0):# | (pl['Year'][iPlot_t1]==pl['Year'][iPlot_t0]):
                flg_PSP=0
            else:
                flg_PSP=1

            # Isolate trees from plot I and inventory year iVisit
            ind_all0=np.where( (tl['ID Plot']==uP[iP]) & (tl['ID Visit']==uVisit[iVisit]) & (tl['Flag WithinPlot']==1) )[0]

            if flg_PSP==1:

                # Isolate trees from plot I and inventory year iVisit+1
                ind_all1=np.where( (tl['ID Plot']==uP[iP]) & (tl['ID Visit']==uVisit[iVisit]+1) & (tl['Flag WithinPlot']==1) )[0]

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
            sobs['ID Source'][cnt]=meta['LUT']['Source'][jurL[iJur]]
            sobs['ID Plot'][cnt]=pl['ID Plot'][iPlot_t0]
            sobs['ID Visit'][cnt]=uVisit[iVisit]
            sobs['Lat'][cnt]=np.round(pl['Lat'][iPlot_t0],decimals=4)
            sobs['Lon'][cnt]=np.round(pl['Lon'][iPlot_t0],decimals=4)
            sobs['X'][cnt]=np.round(pl['X'][iPlot_t0],decimals=2)
            sobs['Y'][cnt]=np.round(pl['Y'][iPlot_t0],decimals=2)
            sobs['Year t0'][cnt]=pl['Year'][iPlot_t0]
            sobs['Month t0'][cnt]=pl['Month'][iPlot_t0]
            sobs['Age VRI t0'][cnt]=pl['Age VRI'][iPlot_t0]
            sobs['Plot Type'][cnt]=pl['Plot Type'][iPlot_t0]
            sobs['Num Plots'][cnt]=pl['Num Plots'][iPlot_t0]
            sobs['Ecozone BC L1'][cnt]=pl['Ecozone BC L1'][iPlot_t0]
            sobs['Ecozone BC L2'][cnt]=pl['Ecozone BC L2'][iPlot_t0]
            #sobs['Area'][cnt]=pl['Area'][iPlot_t0]
            #sobs['Utilization Level'][cnt]=np.nan
            #sobs['Management'][cnt]=pl['Management'][iPlot_t0]
            #sobs['StandOrigin'][cnt]=pl['StandOrigin'][iPlot_t0]
            #sobs['SI'][cnt]=pl['SI_SourceDB'][iPlot_t0]

            if flg_PSP==1:
                sobs['Year t1'][cnt]=pl['Year'][iPlot_t1]
                sobs['Month t1'][cnt]=pl['Month'][iPlot_t1]
                dt=sobs['Year t1'][cnt]-sobs['Year t0'][cnt]
                sobs['Delta t'][cnt]=dt
                sobs['Age VRI t1'][cnt]=pl['Age VRI'][iPlot_t1]
            else:
                dt=0

            # Extract trees

            # All trees at t0
            id0_all=tl['ID Tree'][ind_all0]
            resid0_all=tl['Resid'][ind_all0]
            age0_all=tl['Age'][ind_all0]
            dbh0_all=tl['DBH'][ind_all0]
            ba0_all=np.pi*(tl['DBH'][ind_all0]/2)**2
            h0_all=tl['H'][ind_all0]
            h0_obs_all=tl['H Obs'][ind_all0]
            vws0_all=tl['Vws'][ind_all0]
            vntwb0_all=tl['Vntwb'][ind_all0]
            pft0_all=tl['ID PFT'][ind_all0]
            vital_status0_all=tl['Vital Status'][ind_all0]
            stature0_all=tl['Stature'][ind_all0]
            aef0_all=tl['AEF'][ind_all0]
            spc0_all=tl['ID Species'][ind_all0]
            da0_all=tl['ID DA1'][ind_all0]
            flag_ibm0_all=tl['Flag IBM'][ind_all0]
            flag_idw0_all=tl['Flag IDW'][ind_all0]
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
                age1_all=tl['Age'][ind_all1]
                dbh1_all=tl['DBH'][ind_all1]
                ba1_all=np.pi*(tl['DBH'][ind_all1]/2)**2
                h1_all=tl['H'][ind_all1]
                h1_obs_all=tl['H Obs'][ind_all1]
                vws1_all=tl['Vws'][ind_all1]
                vntwb1_all=tl['Vntwb'][ind_all1]
                pft1_all=tl['ID PFT'][ind_all1]
                vital_status1_all=tl['Vital Status'][ind_all1]
                stature1_all=tl['Stature'][ind_all1]
                aef1_all=tl['AEF'][ind_all1]
                spc1_all=tl['ID Species'][ind_all1]
                da1_all=tl['ID DA1'][ind_all1]
                flag_ibm1_all=tl['Flag IBM'][ind_all1]
                flag_idw1_all=tl['Flag IDW'][ind_all1]
                Cag1_all=tl['Cag'][ind_all1]
                Cbk1_all=tl['Cbk'][ind_all1]
                Cbr1_all=tl['Cbr'][ind_all1]
                Cf1_all=tl['Cf'][ind_all1]
                Cr1_all=tl['Cr'][ind_all1]
                Csw1_all=tl['Csw'][ind_all1]
                Csw125_1_all=tl['Csw'][ind_all1]
                Ctot1_all=tl['Ctot'][ind_all1]
                Cf_Litterfall1_all=tl['Cf Litterfall'][ind_all1]
                Cbk_Litterfall1_all=tl['Cbk Litterfall'][ind_all1]
                Cbr_Litterfall1_all=tl['Cbr Litterfall'][ind_all1]
                Cr_Litterfall1_all=tl['Cr Litterfall'][ind_all1]

            # Define indices to trees with a specific status and fate

            # Index to trees that were alive in t0 (including those that died)
            ind_live0=np.where(vital_status0_all==meta['LUT']['Vital Status']['Live'])[0]
            ind_live_resid0=np.where( (vital_status0_all==meta['LUT']['Vital Status']['Live']) & (resid0_all==1) )[0]
            ind_live_fallen0=np.where( (vital_status0_all==meta['LUT']['Vital Status']['Live']) & (stature0_all==meta['LUT']['Stature']['Fallen']) )[0]
            ind_live_conifer0=np.where( (vital_status0_all==meta['LUT']['Vital Status']['Live']) & (pft0_all==meta['LUT']['PFT']['Coniferous']) )[0]
            ind_live_decid0=np.where( (vital_status0_all==meta['LUT']['Vital Status']['Live']) & (pft0_all==meta['LUT']['PFT']['Deciduous']) )[0]
            ind_live_pine0=np.where( (vital_status0_all==meta['LUT']['Vital Status']['Live']) & (spc0_all==meta['LUT']['Species']['PL']) )[0]
            ind_live_bl0=np.where( (vital_status0_all==meta['LUT']['Vital Status']['Live']) & (spc0_all==meta['LUT']['Species']['BL']) )[0]
            ind_live_fd0=np.where( (vital_status0_all==meta['LUT']['Vital Status']['Live']) & (spc0_all==meta['LUT']['Species']['FD']) )[0]
            ind_live_sw0=np.where( (vital_status0_all==meta['LUT']['Vital Status']['Live']) & (spc0_all==meta['LUT']['Species']['SW']) )[0]
            ind_live_se0=np.where( (vital_status0_all==meta['LUT']['Vital Status']['Live']) & (spc0_all==meta['LUT']['Species']['SE']) )[0]

            # Index to trees that were dead at t0
            ind_dead0=np.where(vital_status0_all==meta['LUT']['Vital Status']['Dead'])[0]
            ind_dead_fallen0=np.where( (vital_status0_all==meta['LUT']['Vital Status']['Dead']) & (stature0_all==meta['LUT']['Stature']['Fallen']) )[0]
            ind_dead_decid0=np.where( (vital_status0_all==meta['LUT']['Vital Status']['Dead']) & (pft0_all==meta['LUT']['PFT']['Deciduous']) )[0]
            ind_dead_pine0=np.where( (vital_status0_all==meta['LUT']['Vital Status']['Dead']) & (spc0_all==meta['LUT']['Species']['PL']) )[0]
            ind_dead_bl0=np.where( (vital_status0_all==meta['LUT']['Vital Status']['Dead']) & (spc0_all==meta['LUT']['Species']['BL']) )[0]
            ind_dead_fd0=np.where( (vital_status0_all==meta['LUT']['Vital Status']['Dead']) & (spc0_all==meta['LUT']['Species']['FD']) )[0]
            ind_dead_sw0=np.where( (vital_status0_all==meta['LUT']['Vital Status']['Dead']) & (spc0_all==meta['LUT']['Species']['SW']) )[0]
            ind_dead_se0=np.where( (vital_status0_all==meta['LUT']['Vital Status']['Dead']) & (spc0_all==meta['LUT']['Species']['SE']) )[0]

            if flg_PSP==1:

                # Index to trees that were alive in t1 (including those that were
                # Recruted)
                ind_live1=np.where(vital_status1_all==meta['LUT']['Vital Status']['Live'])[0]
                ind_live_fallen1=np.where( (vital_status1_all==meta['LUT']['Vital Status']['Live']) & (stature1_all==meta['LUT']['Stature']['Fallen']) )[0]
                ind_live_conifer1=np.where( (vital_status1_all==meta['LUT']['Vital Status']['Live']) & (stature1_all==meta['LUT']['PFT']['Coniferous']) )[0]
                ind_live_decid1=np.where( (vital_status1_all==meta['LUT']['Vital Status']['Live']) & (stature1_all==meta['LUT']['PFT']['Deciduous']) )[0]
                ind_live_pine1=np.where( (vital_status1_all==meta['LUT']['Vital Status']['Live']) & (spc1_all==meta['LUT']['Species']['PL']) )[0]

                # Trees that go uncounted at t1
                ind_live_lost1=np.where( (vital_status0_all==meta['LUT']['Vital Status']['Live']) & (np.isin(id0_all,id1_all)==False) )[0]
                ind_dead_lost1=np.where( (vital_status0_all==meta['LUT']['Vital Status']['Dead']) & (np.isin(id0_all,id1_all)==False) )[0]

                # Index to trees that were dead at t1
                ind_dead1=np.where(vital_status1_all==meta['LUT']['Vital Status']['Dead'])[0]
                ind_dead_fallen1=np.where( (vital_status1_all==meta['LUT']['Vital Status']['Dead']) & (stature1_all==meta['LUT']['Stature']['Fallen']) )[0]
                ind_dead_pine1=np.where( (vital_status1_all==meta['LUT']['Vital Status']['Dead']) & (spc1_all==meta['LUT']['Species']['PL']) )[0]

                # Index to survivors
                ind_surv=np.where( (vital_status0_all[ia1]==meta['LUT']['Vital Status']['Live']) & (vital_status1_all[ib1]==meta['LUT']['Vital Status']['Live']) )[0]

                ind_surv_HGT10=np.where( (vital_status0_all[ia1]==meta['LUT']['Vital Status']['Live']) & (vital_status1_all[ib1]==meta['LUT']['Vital Status']['Live']) & (h0_all[ia1]>10) )[0]

                # Index to trees that were counted as dead and then counted as alive
                ind_reborn=np.where( (vital_status0_all[ia1]==meta['LUT']['Vital Status']['Dead']) & (vital_status1_all[ib1]==meta['LUT']['Vital Status']['Live']) )[0]

                # Index to trees that died during interval
                ind_mort=np.where( (vital_status0_all[ia1]==meta['LUT']['Vital Status']['Live']) & (vital_status1_all[ib1]==meta['LUT']['Vital Status']['Dead']) )[0]

                # Index to trees that died during interval
                ind_mort_da=[]
                for iDA in range(meta['LUT Tables']['Damage Agents']['ID'].size):
                    ind_mort_da.append(np.where( (vital_status0_all[ia1]==meta['LUT']['Vital Status']['Live']) & (vital_status1_all[ib1]==meta['LUT']['Vital Status']['Dead']) & (da0_all[ia1]==meta['LUT Tables']['Damage Agents']['ID'][iDA]) |
                                                 (vital_status0_all[ia1]==meta['LUT']['Vital Status']['Live']) & (vital_status1_all[ib1]==meta['LUT']['Vital Status']['Dead']) & (da1_all[ib1]==meta['LUT Tables']['Damage Agents']['ID'][iDA]) )[0])

                ind_lost_da=[]
                for iDA in range(meta['LUT Tables']['Damage Agents']['ID'].size):
                    ind_lost_da.append(np.where( (vital_status0_all==meta['LUT']['Vital Status']['Live']) &
                                                (np.isin(id0_all,id1_all)==False) &
                                                (da0_all==meta['LUT Tables']['Damage Agents']['ID'][iDA]) )[0])

                # Mortality due to IBM
                #ind_mort_ibm=np.where( (vital_status0_all[ia1]==meta['LUT']['Vital Status']['Live']) & (vital_status1_all[ib1]==meta['LUT']['Vital Status']['Dead']) & (flag_ibm0_all[ia1]==0) &
                #                             (vital_status0_all[ia1]==meta['LUT']['Vital Status']['Live']) & (vital_status1_all[ib1]==meta['LUT']['Vital Status']['Dead']) & (flag_ibm1_all[ia1]==1) )[0]

                # Index to Recruted trees that died during interval
                ind_mort_r=np.where(vital_status1_all[ib2]==meta['LUT']['Vital Status']['Dead'])[0]

                # Index to Recruted trees that remain alive
                ind_rec=np.where(vital_status1_all[ib2]==meta['LUT']['Vital Status']['Live'])[0]

                # Index to harvested trees
                ind_harv=np.where( (vital_status0_all[ia1]==meta['LUT']['Vital Status']['Live']) & (vital_status1_all[ib1]==meta['LUT']['Vital Status']['Removed']) )[0]

            # Cleaning of trees that die
            # *** Sometimes they don't measure the DBH and H of the trees that die
            # or are harvested. They might have DBH <=0 or NaN. In Quebec, they are
            # -999.These were set to the DBH and H of the previous measurement. ***
            if flg_PSP==1:

                # Mortality
                iBadD=np.where( (dbh0_all[ia1[ind_mort]]>0) & (dbh1_all[ib1[ind_mort]]<=0) | (dbh0_all[ia1[ind_mort]]>0) & (np.isnan(dbh1_all[ib1[ind_mort]])==1) )[0]
                if iBadD.size!=0:
                    aef1_all[ib1[ind_mort[iBadD]]]=aef0_all[ia1[ind_mort[iBadD]]]
                    age1_all[ib1[ind_mort[iBadD]]]=age0_all[ia1[ind_mort[iBadD]]]
                    dbh1_all[ib1[ind_mort[iBadD]]]=dbh0_all[ia1[ind_mort[iBadD]]]
                    ba1_all[ib1[ind_mort[iBadD]]]=ba0_all[ia1[ind_mort[iBadD]]]
                    h1_all[ib1[ind_mort[iBadD]]]=h0_all[ia1[ind_mort[iBadD]]]
                    h1_obs_all[ib1[ind_mort[iBadD]]]=h0_obs_all[ia1[ind_mort[iBadD]]]
                    Csw1_all[ib1[ind_mort[iBadD]]]=Csw0_all[ia1[ind_mort[iBadD]]]
                    Csw125_1_all[ib1[ind_mort[iBadD]]]=Csw125_0_all[ia1[ind_mort[iBadD]]]
                    Cag1_all[ib1[ind_mort[iBadD]]]=Cag0_all[ia1[ind_mort[iBadD]]]
                    Ctot1_all[ib1[ind_mort[iBadD]]]=Ctot0_all[ia1[ind_mort[iBadD]]]

                # Harvesting
                iBadD=np.where( (dbh0_all[ia1[ind_harv]]>0) & (dbh1_all[ib1[ind_harv]]<=0) | (dbh0_all[ia1[ind_harv]]>0) & (np.isnan(dbh1_all[ib1[ind_harv]])==1) )[0]
                if iBadD.size!=0:
                    aef1_all[ib1[ind_harv[iBadD]]]=aef0_all[ia1[ind_harv[iBadD]]]
                    age1_all[ib1[ind_harv[iBadD]]]=age0_all[ia1[ind_harv[iBadD]]]
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
            uSpc0=np.unique(spc0_all[ind_live0])

            # Total basal area of stand
            BA_tot=np.nansum((np.pi*(dbh0_all[ind_live0]/200)**2)*aef0_all[ind_live0])
            N_tot=np.sum(aef0_all[ind_live0])

            SpcComp=np.zeros((uSpc0.size,3))
            SpcComp[:,0]=uSpc0
            for iSpc in range(uSpc0.size):
                indS=np.where(spc0_all[ind_live0]==uSpc0[iSpc])[0]

                # Composition by basal area
                BA_spc=np.nansum((np.pi*(dbh0_all[ind_live0[indS]]/200)**2)*aef0_all[ind_live0[indS]])
                SpcComp[iSpc,1]=BA_spc/BA_tot

                # Composition by tree density
                N_spc=np.sum(aef0_all[ind_live0[indS]])
                SpcComp[iSpc,2]=N_spc/N_tot

            ord=np.argsort(SpcComp[:,1])
            SpcComp=np.flip(SpcComp[ord,:],axis=0)

            for iSpc in range( np.minimum(4,uSpc0.size) ):
                sobs['Spc' + str(iSpc+1) + ' L ID t0'][cnt]=SpcComp[iSpc,0]
                sobs['Spc' + str(iSpc+1) + ' L %BA t0'][cnt]=np.round(SpcComp[iSpc,1]*100,decimals=0)
                sobs['Spc' + str(iSpc+1) + ' L %N t0'][cnt]=np.round(SpcComp[iSpc,2]*100,decimals=0)

            # Plant functional type
            indPFT=np.where(pft0_all[ind_live0]==meta['LUT']['PFT']['Coniferous'])[0]

            BA_pft=np.nansum((np.pi*(dbh0_all[indPFT]/200)**2)*aef0_all[indPFT])
            sobs['Conifer L %BA t0'][cnt]=np.round(BA_pft/BA_tot*100,decimals=0)
            sobs['Deciduous L %BA t0'][cnt]=np.round(100-sobs['Conifer L %BA t0'][cnt],decimals=0)

            N_pft=np.sum(aef0_all[ind_live0[indPFT]])
            sobs['Conifer L %N t0'][cnt]=np.round(N_pft/N_tot*100,decimals=0)
            sobs['Deciduous L %N t0'][cnt]=np.round(100-sobs['Conifer L %N t0'][cnt],decimals=0)

            # Age
            sobs['Age Mean t0'][cnt]=np.round(np.nanmean(age0_all[ind_live0]),decimals=0)
            sobs['Age Med t0'][cnt]=np.round(np.nanmedian(age0_all[ind_live0]),decimals=0)
            if ind_live0.size>0:
                sobs['Age Min t0'][cnt]=np.round(np.nanmin(age0_all[ind_live0]),decimals=0)
                sobs['Age Max t0'][cnt]=np.round(np.nanmax(age0_all[ind_live0]),decimals=0)

            #------------------------------------------------------------------
            # Populate state variables
            #------------------------------------------------------------------

            # Stand density (stems ha-1)
            sobs['N L t0'][cnt]=np.round(np.nansum(aef0_all[ind_live0])/sobs['Num Plots'][cnt],decimals=0)
            sobs['N D t0'][cnt]=np.round(np.nansum(aef0_all[ind_dead0])/sobs['Num Plots'][cnt],decimals=0)
            sobs['N L PL t0'][cnt]=np.round(np.nansum(aef0_all[ind_live_pine0])/sobs['Num Plots'][cnt],decimals=0)
            sobs['N D PL t0'][cnt]=np.round(np.nansum(aef0_all[ind_dead_pine0])/sobs['Num Plots'][cnt],decimals=0)

            #sobs['N_R'][cnt]=np.round(sum(aef1_all[ib2[ind_rec]])/dt*100)/100

            # Diameter, arithmetic mean (cm)
            sobs['Dam L t0'][cnt]=np.round(np.nansum(dbh0_all[ind_live0]*aef0_all[ind_live0])/np.nansum(aef0_all[ind_live0]),decimals=1)

            # Diameter, quadratic mean (cm)
            sobs['Dqm L t0'][cnt]=np.round(np.sqrt(np.nanmean(dbh0_all[ind_live0]**2)),decimals=1)

            # Basal area (m2 ha-1)
            sobs['BA L t0'][cnt]=np.round(np.nansum((np.pi*(dbh0_all[ind_live0]/200)**2)*aef0_all[ind_live0])/sobs['Num Plots'][cnt],decimals=1)

            # Height (m)
            sobs['H L t0'][cnt]=np.round(np.nansum(h0_all[ind_live0]*aef0_all[ind_live0])/np.nansum(aef0_all[ind_live0]),decimals=1)

            # Volume whole stem (m3/ha)
            sobs['Vws L t0'][cnt]=np.round(np.nansum(aef0_all[ind_live0]*vws0_all[ind_live0])/sobs['Num Plots'][cnt],decimals=0)

            # Net merch volume (m3/ha)
            sobs['Vntwb L t0'][cnt]=np.round(np.nansum(aef0_all[ind_live0]*vntwb0_all[ind_live0])/sobs['Num Plots'][cnt],decimals=0)
            sobs['Vntwb D t0'][cnt]=np.round(np.nansum(aef0_all[ind_dead0]*vntwb0_all[ind_dead0])/sobs['Num Plots'][cnt],decimals=0)

            # Aboveground carbon, live (Mg C ha-1)
            sobs['Cag L t0'][cnt]=np.round(np.nansum(aef0_all[ind_live0]*0.001*Cag0_all[ind_live0])/sobs['Num Plots'][cnt],decimals=0)
            sobs['Cag D t0'][cnt]=np.round(np.nansum(aef0_all[ind_dead0]*0.001*Cag0_all[ind_dead0])/sobs['Num Plots'][cnt],decimals=0)
            sobs['Cag L Fallen t0'][cnt]=np.round(np.nansum(aef0_all[ind_live_fallen0]*0.001*Cag0_all[ind_live_fallen0])/sobs['Num Plots'][cnt],decimals=0)
            sobs['Cag D Fallen t0'][cnt]=np.round(np.nansum(aef0_all[ind_dead_fallen0]*0.001*Cag0_all[ind_dead_fallen0])/sobs['Num Plots'][cnt],decimals=0)

            # Bark carbon, live (Mg C ha-1)
            sobs['Cbk L t0'][cnt]=np.round(np.nansum(aef0_all[ind_live0]*0.001*Cbk0_all[ind_live0])/sobs['Num Plots'][cnt],decimals=0)

            # Branch carbon, live (Mg C ha-1)
            sobs['Cbr L t0'][cnt]=np.round(np.nansum(aef0_all[ind_live0]*0.001*Cbr0_all[ind_live0])/sobs['Num Plots'][cnt],decimals=0)

            # Foliage carbon, live (Mg C ha-1)
            sobs['Cf L t0'][cnt]=np.round(np.nansum(aef0_all[ind_live0]*0.001*Cf0_all[ind_live0])/sobs['Num Plots'][cnt],decimals=0)

            # Root carbon, live (Mg C ha-1)
            sobs['Cr L t0'][cnt]=np.round(np.nansum(aef0_all[ind_live0]*0.001*Cr0_all[ind_live0])/sobs['Num Plots'][cnt],decimals=0)

            # Stemwood carbon, live (Mg C ha-1)
            sobs['Csw L t0'][cnt]=np.round(np.nansum(aef0_all[ind_live0]*0.001*Csw0_all[ind_live0])/sobs['Num Plots'][cnt],decimals=0)
            sobs['Csw D t0'][cnt]=np.round(np.nansum(aef0_all[ind_dead0]*0.001*Csw0_all[ind_dead0])/sobs['Num Plots'][cnt],decimals=0)

            sobs['Csw Indiv Med t0'][cnt]=np.round(np.nanmedian(Csw0_all[ind_live0]),decimals=0)
            sobs['Csw Indiv Mean t0'][cnt]=np.round(np.nanmean(Csw0_all[ind_live0]),decimals=0)
            try:
                sobs['Csw Indiv Max t0'][cnt]=np.round(np.nanmax(Csw0_all[ind_live0]),decimals=0)
            except:
                pass

            # Stemwood carbon, live (Mg C ha-1)
            sobs['Csw125 L t0'][cnt]=np.round(np.nansum(aef0_all[ind_live0]*0.001*Csw125_0_all[ind_live0])/sobs['Num Plots'][cnt],decimals=0)

            # Dead wood (Mg C ha-1)
            sobs['Cdw t0'][cnt]=np.round(np.nansum(aef0_all[ind_dead0]*0.001*(Csw0_all[ind_dead0]+Cbk0_all[ind_dead0]+Cbr0_all[ind_dead0]))/sobs['Num Plots'][cnt],decimals=0)

            # Total carbon, live (Mg C ha-1)
            sobs['Ctot L t0'][cnt]=np.round(np.nansum(aef0_all[ind_live0]*0.001*Ctot0_all[ind_live0])/sobs['Num Plots'][cnt],decimals=0)
            sobs['Ctot D t0'][cnt]=np.round(np.nansum(aef0_all[ind_dead0]*0.001*Ctot0_all[ind_dead0])/sobs['Num Plots'][cnt],decimals=0)
            sobs['Ctot D Fallen t0'][cnt]=np.round(np.nansum(aef0_all[ind_dead_fallen0]*0.001*Ctot0_all[ind_dead_fallen0])/sobs['Num Plots'][cnt],decimals=0)
            sobs['Ctot L Resid t0'][cnt]=np.round(np.nansum(aef0_all[ind_live_resid0]*0.001*Ctot0_all[ind_live_resid0])/sobs['Num Plots'][cnt],decimals=0)

            sobs['Ctot L PL t0'][cnt]=np.round(np.nansum(aef0_all[ind_live_pine0]*0.001*Ctot0_all[ind_live_pine0])/sobs['Num Plots'][cnt],decimals=0)
            sobs['Ctot D PL t0'][cnt]=np.round(np.nansum(aef0_all[ind_dead_pine0]*0.001*Ctot0_all[ind_dead_pine0])/sobs['Num Plots'][cnt],decimals=0)

            sobs['Ctot L BL t0'][cnt]=np.round(np.nansum(aef0_all[ind_live_bl0]*0.001*Ctot0_all[ind_live_bl0])/sobs['Num Plots'][cnt],decimals=0)
            sobs['Ctot D BL t0'][cnt]=np.round(np.nansum(aef0_all[ind_dead_bl0]*0.001*Ctot0_all[ind_dead_bl0])/sobs['Num Plots'][cnt],decimals=0)

            sobs['Ctot L Decid t0'][cnt]=np.round(np.nansum(aef0_all[ind_live_decid0]*0.001*Ctot0_all[ind_live_decid0])/sobs['Num Plots'][cnt],decimals=0)
            sobs['Ctot D Decid t0'][cnt]=np.round(np.nansum(aef0_all[ind_dead_decid0]*0.001*Ctot0_all[ind_dead_decid0])/sobs['Num Plots'][cnt],decimals=0)

            sobs['Ctot L FD t0'][cnt]=np.round(np.nansum(aef0_all[ind_live_fd0]*0.001*Ctot0_all[ind_live_fd0])/sobs['Num Plots'][cnt],decimals=0)
            sobs['Ctot D FD t0'][cnt]=np.round(np.nansum(aef0_all[ind_dead_fd0]*0.001*Ctot0_all[ind_dead_fd0])/sobs['Num Plots'][cnt],decimals=0)

            sobs['Ctot L SW t0'][cnt]=np.round(np.nansum(aef0_all[ind_live_sw0]*0.001*Ctot0_all[ind_live_sw0])/sobs['Num Plots'][cnt],decimals=0)
            sobs['Ctot D SW t0'][cnt]=np.round(np.nansum(aef0_all[ind_dead_sw0]*0.001*Ctot0_all[ind_dead_sw0])/sobs['Num Plots'][cnt],decimals=0)

            sobs['Ctot L SE t0'][cnt]=np.round(np.nansum(aef0_all[ind_live_se0]*0.001*Ctot0_all[ind_live_se0])/sobs['Num Plots'][cnt],decimals=0)
            sobs['Ctot D SE t0'][cnt]=np.round(np.nansum(aef0_all[ind_dead_se0]*0.001*Ctot0_all[ind_dead_se0])/sobs['Num Plots'][cnt],decimals=0)

            if (flg_PSP==1) & (sobs['Plot Type'][cnt]!=meta['LUT']['Plot Type BC']['VRI']) & (dt>0):

                # Stand density (stems ha-1)
                sobs['N L t1'][cnt]=np.round(np.nansum(aef1_all[ind_live1])/sobs['Num Plots'][cnt],decimals=0)
                sobs['N D t1'][cnt]=np.round(np.nansum(aef1_all[ind_dead1])/sobs['Num Plots'][cnt],decimals=0)

                sobs['N L PL t1'][cnt]=np.round(np.nansum(aef1_all[ind_live_pine1])/sobs['Num Plots'][cnt],decimals=0)
                sobs['N D PL t1'][cnt]=np.round(np.nansum(aef1_all[ind_dead_pine1])/sobs['Num Plots'][cnt],decimals=0)

                # Demographic lost (trees yr-1)
                sobs['N L Lost'][cnt]=np.round(np.nansum(aef0_all[ind_live_lost1])/dt/sobs['Num Plots'][cnt],decimals=1)
                sobs['N D Lost'][cnt]=np.round(np.nansum(aef0_all[ind_dead_lost1])/dt/sobs['Num Plots'][cnt],decimals=1)

                # Diameter, arithmetic mean (cm)
                sobs['Dam L t1'][cnt]=np.round(np.nanmean(dbh1_all[ind_live1]),decimals=1)

                # Diameter, quadratic mean (cm)
                sobs['Dqm L t1'][cnt]=np.round(np.sqrt(np.mean(dbh1_all[ind_live1]**2)),decimals=1)

                # Basal area (m2 ha-1)
                sobs['BA L t1'][cnt]=np.round(np.nansum((np.pi*(dbh1_all[ind_live1]/200)**2)*aef1_all[ind_live1])/sobs['Num Plots'][cnt],decimals=1)

                # Height (m)
                sobs['H L t1'][cnt]=np.round(np.nanmean(h1_all[ind_live1]),decimals=1)

                # Volume whole stem (m3/ha)
                sobs['Vws L t1'][cnt]=np.round(np.nansum(aef1_all[ind_live1]*vws1_all[ind_live1])/sobs['Num Plots'][cnt],decimals=0)

                # Aboveground carbon, live (Mg C ha-1)
                sobs['Cag L t1'][cnt]=np.round(np.nansum(aef1_all[ind_live1]*0.001*Cag1_all[ind_live1])/sobs['Num Plots'][cnt],decimals=0)
                sobs['Cag D t1'][cnt]=np.round(np.nansum(aef1_all[ind_dead1]*0.001*Cag1_all[ind_dead1])/sobs['Num Plots'][cnt],decimals=0)
                sobs['Cag L Fallen t1'][cnt]=np.round(np.nansum(aef1_all[ind_live_fallen1]*0.001*Cag1_all[ind_live_fallen1])/sobs['Num Plots'][cnt],decimals=0)
                sobs['Cag D Fallen t1'][cnt]=np.round(np.nansum(aef1_all[ind_dead_fallen1]*0.001*Cag1_all[ind_dead_fallen1])/sobs['Num Plots'][cnt],decimals=0)

                # Bark carbon, live (Mg C ha-1)
                sobs['Cbk L t1'][cnt]=np.round(np.nansum(aef1_all[ind_live1]*0.001*Cbk1_all[ind_live1])/sobs['Num Plots'][cnt],decimals=0)

                # Branch carbon, live (Mg C ha-1)
                sobs['Cbr L t1'][cnt]=np.round(np.nansum(aef1_all[ind_live1]*0.001*Cbr1_all[ind_live1])/sobs['Num Plots'][cnt],decimals=0)

                # Foliage carbon, live (Mg C ha-1)
                sobs['Cf L t1'][cnt]=np.round(np.nansum(aef1_all[ind_live1]*0.001*Cf1_all[ind_live1])/sobs['Num Plots'][cnt],decimals=0)

                # Root carbon, live (Mg C ha-1)
                sobs['Cr L t1'][cnt]=np.round(np.nansum(aef1_all[ind_live1]*0.001*Cr1_all[ind_live1])/sobs['Num Plots'][cnt],decimals=0)

                # Stemwood carbon, live (Mg C ha-1)
                sobs['Csw L t1'][cnt]=np.round(np.nansum(aef1_all[ind_live1]*0.001*Csw1_all[ind_live1])/sobs['Num Plots'][cnt],decimals=0)
                sobs['Csw D t1'][cnt]=np.round(np.nansum(aef1_all[ind_dead1]*0.001*Csw1_all[ind_dead1])/sobs['Num Plots'][cnt],decimals=0)

                # Stemwood carbon, live (Mg C ha-1)
                sobs['Csw125 L t1'][cnt]=np.round(np.nansum(aef1_all[ind_live1]*0.001*Csw125_1_all[ind_live1])/sobs['Num Plots'][cnt],decimals=0)

                # Total carbon (Mg C ha-1)
                sobs['Ctot L t1'][cnt]=np.round(np.nansum(aef1_all[ind_live1]*0.001*Ctot1_all[ind_live1])/sobs['Num Plots'][cnt],decimals=0)
                sobs['Ctot D t1'][cnt]=np.round(np.nansum(aef1_all[ind_dead1]*0.001*Ctot1_all[ind_dead1])/sobs['Num Plots'][cnt],decimals=0)

                sobs['Ctot L PL t1'][cnt]=np.round(np.nansum(aef1_all[ind_live_pine1]*0.001*Ctot1_all[ind_live_pine1])/sobs['Num Plots'][cnt],decimals=0)
                sobs['Ctot D PL t1'][cnt]=np.round(np.nansum(aef1_all[ind_dead_pine1]*0.001*Ctot1_all[ind_dead_pine1])/sobs['Num Plots'][cnt],decimals=0)

                #--------------------------------------------------------------
                # Growth of survivors
                #--------------------------------------------------------------

                # Stand-level growth:

                # Diameter, arithmetic mean, growth (cm yr-1)
                Dam_t0_surv=dbh0_all[ia1[ind_surv]]
                Dam_t1_surv=dbh1_all[ib1[ind_surv]]
                dDam_surv=(Dam_t1_surv-Dam_t0_surv)/dt
                sobs['Dam G Surv'][cnt]=np.round(np.nanmean(dDam_surv),decimals=2)

                # Basal area growth (m2 ha-1)
                BA_t0_surv=(np.pi*(dbh0_all[ia1[ind_surv]]/200)**2)*aef0_all[ia1[ind_surv]]
                BA_t1_surv=(np.pi*(dbh1_all[ib1[ind_surv]]/200)**2)*aef1_all[ib1[ind_surv]]
                dBA_surv=(BA_t1_surv-BA_t0_surv)/dt
                dBA_surv_sum=np.nansum(dBA_surv)
                sobs['BA G Surv'][cnt]=np.round(np.nansum(dBA_surv)/sobs['Num Plots'][cnt],decimals=2)

                # Height growth (m yr-1)
                H_t0_surv=h0_all[ia1[ind_surv]]
                H_t1_surv=h1_all[ib1[ind_surv]]
                dH_surv=(H_t1_surv-H_t0_surv)/dt
                sobs['H G Surv'][cnt]=np.round(np.nanmean(dH_surv)/sobs['Num Plots'][cnt],decimals=2)

                # Volume growth (m3 ha-1 yr-1)
                Vws_t0_surv=aef0_all[ia1[ind_surv]]*vws0_all[ia1[ind_surv]]
                Vws_t1_surv=aef1_all[ib1[ind_surv]]*vws1_all[ib1[ind_surv]]
                dVws_surv=(Vws_t1_surv-Vws_t0_surv)/dt
                sobs['Vws G Surv'][cnt]=np.round(np.nansum(dVws_surv)/sobs['Num Plots'][cnt],decimals=2)

                # Aboveground carbon growth (Mg C ha-1 yr-1)
                Cag_t0_surv=aef0_all[ia1[ind_surv]]*0.001*Cag0_all[ia1[ind_surv]]
                Cag_t1_surv=aef1_all[ib1[ind_surv]]*0.001*Cag1_all[ib1[ind_surv]]
                dCag_surv=(Cag_t1_surv-Cag_t0_surv)/dt
                sobs['Cag G Surv'][cnt]=np.round(np.nansum(dCag_surv)/sobs['Num Plots'][cnt],decimals=2)

                # Bark carbon growth (Mg C ha-1 yr-1)
                Cbk_t0_surv=aef0_all[ia1[ind_surv]]*0.001*Cbk0_all[ia1[ind_surv]]
                Cbk_t1_surv=aef1_all[ib1[ind_surv]]*0.001*Cbk1_all[ib1[ind_surv]]
                dCbk_surv=(Cbk_t1_surv-Cbk_t0_surv)/dt
                sobs['Cbk G Surv'][cnt]=np.round(np.nansum(dCbk_surv)/sobs['Num Plots'][cnt],decimals=2)

                # Branch carbon growth (Mg C ha-1 yr-1)
                Cbr_t0_surv=aef0_all[ia1[ind_surv]]*0.001*Cbr0_all[ia1[ind_surv]]
                Cbr_t1_surv=aef1_all[ib1[ind_surv]]*0.001*Cbr1_all[ib1[ind_surv]]
                dCbr_surv=(Cbr_t1_surv-Cbr_t0_surv)/dt
                sobs['Cbr G Surv'][cnt]=np.round(np.nansum(dCbr_surv)/sobs['Num Plots'][cnt],decimals=2)

                # Foliage carbon growth (Mg C ha-1 yr-1)
                Cf_t0_surv=aef0_all[ia1[ind_surv]]*0.001*Cf0_all[ia1[ind_surv]]
                Cf_t1_surv=aef1_all[ib1[ind_surv]]*0.001*Cf1_all[ib1[ind_surv]]
                dCf_surv=(Cf_t1_surv-Cf_t0_surv)/dt
                sobs['Cf G Surv'][cnt]=np.round(np.nansum(dCf_surv)/sobs['Num Plots'][cnt],decimals=2)

                # Stemwood carbon growth (Mg C ha-1 yr-1)
                Csw_t0_surv=aef0_all[ia1[ind_surv]]*0.001*Csw0_all[ia1[ind_surv]]
                Csw_t1_surv=aef1_all[ib1[ind_surv]]*0.001*Csw1_all[ib1[ind_surv]]
                dCsw_surv=(Csw_t1_surv-Csw_t0_surv)/dt
                sobs['Csw G Surv'][cnt]=np.round(np.nansum(dCsw_surv)/sobs['Num Plots'][cnt],decimals=2)

                # Total carbon growth (Mg C ha-1 yr-1)
                Ctot_t0_surv=aef0_all[ia1[ind_surv]]*0.001*Ctot0_all[ia1[ind_surv]]
                Ctot_t1_surv=aef1_all[ib1[ind_surv]]*0.001*Ctot1_all[ib1[ind_surv]]
                dCtot_surv=(Ctot_t1_surv-Ctot_t0_surv)/dt
                sobs['Ctot G Surv'][cnt]=np.round(np.nansum(dCtot_surv)/sobs['Num Plots'][cnt],decimals=2)

                # Total carbon growth (Mg C ha-1 yr-1)
                Ctot_t0_surv=aef0_all[ia1[ind_surv_HGT10]]*0.001*Ctot0_all[ia1[ind_surv_HGT10]]
                Ctot_t1_surv=aef1_all[ib1[ind_surv_HGT10]]*0.001*Ctot1_all[ib1[ind_surv_HGT10]]
                dCtot_surv=(Ctot_t1_surv-Ctot_t0_surv)/dt
                sobs['Ctot G Surv HGT10'][cnt]=np.round(np.nansum(dCtot_surv)/sobs['Num Plots'][cnt],decimals=2)

                # Litterfall of survivors:

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
                sobs['Ctot Litf'][cnt]=np.round( (Cf_Litterfall+Cbk_Litterfall+Cbr_Litterfall+Cr_Litterfall)/sobs['Num Plots'][cnt] ,decimals=2)

                # Tree-level growth of survivors:

                # Find where to add tree level data
                iStart=np.where(np.isnan(tobs['ID Tree'])==True)[0]
                iFill=np.arange(iStart[0],iStart[0]+ind_surv.size,1).astype(int)

                tobs['ID Source'][iFill]=sobs['ID Source'][cnt]
                tobs['ID Plot'][iFill]=pl['ID Plot'][iPlot_t0]
                tobs['ID Visit'][iFill]=sobs['ID Visit'][cnt]
                tobs['ID Tree'][iFill]=id0_all[ia1[ind_surv]]
                tobs['ID Species'][iFill]=spc0_all[ia1[ind_surv]]
                tobs['ID PFT'][iFill]=pft0_all[ia1[ind_surv]]
                tobs['Lat'][iFill]=sobs['Lat'][cnt]
                tobs['Lon'][iFill]=sobs['Lon'][cnt]
                tobs['X'][iFill]=sobs['X'][cnt]
                tobs['Y'][iFill]=sobs['Y'][cnt]
                tobs['Year t0'][iFill]=sobs['Year t0'][cnt]
                tobs['Year t1'][iFill]=sobs['Year t1'][cnt]
                tobs['Delta t'][iFill]=sobs['Delta t'][cnt]
                tobs['Plot Type'][iFill]=sobs['Plot Type'][cnt]
                tobs['Ecozone BC L1'][iFill]=sobs['Ecozone BC L1'][cnt]
                tobs['Ecozone BC L2'][iFill]=sobs['Ecozone BC L2'][cnt]

                tobs['Stand Age VRI t0'][iFill]=sobs['Age VRI t0'][cnt]
                tobs['Stand Age VRI t1'][iFill]=sobs['Age VRI t1'][cnt]
                tobs['Stand Age Med t0'][iFill]=sobs['Age Med t0'][cnt]
                tobs['Stand Spc1 ID'][iFill]=sobs['Spc1 L ID t0'][cnt]
                tobs['Stand Spc1 %BA'][iFill]=sobs['Spc1 L %BA t0'][cnt]
                tobs['Stand Spc1 %N'][iFill]=sobs['Spc1 L %N t0'][cnt]
                tobs['Stand N L t0'][iFill]=sobs['N L t0'][cnt]
                tobs['Stand N L t1'][iFill]=sobs['N L t1'][cnt]
                tobs['Stand Cag L t0'][iFill]=sobs['Cag L t0'][cnt]
                tobs['Stand Cag L t1'][iFill]=sobs['Cag L t1'][cnt]

                tobs['AEF'][iFill]=aef0_all[ia1[ind_surv]]
                tobs['Vital Status t0'][iFill]=1*np.ones(ind_surv.size)
                tobs['Vital Status t1'][iFill]=1*np.ones(ind_surv.size)
                tobs['Age t0'][iFill]=age0_all[ia1[ind_surv]]
                tobs['Resid'][iFill]=resid0_all[ia1[ind_surv]]
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

                tobs['Vws t0'][iFill]=vws0_all[ia1[ind_surv]]
                tobs['Vws t1'][iFill]=vws1_all[ib1[ind_surv]]

                tobs['Cag t0'][iFill]=Cag0_all[ia1[ind_surv]]
                tobs['Cag t1'][iFill]=Cag1_all[ib1[ind_surv]]
                tobs['Cag G'][iFill]=(Cag1_all[ib1[ind_surv]]-Cag0_all[ia1[ind_surv]])/dt

                tobs['Cbk t0'][iFill]=Cbk0_all[ia1[ind_surv]]
                tobs['Cbk t1'][iFill]=Cbk1_all[ib1[ind_surv]]

                tobs['Cbr t0'][iFill]=Cbr0_all[ia1[ind_surv]]
                tobs['Cbr t1'][iFill]=Cbr1_all[ib1[ind_surv]]

                tobs['Cf t0'][iFill]=Cf0_all[ia1[ind_surv]]
                tobs['Cf t1'][iFill]=Cf1_all[ib1[ind_surv]]

                tobs['Cr t0'][iFill]=Cr0_all[ia1[ind_surv]]
                tobs['Cr t1'][iFill]=Cr1_all[ib1[ind_surv]]

                tobs['Csw t0'][iFill]=Csw0_all[ia1[ind_surv]]
                tobs['Csw t1'][iFill]=Csw1_all[ib1[ind_surv]]
                tobs['Csw G'][iFill]=(Csw1_all[ib1[ind_surv]]-Csw0_all[ia1[ind_surv]])/dt

                # Competition indices
                for i in range(iFill.size):
                    iLarger=np.where( tobs['Cag t0'][iFill]>tobs['Cag t0'][iFill[i]] )[0]
                    tobs['Stand Cag L Larger t0'][iFill[i]]=np.sum(tobs['AEF'][iFill[iLarger]]*tobs['Cag t0'][iFill[iLarger]])*0.001
                    tobs['Stand BA L Larger t0'][iFill[i]]=np.sum(tobs['AEF'][iFill[iLarger]]*(np.pi*(tobs['DBH t0'][iFill[iLarger]]/200)**2))

                #--------------------------------------------------------------
                # Mortality
                #--------------------------------------------------------------

                if ind_mort.size==0:

                    # Ensure that an absence of dead trees translates into a rate of 0 rather
                    # than NaN (which implies an absence of monitoring)
                    sobs['N Mort'][cnt]=0
                    sobs['Vws Mort'][cnt]=0
                    sobs['Cag Mort'][cnt]=0
                    sobs['Cbk Mort'][cnt]=0
                    sobs['Cbr Mort'][cnt]=0
                    sobs['Cf Mort'][cnt]=0
                    sobs['Csw Mort'][cnt]=0
                    sobs['Ctot Mort'][cnt]=0

                else:

                    # Demographic mortality (trees yr-1)
                    sobs['N Mort'][cnt]=np.round(np.nansum(aef0_all[ia1[ind_mort]])/dt/sobs['Num Plots'][cnt],decimals=1)

                    # Volume (m3 ha-1 yr-1)
                    Vws0_dying=np.nansum(aef0_all[ia1[ind_mort]]*vws0_all[ia1[ind_mort]])
                    Vws1_dying=np.nansum(aef1_all[ib1[ind_mort]]*vws1_all[ib1[ind_mort]])
                    sobs['Vws Mort'][cnt]=np.round(Vws1_dying/dt/sobs['Num Plots'][cnt],decimals=2)

                    # Aboveground carbon mortality (Mg C ha-1 yr-1)
                    Cag0_dying=np.nansum(aef0_all[ia1[ind_mort]]*0.001*Cag0_all[ia1[ind_mort]])
                    Cag1_dying=np.nansum(aef1_all[ib1[ind_mort]]*0.001*Cag1_all[ib1[ind_mort]])
                    sobs['Cag Mort'][cnt]=np.round(Cag1_dying/dt/sobs['Num Plots'][cnt],decimals=2)

                    # Bark carbon mortality (Mg C ha-1 yr-1)
                    Cbk0_dying=np.nansum(aef0_all[ia1[ind_mort]]*0.001*Cbk0_all[ia1[ind_mort]])
                    Cbk1_dying=np.nansum(aef1_all[ib1[ind_mort]]*0.001*Cbk1_all[ib1[ind_mort]])
                    sobs['Cbk Mort'][cnt]=np.round(Cbk1_dying/dt/sobs['Num Plots'][cnt],decimals=2)

                    # Branch carbon mortality (Mg C ha-1 yr-1)
                    Cbr0_dying=np.nansum(aef0_all[ia1[ind_mort]]*0.001*Cbr0_all[ia1[ind_mort]])
                    Cbr1_dying=np.nansum(aef1_all[ib1[ind_mort]]*0.001*Cbr1_all[ib1[ind_mort]])
                    sobs['Cbr Mort'][cnt]=np.round(Cbr1_dying/dt/sobs['Num Plots'][cnt],decimals=2)

                    # Foliage carbon mortality (Mg C ha-1 yr-1)
                    Cf0_dying=np.nansum(aef0_all[ia1[ind_mort]]*0.001*Cf0_all[ia1[ind_mort]])
                    Cf1_dying=np.nansum(aef1_all[ib1[ind_mort]]*0.001*Cf1_all[ib1[ind_mort]])
                    sobs['Cf Mort'][cnt]=np.round(Cf1_dying/dt/sobs['Num Plots'][cnt],decimals=2)

                    # Stemwood carbon mortality (Mg C ha-1 yr-1)
                    Csw0_dying=np.nansum(aef0_all[ia1[ind_mort]]*0.001*Csw0_all[ia1[ind_mort]])
                    Csw1_dying=np.nansum(aef1_all[ib1[ind_mort]]*0.001*Csw1_all[ib1[ind_mort]])
                    sobs['Csw Mort'][cnt]=np.round(Csw1_dying/dt/sobs['Num Plots'][cnt],decimals=2)

                    # Total carbon mortality (Mg C ha-1 yr-1)
                    Ctot0_dying=np.nansum(aef0_all[ia1[ind_mort]]*0.001*Ctot0_all[ia1[ind_mort]])
                    Ctot1_dying=np.nansum(aef1_all[ib1[ind_mort]]*0.001*Ctot1_all[ib1[ind_mort]])
                    sobs['Ctot Mort'][cnt]=np.round(Ctot1_dying/dt/sobs['Num Plots'][cnt],decimals=2)

                # Total mortality due to damage agents (Mg C ha-1 yr-1)
                for iDA in range(meta['LUT Tables']['Damage Agents']['ID'].size):

                    nam=meta['LUT Tables']['Damage Agents']['Value'][iDA]
                    ind_mort_da0=ind_mort_da[iDA]

                    if ind_mort_da0.size==0:
                        sobs['Ctot Mort+Lost ' + nam][cnt]=np.round(0.0,decimals=2)
                    else:
                        nam=meta['LUT Tables']['Damage Agents']['Value'][iDA]
                        Ctot0_dying=np.nansum(aef0_all[ia1[ind_mort_da0]]*0.001*Ctot0_all[ia1[ind_mort_da0]])
                        Ctot1_dying=np.nansum(aef1_all[ib1[ind_mort_da0]]*0.001*Ctot1_all[ib1[ind_mort_da0]])
                        sobs['Ctot Mort+Lost ' + nam][cnt]=np.round(Ctot1_dying/dt/sobs['Num Plots'][cnt],decimals=2)

                # Tree-level natural mortality

                # Find where to add to tree-level data
                iStart=np.where(np.isnan(tobs['ID Tree'])==True)[0]
                iFillD=np.arange(iStart[0],iStart[0]+ind_mort.size,1).astype(int)

                tobs['ID Source'][iFillD]=sobs['ID Source'][cnt]
                tobs['ID Plot'][iFillD]=pl['ID Plot'][iPlot_t0]
                tobs['ID Visit'][iFillD]=sobs['ID Visit'][cnt]
                tobs['ID Tree'][iFillD]=id0_all[ia1[ind_mort]]
                tobs['ID Species'][iFillD]=spc0_all[ia1[ind_mort]]
                tobs['ID PFT'][iFillD]=pft0_all[ia1[ind_mort]]
                tobs['Lat'][iFillD]=sobs['Lat'][cnt]
                tobs['Lon'][iFillD]=sobs['Lon'][cnt]
                tobs['X'][iFill]=sobs['X'][cnt]
                tobs['Y'][iFill]=sobs['Y'][cnt]

                tobs['Year t0'][iFillD]=sobs['Year t0'][cnt]
                tobs['Year t1'][iFillD]=sobs['Year t1'][cnt]
                tobs['Delta t'][iFillD]=sobs['Delta t'][cnt]
                tobs['Plot Type'][iFillD]=sobs['Plot Type'][cnt]
                tobs['Ecozone BC L1'][iFillD]=sobs['Ecozone BC L1'][cnt]
                tobs['Ecozone BC L2'][iFillD]=sobs['Ecozone BC L2'][cnt]

                tobs['Stand Age VRI t0'][iFillD]=sobs['Age VRI t0'][cnt]
                tobs['Stand Age VRI t1'][iFillD]=sobs['Age VRI t1'][cnt]
                tobs['Stand Age Med t0'][iFillD]=sobs['Age Med t0'][cnt]
                tobs['Stand Spc1 ID'][iFillD]=sobs['Spc1 L ID t0'][cnt]
                tobs['Stand Spc1 %BA'][iFillD]=sobs['Spc1 L %BA t0'][cnt]
                tobs['Stand Spc1 %N'][iFillD]=sobs['Spc1 L %N t0'][cnt]
                tobs['Stand N L t0'][iFillD]=sobs['N L t0'][cnt]
                tobs['Stand N L t1'][iFillD]=sobs['N L t1'][cnt]
                tobs['Stand Cag L t0'][iFillD]=sobs['Cag L t0'][cnt]
                tobs['Stand Cag L t1'][iFillD]=sobs['Cag L t1'][cnt]
                tobs['Stand Management'][iFillD]=sobs['Management'][cnt]

                tobs['AEF'][iFillD]=aef0_all[ia1[ind_mort]]
                tobs['Vital Status t0'][iFillD]=1*np.ones(ind_mort.size)
                tobs['Vital Status t1'][iFillD]=0*np.ones(ind_mort.size)
                tobs['Age t0'][iFillD]=age0_all[ia1[ind_mort]]
                tobs['Resid'][iFillD]=resid0_all[ia1[ind_mort]]
                tobs['Mortality'][iFillD]=np.ones(ind_mort.size)
                tobs['Recruitment'][iFillD]=0*np.ones(ind_mort.size)
                tobs['DA t0'][iFillD]=da0_all[ia1[ind_mort]]
                tobs['DA t1'][iFillD]=da1_all[ib1[ind_mort]]

                tobs['DBH t0'][iFillD]=dbh0_all[ia1[ind_mort]]
                tobs['DBH t1'][iFillD]=dbh1_all[ib1[ind_mort]]
                tobs['BA t0'][iFillD]=ba0_all[ia1[ind_mort]]
                tobs['BA t1'][iFillD]=ba1_all[ib1[ind_mort]]
                tobs['H t0'][iFillD]=h0_all[ia1[ind_mort]]
                tobs['H t1'][iFillD]=h1_all[ib1[ind_mort]]
                tobs['H Obs t0'][iFillD]=h0_all[ia1[ind_mort]]
                tobs['H Obs t1'][iFillD]=h1_all[ib1[ind_mort]]
                tobs['Vws t0'][iFillD]=vws0_all[ia1[ind_mort]]
                tobs['Vws t1'][iFillD]=vws1_all[ib1[ind_mort]]
                tobs['Cag t0'][iFillD]=Cag0_all[ia1[ind_mort]]
                tobs['Cag t1'][iFillD]=Cag1_all[ib1[ind_mort]]
                tobs['Cbk t0'][iFillD]=Cbk0_all[ia1[ind_mort]]
                tobs['Cbk t1'][iFillD]=Cbk1_all[ib1[ind_mort]]
                tobs['Cbr t0'][iFillD]=Cbr0_all[ia1[ind_mort]]
                tobs['Cbr t1'][iFillD]=Cbr1_all[ib1[ind_mort]]
                tobs['Cf t0'][iFillD]=Cf0_all[ia1[ind_mort]]
                tobs['Cf t1'][iFillD]=Cf1_all[ib1[ind_mort]]
                tobs['Cr t0'][iFillD]=Cr0_all[ia1[ind_mort]]
                tobs['Cr t1'][iFillD]=Cr1_all[ib1[ind_mort]]
                tobs['Csw t0'][iFillD]=Csw0_all[ia1[ind_mort]]
                tobs['Csw t1'][iFillD]=Csw1_all[ib1[ind_mort]]

                # Competition indices
                for i in range(iFillD.size):
                    iLarger=np.where(tobs['Cag t0'][iFill]>tobs['Cag t0'][iFillD[i]])[0]
                    tobs['Stand Cag L Larger t0'][iFillD[i]]=np.sum(tobs['AEF'][iFill[iLarger]]*tobs['Cag t0'][iFill[iLarger]])*0.001
                    tobs['Stand BA L Larger t0'][iFillD[i]]=np.sum(tobs['AEF'][iFill[iLarger]]*(np.pi*(tobs['DBH t0'][iFill[iLarger]]/200)**2))

                #--------------------------------------------------------------
                # Lost trees
                #--------------------------------------------------------------

                # Volume lost
                sobs['Vws L Lost'][cnt]=np.round(np.nansum(vws0_all[ind_live_lost1]*vws0_all[ind_live_lost1])/dt/sobs['Num Plots'][cnt],decimals=1)

                # Total carbon lost
                sobs['Ctot L Lost'][cnt]=np.round(np.nansum(aef0_all[ind_live_lost1]*0.001*Ctot0_all[ind_live_lost1])/dt/sobs['Num Plots'][cnt],decimals=1)
                sobs['Ctot D Lost'][cnt]=np.round(np.nansum(aef0_all[ind_dead_lost1]*0.001*Ctot0_all[ind_dead_lost1])/dt/sobs['Num Plots'][cnt],decimals=1)

                sobs['Csw L Lost'][cnt]=np.round(np.nansum(aef0_all[ind_live_lost1]*0.001*Csw0_all[ind_live_lost1])/dt/sobs['Num Plots'][cnt],decimals=1)

                nam=meta['LUT Tables']['Damage Agents']['Value'][iDA]
                ind_mort_da0=ind_mort_da[iDA]

                # Add lost carbon to mortality
                sobs['Ctot Mort+Lost'][cnt]=np.round(sobs['Ctot Mort'][cnt]+np.nansum(aef0_all[ind_live_lost1]*0.001*Ctot0_all[ind_live_lost1])/dt/sobs['Num Plots'][cnt],decimals=1)
                for iDA in range(meta['LUT Tables']['Damage Agents']['ID'].size):
                    ind_lost_da0=ind_lost_da[iDA]
                    nam=meta['LUT Tables']['Damage Agents']['Value'][iDA]
                    sobs['Ctot Mort+Lost ' + nam][cnt]=np.round(sobs['Ctot Mort+Lost ' + nam][cnt]+np.nansum(aef0_all[ind_lost_da0]*0.001*Ctot0_all[ind_lost_da0])/dt/sobs['Num Plots'][cnt],decimals=1)
                sobs['Csw Mort+Lost'][cnt]=np.round(sobs['Csw Mort'][cnt]+np.nansum(aef0_all[ind_live_lost1]*0.001*Csw0_all[ind_live_lost1])/dt/sobs['Num Plots'][cnt],decimals=1)
                sobs['Vws Mort+Lost'][cnt]=np.round(sobs['Vws Mort'][cnt]+np.nansum(aef0_all[ind_live_lost1]*vws0_all[ind_live_lost1])/dt/sobs['Num Plots'][cnt],decimals=1)

                # Tree-level mortality (lost trees)

                # Find where to add to tree-level data
                iStart=np.where(np.isnan(tobs['ID Tree'])==True)[0]
                iFillD=np.arange(iStart[0],iStart[0]+ind_live_lost1.size,1).astype(int)

                tobs['ID Source'][iFillD]=sobs['ID Source'][cnt]
                tobs['ID Plot'][iFillD]=pl['ID Plot'][iPlot_t0]
                tobs['ID Visit'][iFillD]=sobs['ID Visit'][cnt]
                tobs['ID Tree'][iFillD]=id0_all[ind_live_lost1]
                tobs['ID Species'][iFillD]=spc0_all[ind_live_lost1]
                tobs['ID PFT'][iFillD]=pft0_all[ind_live_lost1]
                tobs['Lat'][iFillD]=sobs['Lat'][cnt]
                tobs['Lon'][iFillD]=sobs['Lon'][cnt]
                tobs['X'][iFill]=sobs['X'][cnt]
                tobs['Y'][iFill]=sobs['Y'][cnt]

                tobs['Year t0'][iFillD]=sobs['Year t0'][cnt]
                tobs['Year t1'][iFillD]=sobs['Year t1'][cnt]
                tobs['Delta t'][iFillD]=sobs['Delta t'][cnt]
                tobs['Plot Type'][iFillD]=sobs['Plot Type'][cnt]
                tobs['Ecozone BC L1'][iFillD]=sobs['Ecozone BC L1'][cnt]
                tobs['Ecozone BC L2'][iFillD]=sobs['Ecozone BC L2'][cnt]

                tobs['Stand Age VRI t0'][iFillD]=sobs['Age VRI t0'][cnt]
                tobs['Stand Age VRI t1'][iFillD]=sobs['Age VRI t1'][cnt]
                tobs['Stand Age Med t0'][iFillD]=sobs['Age Med t0'][cnt]
                tobs['Stand Spc1 ID'][iFillD]=sobs['Spc1 L ID t0'][cnt]
                tobs['Stand Spc1 %BA'][iFillD]=sobs['Spc1 L %BA t0'][cnt]
                tobs['Stand Spc1 %N'][iFillD]=sobs['Spc1 L %N t0'][cnt]
                tobs['Stand N L t0'][iFillD]=sobs['N L t0'][cnt]
                tobs['Stand N L t1'][iFillD]=sobs['N L t1'][cnt]
                tobs['Stand Cag L t0'][iFillD]=sobs['Cag L t0'][cnt]
                tobs['Stand Cag L t1'][iFillD]=sobs['Cag L t1'][cnt]
                tobs['Stand Management'][iFillD]=sobs['Management'][cnt]

                tobs['AEF'][iFillD]=aef0_all[ind_live_lost1]
                tobs['Vital Status t0'][iFillD]=1*np.ones(ind_live_lost1.size)
                tobs['Vital Status t1'][iFillD]=0*np.ones(ind_live_lost1.size)
                tobs['Age t0'][iFillD]=age0_all[ind_live_lost1]
                tobs['Mortality'][iFillD]=np.ones(ind_live_lost1.size)
                tobs['Recruitment'][iFillD]=0*np.ones(ind_live_lost1.size)
                tobs['DA t0'][iFillD]=da0_all[ind_live_lost1]
                tobs['DA t1'][iFillD]=np.nan

                tobs['DBH t0'][iFillD]=dbh0_all[ind_live_lost1]
                tobs['DBH t1'][iFillD]=np.nan
                tobs['BA t0'][iFillD]=ba0_all[ind_live_lost1]
                tobs['BA t1'][iFillD]=np.nan
                tobs['H t0'][iFillD]=h0_all[ind_live_lost1]
                tobs['H t1'][iFillD]=np.nan
                tobs['H Obs t0'][iFillD]=h0_all[ind_live_lost1]
                tobs['H Obs t1'][iFillD]=np.nan
                tobs['Vws t0'][iFillD]=vws0_all[ind_live_lost1]
                tobs['Vws t1'][iFillD]=np.nan
                tobs['Cag t0'][iFillD]=Cag0_all[ind_live_lost1]
                tobs['Cag t1'][iFillD]=np.nan
                tobs['Cbk t0'][iFillD]=Cbk0_all[ind_live_lost1]
                tobs['Cbk t1'][iFillD]=np.nan
                tobs['Cbr t0'][iFillD]=Cbr0_all[ind_live_lost1]
                tobs['Cbr t1'][iFillD]=np.nan
                tobs['Cf t0'][iFillD]=Cf0_all[ind_live_lost1]
                tobs['Cf t1'][iFillD]=np.nan
                tobs['Cr t0'][iFillD]=Cr0_all[ind_live_lost1]
                tobs['Cr t1'][iFillD]=np.nan
                tobs['Csw t0'][iFillD]=Csw0_all[ind_live_lost1]
                tobs['Csw t1'][iFillD]=np.nan

                # Competition indices
                for i in range(iFillD.size):
                    iLarger=np.where(tobs['Cag t0'][iFill]>tobs['Cag t0'][iFillD[i]])[0]
                    tobs['Stand Cag L Larger t0'][iFillD[i]]=np.sum(tobs['AEF'][iFill[iLarger]]*tobs['Cag t0'][iFill[iLarger]])*0.001
                    tobs['Stand BA L Larger t0'][iFillD[i]]=np.sum(tobs['AEF'][iFill[iLarger]]*(np.pi*(tobs['DBH t0'][iFill[iLarger]]/200)**2))

                #------------------------------------------------------------------
                # Recruitment (ingrowth)
                #------------------------------------------------------------------

                # Stand-level recruitment:

                if (ib2.size==0) | (ind_rec.size==0):

                    # Ensure that an absence of Recruitment translates into a rate of 0
                    # rather than NaN (which implies an absence of monitoring)
                    sobs['Vws G Recr'][cnt]=0
                    sobs['Cag G Recr'][cnt]=0
                    sobs['Csw G Recr'][cnt]=0
                    sobs['Cag G Recr'][cnt]=0
                    sobs['Ctot G Recr'][cnt]=0
                    sobs['N Recr'][cnt]=0

                else:

                    # Demographic Recruitment (trees yr-1)
                    sobs['N Recr'][cnt]=np.round(np.nansum(aef1_all[ib2[ind_rec]])/dt/sobs['Num Plots'][cnt],decimals=1)

                    # Volume growth of Recruitment (m3 ha yr-1)
                    Vws_r=np.nansum(aef1_all[ib2[ind_rec]]*vws1_all[ib2[ind_rec]])
                    sobs['Vws G Recr'][cnt]=np.round(Vws_r/dt/sobs['Num Plots'][cnt],decimals=2)

                    # Aboveground biomass growth of Recruitment
                    Cag_r=np.nansum(aef1_all[ib2[ind_rec]]*0.001*Cag1_all[ib2[ind_rec]])
                    sobs['Cag G Recr'][cnt]=np.round(Cag_r/dt/sobs['Num Plots'][cnt],decimals=2)

                    # Bark biomass growth of Recruitment
                    Cbk_r=np.nansum(aef1_all[ib2[ind_rec]]*0.001*Cbk1_all[ib2[ind_rec]])
                    sobs['Cbk G Recr'][cnt]=np.round(Cbk_r/dt/sobs['Num Plots'][cnt],decimals=2)

                    # Branch biomass growth of Recruitment
                    Cbr_r=np.nansum(aef1_all[ib2[ind_rec]]*0.001*Cbr1_all[ib2[ind_rec]])
                    sobs['Cbr G Recr'][cnt]=np.round(Cbr_r/dt/sobs['Num Plots'][cnt],decimals=2)

                    # Foliage biomass growth of Recruitment
                    Cf_r=np.nansum(aef1_all[ib2[ind_rec]]*0.001*Cf1_all[ib2[ind_rec]])
                    sobs['Cf G Recr'][cnt]=np.round(Cf_r/dt/sobs['Num Plots'][cnt],decimals=2)

                    # Stemwood biomass growth of Recruitment
                    Csw_r=np.nansum(aef1_all[ib2[ind_rec]]*0.001*Csw1_all[ib2[ind_rec]])
                    sobs['Csw G Recr'][cnt]=np.round(Csw_r/dt/sobs['Num Plots'][cnt],decimals=2)

                    # Total biomass growth of Recruitment
                    Ctot_r=np.nansum(aef1_all[ib2[ind_rec]]*0.001*Ctot1_all[ib2[ind_rec]])
                    sobs['Ctot G Recr'][cnt]=np.round(Ctot_r/dt/sobs['Num Plots'][cnt],decimals=2)

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

                    # Total Litterfall of recruits
                    Ctot_Litterfall=Cf_Litterfall+Cbk_Litterfall+Cbr_Litterfall+Cr_Litterfall
                    sobs['Ctot Litf'][cnt]=np.round((sobs['Ctot Litf'][cnt]+Ctot_Litterfall)/sobs['Num Plots'][cnt],decimals=2)

                # Tree-level recruitment

                # Find where to add to tree-level data
                iStart=np.where(np.isnan(tobs['ID Tree'])==True)[0]
                iFill=np.arange(iStart[0],iStart[0]+ind_rec.size,1).astype(int)

                tobs['ID Source'][iFill]=sobs['ID Source'][cnt]
                tobs['ID Plot'][iFill]=pl['ID Plot'][iPlot_t0]
                tobs['ID Visit'][iFill]=sobs['ID Visit'][cnt]
                tobs['ID Tree'][iFill]=id1_all[ib2[ind_rec]]
                tobs['ID Species'][iFill]=spc1_all[ib2[ind_rec]]
                tobs['ID PFT'][iFill]=pft1_all[ib2[ind_rec]]
                tobs['Lat'][iFill]=sobs['Lat'][cnt]
                tobs['Lon'][iFill]=sobs['Lon'][cnt]
                tobs['X'][iFill]=sobs['X'][cnt]
                tobs['Y'][iFill]=sobs['Y'][cnt]
                tobs['Year t0'][iFill]=sobs['Year t0'][cnt]
                tobs['Year t1'][iFill]=sobs['Year t1'][cnt]
                tobs['Delta t'][iFill]=sobs['Delta t'][cnt]
                tobs['Plot Type'][iFill]=sobs['Plot Type'][cnt]
                tobs['Ecozone BC L1'][iFill]=sobs['Ecozone BC L1'][cnt]
                tobs['Ecozone BC L2'][iFill]=sobs['Ecozone BC L2'][cnt]

                tobs['Stand Age VRI t0'][iFill]=sobs['Age VRI t0'][cnt]
                tobs['Stand Age VRI t1'][iFill]=sobs['Age VRI t1'][cnt]
                tobs['Stand Age Med t0'][iFill]=sobs['Age Med t0'][cnt]
                tobs['Stand Spc1 ID'][iFill]=sobs['Spc1 L ID t0'][cnt]
                tobs['Stand Spc1 %BA'][iFill]=sobs['Spc1 L %BA t0'][cnt]
                tobs['Stand Spc1 %N'][iFill]=sobs['Spc1 L %N t0'][cnt]
                tobs['Stand N L t0'][iFill]=sobs['N L t0'][cnt]
                tobs['Stand N L t1'][iFill]=sobs['N L t1'][cnt]
                tobs['Stand Cag L t0'][iFill]=sobs['Cag L t0'][cnt]
                tobs['Stand Cag L t1'][iFill]=sobs['Cag L t1'][cnt]
                tobs['Stand Management'][iFill]=sobs['Management'][cnt]

                tobs['AEF'][iFill]=aef1_all[ib2[ind_rec]]
                tobs['Age t0'][iFill]=age1_all[ib2[ind_rec]]
                tobs['Mortality'][iFill]=np.zeros(ind_rec.size)
                tobs['Vital Status t0'][iFill]=np.ones(ind_rec.size)
                tobs['Vital Status t1'][iFill]=np.ones(ind_rec.size)
                tobs['Recruitment'][iFill]=np.ones(ind_rec.size)
                tobs['DA t0'][iFill]=np.nan
                tobs['DA t1'][iFill]=da1_all[ib2[ind_rec]]
                tobs['DBH t0'][iFill]=np.nan
                tobs['DBH t1'][iFill]=dbh1_all[ib2[ind_rec]]
                tobs['BA t0'][iFill]=np.nan
                tobs['BA t1'][iFill]=ba1_all[ib2[ind_rec]]
                tobs['H t0'][iFill]=np.nan
                tobs['H t1'][iFill]=h1_all[ib2[ind_rec]]
                tobs['H Obs t0'][iFill]=np.nan
                tobs['H Obs t1'][iFill]=h1_all[ib2[ind_rec]]
                tobs['Vws t0'][iFill]=np.nan
                tobs['Vws t1'][iFill]=vws1_all[ib2[ind_rec]]
                tobs['Cag t0'][iFill]=np.nan
                tobs['Cag t1'][iFill]=Cag1_all[ib2[ind_rec]]
                tobs['Cbk t0'][iFill]=np.nan
                tobs['Cbk t1'][iFill]=Cbk1_all[ib2[ind_rec]]
                tobs['Cbr t0'][iFill]=np.nan
                tobs['Cbr t1'][iFill]=Cbr1_all[ib2[ind_rec]]
                tobs['Cf t0'][iFill]=np.nan
                tobs['Cf t1'][iFill]=Cf1_all[ib2[ind_rec]]
                tobs['Cr t0'][iFill]=np.nan
                tobs['Cr t1'][iFill]=Cr1_all[ib2[ind_rec]]
                tobs['Csw t0'][iFill]=np.nan
                tobs['Csw t1'][iFill]=Csw1_all[ib2[ind_rec]]

                # Stand-level harvest removals

                if ind_harv.size==0:

                    # Ensure that an absence of harvested trees translates into a rate of 0 rather
                    # than NaN (which implies an absence of monitoring)
                    #sobs['N Lost'][cnt]=0
                    sobs['Vws Harv'][cnt]=0
                    sobs['Csw Harv'][cnt]=0
                    sobs['Cag Harv'][cnt]=0
                    #sobs['Ctot Harv'][cnt]=0

                else:

                    # Demographic harvesting (trees yr-1).
                    sobs['N Harv'][cnt]=np.round(np.nansum(aef0_all[ia1[ind_harv]])/dt/sobs['Num Plots'][cnt],decimals=2)

                    # Volume harvesting (m3 ha-1 yr-1)
                    Vws_h=np.nansum(aef0_all[ia1[ind_harv]]*vws0_all[ia1[ind_harv]])
                    sobs['Vws Harv'][cnt]=np.round(Vws_h/dt/sobs['Num Plots'][cnt],decimals=2)

                    # Aboveground carbon harvesting (Mg C ha-1 yr-1)
                    Cag_h=np.nansum(aef0_all[ia1[ind_harv]]*0.001*Cag0_all[ia1[ind_harv]])
                    sobs['Cag Harv'][cnt]=np.round(Cag_h/dt/sobs['Num Plots'][cnt],decimals=2)

                    # Stemwood carbon harvesting (Mg C ha-1 yr-1)
                    Csw_h=np.nansum(aef0_all[ia1[ind_harv]]*0.001*Csw0_all[ia1[ind_harv]])
                    sobs['Csw Harv'][cnt]=np.round(Csw_h/dt/sobs['Num Plots'][cnt],decimals=2)

                    # Total carbon harvesting (Mg C ha-1 yr-1)
                    Ctot_h=np.nansum(aef0_all[ia1[ind_harv]]*0.001*Ctot0_all[ia1[ind_harv]])
                    #sobs['Ctot Harv'][cnt]=np.round(Ctot_h/dt/sobs['Num Plots'][cnt],decimals=2)

            # Update stand level counter
            cnt=cnt+1

    # Summary variables
    sobs['Ctot Net']=sobs['Ctot G Surv']+sobs['Ctot G Recr']-sobs['Ctot Mort']-sobs['Ctot L Lost']
    sobs['Csw Net']=sobs['Csw G Surv']+sobs['Csw G Recr']-sobs['Csw Mort']-sobs['Csw L Lost']
    sobs['Vws Net']=sobs['Vws G Surv']+sobs['Vws G Recr']-sobs['Vws Mort']-sobs['Vws L Lost']

    sobs['N Net']=sobs['N Recr']-sobs['N Mort']-sobs['N L Lost']

    ind=np.where(sobs['N L t0']>0)[0]
    sobs['N Recr (%)'][ind]=sobs['N Recr'][ind]/sobs['N L t0'][ind]*100
    sobs['N Mort (%)'][ind]=sobs['N Mort'][ind]/sobs['N L t0'][ind]*100
    sobs['N L Lost (%)'][ind]=sobs['N L Lost'][ind]/sobs['N L t0'][ind]*100
    sobs['N Mort+Lost (%)'][ind]=(sobs['N Mort'][ind]+sobs['N L Lost'][ind])/sobs['N L t0'][ind]*100
    sobs['N Net (%)'][ind]=sobs['N Net'][ind]/sobs['N L t0'][ind]*100

    # Convert some variable data types
    vL=['ID Source','ID Plot','ID Visit','Jourisiction','Ecozone CA L1','Ecozone BC L1','Ecozone BC L2','Plot Type','Management']
    for v in vL:
        ind=np.where(np.isnan(sobs[v])==False)[0]
        sobs[v][ind]=sobs[v][ind].astype(int)

    # Save
    d1={'tobs':tobs,'sobs':sobs}
    gu.opickle(meta['Paths']['DB'] + '\\Processed\\L2\\L2_' + jurL[iJur] + '.pkl',d1)

    # Export to spreadsheet
    flg=0
    if flg==1:
        df=pd.DataFrame(sobs)
        df.to_excel(meta['Paths']['DB'] + '\\Processed\\L2\\L2_sobs_' + jurL[iJur] + '.xlsx',index=False)

#%% Plot type filter (CMI, NFI, and VRI)

# Stand level

sobs['PTF CNV']=np.zeros(sobs['ID Plot'].size)
ind=np.where( (sobs['Plot Type']==meta['LUT']['Plot Type BC']['CMI']) |
             (sobs['Plot Type']==meta['LUT']['Plot Type BC']['NFI']) |
             (sobs['Plot Type']==meta['LUT']['Plot Type BC']['VRI']) )[0]
sobs['PTF CNV'][ind]=1
np.sum(sobs['PTF CNV'])

sobs['PTF CN']=np.zeros(sobs['ID Plot'].size)
ind=np.where( (sobs['Plot Type']==meta['LUT']['Plot Type BC']['CMI']) | (sobs['Plot Type']==meta['LUT']['Plot Type BC']['NFI']) )[0]
sobs['PTF CN'][ind]=1

sobs['PTF CNY']=np.zeros(sobs['ID Plot'].size)
ind=np.where( (sobs['Plot Type']==meta['LUT']['Plot Type BC']['CMI']) |
             (sobs['Plot Type']==meta['LUT']['Plot Type BC']['NFI']) |
             (sobs['Plot Type']==meta['LUT']['Plot Type BC']['YSM']) )[0]
sobs['PTF CNY'][ind]=1

sobs['PTF YSM']=np.zeros(sobs['ID Plot'].size)
ind=np.where( (sobs['Plot Type']==meta['LUT']['Plot Type BC']['YSM']) )[0]
sobs['PTF YSM'][ind]=1

# Tree level

tobs['PTF CNV']=np.zeros(tobs['ID Plot'].size)
ind=np.where( (tobs['Plot Type']==meta['LUT']['Plot Type BC']['CMI']) |
             (tobs['Plot Type']==meta['LUT']['Plot Type BC']['NFI']) |
             (tobs['Plot Type']==meta['LUT']['Plot Type BC']['VRI']) )[0]
tobs['PTF CNV'][ind]=1
np.sum(tobs['PTF CNV'])

tobs['PTF CN']=np.zeros(tobs['ID Plot'].size)
ind=np.where( (tobs['Plot Type']==meta['LUT']['Plot Type BC']['CMI']) | (tobs['Plot Type']==meta['LUT']['Plot Type BC']['NFI']) )[0]
tobs['PTF CN'][ind]=1

#%% Harvest indicator

flg=1
if flg==1:
    # Import reference crs
    gdf_bm=gpd.read_file(r'C:\Users\rhember\Documents\Data\Basemaps\Basemaps.gdb',layer='NRC_POLITICAL_BOUNDARIES_1M_SP')

    # Convert ground plot database to geodataframe
    points=[]
    for i in range(sobs['X'].size):
        points.append( Point(sobs['X'][i],sobs['Y'][i]) )
    d={'geometry':points}
    for k in sobs.keys():
        d[k]=sobs[k]
    gdf_sobs=gpd.GeoDataFrame(d)
    gdf_sobs.crs=gdf_bm.crs

    # Import opening layer
    flg=1
    if flg==1:
        op={}
        op['Path']=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20220422\Results.gdb'
        op['Layer']='RSLT_OPENING_SVW'; # fiona.listlayers(op['Path'])
        op['crs']=gdf_bm.crs
        op['Keep Geom']='On'
        op['Select Openings']=np.array([])
        op['SBC']=np.array([])
        op['FSC']=np.array([])
        op['ROI']=[]
        op['gdf']=qgdb.Query_Openings(op,[])
        #op['gdf']['crs']=gdf_bm.crs

    # Overlay
    gdf_j=gpd.overlay(gdf_sobs,op['gdf'],how='intersection')

    # Get decimal date of harvest
    Year_Harv=np.nan*np.ones(sobs['ID Plot'].size)
    Mon_Harv=np.nan*np.ones(sobs['ID Plot'].size)
    for i in range(len(gdf_j)):
        ind0=np.where((sobs['ID Plot']==gdf_j['ID Plot'][i]))[0]
        if (gdf_j.loc[i,'DENUDATION_1_DISTURBANCE_CODE']=='L') | (gdf_j.loc[i,'DENUDATION_1_DISTURBANCE_CODE']=='S'):
            if gdf_j.loc[i,'DENUDATION_1_COMPLETION_DATE'] is not None:
                try:
                    Year_Harv[ind0]=int(gdf_j.loc[i,'DENUDATION_1_COMPLETION_DATE'][0:4])
                    Mon_Harv[ind0]=int(gdf_j.loc[i,'DENUDATION_1_COMPLETION_DATE'][5:7])
                except:
                    pass
        if (gdf_j.loc[i,'DENUDATION_2_DISTURBANCE_CODE']=='L') | (gdf_j.loc[i,'DENUDATION_2_DISTURBANCE_CODE']=='S'):
            if gdf_j.loc[i,'DENUDATION_2_COMPLETION_DATE'] is not None:
                try:
                    Year_Harv[ind0]=int(gdf_j.loc[i,'DENUDATION_2_COMPLETION_DATE'][0:4])
                    Mon_Harv[ind0]=int(gdf_j.loc[i,'DENUDATION_2_COMPLETION_DATE'][5:7])
                except:
                    pass
    DY_Harv=Year_Harv+Mon_Harv/12
    sobs['Year Harvest']=DY_Harv

    DY_Plot_t0=sobs['Year t0']+sobs['Month t0']/12
    DY_Plot_t1=sobs['Year t1']+sobs['Month t1']/12

    #sobs['Ctot Mort Harv']=np.nan*sobs['Ctot Mort']
    #sobs['Vws Mort Harv']=np.nan*sobs['Vws Mort']

    #ind=np.where(sobs['Ctot Mort']>=0)[0]
    #sobs['Ctot Mort Harv'][ind]=0

    #ind=np.where(sobs['Vws Mort']>=0)[0]
    #sobs['Vws Mort Harv'][ind]=0

    sobs['Harv Occurrence']=np.zeros(sobs['Year t0'].size)
    ind_H=np.where( (DY_Harv>DY_Plot_t0) & (DY_Harv<=DY_Plot_t1) )[0]
    sobs['Harv Occurrence'][ind_H]=1

    sobs['Ctot Mort+Lost Harv'][ind_H]=sobs['Harv Occurrence'][ind_H]*(sobs['Ctot Mort+Lost Harv'][ind_H]-np.nan_to_num(sobs['Ctot Mort+Lost Fire'][ind_H]))
    #sobs['Vws Mort+Lost Harv'][ind_H]=sobs['Harv Occurrence'][ind_H]*(sobs['Vws Mort+Lost'][ind_H]-np.nan_to_num(sobs['Vws Mort+Lost Fire'][ind_H]))

#%% Forest health indicator

ind=gis.GetGridIndexToPoints(zRef,sobs['X'],sobs['Y'])
tv=np.arange(1990,2022,1)

vL=['IBM','IBB','IBD','IBS','IDW']
for v in vL:
    z1=np.zeros( (tv.size,ind[0].size) ,dtype='int8')
    for iT in range(tv.size):
        z=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\PEST_INFESTATION_POLY_' + v + '_SeverityClass_' + str(tv[iT]) + '.tif')
        z1[iT,:]=z['Data'][ind].copy()

    #sobs[v + ' Severity Max']=np.zeros(sobs['Year t0'].size)
    sobs[v + ' Occurrence']=np.zeros(sobs['Year t0'].size)
    for i in range(sobs['Year t0'].size):
        ind1=np.where( (z1[:,i]>=2) & (tv>=sobs['Year t0'][i]) & (tv<sobs['Year t1'][i]) )[0]
        if ind1.size==0:
            continue
        #sobs[v + ' Severity Max'][i]=np.max(aos[ind,i])
        sobs[v + ' Occurrence'][i]=1

#%% Wildfire indicator

ind=gis.GetGridIndexToPoints(zRef,sobs['X'],sobs['Y'])
tv=np.arange(1990,2022,1)

z1=np.zeros( (tv.size,ind[0].size) ,dtype='int8')
for iT in range(tv.size):
    z=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\PROT_HISTORICAL_FIRE_POLYS_SP_' + str(tv[iT]) + '.tif')
    z1[iT,:]=z['Data'][ind].copy()

sobs['Wildfire Occurrence']=np.zeros(sobs['Year t0'].size)
for i in range(sobs['Year t0'].size):
    ind1=np.where( (z1[:,i]==1) & (tv>=sobs['Year t0'][i]) & (tv<sobs['Year t1'][i]) )[0]
    if ind1.size==0:
        continue
    sobs['Wildfire Occurrence'][i]=1

#%% Add VRI variables

# Crown closure
z=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\crownc.tif')
ind=gis.GetGridIndexToPoints(z,sobs['X'],sobs['Y'])
z0=z['Data'][ind]
sobs['Crown Closure (VRI)']=z0

# Tree Density Class
z=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\LandUseLandCover\ForestDensityClass.tif')
ind=gis.GetGridIndexToPoints(z,sobs['X'],sobs['Y'])
z0=z['Data'][ind]
sobs['Tree Density Class (VRI)']=z0

#%% Add climate data

# Thornthwaite's climate classification
z=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\BC1ha_ThornthwaiteClimateClassCondensed_norm_1971to2000.tif')
ind=gis.GetGridIndexToPoints(z,sobs['X'],sobs['Y'])
z0=z['Data'][ind]
sobs['Climate Class']=z0

# Temperature
z=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\Seasonal\BC1ha_mat_norm_1971to2000_si_hist_v1.tif')
ind=gis.GetGridIndexToPoints(z,sobs['X'],sobs['Y'])
z0=z['Data'][ind]
sobs['ta_ann_n']=z0

# Climatic water deficit
z=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\Seasonal\wbm_cwd_gs_norm_1971to2000.tif')
ind=gis.GetGridIndexToPoints(z,sobs['X'],sobs['Y'])
z0=z['Data'][ind]
sobs['cwd_gs_n']=z0

# Potential evapotranspiration
z=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\Seasonal\wbm_etp_gs_norm_1971to2000.tif')
ind=gis.GetGridIndexToPoints(z,sobs['X'],sobs['Y'])
z0=z['Data'][ind]
sobs['etp_gs_n']=z0

# Soil water content
z=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Climate\Seasonal\wbm_ws_gs_norm_1971to2000.tif')
ind=gis.GetGridIndexToPoints(z,sobs['X'],sobs['Y'])
z0=z['Data'][ind]
sobs['ws_gs_n']=z0

#%% Remove plots with bad coordinates

#ind=np.where( (sobs['Lat']>30) & (sobs['Lon']!=0) & (sobs['Lat']!=0) )[0]
#for k in sobs.keys():
#    sobs[k]=sobs[k][ind]

#%% Save

df=pd.DataFrame(sobs)
df.to_excel(meta['Paths']['DB'] + '\\Processed\\L2\\L2_sobs_' + jurL[iJur] + '.xlsx',index=False)

d1={'tobs':tobs,'sobs':sobs}
gu.opickle(meta['Paths']['DB'] + '\\Processed\\L2\\L2_' + jurL[iJur] + '.pkl',d1)
