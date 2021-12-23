
'''
ROLLUP FCI - FROM INVENTORY - PREPARE INPUTS FOR TIPSY AND CBRUNNER

'''

#%% Import Python modules

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import geopandas as gpd
import openpyxl
import time
import glob
import gc as garc
from scipy import stats
from numpy import matlib as mb
from shapely.geometry import Point
from fcgadgets.macgyver import utilities_general as gu
from fcgadgets.macgyver import utilities_gis as gis
from fcgadgets.macgyver import utilities_inventory as invu
from fcgadgets.cbrunner import cbrun_utilities as cbu
from fcgadgets.cbrunner import cbrun as cbr

#%% Import paths

meta=gu.ipickle(r'D:\Data\FCI_Projects\SummaryReforestation\Inputs\Metadata.pkl')
meta['Paths']['Model Code']=r'C:\Users\rhember\Documents\Code_Python\fcgadgets\cbrunner'

#%% Import fcgadgets parameters

par=invu.Load_Params(meta)

#%% Import FCI-funded areas (from RESULTS ATU layer)

atu_multipolygons=gu.ipickle(meta['Paths']['Geospatial'] + '\\atu_multipolygons.pkl')

#%% Open invnetory files

sxy,atu,op,burnsev,vri,fcinv,fcsilv,fcres,pl,fire,pest,cut,lul,park,ogmal=invu.LoadSparseGeospatialInputs(meta)

#%% Import indices

idx=invu.LoadSparseGridIndex(meta)

#%% Configure project and scenario parameters

meta['N Stand']=sxy['x'].size
meta['N Scenario']=2
#meta['N Scenario']=5

# Index to scenario types
iB1=0 # Baseline
iP1=1 # For Reforestation underplanting
iP2=2 # For Dwarf Mistletoe
iP3=3 # For Site Preparation
iP4=4 # For Class A seed

# Import project parameters from spreadsheet
meta=cbu.ImportProjectConfig(meta)

#%% Specify the scenarios are affected by simulated disturbances on the fly

meta['Scenario Switch']['Dist on Fly']={}
meta['Scenario Switch']['Dist on Fly']['Harvesting (historical)']=[0,0]#,0,0,0]
meta['Scenario Switch']['Dist on Fly']['Harvesting (future)']=[0,0]#,0,0,0]
meta['Scenario Switch']['Dist on Fly']['Breakup']=[0,0]#,0,0,0]

meta['Scenario Switch']['Dist on Fly']['Harvesting (future) Year Start']=[2022,2022]#,2020,2020,2020]

#%% Extract disturbance/management history from inventory layers
# Takes < 2 min

t0=time.time()
dmec=invu.PrepDMEC(idx,meta,par,atu,pl,op,fcinv,vri,cut,fire,burnsev,pest)
t1=time.time()
print((t1-t0)/60)

#%% Exclude duplicate events
# Takes < 1 min

if meta['Exclude duplicate events']=='On':
    dmec=invu.Exclude_Duplicate_Events(meta,dmec)

#%% Exclude unidentified activities or disturbances
# Takes < 1 min

if meta['Exclude unidentified events']=='On':
    dmec=invu.Exclude_Unidentified_Events(meta,dmec)

#%% Add FCI Funded indicator

uPP=np.array(['Dummy Data'])

t0=time.time()
dmec=invu.IsFCIFunding(meta,dmec,uPP)
t1=time.time()
print((t1-t0)/60)

#%% Define strings to fill when adding events with only basic info

StringsToFill=['Month','Day','FCI Funded','FIA_PROJECT_ID','SILV_FUND_SOURCE_CODE','OPENING_ID','ACTUAL_TREATMENT_AREA','ACTUAL_PLANTED_NUMBER', \
            'PL_SPECIES_CD1','PL_SPECIES_PCT1','PL_SPECIES_GW1','PL_SPECIES_CD2','PL_SPECIES_PCT2','PL_SPECIES_GW2', \
            'PL_SPECIES_CD3','PL_SPECIES_PCT3','PL_SPECIES_GW3','PL_SPECIES_CD4','PL_SPECIES_PCT4','PL_SPECIES_GW4', \
            'PL_SPECIES_CD5','PL_SPECIES_PCT5','PL_SPECIES_GW5']

#%% Don't allow stand establishment milestones to be the first event in the modern era

# There needs to be an event preceding planting or direct seeding for the
# pre-project growth curve transitions to the baseline growth curve.
# This appears to be happening very rarely.

# Assumptions
AgeSincePrev=2
ID_Type_Prev=meta['LUT']['Dist']['Wildfire']
Sev_Prev=100

Ntot=np.array([])
Nfix=np.array([])
for iStand in range(meta['N Stand Full']):

    iA=np.where( (dmec[iStand]['ID_Type']==meta['LUT']['Dist']['Planting']) | \
                (dmec[iStand]['ID_Type']==meta['LUT']['Dist']['Direct Seeding']) )[0]
    
    if iA.size==0: 
        continue
    
    Ntot=np.append(Ntot,1)
    
    if (dmec[iStand]['Year'][iA[0]]==np.min(dmec[iStand]['Year'])):
        
        Nfix=np.append(Nfix,1) 
        
        dmec[iStand]['Year']=np.append(dmec[iStand]['Year'],dmec[iStand]['Year'][iA[0]]-AgeSincePrev)
        dmec[iStand]['ID_Type']=np.append(dmec[iStand]['ID_Type'],ID_Type_Prev)
        dmec[iStand]['MortalityFactor']=np.append(dmec[iStand]['MortalityFactor'],Sev_Prev)
        dmec[iStand]['GrowthFactor']=np.append(dmec[iStand]['GrowthFactor'],np.array(0,dtype='int16'))
        for v in StringsToFill:
           dmec[iStand][v]=np.append(dmec[iStand][v],-999)

# What fraction of planting events needed fixing (ie adding a first event in the record)
print([str(np.round(100*np.sum(Nfix)/np.sum(Ntot)*100)/100) + ' % of stand establishment projects required addition of previous disturbance (wildfire).'])

# Put events in order of calendar date
dmec=invu.PutEventsInOrder(dmec,meta)

#%% Ensure every stand has a modern disturbance

# So that there is at least one event - this is required because there are about 30
# locations with no reognized activities and no natural disturbances (LB, Layout, SP)
# Takes < 1 min

if meta['Ensure every stand has a modern disturbance']=='On':
    name_dist='Wildfire'
    severity=100
    dmec=invu.Ensure_Every_Stand_Has_Modern_Disturbance(meta,dmec,name_dist,severity,StringsToFill)

#%% IDW - Western Spruce Budworm - Fix severity
# The DMEC was populated with the numeric severity ID. Mortality only occurs 
# following repeated outrbreak years. 

# Takes < 1 min
t0=time.time()

if meta['Fix severity of western spruce budworm']=='On':
    dmec=invu.IDW_Fix_Severity(meta,dmec,par)
 
print((time.time()-t0)/60)
      
#%% QA: Check variable lengths

Check=np.zeros(meta['N Stand Full'])
for iStand in range(meta['N Stand Full']):
    n=[]
    for key in dmec[iStand].keys():
        n.append(dmec[iStand][key].size)
    Check[iStand]=np.unique(n).size
print('Inconsistent grid cells = ' + str(np.sum(Check!=1)/meta['N Stand Full']*100) + '%')

#%% Create index to the impacting disturbance for underplanting projects:

# It is critical that the above steps (ensuring there is a previous 
# stand-replacing disturbance) work before this will work properly.

for iStand in range(meta['N Stand Full']):
    
    # Initialize variable
    dmec[iStand]['IndPrevDistForUnderplanting']=-999*np.ones(dmec[iStand]['Year'].size)
    
    # This only applies to planting or direct seeding, only continue if there is 
    # such an event
    iPlant=np.where( (dmec[iStand]['ID_Type']==meta['LUT']['Dist']['Planting']) & (dmec[iStand]['FCI Funded']==1) | \
                     (dmec[iStand]['ID_Type']==meta['LUT']['Dist']['Direct Seeding']) & (dmec[iStand]['FCI Funded']==1) )[0]
    if iPlant.size==0: 
        continue
    
    # We don't want to apply this to projects with delayed regen (eg HB0000076)
    iDelayedRegen=np.where( (dmec[iStand]['GrowthFactor']==-100) )[0]
    if iDelayedRegen.size>0:
        continue
    
    # Find the previous stand-replacing disturbance
    iDist=np.where( (dmec[iStand]['ID_Type']==meta['LUT']['Dist']['Wildfire']) & (np.abs(dmec[iStand]['Year']-dmec[iStand]['Year'][iPlant[0]])<10) )[0]
    if iDist.size!=0:
        dmec[iStand]['IndPrevDistForUnderplanting'][iPlant]=iDist[-1]

#%% Indicator for which scenario is affected by events

# Notes: 
# I think this could be based on FCI Funded, only, but leave as is for now

for iStand in range(meta['N Stand Full']):    
    
    ScnAffected=[]
    ScnAffected.append(np.zeros(dmec[iStand]['Year'].size)) # baseline
    ScnAffected.append(np.zeros(dmec[iStand]['Year'].size)) # P1
    #ScnAffected.append(np.zeros(dmec[iStand]['Year'].size)) # P2
    #ScnAffected.append(np.zeros(dmec[iStand]['Year'].size)) # P3
    #ScnAffected.append(np.zeros(dmec[iStand]['Year'].size)) # P4
    for iYr in range(dmec[iStand]['Year'].size):
        
        ID_Type=dmec[iStand]['ID_Type'][iYr]
        FCI_Funded=dmec[iStand]['FCI Funded'][iYr]        
        
        if FCI_Funded!=1:
            # All events impact all scenarios
            ScnAffected[iB1][iYr]=1
            ScnAffected[iP1][iYr]=1
            #ScnAffected[iP2][iYr]=1
            #ScnAffected[iP3][iYr]=1
            #ScnAffected[iP4][iYr]=1
        else:
            # Only some events impact some scenarios
            if (ID_Type==meta['LUT']['Dist']['Fertilization Aerial']):
                ScnAffected[iP1][iYr]=1
            elif (ID_Type==meta['LUT']['Dist']['Fertilization Hand']):
                ScnAffected[iP1][iYr]=1
            elif (ID_Type==meta['LUT']['Dist']['Planting']):
                ScnAffected[iP1][iYr]=1
                #ScnAffected[iP2][iYr]=1
                #ScnAffected[iP3][iYr]=1
                #ScnAffected[iP4][iYr]=1
            elif (ID_Type==meta['LUT']['Dist']['Direct Seeding']):
                ScnAffected[iP1][iYr]=1
            #elif (ID_Type==meta['LUT']['Dist']['Dwarf Mistletoe Control']):
                #ScnAffected[iP2][iYr]=1
            #elif (ID_Type==meta['LUT']['Dist']['Ripping']):
                #ScnAffected[iP3][iYr]=1
            #elif (ID_Type==meta['LUT']['Dist']['Disc Trenching']):
                #ScnAffected[iP3][iYr]=1                
            else:
                # All other events impact both baseline and project
                ScnAffected[iB1][iYr]=1
                ScnAffected[iP1][iYr]=1
                #ScnAffected[iP2][iYr]=1
                #ScnAffected[iP3][iYr]=1
                #ScnAffected[iP4][iYr]=1
            
    dmec[iStand]['ScnAffected']=[]
    dmec[iStand]['ScnAffected'].append(ScnAffected[iB1])
    dmec[iStand]['ScnAffected'].append(ScnAffected[iP1])
    #dmec[iStand]['ScnAffected'].append(ScnAffected[iP2])
    #dmec[iStand]['ScnAffected'].append(ScnAffected[iP3])
    #dmec[iStand]['ScnAffected'].append(ScnAffected[iP4])

#%% Clean species composition - TIPSY will not recognize certain codes

meta,dmec,vri,fcinv=invu.Clean_Species_Composition(meta,dmec,vri,fcinv)
                
#%% Define best-available (gap-filled) inventory
# Takes < 1 min

t0=time.time()
flag_projects=0
ba,ba_source=invu.CreateBestAvailableInventory(meta,vri,fcinv,flag_projects,idx,sxy)
print((time.time()-t0)/60)

# Save for use in analysis of results
gu.opickle(meta['Paths']['Project'] + '\\Inputs\\ba.pkl',ba)

#%% Reduce number of growth curves by adjusting stand conditions

# Site index
flg=0
if flg==1:
    trig=0
    for i in range(1,55,2):
        ind=np.where(ba['SI']==i)[0]
        if trig==0:
            ba['SI'][ind]=ba['SI'][ind]+1
            trig=1
        else:
            ba['SI'][ind]=ba['SI'][ind]-1
            trig=0

#%% Reduce species composition

flg=0
if flg==1:
    for iStand in range(meta['N Stand']):
        dN=100-(ba['Spc_Pct1'][iStand]+ba['Spc_Pct2'][iStand])
        ba['Spc_Pct1'][iStand]=ba['Spc_Pct1'][iStand]+dN        
        ba['Spc_CD3'][iStand]=-999
        ba['Spc_Pct3'][iStand]=-999
        ba['Spc_CD4'][iStand]=-999
        ba['Spc_Pct4'][iStand]=-999
        ba['Spc_CD5'][iStand]=-999
        ba['Spc_Pct5'][iStand]=-999

#%% Parameterize growth curves

# Takes 4 min
t0=time.time()

# Parameter list
GC_Variable_List=['ID_Stand','ID_Scn','ID_GC','regeneration_method','s1','p1','i1','s2', \
   'p2','s3','p3','s4','p4','s5','p5','gain1','selage1','gain2','selage2', \
   'gain3','selage3','gain4','selage4','gain5','selage5', \
	'init_density','regen_delay','oaf1','oaf2','bec_zone','FIZ', \
	'fert_age1','fert_age2','fert_age3','fert_age4','fert_age5']

# Initialize list
gc=[None]*meta['N Stand Full']

# Loop through stands
for iStand in range(meta['N Stand Full']):
    
    # Initialize growth curve identifiers
    ID_GC=[None]*meta['N Scenario']
    ID_GC[iB1]=1*np.ones(dmec[iStand]['Year'].size) # Baseline
    ID_GC[iP1]=1*np.ones(dmec[iStand]['Year'].size) # P1
    #ID_GC[iP2]=1*np.ones(dmec[iStand]['Year'].size) # P2
    #ID_GC[iP3]=1*np.ones(dmec[iStand]['Year'].size) # P3
    #ID_GC[iP4]=1*np.ones(dmec[iStand]['Year'].size) # P3
    
    # Initialize variables within growth curve dictionary
    gc0=[None]*meta['N Scenario']
    Counter=[None]*meta['N Scenario']
    for iScn in range(len(gc0)):
        gc0[iScn]={}
        for key in GC_Variable_List:
            gc0[iScn][key]=-999*np.ones(12)
        Counter[iScn]=0
    
    #--------------------------------------------------------------------------
    # Add pre-modern inventory curve
    #--------------------------------------------------------------------------
    
    for iScn in range(meta['N Scenario']):
        c=Counter[iScn]
        gc0[iScn]['ID_Stand'][c]=iStand
        gc0[iScn]['ID_Scn'][c]=iScn
        gc0[iScn]['ID_GC'][c]=1   
        gc0[iScn]['regeneration_method'][c]=meta['LUT']['TIPSY']['regeneration_method']['N']
        gc0[iScn]['s1'][c]=ba['Spc_CD1'][iStand]
        gc0[iScn]['p1'][c]=ba['Spc_Pct1'][iStand]
        gc0[iScn]['i1'][c]=ba['SI'][iStand]
        gc0[iScn]['s2'][c]=ba['Spc_CD2'][iStand]
        gc0[iScn]['p2'][c]=ba['Spc_Pct2'][iStand]
        gc0[iScn]['s3'][c]=ba['Spc_CD3'][iStand]
        gc0[iScn]['p3'][c]=ba['Spc_Pct3'][iStand]
        gc0[iScn]['s4'][c]=ba['Spc_CD4'][iStand]
        gc0[iScn]['p4'][c]=ba['Spc_Pct4'][iStand]
        gc0[iScn]['s5'][c]=ba['Spc_CD5'][iStand]
        gc0[iScn]['p5'][c]=ba['Spc_Pct5'][iStand]
        gc0[iScn]['init_density'][c]=2000
        gc0[iScn]['regen_delay'][c]=2
        gc0[iScn]['oaf1'][c]=0.85
        gc0[iScn]['oaf2'][c]=0.95
        gc0[iScn]['bec_zone'][c]=ba['BEC_ZONE_CODE'][iStand]
        gc0[iScn]['FIZ'][c]=ba['FIZ'][iStand]
        Counter[iScn]=Counter[iScn]+1
    
    #--------------------------------------------------------------------------
    # Add events from disturbance/management event history
    #--------------------------------------------------------------------------
    
    for iYr in range(dmec[iStand]['Year'].size):
        
        # Index to previous disturbance for fertilization
        #IndPrevDistForFert=int(dmec[iStand]['IndPrevDistForFert'][iYr])
        
        # Index to previous disturbance for underplanting
        IndPrevDistForUnderplanting=int(dmec[iStand]['IndPrevDistForUnderplanting'][iYr])
        
        # Calculate planting density
        PlantingDensity=dmec[iStand]['ACTUAL_PLANTED_NUMBER'][iYr]/dmec[iStand]['ACTUAL_TREATMENT_AREA'][iYr]
        PlantingDensity=int(PlantingDensity)
        
        # Create a flag that indicates whether there are back-to-back planting
        # Back to back planting I think occurs in some cases because they go back
        # and add a bit. You can tell by looking at the treatment area - if the
        # second planting treatment area is tiny compared to the first, you could
        # ignore it I guess. Certainly not ideal, but I don't see a work around.
        # We are going to ignore the second planting for now.
        Flag_PlantingBackToBack=0
        if iYr>0:
            if (dmec[iStand]['ID_Type'][iYr]==meta['LUT']['Dist']['Planting']):
                Yr=dmec[iStand]['Year'][iYr]
                indPrevPL=np.where( (dmec[iStand]['Year']==Yr-1) & (dmec[iStand]['ID_Type']==meta['LUT']['Dist']['Planting']) )[0]
                if indPrevPL.size>0:
                    Flag_PlantingBackToBack=1
        
        #----------------------------------------------------------------------
        # Planting, not FCI funded
        #----------------------------------------------------------------------
        
        if (dmec[iStand]['ID_Type'][iYr]==meta['LUT']['Dist']['Planting']) & (dmec[iStand]['FCI Funded'][iYr]==0) & (Flag_PlantingBackToBack==0):
                
            ID_GC[iB1][iYr:]=ID_GC[iB1][iYr]+1
            ID_GC[iP1][iYr:]=ID_GC[iP1][iYr]+1
            #ID_GC[iP2][iYr:]=ID_GC[iP2][iYr]+1
            #ID_GC[iP3][iYr:]=ID_GC[iP3][iYr]+1
            #ID_GC[iP4][iYr:]=ID_GC[iP4][iYr]+1
            
            ScnToInclude=[iB1,iP1]
            #ScnToInclude=[iB1,iP1,iP2,iP3,iP4]
            for i in range(len(ScnToInclude)):
                iScn=ScnToInclude[i]
                c=Counter[iScn]
                gc0[iScn]['ID_Stand'][c]=iStand
                gc0[iScn]['ID_Scn'][c]=iScn
                gc0[iScn]['ID_GC'][c]=ID_GC[iScn][iYr]
                gc0[iScn]['regeneration_method'][c]=meta['LUT']['TIPSY']['regeneration_method']['P']
                gc0[iScn]['i1'][c]=ba['SI'][iStand]
                if (dmec[iStand]['PL_SPECIES_CD1'][iYr]>0) & (dmec[iStand]['PL_SPECIES_PCT1'][iYr]>0):
                    gc0[iScn]['s1'][c]=dmec[iStand]['PL_SPECIES_CD1'][iYr]
                    gc0[iScn]['p1'][c]=dmec[iStand]['PL_SPECIES_PCT1'][iYr]
                    gc0[iScn]['gain1'][c]=dmec[iStand]['PL_SPECIES_GW1'][iYr]
                    gc0[iScn]['selage1'][c]=10
                    
                    gc0[iScn]['s2'][c]=dmec[iStand]['PL_SPECIES_CD2'][iYr]
                    gc0[iScn]['p2'][c]=dmec[iStand]['PL_SPECIES_PCT2'][iYr]
                    gc0[iScn]['gain2'][c]=dmec[iStand]['PL_SPECIES_GW2'][iYr]
                    gc0[iScn]['selage2'][c]=10
                    
                    gc0[iScn]['s3'][c]=dmec[iStand]['PL_SPECIES_CD3'][iYr]
                    gc0[iScn]['p3'][c]=dmec[iStand]['PL_SPECIES_PCT3'][iYr]
                    gc0[iScn]['gain3'][c]=dmec[iStand]['PL_SPECIES_GW3'][iYr]
                    gc0[iScn]['selage3'][c]=10
                    
                    gc0[iScn]['s4'][c]=dmec[iStand]['PL_SPECIES_CD4'][iYr]
                    gc0[iScn]['p4'][c]=dmec[iStand]['PL_SPECIES_PCT4'][iYr]
                    gc0[iScn]['gain4'][c]=dmec[iStand]['PL_SPECIES_GW4'][iYr]
                    gc0[iScn]['selage4'][c]=10
                    
                    gc0[iScn]['s5'][c]=dmec[iStand]['PL_SPECIES_CD5'][iYr]
                    gc0[iScn]['p5'][c]=dmec[iStand]['PL_SPECIES_PCT5'][iYr]                    
                    gc0[iScn]['gain5'][c]=dmec[iStand]['PL_SPECIES_GW5'][iYr]
                    gc0[iScn]['selage5'][c]=10
                else:
                    gc0[iScn]['s1'][c]=ba['Spc_CD1'][iStand]
                    gc0[iScn]['p1'][c]=ba['Spc_Pct1'][iStand]
  
                    gc0[iScn]['s2'][c]=ba['Spc_CD2'][iStand]
                    gc0[iScn]['p2'][c]=ba['Spc_Pct2'][iStand]
                    
                    gc0[iScn]['s3'][c]=ba['Spc_CD3'][iStand]
                    gc0[iScn]['p3'][c]=ba['Spc_Pct3'][iStand]
                    
                    gc0[iScn]['s4'][c]=ba['Spc_CD4'][iStand]
                    gc0[iScn]['p4'][c]=ba['Spc_Pct4'][iStand]
                    
                    gc0[iScn]['s5'][c]=ba['Spc_CD5'][iStand]
                    gc0[iScn]['p5'][c]=ba['Spc_Pct5'][iStand]
                
                gc0[iScn]['init_density'][c]=int(PlantingDensity)
                gc0[iScn]['regen_delay'][c]=0
                gc0[iScn]['oaf1'][c]=0.85
                gc0[iScn]['oaf2'][c]=0.95
                gc0[iScn]['bec_zone'][c]=ba['BEC_ZONE_CODE'][iStand]
                gc0[iScn]['FIZ'][c]=ba['FIZ'][iStand]
                Counter[iScn]=Counter[iScn]+1
        
        #----------------------------------------------------------------------
        # Planting, FCI funded
        #----------------------------------------------------------------------
        
        if (dmec[iStand]['ID_Type'][iYr]==meta['LUT']['Dist']['Planting']) & \
            (dmec[iStand]['FCI Funded'][iYr]==1) & \
            (Flag_PlantingBackToBack==0):
            
            # Adjust the baseline growth curve at the time of the preceding disturbance
            # *** This will only work with previous stand-replacing disturbances ***
            if IndPrevDistForUnderplanting!=-999:
                ID_GC[iB1][IndPrevDistForUnderplanting:]=ID_GC[iB1][IndPrevDistForUnderplanting]+1
            
            # Adjust the project growth curve at the time of planting
            ID_GC[iP1][iYr:]=ID_GC[iP1][iYr]+1
            #ID_GC[iP2][iYr:]=ID_GC[iP2][iYr]+1
            #ID_GC[iP3][iYr:]=ID_GC[iP3][iYr]+1
            #ID_GC[iP4][iYr:]=ID_GC[iP4][iYr]+1
            
            # Determine whether there was Dwarf Mistletoe control or site prep
            iDMC=np.where( (dmec[iStand]['ID_Type']==meta['LUT']['Dist']['Dwarf Mistletoe Control']) )[0]
            iRip=np.where( (dmec[iStand]['ID_Type']==meta['LUT']['Dist']['Ripping']) )[0]
            iDiscT=np.where( (dmec[iStand]['ID_Type']==meta['LUT']['Dist']['Disc Trenching']) )[0]            
            
            ScnToInclude=[iB1,iP1]
            #ScnToInclude=[iB1,iP1,iP2,iP3,iP4]
            for i in range(len(ScnToInclude)):
                iScn=ScnToInclude[i]
                c=Counter[iScn]
                gc0[iScn]['ID_Stand'][c]=iStand
                gc0[iScn]['ID_Scn'][c]=iScn
                gc0[iScn]['ID_GC'][c]=ID_GC[iScn][iYr]
                
                # Regeneration method
                if (iScn==iB1):
                    gc0[iScn]['regeneration_method'][c]=meta['LUT']['TIPSY']['regeneration_method']['N']
                else:
                    gc0[iScn]['regeneration_method'][c]=meta['LUT']['TIPSY']['regeneration_method']['P']
                
                # Site index
                gc0[iScn]['i1'][c]=ba['SI'][iStand]
                
                if (iScn==iP2) & (iDMC.size>0):
                    SI_add_For_DMC=1
                    gc0[iScn]['i1'][c]=ba['SI'][iStand]+SI_add_For_DMC
                
                if (iScn==iP3) & (iRip.size>0) | (iScn==iP3) & (iDiscT.size>0):
                    SI_add_For_SP=1
                    gc0[iScn]['i1'][c]=ba['SI'][iStand]+SI_add_For_SP                
                    
                # Species composition and genetic worth
                if (iScn==iP4) & (dmec[iStand]['PL_SPECIES_CD1'][iYr]>0) & (dmec[iStand]['PL_SPECIES_PCT1'][iYr]>0):
                    
                    # Include genetic gain
                    gc0[iScn]['s1'][c]=dmec[iStand]['PL_SPECIES_CD1'][iYr]
                    gc0[iScn]['p1'][c]=dmec[iStand]['PL_SPECIES_PCT1'][iYr]
                    gc0[iScn]['gain1'][c]=int(dmec[iStand]['PL_SPECIES_GW1'][iYr])
                    gc0[iScn]['selage1'][c]=10
                    
                    gc0[iScn]['s2'][c]=dmec[iStand]['PL_SPECIES_CD2'][iYr]
                    gc0[iScn]['p2'][c]=dmec[iStand]['PL_SPECIES_PCT2'][iYr]
                    gc0[iScn]['gain2'][c]=int(dmec[iStand]['PL_SPECIES_GW2'][iYr])
                    gc0[iScn]['selage2'][c]=10
                    
                    gc0[iScn]['s3'][c]=dmec[iStand]['PL_SPECIES_CD3'][iYr]
                    gc0[iScn]['p3'][c]=dmec[iStand]['PL_SPECIES_PCT3'][iYr]
                    gc0[iScn]['gain3'][c]=int(dmec[iStand]['PL_SPECIES_GW3'][iYr])
                    gc0[iScn]['selage3'][c]=10
                    
                    gc0[iScn]['s4'][c]=dmec[iStand]['PL_SPECIES_CD4'][iYr]
                    gc0[iScn]['p4'][c]=dmec[iStand]['PL_SPECIES_PCT4'][iYr]
                    gc0[iScn]['gain4'][c]=int(dmec[iStand]['PL_SPECIES_GW4'][iYr])
                    gc0[iScn]['selage4'][c]=10
                    
                    gc0[iScn]['s5'][c]=dmec[iStand]['PL_SPECIES_CD5'][iYr]
                    gc0[iScn]['p5'][c]=dmec[iStand]['PL_SPECIES_PCT5'][iYr]                    
                    gc0[iScn]['gain5'][c]=int(dmec[iStand]['PL_SPECIES_GW5'][iYr])
                    gc0[iScn]['selage5'][c]=10
                    
                elif (iScn!=iB1) & (iScn!=iP4) & (dmec[iStand]['PL_SPECIES_CD1'][iYr]>0) & (dmec[iStand]['PL_SPECIES_PCT1'][iYr]>0):
                    
                    # Use planted species, but exclude genetic gain
                    gc0[iScn]['s1'][c]=dmec[iStand]['PL_SPECIES_CD1'][iYr]
                    gc0[iScn]['p1'][c]=dmec[iStand]['PL_SPECIES_PCT1'][iYr]
                    
                    gc0[iScn]['s2'][c]=dmec[iStand]['PL_SPECIES_CD2'][iYr]
                    gc0[iScn]['p2'][c]=dmec[iStand]['PL_SPECIES_PCT2'][iYr]
                    
                    gc0[iScn]['s3'][c]=dmec[iStand]['PL_SPECIES_CD3'][iYr]
                    gc0[iScn]['p3'][c]=dmec[iStand]['PL_SPECIES_PCT3'][iYr]
                    
                    gc0[iScn]['s4'][c]=dmec[iStand]['PL_SPECIES_CD4'][iYr]
                    gc0[iScn]['p4'][c]=dmec[iStand]['PL_SPECIES_PCT4'][iYr]
                    
                    gc0[iScn]['s5'][c]=dmec[iStand]['PL_SPECIES_CD5'][iYr]
                    gc0[iScn]['p5'][c]=dmec[iStand]['PL_SPECIES_PCT5'][iYr]                    
                
                elif (iScn!=iB1) & (iScn!=iP4) & (dmec[iStand]['PL_SPECIES_CD1'][iYr]<0):
                    
                    # Planting, but planting information missing (this happens for
                    # a disc trenching project by BCTS)
                    # Assume 100% pine
                    gc0[iScn]['s1'][c]=58
                    gc0[iScn]['p1'][c]=100
                
                else:
                    
                    # Baseline
                    gc0[iScn]['s1'][c]=ba['Spc_CD1'][iStand]
                    gc0[iScn]['p1'][c]=ba['Spc_Pct1'][iStand]                
                     
                    gc0[iScn]['s2'][c]=ba['Spc_CD2'][iStand]
                    gc0[iScn]['p2'][c]=ba['Spc_Pct2'][iStand]
                    
                    gc0[iScn]['s3'][c]=ba['Spc_CD3'][iStand]
                    gc0[iScn]['p3'][c]=ba['Spc_Pct3'][iStand]
                    
                    gc0[iScn]['s4'][c]=ba['Spc_CD4'][iStand]
                    gc0[iScn]['p4'][c]=ba['Spc_Pct4'][iStand]
                    
                    gc0[iScn]['s5'][c]=ba['Spc_CD5'][iStand]
                    gc0[iScn]['p5'][c]=ba['Spc_Pct5'][iStand]                    
                
                # Regeneration delay
                if iScn==iB1:
                    gc0[iScn]['init_density'][c]=200
                    gc0[iScn]['regen_delay'][c]=2
                else:
                    gc0[iScn]['init_density'][c]=int(PlantingDensity)
                    gc0[iScn]['regen_delay'][c]=0
                
                # OAF1
                gc0[iScn]['oaf1'][c]=0.85
                
                # OAF2
                gc0[iScn]['oaf2'][c]=0.95                   
                
                gc0[iScn]['bec_zone'][c]=ba['BEC_ZONE_CODE'][iStand]
                gc0[iScn]['FIZ'][c]=ba['FIZ'][iStand]
                gc0[iScn]['fert_age1'][c]=-999
                Counter[iScn]=Counter[iScn]+1
    
        #----------------------------------------------------------------------
        # Direct seeding, not FCI funded
        #----------------------------------------------------------------------
        
        if (dmec[iStand]['ID_Type'][iYr]==meta['LUT']['Dist']['Direct Seeding']) & \
            (dmec[iStand]['FCI Funded'][iYr]==0):
            
            ID_GC[iB1][iYr:]=ID_GC[iB1][iYr]+1
            ID_GC[iP1][iYr:]=ID_GC[iP1][iYr]+1
            #ID_GC[iP2][iYr:]=ID_GC[iP2][iYr]+1
            #ID_GC[iP3][iYr:]=ID_GC[iP3][iYr]+1
            #ID_GC[iP4][iYr:]=ID_GC[iP4][iYr]+1
            
            ScnToInclude=[iB1,iP1]
            #ScnToInclude=[iB1,iP1,iP2,iP3,iP4]
            for i in range(len(ScnToInclude)):
                iScn=ScnToInclude[i]
                c=Counter[iScn]
                gc0[iScn]['ID_Stand'][c]=iStand
                gc0[iScn]['ID_Scn'][c]=iScn
                gc0[iScn]['ID_GC'][c]=ID_GC[iScn][iYr]
                gc0[iScn]['regeneration_method'][c]=meta['LUT']['TIPSY']['regeneration_method']['N']
                gc0[iScn]['i1'][c]=ba['SI'][iStand]
                gc0[iScn]['s1'][c]=ba['Spc_CD1'][iStand]
                gc0[iScn]['p1'][c]=ba['Spc_Pct1'][iStand]                
                gc0[iScn]['s2'][c]=ba['Spc_CD2'][iStand]
                gc0[iScn]['p2'][c]=ba['Spc_Pct2'][iStand]
                gc0[iScn]['s3'][c]=ba['Spc_CD3'][iStand]
                gc0[iScn]['p3'][c]=ba['Spc_Pct3'][iStand]
                gc0[iScn]['s4'][c]=ba['Spc_CD4'][iStand]
                gc0[iScn]['p4'][c]=ba['Spc_Pct4'][iStand]
                gc0[iScn]['s5'][c]=ba['Spc_CD5'][iStand]
                gc0[iScn]['p5'][c]=ba['Spc_Pct5'][iStand]
                gc0[iScn]['init_density'][c]=2000
                gc0[iScn]['regen_delay'][c]=0                
                gc0[iScn]['oaf1'][c]=0.85
                gc0[iScn]['oaf2'][c]=0.95
                gc0[iScn]['bec_zone'][c]=ba['BEC_ZONE_CODE'][iStand]
                gc0[iScn]['FIZ'][c]=ba['FIZ'][iStand]
                gc0[iScn]['fert_age1'][c]=-999
                Counter[iScn]=Counter[iScn]+1
                
        #----------------------------------------------------------------------
        # Direct seeding, FCI funded
        #----------------------------------------------------------------------
        
        if (dmec[iStand]['ID_Type'][iYr]==meta['LUT']['Dist']['Direct Seeding']) & \
            (dmec[iStand]['FCI Funded'][iYr]==1):
            
            # Adjust baseline growth curve at time of previous disturbance
            if IndPrevDistForUnderplanting!=-999:
                ID_GC[iB1][IndPrevDistForUnderplanting:]=ID_GC[iB1][IndPrevDistForUnderplanting]+1
            
            # Adjust project growth curve at time of seeding
            ID_GC[iP1][iYr:]=ID_GC[iP1][iYr]+1
            #ID_GC[iP2][iYr:]=ID_GC[iP2][iYr]+1
            #ID_GC[iP3][iYr:]=ID_GC[iP3][iYr]+1
            #ID_GC[iP4][iYr:]=ID_GC[iP4][iYr]+1
            
            ScnToInclude=[iB1,iP1]
            #ScnToInclude=[iB1,iP1,iP2,iP3,iP4]
            for i in range(len(ScnToInclude)):
                iScn=ScnToInclude[i]
                c=Counter[iScn]
                gc0[iScn]['ID_Stand'][c]=iStand
                gc0[iScn]['ID_Scn'][c]=iScn
                gc0[iScn]['ID_GC'][c]=ID_GC[iScn][iYr]
                gc0[iScn]['regeneration_method'][c]=meta['LUT']['TIPSY']['regeneration_method']['N']
                gc0[iScn]['i1'][c]=ba['SI'][iStand]
                gc0[iScn]['s1'][c]=ba['Spc_CD1'][iStand]
                gc0[iScn]['p1'][c]=ba['Spc_Pct1'][iStand]                
                gc0[iScn]['s2'][c]=ba['Spc_CD2'][iStand]
                gc0[iScn]['p2'][c]=ba['Spc_Pct2'][iStand]
                gc0[iScn]['s3'][c]=ba['Spc_CD3'][iStand]
                gc0[iScn]['p3'][c]=ba['Spc_Pct3'][iStand]
                gc0[iScn]['s4'][c]=ba['Spc_CD4'][iStand]
                gc0[iScn]['p4'][c]=ba['Spc_Pct4'][iStand]
                gc0[iScn]['s5'][c]=ba['Spc_CD5'][iStand]
                gc0[iScn]['p5'][c]=ba['Spc_Pct5'][iStand]
              
                if iScn==iB1:
                    # This should be irrelevent - growth factor will force it
                    # to be zero.
                    gc0[iScn]['init_density'][c]=200
                    gc0[iScn]['regen_delay'][c]=2
                elif iScn!=iB1:
                    gc0[iScn]['init_density'][c]=2000
                    gc0[iScn]['regen_delay'][c]=0
                
                gc0[iScn]['oaf1'][c]=0.85
                gc0[iScn]['oaf2'][c]=0.95
                gc0[iScn]['bec_zone'][c]=ba['BEC_ZONE_CODE'][iStand]
                gc0[iScn]['FIZ'][c]=ba['FIZ'][iStand]                
                Counter[iScn]=Counter[iScn]+1
    
    #--------------------------------------------------------------------------
    # Add growth curve identifiers to disturbance/management dictionary
    #--------------------------------------------------------------------------
    
    gc[iStand]=gc0
    dmec[iStand]['ID_GC']=ID_GC

#------------------------------------------------------------------------------
# Get rid of rows with no info
#------------------------------------------------------------------------------

for iStand in range(meta['N Stand Full']):
    for iScn in range(meta['N Scenario']):
        ind=np.where(gc[iStand][iScn]['ID_Stand']!=-999)[0]
        for key in GC_Variable_List:
            gc[iStand][iScn][key]=gc[iStand][iScn][key][ind]

print((time.time()-t0)/60)

#------------------------------------------------------------------------------
# Save
#------------------------------------------------------------------------------

gu.opickle(meta['Paths']['Project'] + '\\Inputs\\gc.pkl',gc)


#%% Adjust mortality factors that only affect specific tree species
# Takes < 1 min

t0=time.time()
if meta['Adjust species-specific mortality']=='On':
    dmec=invu.AdjustSpeciesSpecificMortality(meta,dmec,par,gc,iB1)
print((time.time()-t0)/60)  

#%% Save dmec

gu.opickle(meta['Paths']['Project'] + '\\Inputs\\dmec.pkl',dmec)      

#%% Extract a set of unique growth curves
# Decompose the full set of stands into a subset of unique stand types.
# Exclude the first three columns, as they are all different.

# Takes < 1 min
t0=time.time()

ugc={}
ugc['GC_Variable_List']=np.array(GC_Variable_List)[3:]

# Calculate unique stand types
#ugc['Full']=np.zeros((4e6,len(GC_Variable_List)-3))
ugc['Full']=np.zeros((int(4e6),len(GC_Variable_List)))

cnt=0
for iStand in range(meta['N Stand']):
    for iScn in range(meta['N Scenario']):
        for iGC in range(gc[iStand][iScn]['ID_GC'].size):
            for k in range(len(GC_Variable_List)):
                key=GC_Variable_List[k]
                ugc['Full'][cnt,k]=gc[iStand][iScn][key][iGC]
            cnt=cnt+1
ugc['Full']=ugc['Full'][0:cnt,:]
    
# Unique growth curves
# The 'Inverse' variable acts as the crosswalk between the full and unique gc arrays
ugc['Unique'],ugc['Index'],ugc['Inverse']=np.unique(ugc['Full'][:,3:],return_index=True,return_inverse=True,axis=0)

# Save
flg=0
if flg==1:
    gu.opickle(meta['Path Project'] + '\\Inputs\\ugc.pkl',ugc)

## Add the crosswalk to the gc dictionary
#flg=0
#if flg==1:
#    gcFull=np.column_stack((ugc['Full'][:,0:3],ugc['Inverse']))
#    df=pd.DataFrame(data=gcFull,columns={'ID_Stand','ID_Scn','ID_GC','ID_GC_Unique'})
#    gu.opickle(r'C:\Users\rhember\Documents\Data\FCI_Projects\FertilizationSummary\ForOthers\gc_full.pkl',gcFull)
#    df.to_csv(r'C:\Users\rhember\Documents\Data\FCI_Projects\FertilizationSummary\ForOthers\gc_full.csv')

print((time.time()-t0)/60)

#%% Export to BatchTIPSY parameter spreadsheet
# Notes: Make sure all lines are deleted from preivous run, or they will not be
# overwritten.
# Takes 2 min

t0=time.time()
cbu.Write_BatchTIPSY_Input_Spreadsheet(meta,ugc)
print((time.time()-t0)/60)

#%% Populate BatchTIPSY.exe input variable (.dat) file 

# Takes 1 minutes!

# Collect garbage
garc.collect()

t0=time.time()
cbu.Write_BatchTIPSY_Input_File(meta)
print((time.time()-t0)/60)

# ************************* MANUAL OPERATION **********************************
#------------------------------------------------------------------------------
# Ensure BatchTIPSY.exe config (.BPS) file has up-to-date paths 
# Run BatchTIPSY.exe
#------------------------------------------------------------------------------
# ************************* MANUAL OPERATION **********************************

#%% Timber harvesting land base

# Define THLB
thlb_flag_Actual,thlb_flag_Baseline=invu.DefineTHLB(meta,ba,dmec,fcres,lul,ogmal,park)

# Define differences in THLB by scenario
thlb_by_scenario=[None]*meta['N Scenario']
thlb_by_scenario[0]=thlb_flag_Actual
thlb_by_scenario[1]=thlb_flag_Actual
#thlb_by_scenario[2]=thlb_flag_Actual
#thlb_by_scenario[3]=thlb_flag_Actual
#thlb_by_scenario[4]=thlb_flag_Actual

#%% Prepare inventory

# Loop through batches, saving inventory to file
for iScn in range(meta['N Scenario']):
    
    # Loop through batches, saving inventory to file
    for iBat in range(meta['N Batch']):
      
        inv={}
    
        # Index to batch
        indBat=cbu.IndexToBatch(meta,iBat)
        
        N_StandsInBatch=len(indBat)
    
        # Initialize inventory variables
        inv['Lat']=np.zeros((1,N_StandsInBatch))
        inv['Lon']=np.zeros((1,N_StandsInBatch))
        inv['X']=inv['Lat']
        inv['Y']=inv['Lon']
        
        # BEC zone
        inv['ID_BECZ']=np.zeros((1,N_StandsInBatch),dtype=np.int)    
        for i in range(inv['ID_BECZ'].size):
            try:
                inv['ID_BECZ'][0,i]=ba['BEC_ZONE_CODE'][indBat[i]]
            except:
                inv['ID_BECZ'][0,i]=meta['LUT BGC Zone']['SBS'] 
    
        # Timber harvesting landbase (1=yes, 0=no)
        # Conditional upon scenario switch
        inv['THLB']=thlb_by_scenario[iScn][:,indBat]
    
        # Temperature will be updated automatically
        inv['MAT']=4*np.ones((1,N_StandsInBatch))
    
        if meta['Biomass Module']=='Sawtooth':
            inv['Srs1_ID']=meta['LUT Spc'][meta['Scenario'][iScn]['SRS1_CD']]*np.ones((1,N_StandsInBatch),dtype=np.int)
        else:
            inv['Srs1_ID']=9999*np.ones((1,N_StandsInBatch),dtype=np.int)
        #inv['Srs1_Pct']=100*np.ones((1,N_StandsInBatch),dtype=np.int)
        #inv['Srs2_ID']=0*np.ones((1,N_StandsInBatch),dtype=np.int)
        #inv['Srs2_Pct']=0*np.ones((1,N_StandsInBatch),dtype=np.int)
        #inv['Srs3_ID']=0*np.ones((1,N_StandsInBatch),dtype=np.int)
        #inv['Srs3_Pct']=0*np.ones((1,N_StandsInBatch),dtype=np.int)
        
        inv['Spc1_ID']=0*np.ones((1,N_StandsInBatch),dtype=np.int)
        inv['Spc1_Pct']=0*np.ones((1,N_StandsInBatch),dtype=np.int)
        inv['Spc2_ID']=0*np.ones((1,N_StandsInBatch),dtype=np.int)
        inv['Spc2_Pct']=0*np.ones((1,N_StandsInBatch),dtype=np.int)
        inv['Spc3_ID']=0*np.ones((1,N_StandsInBatch),dtype=np.int)
        inv['Spc3_Pct']=0*np.ones((1,N_StandsInBatch),dtype=np.int)
        inv['Spc3_ID']=0*np.ones((1,N_StandsInBatch),dtype=np.int)
        inv['Spc3_Pct']=0*np.ones((1,N_StandsInBatch),dtype=np.int)

        # Save
        gu.opickle(meta['Paths']['Input Scenario'][iScn] + '\\Inventory_Bat' + cbu.FixFileNum(iBat) + '.pkl',inv)

#%% Prepare disturbance/management event history

tv=np.arange(meta['Year Start'],meta['Year End']+1,1)

# Loop through scenarios, ensembles and batches
for iScn in range(meta['N Scenario']):    
    for iEns in range(meta['N Ensemble']):        
        for iBat in range(meta['N Batch']):
    
            # Index to batch
            indBat=cbu.IndexToBatch(meta,iBat)
            
            # Initialize dictionary
            ec={}
            ec['ID_Type']=np.zeros((meta['Year'].size,indBat.size,meta['Max Events Per Year']),dtype='int16')
            ec['MortalityFactor']=np.zeros((meta['Year'].size,indBat.size,meta['Max Events Per Year']),dtype='int16')
            ec['GrowthFactor']=np.zeros((meta['Year'].size,indBat.size,meta['Max Events Per Year']),dtype='int16')
            ec['ID_GrowthCurve']=np.zeros((meta['Year'].size,indBat.size,meta['Max Events Per Year']),dtype='int16')
            
            for iS in range(len(indBat)):
                
                # Index to stand
                iStandFull=indBat[iS]
            
                #--------------------------------------------------------------
                # Pre-industrial disturbances
                #--------------------------------------------------------------
    
                # Pre-industrial disturbance interval
                if ba['FIZ'][iStandFull]==meta['LUT']['TIPSY']['FIZ']['C']:
                    ivl_pi=250
                else:
                    ivl_pi=125
    
                # Timing of transition between pre-industrial and modern eras
                YearRef=dmec[iStandFull]['Year'][0]
                AgeRef=120
        
                YrRegCyc=np.arange(YearRef-AgeRef-100*ivl_pi,YearRef-AgeRef+ivl_pi,ivl_pi)
                Year=YrRegCyc[np.where(YrRegCyc>=meta['Year'][0])[0]]
                ID_Type=meta['LUT']['Dist']['Wildfire']*np.ones(Year.size)
                MortF=100*np.ones(Year.size)
                GrowthF=0*np.ones(Year.size)
                ID_GrowthCurve=1*np.ones(Year.size)
                ec=cbu.CompileEvents(ec,tv,iS,ID_Type,Year,MortF,GrowthF,ID_GrowthCurve)
          
                #--------------------------------------------------------------
                # Add modern era events
                #--------------------------------------------------------------
                
                ind=np.where(dmec[iStandFull]['ScnAffected'][iScn]==1)[0]
                if ind.size>0:
                    ID_Type=dmec[iStandFull]['ID_Type'][ind]
                    Year=dmec[iStandFull]['Year'][ind]
                    MortF=dmec[iStandFull]['MortalityFactor'][ind]
                    GrowthF=dmec[iStandFull]['GrowthFactor'][ind]
                    ID_GrowthCurve=dmec[iStandFull]['ID_GC'][iScn][ind]
                    ec=cbu.CompileEvents(ec,tv,iS,ID_Type,Year,MortF,GrowthF,ID_GrowthCurve)
            
            #------------------------------------------------------------------
            # Compress by indexing into the elements with information
            #------------------------------------------------------------------
            
            ec['idx']=np.where(ec['ID_Type']>0)
            ec['ID_Type']=ec['ID_Type'][ec['idx']]
            ec['MortalityFactor']=ec['MortalityFactor'][ec['idx']]
            ec['GrowthFactor']=ec['GrowthFactor'][ec['idx']]
            ec['ID_GrowthCurve']=ec['ID_GrowthCurve'][ec['idx']]
            
            #------------------------------------------------------------------    
            # Save to file
            #------------------------------------------------------------------
            
            gu.opickle(meta['Paths']['Input Scenario'][iScn] + '\\Events_Ens' + cbu.FixFileNum(iEns) + '_Bat' + cbu.FixFileNum(iBat) + '.pkl',ec)

#%% Prepare growth curves

t0=time.time()
cbu.PrepGrowthCurvesUniqueForCBR(meta,ugc)
print((time.time()-t0)/60)

#%% Save metadata

gu.opickle(meta['Paths']['Project'] + '\\Inputs\\Metadata.pkl',meta)

#%% Delete output files

for iScn in range(meta['N Scenario']):
    files=glob.glob(meta['Paths']['Output Scenario'][iScn] + '\\*')
    for f in files:
        os.remove(f)

#%% Run simulation

cbr.MeepMeep(meta)

#%% Save everything

flg=0
if flg==1:
    gu.opickle(meta['Paths']['Project'] + '\\Inputs\\dmec.pkl',dmec)
    gu.opickle(meta['Paths']['Project'] + '\\Inputs\\gcd.pkl',gcd)
    


