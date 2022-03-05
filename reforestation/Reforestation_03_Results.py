
''' REFORESTATION (ALL) - RESULTS '''

#%% Prepare session

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd
import openpyxl
import statistics
from scipy import stats
from numpy import matlib as mb
import gc
import fiona
import time

from fcgadgets.macgyver import utilities_general as gu
from fcgadgets.macgyver import utilities_gis as gis
from fcgadgets.macgyver import utilities_inventory as invu
from fcgadgets.cbrunner import cbrun_utilities as cbu
from fcgadgets.cbrunner import cbrun as cbr

#%% Set figure properties

fs=8
params={'font.sans-serif':'Arial','font.size':fs,'axes.edgecolor':'black','axes.labelsize':fs,'axes.labelcolor':'black','axes.titlesize':fs,'axes.linewidth':0.5,'lines.linewidth':0.5,
        'text.color':'black','xtick.color':'black','xtick.labelsize':fs,'xtick.major.width':0.5,'xtick.major.size':3,'xtick.direction':'in','ytick.color':'black','ytick.labelsize':fs,
        'ytick.major.width':0.5,'ytick.major.size':3,'ytick.direction':'in','legend.fontsize':fs,'savefig.dpi':300,'savefig.transparent':True}
plt.rcParams.update(params)

#%% Import project configuration and inputs

# Configuration
meta=gu.ipickle(r'D:\Data\FCI_Projects\SummaryReforestation\Inputs\Metadata.pkl')

# This project is not actually running the model so it never populates the model code path - do it manually
meta['Paths']['Model Code']=r'C:\Users\rhember\Documents\Code_Python\fcgadgets\cbrunner'

# Define path for figures
#meta['Paths']['Figures']=r'C:\Users\rhember\OneDrive - Government of BC\Figures\FCI Project Summary'

# Time series of saved results
#tv=np.arange(meta['Year Start Saving'],meta['Year End']+1,1)

# Import multipolygons
atu_multipolygons=gu.ipickle(meta['Paths']['Geospatial'] + '\\atu_multipolygons.pkl')

# Import sxy
sxy=gu.ipickle(meta['Paths']['Geospatial'] + '\\sxy.pkl')

# Unique multipolygons in the modelling
uMP=np.unique(sxy['ID_atu_multipolygons'])

# Inventory
sxy,atu,op,burnsev,vri,fcinv,fcsilv,fcres,pl,fire,pest,cut,lul,park,ogmal=invu.LoadSparseGeospatialInputs(meta)

# Age at fertilization
#AgeAtFert=gu.ipickle(meta['Paths']['Project'] + '\\Inputs\\AgeAtFert.pkl')

# Import best-available inventory
ba=gu.ipickle(meta['Paths']['Project'] + '\\Inputs\\ba.pkl')

# Import best-available inventory
dmec=gu.ipickle(meta['Paths']['Project'] + '\\Inputs\\dmec.pkl')

# Added because Im running results before running model
meta=invu.Load_LUTs(meta)

# Import fcgadgets parameters
par=invu.Load_Params(meta)

# Import MP subsampling parameters
ail=gu.ipickle(r'D:\Data\FCI_Projects\SummaryReforestation\Inputs\AnnualImplementationLevel.pkl')

# Area expansion factor
AEF=1.0/ail['samp_rate_mp']

# Generate a list of project IDs for each unique multipolygon
uP_ByMP=[None]*uMP.size
for iMP in range(uMP.size):
    uP_ByMP[iMP]=atu_multipolygons[uMP[iMP]]['FIA_PROJECT_ID']
uP_ByMP=np.array(uP_ByMP)

# Generate a list of project IDs for each unique openings
uO_ByMP=[None]*uMP.size
for iMP in range(uMP.size):
    uO_ByMP[iMP]=atu_multipolygons[uMP[iMP]]['OPENING_ID']
uO_ByMP=np.array(uO_ByMP)

#%% Import MOS by whole project
# Takes 8 min

flg=0
if flg==1:
    t0=time.time()
    mos=cbu.ModelOutputStats(meta,'On')
    print((time.time()-t0)/60)

#mos=gu.ipickle(meta['Paths']['Project'] + '\\Outputs\\MOS.pkl')

#%% Import MOS by multipolygon
# Takes 11 min

flg=0
if flg==1:
    t0=time.time()
    cbu.MosByMultipolygon(meta,include_area='On')
    print((time.time()-t0)/60)

# Import MOS by multipolygon
#MosByMP=gu.ipickle(meta['Paths']['Project'] + '\\Outputs\MosByMultipolygon.pkl')

#%% Export tabular summaries

flg=0
if flg==1:
    # Export summary of activities by sparse grid cells
    dS=invu.ExportSummaryActivities_BySXY(meta,par,atu_multipolygons,sxy,atu,op,burnsev,vri,fcinv,fcsilv,fcres,pl,fire,pest,cut,lul,park,ogmal)

flg=0
if flg==1:
    dMP=invu.ExportSummaryAttributes_ByMP()

flg=0
if flg==1:
    d=invu.ExportSummaryAttributes_BySXY()

flg=0
if flg==1:
    invu.ExportSummaryEvents_ByMP()

flg=0
if flg==1:
    invu.ExportSummaryAttributes_AreaWeighted_ByActivityType()

#%% Plot time series
    
flg=0
if flg==1:
    iB=0
    iP=1
    it=np.where( (tv>=1901) & (tv<=2100) )[0]
    invu.PlotTimeSeries_ByPP()

#%% How much of the 17-18 wildfire areas have been harvested and/or planted

# I ran a query of the burn severity rating layer using the utilities_explore_gdb.py to find that
A_burned=1924393 # ha burned in 17-18.

# Indicator for activities that occur after July 2017
atu['Time Condition ATU']=np.zeros(atu['Year'].size)
ind=np.where( (atu['Year']==2017) & (atu['Month']>=7) | (atu['Year']>=2018) )[0]
atu['Time Condition ATU'][ind]=1

op['Time Condition 1']=np.zeros(op['Year_Denu1_Comp'].size)
ind=np.where( (op['Year_Denu1_Comp']==2017) & (op['Month_Denu1_Comp']>=7) | (op['Year_Denu1_Comp']>=2018) )[0]
op['Time Condition 1'][ind]=1

op['Time Condition 2']=np.zeros(op['Year_Denu2_Comp'].size)
ind=np.where( (op['Year_Denu2_Comp']==2017) & (op['Month_Denu2_Comp']>=7) | (op['Year_Denu2_Comp']>=2018) )[0]
op['Time Condition 2'][ind]=1

# Loop through unique openings
uOP=np.unique(op['OPENING_ID'])
d={}
d['Area PL']=np.zeros(uOP.size)
d['Area SP']=np.zeros(uOP.size)
d['Area Harvest']=np.zeros(uOP.size)
d['Area Harvest RSLT']=np.zeros(uOP.size)
d['FSC PL']=np.array(['' for _ in range(uOP.size)],dtype=object)
d['Flg Dist L RSLT']=np.zeros(uOP.size)
d['Flg Dist R RSLT']=np.zeros(uOP.size)
d['Flg Dist S RSLT']=np.zeros(uOP.size)
for iOP in range(uOP.size):
    
    # Index to unique ATUs
    ind0=np.where( (op['OPENING_ID']==uOP[iOP]) )[0]
    idx0=op['IdxToSXY'][ind0]

    # Did it burn?
    c,ia,ib=np.intersect1d(idx0,fire['IdxToSXY'],return_indices=True)
    if ib.size==0:
        continue
    ind1=np.where( (fire['FIRE_YEAR'][ib]>=2017) & (fire['FIRE_YEAR'][ib]<=2018) )[0]
    if ind1.size==0:
        continue
 
    # Planting
    ind1=np.where( (atu['OPENING_ID']==uOP[iOP]) & \
                  (atu['SILV_BASE_CODE']==meta['LUT']['ATU']['SILV_BASE_CODE']['PL'][0]) & \
                  (atu['SILV_TECHNIQUE_CODE']!=meta['LUT']['ATU']['SILV_TECHNIQUE_CODE']['SE'][0]) & \
                  (atu['SILV_TECHNIQUE_CODE']!=meta['LUT']['ATU']['SILV_TECHNIQUE_CODE']['CG'][0]) & \
                  (atu['SILV_METHOD_CODE']!=meta['LUT']['ATU']['SILV_METHOD_CODE']['LAYOT'][0]) & \
                  (atu['RESULTS_IND']=='Y') & \
                  (atu['ACTUAL_TREATMENT_AREA']!=None) & \
                  (atu['Year']!=None) & \
                  (atu['Time Condition']==1) )[0]
    if ind1.size>0:
        ind1=ind1[0]
        #print(ind1.size)
        d['Area PL'][iOP]=d['Area PL'][iOP]+np.sum(atu['ACTUAL_TREATMENT_AREA'][ind1])
        #d['Area PL'][iOP]=np.sum(atu['ACTUAL_TREATMENT_AREA'][ind1])
        d['FSC PL'][iOP]=cbu.lut_n2s(meta['LUT']['ATU']['SILV_FUND_SOURCE_CODE'],atu['SILV_FUND_SOURCE_CODE'][ind1])[0]
    
    # Site prep
    ind1=np.where( (atu['OPENING_ID']==uOP[iOP]) & (atu['SILV_BASE_CODE']==meta['LUT']['ATU']['SILV_BASE_CODE']['SP'][0]) & (atu['Time Condition']==1) )[0]
    if ind1.size>0:
        d['Area SP'][iOP]=d['Area SP'][iOP]+np.sum(atu['ACTUAL_TREATMENT_AREA'][ind1])
    
    # Harvest from results
    ind=np.where( (op['DENUDATION_1_DISTURBANCE_CODE'][ind0]==meta['LUT']['OP']['DENUDATION_1_DISTURBANCE_CODE']['L']) & (op['Time Condition 1'][ind0]==1) | (op['DENUDATION_2_DISTURBANCE_CODE'][ind0]==meta['LUT']['OP']['DENUDATION_2_DISTURBANCE_CODE']['L']) & (op['Time Condition 2'][ind0]==1) )[0]
    if ind.size>0:
        d['Flg Dist L RSLT'][iOP]=1
        d['Area Harvest RSLT'][iOP]=op['OPENING_GROSS_AREA'][ind0[ind[0]]]
    ind=np.where( (op['DENUDATION_1_DISTURBANCE_CODE'][ind0]==meta['LUT']['OP']['DENUDATION_1_DISTURBANCE_CODE']['S']) & (op['Time Condition 1'][ind0]==1) | (op['DENUDATION_2_DISTURBANCE_CODE'][ind0]==meta['LUT']['OP']['DENUDATION_2_DISTURBANCE_CODE']['S']) & (op['Time Condition 2'][ind0]==1) )[0]
    if ind.size>0:
        d['Flg Dist S RSLT'][iOP]=1
        d['Area Harvest RSLT'][iOP]=op['OPENING_GROSS_AREA'][ind0[ind[0]]]
    ind=np.where( (op['DENUDATION_1_DISTURBANCE_CODE'][ind0]==meta['LUT']['OP']['DENUDATION_1_DISTURBANCE_CODE']['R']) & (op['Time Condition 1'][ind0]==1) | (op['DENUDATION_2_DISTURBANCE_CODE'][ind0]==meta['LUT']['OP']['DENUDATION_2_DISTURBANCE_CODE']['R']) & (op['Time Condition 2'][ind0]==1) )[0]
    if ind.size>0:
        d['Flg Dist R RSLT'][iOP]=1
        d['Area Harvest RSLT'][iOP]=op['OPENING_GROSS_AREA'][ind0[ind[0]]]
    
    # Harvest area
    indH=np.where( (cut['HARVEST_YEAR']>=2017) & (cut['OPENING_ID']==uOP[iOP]) )[0]
    if indH.size>0:
        d['Area Harvest'][iOP]=np.sum(np.unique(cut['AREA_HA'][indH]))

# Summarize in dictionary
dS={}
dS['FSC']=np.unique(d['FSC PL'])
dS['A PL']=np.zeros(dS['FSC'].size)
dS['A H']=np.zeros(dS['FSC'].size)
dS['A H RSLT']=np.zeros(dS['FSC'].size)
for i in range(dS['FSC'].size):
    ind=np.where(d['FSC PL']==dS['FSC'][i])[0]
    dS['A PL'][i]=AEF*np.sum(d['Area PL'][ind])
    dS['A H'][i]=AEF*np.sum(d['Area Harvest'][ind])
    dS['A H RSLT'][i]=AEF*np.sum(d['Area Harvest RSLT'][ind])

# Export to dataframe
df=pd.DataFrame(dS)
df

#%% Inspect growth curve parameters for multipolygon 
# The resulting array will show growth curves 1, 2, ...

flg=0
if flg==1:

    gc=gu.ipickle(meta['Paths']['Project'] + '\\Inputs\\gc.pkl')

    iMP=551
    
    indMP=np.where(sxy['ID_atu_multipolygons']==iMP)[0]

    gc[indMP[0]][iP]


#%% Plot time series for each project by multipolygon
# Takes a long time
    
flg=0
if flg==1:
    iB=0
    iP=1
    ivlMP=1
    iScnForArea=1
    ivlT=5
    iScn=0
    it=np.where( (tv>=1901) & (tv<=2100) )[0]

    cbu.QA_Plot_ByMultiPolygon(meta,uMP,ivlMP,iScnForArea,ivlT,tv,it,MosByMP,iB,iP)

#%% Export summary by grid cell

flg=0
if flg==1:
    project_name='FCI'
    include_planting='Off'
    invu.ExportSummaryByGridCell(meta,atu_multipolygons,dAdmin,sxy,atu,fcinv,vri,pl,op,include_planting,project_name)
    
#%% Summary by BGC
# 68% of completed LCELF funds (excluding FFT) are CWH

def Summary_ByBGC():
    
    uAT=np.unique(dMP['Project Type'])
    
    for iAT in range(uAT.size):
        
        indAT=np.where(dMP['Project Type']==uAT[iAT])[0]
        
        # Unique BGC classes
        uC=np.unique(np.array([dMP['BGCz'][indAT],dMP['BGCsz'][indAT],dS['BGCv'][indAT]]).T,axis=0)
        
        for iC in range(uC.size):
            ind=np.where( (dMP['Project Type']==uAT[iAT]) & (dMP['BGCz'][indAT],dMP['BGCsz'][indAT],dS['BGCv'][indAT]) )[0]
    
    ind=np.where( (dMP['SILV_BASE_CODE']=='FE') )[0]
    uBGC=np.unique(dMP['BGC Zone Mode'][ind])

    A=np.zeros(uBGC.size)
    for i in range(uBGC.size):
        ind=np.where( (dMP['BGC Zone Mode']==uBGC[i]) & (dMP['SILV_BASE_CODE']=='FE') )[0]
        A[i]=np.sum(dMP['ACTUAL_TREATMENT_AREA'][ind])

    ind=np.where( (uBGC=='CWH') )[0]
    A_cwh=np.sum(A[ind])
    A_tot=np.sum(A)
    A_cwh/A_tot*100

    # Area-weighted site index = 26 (May 3, 2021)
    ind=np.where( (dMP['SILV_BASE_CODE']=='FE') )[0]
    SI_w=np.sum(dMP['ACTUAL_TREATMENT_AREA'][ind]*dMP['SI_ba_mean'][ind])/np.sum(dMP['ACTUAL_TREATMENT_AREA'][ind])
    SI_w
    
    return




#%% Query planting projects by BGC
#
#dS=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\FCI_Projects\FCI_RollupFCI_Inv\Inputs\SummaryBySparseGridCell.xlsx')
#
#ind1=np.where( (dS['FSC']=='FCE') & (dS['SBC']=='PL') | (dS['FSC']=='FCM') & (dS['SBC']=='PL') )[0]
#
#u=np.unique(np.array([dS['BGCz'][ind1],dS['BGCsz'][ind1],dS['BGCv'][ind1]]).T,axis=0)
#
#N=np.zeros(u.shape[0])
#for i in range(u.shape[0]):
#    ind2=np.where( (dS['BGCz'][ind1]==u[i,0]) & (dS['BGCsz'][ind1]==u[i,1]) & (dS['BGCv'][ind1]==u[i,2]) )[0]    
#    N[i]=np.round(np.sum(dS['AEF_ATU'][ind1[ind2]]),2)
#sts=np.column_stack([u,N])
#df=pd.DataFrame(sts)
#df.to_excel(r'C:\Users\rhember\Documents\Data\FCI_Projects\FCI_RollupFCI_Inv\Inputs\SummaryByBGCZone.xlsx')
#
## Nutrient management
#
#ind1=np.where( (dS['FSC']=='FCE') & (dS['SBC']=='FE') | (dS['FSC']=='FCM') & (dS['SBC']=='FE') )[0]
#
#u=np.unique(np.array([dS['BGCz'][ind1],dS['BGCsz'][ind1],dS['BGCv'][ind1]]).T,axis=0)
#
#N=np.zeros(u.shape[0])
#for i in range(u.shape[0]):
#    ind2=np.where( (dS['BGCz'][ind1]==u[i,0]) & (dS['BGCsz'][ind1]==u[i,1]) & (dS['BGCv'][ind1]==u[i,2]) )[0]    
#    N[i]=np.round(np.sum(dS['AEF_ATU'][ind1[ind2]]),2)
#sts=np.column_stack([u,N])
#df=pd.DataFrame(sts)
#df.to_excel(r'C:\Users\rhember\Documents\Data\FCI_Projects\FCI_RollupFCI_Inv\Inputs\SummaryByBGCZone_NutrientManagement.xlsx')
