
'''
BC FOREST CARBON SUMMARY - NUTRIENT MANAGEMENT FUTURE
'''

#%% Import Python modules

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
import warnings
import fcgadgets.macgyver.utilities_general as gu
import fcgadgets.macgyver.utilities_gis as gis
import fcgadgets.macgyver.utilities_inventory as invu
import fcgadgets.cbrunner.cbrun_utilities as cbu
import fcgadgets.cbrunner.cbrun as cbr
warnings.filterwarnings("ignore")

#%% Initiate metadata structure and import paths

meta=gu.ipickle(r'D:\Data\FCI_Projects\BCFCS_NM_20k_Future\Inputs\Metadata.pkl')

flg=0
if flg==1:
    meta['Paths']['Project']=r'D:\Data\FCI_Projects\BCFCS_NM_20k_Future'
    gu.opickle(r'D:\Data\FCI_Projects\BCFCS_NM_20k_Future\Inputs\Metadata.pkl',meta)

#%% Project adjustments

#meta['Project']['Flag Tracking Projects']=0

#%% Import inventory

geos,atu,op,burnsev,vri,fcinv,fcsilv,fcres,pl,fire,pest,cut,lul,park,ogmal,ogsr=invu.LoadSparseGeospatialInputs(meta)
idx=invu.LoadSparseGridIndex(meta)

#%% Import project configuration

meta=cbu.ImportProjectConfig(meta,geos=geos)

#%% Index to baseline scenarios

meta['Project']['Baseline Indices']=np.array([0,2])
meta['Project']['Actual Indices']=np.array([1,3])

#%% Process project 1

meta,dmec,ba,lsc=invu.ProcessProjectInputs1(meta,geos,atu,op,burnsev,vri,fcinv,fcsilv,fcres,pl,fire,pest,cut,lul,park,ogmal,ogsr,idx)

#%% Process project 2

gc,ugc,dmec,ba=invu.ProcessProjectInputs2(meta,ba,dmec)

#%% Run BatchTIPSY

# ************************* MANUAL OPERATION **********************************
#------------------------------------------------------------------------------
# Ensure BatchTIPSY.exe config (.BPS) file has up-to-date paths
# Run BatchTIPSY.exe
#------------------------------------------------------------------------------
# ************************* MANUAL OPERATION **********************************

#%% Define strata (for analyzing results)

meta['Project']['Strata']={}
meta['Project']['Strata']['Project']={}
meta['Project']['Strata']['Project']['Unique CD']=np.array(['All'],dtype=object)
meta['Project']['Strata']['Project']['Unique ID']=np.ones(1)
meta['Project']['Strata']['Project']['ID']=np.ones(meta['Project']['N Stand'],dtype=int)
meta['Project']['Strata']['Spatial']={}
meta['Project']['Strata']['Spatial']['Unique CD']=np.array(['All'],dtype=object)
meta['Project']['Strata']['Spatial']['Unique ID']=np.ones(1)
meta['Project']['Strata']['Spatial']['ID']=np.ones(meta['Project']['N Stand'],dtype=int)

flg=0
if flg==1:
    # By last year of implementation
    t=np.arange(1970,2022,1)
    meta['Project']['Strata']['Project']['Unique ID']=np.append(meta['Project']['Strata']['Project']['Unique ID'],t)
    meta['Project']['Strata']['Project']['Unique CD']=np.append(meta['Project']['Strata']['Project']['Unique CD'],t.astype(str))

    iScn=1
    for iStand in range(meta['Project']['N Stand']):
        ind=np.where(dmec[iScn][iStand]['ID_Type']==meta['LUT']['Dist']['Fertilization Aerial'])[0]
        if ind.size>0:
            Year0=np.floor(dmec[iScn][iStand]['Year'][ind[-1]])
            iT=np.where(t==Year0)[0]
            if iT.size==0:
                continue
            meta['Project']['Strata']['Project']['ID'][iStand]=t[iT]

#%% Prepare project inputs 3

meta,dmec,ba,thlb=invu.ProcessProjectInputs3(meta,geos,atu,op,burnsev,vri,fcinv,fcsilv,fcres,pl,fire,pest,cut,lul,park,ogmal,ogsr,idx,ba,dmec,lsc,gc,ugc)

#%% Run simulation

t0=time.time()
meta=cbr.MeepMeep(meta)
print((time.time()-t0)/60)

meta['Project']['Run Time Summary']

#%% If running multiple instances:

flg=0
if flg==1:

    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import time
    import os
    import glob
    import fcgadgets.macgyver.utilities_general as gu
    import fcgadgets.macgyver.utilities_gis as gis
    import fcgadgets.macgyver.utilities_inventory as invu
    import fcgadgets.cbrunner.cbrun_utilities as cbu
    import fcgadgets.cbrunner.cbrun as cbr

    meta=gu.ipickle(r'D:\Data\FCI_Projects\BCFCS_NM_20k_Future\Inputs\Metadata.pkl')
    # cbu.DeleteAllOutputFiles(meta)
    cbr.MeepMeep(meta)

#%% Calculate summaries for future simulations
# Econ takes 33 min, the others are around 10 min

flg=1
if flg==1:
    cbu.Calc_MOS_FromPoints_GHG(meta)
    cbu.Calc_MOS_FromPoints_Econ(meta)
    cbu.Calc_MOS_FromPoints_Area(meta)
    cbu.Calc_MOS_FromPoints_MortalityByAgent(meta)

