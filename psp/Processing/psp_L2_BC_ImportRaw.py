
'''
Site ID 6004993 is missing from header
No site series

'''

#%% Import modules

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import gc as garc
import time
import matplotlib.pyplot as plt
from fcgadgets.macgyver import utilities_general as gu
from fcexplore.psp.Processing import psp_utilities as utl

#%% Import project info

meta={}
meta['Paths']={}
meta['Paths']['DB']=r'C:\Users\rhember\Documents\Data\GroundPlots\PSP-NADB2'
meta['Paths']['Raw Data']=meta['Paths']['DB'] + '\\Given\BC\Received 2022-10-13'

meta=utl.ImportParameters(meta)

#%% Import data

hd0=gu.ReadExcel(meta['Paths']['Raw Data'] + '\\faib_header.xlsx')
pl0=gu.ReadExcel(meta['Paths']['Raw Data'] + '\\faib_sample_byvisit.xlsx')
tl0=gu.ReadExcel(meta['Paths']['Raw Data'] + '\\faib_tree_detail.xlsx')

#%% Processing of input variables

# Get a unique list of BEC zones to create LUT
flg=0
if flg==1:
    BEC_ZONE=hd0['BEC_ZONE'].copy().astype(object)
    for i in range(BEC_ZONE.size):
        BEC_ZONE[i]=BEC_ZONE[i] + hd0['BEC_SBZ'][i]
    LUT_BEC_ZONE=np.unique(BEC_ZONE)

#%% Compile plot level info
# Each entry is a visit

N=pl0['site_identifier'].size

pl={}
pl['ID Plot']=pl0['site_identifier']
pl['ID Visit']=pl0['visit_number']
pl['Area (ha)']=np.zeros(N,dtype=float)
pl['Year']=np.zeros(N,dtype=int)
pl['Month']=np.zeros(N,dtype=int)
pl['Lat']=np.zeros(N,dtype=float)
pl['Lon']=np.zeros(N,dtype=float)
pl['X']=np.zeros(N,dtype=float)
pl['Y']=np.zeros(N,dtype=float)
pl['Ecozone BC L1']=np.zeros(N,dtype=int)
pl['Ecozone BC L2']=np.zeros(N,dtype=int)
pl['Site Series']=np.zeros(N,dtype=int)
pl['Plot Type']=np.zeros(N,dtype=int)
pl['Age']=np.zeros(N,dtype=float)
pl['SI']=np.zeros(N,dtype=float)
pl['Elev']=np.zeros(N,dtype=float)
pl['Slope']=np.zeros(N,dtype=float)
pl['Aspect']=np.zeros(N,dtype=float)

# u=np.unique(pl0['site_identifier'])
# for i in range(u.size):
#     ind=np.where(pl0['site_identifier']==u[i])[0]
#     pl['ID Plot'][ind]=i

for i in range(pl0['MEAS_DT'].size):
    pl['Year'][i]=pl0['MEAS_DT'][i].year
    pl['Month'][i]=pl0['MEAS_DT'][i].month

for i in range(N):
    ind0=np.where(hd0['site_identifier']==pl['ID Plot'][i])[0]
    if ind0.size==0:
        continue
    pl['Lat'][i]=hd0['Latitude'][ind0]
    pl['Lon'][i]=hd0['Longitude'][ind0]
    pl['X'][i]=hd0['BC_ALBERS_X'][ind0]
    pl['Y'][i]=hd0['BC_ALBERS_Y'][ind0]

    pl['Ecozone BC L1'][i]=meta['LUT']['Ecozone BC L1'][ hd0['BEC_ZONE'][ind0[0]] ]

    pl['Ecozone BC L2'][i]=meta['LUT']['Ecozone BC L2'][ hd0['BEC_ZONE'][ind0[0]] + hd0['BEC_SBZ'][ind0[0]]  ]

u=np.unique(pl0['sampletype'])
for i in range(u.size):
    ind0=np.where( pl0['sampletype'] )[0]
    ind1=np.where( (meta['LUT Tables']['Plot Type']['Given']==u[i]) & (meta['LUT Tables']['Plot Type']['Jurisdiction']=='BC') )[0]
    pl['Plot Type'][ind0]=meta['LUT Tables']['Plot Type']['Value'][ind1]

for i in range(pl0['MEAS_DT'].size):
    pl['Age'][i]=pl0['PROJ_AGE_1'][i]

#%% Compile tree-level data

N_tl=tl0['site_identifier'].size

tl={}
tl['ID Plot']=tl0['site_identifier'].astype('int32')
tl['ID Visit']=tl0['visit_number'].astype(int)
tl['ID Tree']=tl0['TREE_NO'].astype(int)
tl['ID Tree Unique To Jurisdiction']=np.zeros(N_tl,dtype=int)
tl['ID Species']=np.zeros(N_tl,dtype=int)
tl['Vital Status']=np.zeros(N_tl,dtype=int)
tl['DBH']=tl0['DBH']
tl['H']=tl0['HEIGHT']
tl['H Obs']=tl0['HEIGHT']
tl['AEF']=tl0['PHF_TREE']

u=np.unique( np.column_stack( (tl['ID Plot'],tl['ID Tree']) ),axis=0 )
for i in range(u.shape[0]):
    ind=np.where( (tl['ID Plot']==u[i,0]) & (tl['ID Tree']==u[i,1]) )[0]
    tl['ID Tree Unique To Jurisdiction'][i]=i

ind=np.where(tl0['LV_D']=='L')[0]
tl['Vital Status'][ind]=meta['LUT']['Vital Status']['Live']
ind=np.where(tl0['LV_D']=='D')[0]
tl['Vital Status'][ind]=meta['LUT']['Vital Status']['Dead']

u=np.unique(tl0['SPECIES'])
for i in range(u.size):
    ind0=np.where(tl0['SPECIES']==u[i])[0]
    ind1=np.where(meta['Species']['BC']['Code Given']==u[i])[0][0]
    ind2=np.where(meta['Allo B']['Code']==meta['Species']['BC']['Code Final'][ind1] )[0]
    tl['ID Species'][ind0]=meta['Allo B']['ID'][ind2]

tl['Crown Class']=np.zeros(N_tl,dtype=int)
for k in meta['LUT']['Crown Class'].keys():
    ind=np.where(tl0['CR_CL']==k)[0]
    tl['Crown Class'][ind]=meta['LUT']['Crown Class'][k]

tl['ID DA1']=np.zeros(N_tl,dtype=int)
tl['ID DA2']=np.zeros(N_tl,dtype=int)
tl['ID DA3']=np.zeros(N_tl,dtype=int)
u=np.unique(tl0['DAM_AGNA'])
for i in range(u.size):
    if u[i]=='nan':
        continue
    ind0=np.where(tl0['DAM_AGNA']==u[i])[0]
    ind1=np.where(meta['LUT Tables']['DA BC']['Given']==u[i])[0][0]
    tl['ID DA1'][ind0]=meta['LUT Tables']['DA BC']['Final'][ind1]
u=np.unique(tl0['DAM_AGNB'])
for i in range(u.size):
    if u[i]=='nan':
        continue
    ind0=np.where(tl0['DAM_AGNB']==u[i])[0]
    ind1=np.where(meta['LUT Tables']['DA BC']['Given']==u[i])[0]
    if ind1.size==0:
        continue
    tl['ID DA2'][ind0]=meta['LUT Tables']['DA BC']['Final'][ind1[0]]
u=np.unique(tl0['DAM_AGNC'])
for i in range(u.size):
    if u[i]=='nan':
        continue
    ind0=np.where(tl0['DAM_AGNC']==u[i])[0]
    ind1=np.where(meta['LUT Tables']['DA BC']['Given']==u[i])[0]
    if ind1.size==0:
        continue
    tl['ID DA3'][ind0]=meta['LUT Tables']['DA BC']['Final'][ind1[0]]

tl['Age']=-999*np.ones(N_tl,dtype=float)
ind=np.where(tl0['AGE_TOT']>0)[0]
tl['Age'][ind]=tl0['AGE_TOT'][ind]

#%% Save

d={'pl':pl,'tl':tl}
gu.opickle(r'C:\Users\rhember\Documents\Data\GroundPlots\PSP-NADB2\Processed\L1\L1_BC.pkl',d)

#%%
