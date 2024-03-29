
'''
Site ID 6004993 is missing from header
No site series
Where does the whole-stem volume calculation come from?
VRI remeasurements:
    1) Why are there remeasurements
    2) Tree ID does not align in the Integrated Plot Centre, which leads to bogus p
    growth and mortlaity. I assume they are not meant to be used for growth, but
    why do they get remeasured?

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
from scipy.optimize import curve_fit

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
pl['Num Plots']=np.zeros(N,dtype=int)
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
    ind0=np.where( pl0['sampletype']==u[i] )[0]
    ind1=np.where( (meta['LUT Tables']['Plot Type']['Jurisdiction']=='BC') & (meta['LUT Tables']['Plot Type']['Given']==u[i]) )[0]
    pl['Plot Type'][ind0]=meta['LUT Tables']['Plot Type']['Value'][ind1[0]]

for i in range(pl0['MEAS_DT'].size):
    pl['Num Plots'][i]=pl0['NO_PLOTS'][i]
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
tl['Vws']=tl0['VOL_WSV']
tl['Vntwb']=tl0['VOL_NTWB']
tl['AEF']=tl0['PHF_TREE']
tl['Flag WithinPlot']=np.ones(N_tl,dtype=int)

ind=np.where(tl0['MEAS_INTENSE']=='OUT_OF_PLOT')[0]
tl['Flag WithinPlot'][ind]=0

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

# Mountain pine beetle
tl['Flag IBM']=np.zeros(N_tl,dtype=int)
ind0=np.where( (tl0['DAM_AGNA']=='IBM') | (tl0['DAM_AGNB']=='IBM') | (tl0['DAM_AGNC']=='IBM') )[0]
tl['Flag IBM'][ind0]=1

# Western spruce budworm
tl['Flag IDW']=np.zeros(N_tl,dtype=int)
ind0=np.where( (tl0['DAM_AGNA']=='IDW') | (tl0['DAM_AGNB']=='IDW') | (tl0['DAM_AGNC']=='IDW') )[0]
tl['Flag IDW'][ind0]=1

tl['Age']=-999*np.ones(N_tl,dtype=float)
ind=np.where(tl0['AGE_TOT']>0)[0]
tl['Age'][ind]=tl0['AGE_TOT'][ind]

tl['Stature']=np.zeros(N_tl,dtype=int)
ind=np.where(tl0['s_f']=='S')[0]
tl['Stature'][ind]=meta['LUT']['Stature']['Standing']
ind=np.where(tl0['s_f']=='F')[0]
tl['Stature'][ind]=meta['LUT']['Stature']['Fallen']

#%% Gap-fill heights

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

uS=np.unique(tl['ID Species'])
rs=[None]*uS.size
for iS in range(uS.size):
    ikp=np.where( (tl['ID Species']==uS[iS]) & (tl['DBH']>0) & (tl['H']>0) & (tl['Vital Status']==meta['LUT']['Vital Status']['Live']) & (tl['ID DA1']==meta['LUT']['Damage Agents']['None']) & (tl['Stature']==meta['LUT']['Stature']['Standing']) )[0]
    x=tl['DBH'][ikp]
    y=tl['H'][ikp]
    #plt.plot(x,y,'b.')

    indGF=np.where( (tl['ID Species']==uS[iS]) & (tl['DBH']>0) & (np.isnan(tl['H'])==True) |
                   (tl['ID Species']==uS[iS]) & (tl['DBH']>0) & (tl['H']==0) )[0]

    try:
        popt0=[26,0.1,0.66,2]
        popt,pcov=curve_fit(fun,x,y,popt0)
        yhat=fun(xhat,popt[0],popt[1],popt[2],popt[3])
        plt.plot(xhat,yhat,'g-')
        rs0,txt=gu.GetRegStats(x,y)
        rs[iS]=rs0
        tl['H'][indGF]=np.maximum(0.0,fun(tl['DBH'][indGF],popt[0],popt[1],popt[2],popt[3]))
    except:
        tl['H'][indGF]=np.maximum(0.0,fun(tl['DBH'][indGF],popt_glob[0],popt_glob[1],popt_glob[2],popt_glob[3]))


#%% Save

d={'pl':pl,'tl':tl}
gu.opickle(r'C:\Users\rhember\Documents\Data\GroundPlots\PSP-NADB2\Processed\L1\L1_BC.pkl',d)

#%%
