
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

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcexplore.field_plots.Processing.psp_util as ugp
import fcgadgets.bc1ha.bc1ha_util as u1ha

#%% Import project info
meta=u1ha.Init()
meta=ugp.ImportParameters(meta)
srs=gis.ImportSRSs()
# meta={}
# meta['Paths']={'GP':{}}
# meta['Paths']['GP']['DB']=r'C:\Users\rhember\Documents\Data\GroundPlots\PSP-NADB2'
# meta['Paths']['GP']['Raw Data']=meta['Paths']['GP']['DB'] + '\\Given\BC\Received 2023-03-02'
# meta=ugp.ImportParameters(meta)

#%% Import data
hd0_unfuzzed=gu.ReadExcel(meta['Paths']['GP']['DB'] + '\\Given\\BC\\Received 2022-10-13\\faib_header.xlsx')
hd0=gu.ReadExcel(meta['Paths']['GP']['Raw Data']['BC'] + '\\faib_header.xlsx')
pl0=gu.ReadExcel(meta['Paths']['GP']['Raw Data']['BC'] + '\\faib_sample_byvisit.xlsx')
tl0=gu.ReadExcel(meta['Paths']['GP']['Raw Data']['BC'] + '\\faib_tree_detail.xlsx')

#%% Unique species counts
# u,cnts=np.unique(tl0['SPECIES'],return_counts=True)

#%% Add felled indicator to tree level data
# Received this LUT - replace when felled indicator is added to new DB
# Most harvested trees not in the tree file -> asking Dan to look into fixing
fel=gu.ReadExcel(meta['Paths']['GP']['DB'] + '\\Given\\BC\\Received 2023-07-12\\Compilation_Out_harv_to_Robbie_2023Jul12.xlsx','Query1')

tl0['Felled']=np.zeros(tl0['PLOT'].size,dtype='int8')
N_Missing=0
for i in range(fel['PLOT'].size):
    ind=np.where( (tl0['CLSTR_ID']==fel['CLSTR_ID'][i]) & (tl0['PLOT']==fel['PLOT'][i]) & (tl0['TREE_NO']==fel['TREE_NO'][i]) )[0]
    if ind.size==0:
        #print('Missing')
        print(fel['CLSTR_ID'][i])
        #break
        N_Missing=N_Missing+1
    tl0['Felled'][ind]=1

ind=np.where(tl0['Felled']==1)[0]
np.sum(tl0['Felled'])

#%% Compile plot level info
# Each entry is a visit

N=pl0['site_identifier'].size

pl={}
pl['ID Plot']=pl0['site_identifier']
pl['ID Visit']=pl0['visit_number']
pl['Area (ha)']=np.zeros(N,dtype=float)
pl['Year']=np.zeros(N,dtype='int16')
pl['Month']=np.zeros(N,dtype='int16')
pl['Lat']=np.zeros(N,dtype=float)
pl['Lon']=np.zeros(N,dtype=float)
pl['X']=np.zeros(N,dtype=float)
pl['Y']=np.zeros(N,dtype=float)
pl['Ecozone BC L1']=np.zeros(N,dtype='int16')
pl['Ecozone BC L2']=np.zeros(N,dtype='int16')
pl['Site Series']=np.zeros(N,dtype='int16')
pl['Plot Type']=np.zeros(N,dtype='int16')
pl['Num Plots']=np.zeros(N,dtype='int16')
pl['Age VRI']=np.zeros(N,dtype=float)
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
    pl['X'][i]=hd0['BC_ALBERS_X'][ind0]
    pl['Y'][i]=hd0['BC_ALBERS_Y'][ind0]

    pl['Ecozone BC L1'][i]=meta['LUT']['GP']['Ecozone BC L1'][ hd0['BEC_ZONE'][ind0[0]] ]
    pl['Ecozone BC L2'][i]=meta['LUT']['GP']['Ecozone BC L2'][ hd0['BEC_ZONE'][ind0[0]] + hd0['BEC_SBZ'][ind0[0]]  ]

# Get coordinates from unfuzzed version where possible
for i in range(N):
    ind0=np.where(hd0_unfuzzed['site_identifier']==pl['ID Plot'][i])[0]
    if ind0.size==0:
        print('Missing')
        continue
    pl['X'][i]=hd0_unfuzzed['BC_ALBERS_X'][ind0]
    pl['Y'][i]=hd0_unfuzzed['BC_ALBERS_Y'][ind0]

# Missing lat and long so back it out of ALBERS projected values
pl['Lon'],pl['Lat']=gis.ReprojectCoordinates(srs['String']['BC1ha'],srs['String']['Geographic'],pl['X'],pl['Y'])

u=np.unique(pl0['sampletype'])
for i in range(u.size):
    ind0=np.where( pl0['sampletype']==u[i] )[0]
    ind1=np.where( (meta['Param']['GP']['Raw Tables']['Plot Type']['Territory']=='BC') & (meta['Param']['GP']['Raw Tables']['Plot Type']['Code Given']==u[i]) )[0]
    pl['Plot Type'][ind0]=meta['Param']['GP']['Raw Tables']['Plot Type']['ID'][ind1[0]]

for i in range(pl0['MEAS_DT'].size):
    pl['Num Plots'][i]=pl0['NO_PLOTS'][i]
    pl['Age VRI'][i]=pl0['proj_age_adj'][i]

#%% Compile tree-level data

N_tl=tl0['site_identifier'].size

tl={}
tl['ID Plot']=tl0['site_identifier'].astype('int32')
tl['ID Visit']=tl0['visit_number'].astype('int16')
tl['ID Tree']=tl0['TREE_NO'].astype('int16')
tl['ID Tree Unique To Source']=np.zeros(N_tl,dtype='int16')
tl['ID Species']=np.zeros(N_tl,dtype='int16')
tl['Vital Status']=np.zeros(N_tl,dtype='int16')
tl['Felled']=tl0['Felled']
tl['Age']=tl0['AGE_TOT']
tl['DBH']=tl0['DBH']
tl['H']=tl0['HEIGHT']
tl['H Obs']=tl0['HEIGHT']
tl['Vws']=tl0['VOL_WSV']
tl['Vntwb']=tl0['VOL_NTWB']
tl['AEF']=tl0['PHF_TREE']
tl['Resid']=np.zeros(N_tl,dtype='int16')
tl['Flag WithinPlot']=np.ones(N_tl,dtype='int16')

ind=np.where(tl0['MEAS_INTENSE']=='OUT_OF_PLOT')[0]
tl['Flag WithinPlot'][ind]=0

u=np.unique( np.column_stack( (tl['ID Plot'],tl['ID Tree']) ),axis=0 )
for i in range(u.shape[0]):
    ind=np.where( (tl['ID Plot']==u[i,0]) & (tl['ID Tree']==u[i,1]) )[0]
    tl['ID Tree Unique To Source'][i]=i

ind=np.where(tl0['LV_D']=='L')[0]
tl['Vital Status'][ind]=meta['LUT']['GP']['Vital Status']['Live']
ind=np.where(tl0['LV_D']=='D')[0]
tl['Vital Status'][ind]=meta['LUT']['GP']['Vital Status']['Dead']

#np.unique(tl0['RESIDUAL'])
ind=np.where(tl0['RESIDUAL']=='Y')[0]
tl['Resid'][ind]=1

u=np.unique(tl0['SPECIES'])
for i in range(u.size):
    ind0=np.where(tl0['SPECIES']==u[i])[0]
    ind1=np.where( (meta['Param']['GP']['Raw Tables']['Species Crosswalk']['Source Code']=='BC') & (meta['Param']['GP']['Raw Tables']['Species Crosswalk']['Species Code Given']==u[i]) )[0]
    #meta['LUT']['GP']['Species Given']['BC']['Code Given']==u[i]
    ind2=np.where(meta['Param']['GP']['Raw Tables']['Species']['Code']==meta['Param']['GP']['Raw Tables']['Species Crosswalk']['Species Code Final'][ind1] )[0]
    tl['ID Species'][ind0]=meta['Param']['GP']['Raw Tables']['Species']['ID'][ind2]

tl['Crown Class']=meta['LUT']['GP']['Crown Class']['Unknown']*np.ones(N_tl,dtype='int16')
for i in range(meta['Param']['GP']['Raw Tables']['Crown Class Crosswalk']['ID'].size):
    ind=np.where( tl0['CR_CL']==meta['Param']['GP']['Raw Tables']['Crown Class Crosswalk']['Code Given'][i] )[0]
    tl['Crown Class'][ind]=meta['Param']['GP']['Raw Tables']['Crown Class Crosswalk']['ID'][i]

tl['ID DA1']=np.zeros(N_tl,dtype='int16')
tl['ID DA2']=np.zeros(N_tl,dtype='int16')
tl['ID DA3']=np.zeros(N_tl,dtype='int16')
u=np.unique(tl0['DAM_AGNA'])
for i in range(u.size):
    if u[i]=='nan':
        continue
    ind0=np.where(tl0['DAM_AGNA']==u[i])[0]
    ind1=np.where( (meta['Param']['GP']['Raw Tables']['DA Crosswalk']['Source Code']=='BC') & (meta['Param']['GP']['Raw Tables']['DA Crosswalk']['Code Given']==u[i]) )[0]
    tl['ID DA1'][ind0]=meta['Param']['GP']['Raw Tables']['DA Crosswalk']['ID Final'][ind1]

u=np.unique(tl0['DAM_AGNB'])
for i in range(u.size):
    if u[i]=='nan':
        continue
    ind0=np.where(tl0['DAM_AGNB']==u[i])[0]
    ind1=np.where( (meta['Param']['GP']['Raw Tables']['DA Crosswalk']['Source Code']=='BC') & (meta['Param']['GP']['Raw Tables']['DA Crosswalk']['Code Given']==u[i]) )[0]
    if ind1.size==0:
        continue
    tl['ID DA2'][ind0]=meta['Param']['GP']['Raw Tables']['DA Crosswalk']['ID Final'][ind1]
    
u=np.unique(tl0['DAM_AGNC'])
for i in range(u.size):
    if u[i]=='nan':
        continue
    ind0=np.where(tl0['DAM_AGNC']==u[i])[0]
    ind1=np.where( (meta['Param']['GP']['Raw Tables']['DA Crosswalk']['Source Code']=='BC') & (meta['Param']['GP']['Raw Tables']['DA Crosswalk']['Code Given']==u[i]) )[0]
    if ind1.size==0:
        continue
    tl['ID DA3'][ind0]=meta['Param']['GP']['Raw Tables']['DA Crosswalk']['ID Final'][ind1]

# Mountain pine beetle
tl['Flag IBM']=np.zeros(N_tl,dtype='int16')
ind0=np.where( (tl0['DAM_AGNA']=='IBM') | (tl0['DAM_AGNB']=='IBM') | (tl0['DAM_AGNC']=='IBM') )[0]
tl['Flag IBM'][ind0]=1

# Western spruce budworm
tl['Flag IDW']=np.zeros(N_tl,dtype='int16')
ind0=np.where( (tl0['DAM_AGNA']=='IDW') | (tl0['DAM_AGNB']=='IDW') | (tl0['DAM_AGNC']=='IDW') )[0]
tl['Flag IDW'][ind0]=1

tl['Age']=np.nan*np.ones(N_tl,dtype=float)
ind=np.where(tl0['AGE_TOT']>=0)[0]
tl['Age'][ind]=tl0['AGE_TOT'][ind]

tl['Stature']=meta['LUT']['GP']['Stature']['Unknown']*np.ones(N_tl,dtype='int16')
ind=np.where(tl0['s_f']=='S')[0]
tl['Stature'][ind]=meta['LUT']['GP']['Stature']['Standing']
ind=np.where(tl0['s_f']=='F')[0]
tl['Stature'][ind]=meta['LUT']['GP']['Stature']['Fallen']

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
ikp=np.where( (tl['DBH']>0) & (tl['H']>0) & (tl['Vital Status']==meta['LUT']['GP']['Vital Status']['Live']) & (tl['ID DA1']==meta['LUT']['GP']['Damage Agent']['No damage']) & (tl['Stature']==meta['LUT']['GP']['Stature']['Standing']) )[0]
x=tl['DBH'][ikp]
y=tl['H'][ikp]
popt_glob0=[26,0.1,0.66,2]
popt_glob,pcov=curve_fit(fun,x,y,popt_glob0)
#yhat=fun(xhat,popt_glob[0],popt_glob[1],popt_glob[2],popt_glob[3])
rs_glob,txt=gu.GetRegStats(x,y)

uS=np.unique(tl['ID Species'])
rs=[None]*uS.size
for iS in range(uS.size):
    ikp=np.where( (tl['ID Species']==uS[iS]) & (tl['DBH']>0) & (tl['H']>0) & (tl['Vital Status']==meta['LUT']['GP']['Vital Status']['Live']) & (tl['ID DA1']==meta['LUT']['GP']['Damage Agent']['No damage']) & (tl['Stature']==meta['LUT']['GP']['Stature']['Standing']) )[0]
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
gu.opickle(meta['Paths']['GP']['DB'] + '\\Processed\\L1\\L1_BC.pkl',d)

#%%
