#%% Import modules

from os import listdir
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from fcgadgets.macgyver import utilities_general as gu
import scipy.io
import copy
import geopandas as gpd
import fiona
import scipy.io as spio
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats, linalg
from fcgadgets.macgyver import utilities_inventory as invu
from fcgadgets.macgyver import utilities_query_gdb as qgdb
from fcgadgets.cbrunner import cbrun_utilities as cbu

#%% Graphics parameters

# fs=7
# params={'font.sans-serif':'Arial','font.size':fs,'axes.edgecolor':'black','axes.labelsize':fs,'axes.labelcolor':'black','axes.titlesize':fs,'axes.linewidth':0.5,'lines.linewidth':0.5,
#         'text.color':'black','xtick.color':'black','xtick.labelsize':fs,'xtick.major.width':0.5,'xtick.major.size':3,'xtick.direction':'in','ytick.color':'black','ytick.labelsize':fs,
#         'ytick.major.width':0.5,'ytick.major.size':3,'ytick.direction':'in','legend.fontsize':fs,'savefig.dpi':300,'savefig.transparent':True}
# plt.rcParams.update(params)

gp=gu.SetGraphics('Manuscript')

#%% Prepare

tv=np.arange(1950,2022,1)

# For CRS
gdf_bm=gpd.read_file(r'C:\Users\rhember\Documents\Data\Basemaps\Basemaps.gdb',layer='NRC_POLITICAL_BOUNDARIES_1M_SP')

#%% National Forest Database

d=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\Harvest\Harvest Area\NFD - Area harvested by ownership and harvesting method - EN FR.xlsx')
dNFD={}
dNFD['Area']=np.nan*np.ones(tv.size)
for iT in range(tv.size):
    ind=np.where( (d['Jurisdiction']=='British Columbia') & (d['Year']==tv[iT]) )[0]
    if ind.size==0:
        continue
    dNFD['Area'][iT]=np.nansum(d['Area (hectares)'][ind])

#%% Annual harvest area from consolidated cutblock DB

fin=r'C:\Users\rhember\Documents\Data\ForestInventory\Disturbances\20210401\Disturbances.gdb'

dCC={}
dCC['Year']=tv
dCC['Area Harvested']=np.zeros(tv.size)
with fiona.open(fin,layer='VEG_CONSOLIDATED_CUT_BLOCKS_SP') as source:
    for feat in source:
        if feat['geometry']==None:
            continue
        iT=np.where(tv==feat['properties']['HARVEST_YEAR'])[0]
        if iT.size==0:
            continue
        dCC['Area Harvested'][iT]=dCC['Area Harvested'][iT]+feat['properties']['AREA_HA']
#gu.opickle(r'C:\Users\rhember\Documents\Data\Harvest\Harvest Area\AnnualHarvestAreaFromConCutblocksDB.pkl',d)

#%% Harvest from opening layer

# Import opening layer
metaOP={}
metaOP['Path']=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20220422\Results.gdb'
metaOP['Layer']='RSLT_OPENING_SVW'; # fiona.listlayers(op['Path'])
metaOP['crs']=gdf_bm.crs
metaOP['Keep Geom']='Off'
metaOP['Select Openings']=np.array([])
metaOP['SBC']=np.array([])
metaOP['FSC']=np.array([])
metaOP['ROI']=[]
metaOP['gdf']=qgdb.Query_Openings(metaOP,[])

# Import forest Cover Reserves
metaFCR={}
metaFCR['Path']=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20220422\Results.gdb'
metaFCR['Layer']='RSLT_FOREST_COVER_RESERVE_SVW'; # fiona.listlayers(metaFCR['Path'])
metaFCR['crs']=gdf_bm.crs
metaFCR['Keep Geom']='Off'
metaFCR['Select Openings']=np.array([])
metaFCR['SBC']=np.array([])
metaFCR['FSC']=np.array([])
metaFCR['ROI']=[]
metaFCR['gdf']=qgdb.Query_Openings(metaFCR,[])

# Add reserve areas to opening layer
metaOP['gdf']['Area Reserves']=np.zeros(metaOP['gdf']['OPENING_ID'].size)
for i in range(metaFCR['gdf']['OPENING_ID'].size):
    ind=np.where(metaOP['gdf']['OPENING_ID']==metaFCR['gdf']['OPENING_ID'][i])[0]
    if ind.size==0:
        print('Could not find a match')
        continue
    metaOP['gdf']['Area Reserves'][ind]=metaOP['gdf']['Area Reserves'][ind]+metaFCR['gdf']['GEOMETRY_Area'][i]

# Extract year from date
metaOP['gdf']['DENUDATION_1_COMPLETION_Year']=np.zeros(metaOP['gdf']['OPENING_ID'].size)
metaOP['gdf']['DENUDATION_2_COMPLETION_Year']=np.zeros(metaOP['gdf']['OPENING_ID'].size)
for i in range(metaOP['gdf']['OPENING_ID'].size):
    try:
        metaOP['gdf']['DENUDATION_1_COMPLETION_Year'][i]=np.array(metaOP['gdf']['DENUDATION_1_COMPLETION_DATE'][i][0:4]).astype(float)
    except:
        pass
    try:
        metaOP['gdf']['DENUDATION_2_COMPLETION_Year'][i]=np.array(metaOP['gdf']['DENUDATION_2_COMPLETION_DATE'][i][0:4]).astype(float)
    except:
        pass
#u=np.unique(metaOP['gdf']['DENUDATION_1_DISTURBANCE_CODE'])

metaOP['ts']={}
metaOP['ts']['OPENING_GROSS_AREA']=np.zeros(tv.size)
metaOP['ts']['GEOMETRY_Area']=np.zeros(tv.size)
metaOP['ts']['Area Reserves']=np.zeros(tv.size)
for iT in range(tv.size):
    ind=np.where( (metaOP['gdf']['DENUDATION_1_COMPLETION_Year']==tv[iT]) & (metaOP['gdf']['DENUDATION_1_DISTURBANCE_CODE']=='L') | \
                 (metaOP['gdf']['DENUDATION_1_COMPLETION_Year']==tv[iT]) & (metaOP['gdf']['DENUDATION_1_DISTURBANCE_CODE']=='S') | \
                 (metaOP['gdf']['DENUDATION_1_COMPLETION_Year']==tv[iT]) & (metaOP['gdf']['DENUDATION_1_DISTURBANCE_CODE']=='R') | \
                 (metaOP['gdf']['DENUDATION_2_COMPLETION_Year']==tv[iT]) & (metaOP['gdf']['DENUDATION_2_DISTURBANCE_CODE']=='L') | \
                 (metaOP['gdf']['DENUDATION_2_COMPLETION_Year']==tv[iT]) & (metaOP['gdf']['DENUDATION_2_DISTURBANCE_CODE']=='S') | \
                 (metaOP['gdf']['DENUDATION_2_COMPLETION_Year']==tv[iT]) & (metaOP['gdf']['DENUDATION_2_DISTURBANCE_CODE']=='R') )[0]
    metaOP['ts']['OPENING_GROSS_AREA'][iT]=np.nansum(metaOP['gdf']['OPENING_GROSS_AREA'][ind])
    metaOP['ts']['GEOMETRY_Area'][iT]=np.nansum(metaOP['gdf']['GEOMETRY_Area'][ind])
    metaOP['ts']['Area Reserves'][iT]=np.nansum(metaOP['gdf']['Area Reserves'][ind]/10000)

#%% Annual harvest area from RESULTS

metaAT={}
metaAT['Path']=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20220422\Results.gdb'
metaAT['Layer']='RSLT_ACTIVITY_TREATMENT_SVW'; # fiona.listlayers(op['Path'])
metaAT['crs']=gdf_bm.crs
metaAT['Keep Geom']='Off'
metaAT['Select Openings']=np.array([])
metaAT['SBC']=np.array([])
metaAT['FSC']=np.array([])
metaAT['ROI']=[]
metaAT['gdf']=qgdb.Query_Openings(metaAT,[])

metaAT['gdf']['Year']=np.zeros(metaAT['gdf']['OPENING_ID'].size)
for i in range(metaAT['gdf']['OPENING_ID'].size):
    try:
        metaAT['gdf']['Year'][i]=np.array(metaAT['gdf']['ATU_COMPLETION_DATE'][i][0:4]).astype(float)
    except:
        pass

# There is veritually no direct seeding so I stopped doing planting seperately
metaAT['ts']={}
metaAT['ts']['Area Regen Total']=np.zeros(tv.size)
metaAT['ts']['% Missing']=np.zeros(tv.size)
#metaAT['ts']['Area PL']=np.zeros(tv.size)
#metaAT['ts']['Area DS']=np.zeros(tv.size)
for iT in range(tv.size):
    #ind=np.where( (metaAT['gdf']['Year']==tv[iT]) & (metaAT['gdf']['SILV_BASE_CODE']=='PL') & (metaAT['gdf']['SILV_METHOD_CODE']!='LAYOT') )[0]
    #metaAT['ts']['Area PL'][iT]=np.sum(metaAT['gdf']['ACTUAL_TREATMENT_AREA'][ind])
    #ind=np.where( (metaAT['gdf']['Year']==tv[iT]) & (metaAT['gdf']['SILV_BASE_CODE']=='DS') & (metaAT['gdf']['SILV_METHOD_CODE']!='GS') )[0]
    #metaAT['ts']['Area DS'][iT]=np.sum(metaAT['gdf']['ACTUAL_TREATMENT_AREA'][ind])
    ind=np.where( (metaAT['gdf']['Year']==tv[iT]) & (metaAT['gdf']['RESULTS_IND']=='Y') & (metaAT['gdf']['SILV_BASE_CODE']=='PL') & (metaAT['gdf']['SILV_METHOD_CODE']!='LAYOT') | \
        (metaAT['gdf']['Year']==tv[iT]) & (metaAT['gdf']['RESULTS_IND']=='Y') & (metaAT['gdf']['SILV_BASE_CODE']=='DS') & (metaAT['gdf']['SILV_TECHNIQUE_CODE']!='GS'))[0]
    if ind.size==0:
        continue
    metaAT['ts']['Area Regen Total'][iT]=np.nansum(metaAT['gdf']['ACTUAL_TREATMENT_AREA'][ind])
    ind2=np.where(np.isnan(metaAT['gdf']['ACTUAL_TREATMENT_AREA'][ind])==True)[0]
    metaAT['ts']['% Missing'][iT]=ind2.size/ind.size*100

# Check nans - very little!
#plt.plot(tv,metaAT['ts']['% Missing'],'-ko')

#gu.opickle(r'C:\Users\rhember\Documents\Data\Harvest\Harvest Area\AnnualHarvestAreaFromRESULTS.pkl',metaAT)

#%% Plot

# Stats
iT=np.where( (tv>1989) & (tv<=2020) & (np.isnan(metaAT['ts']['Area Regen Total'])==False) )[0]
A_harv=np.mean(dNFD['Area'][iT]/1e3)
A_regen=np.mean(metaAT['ts']['Area Regen Total'][iT]/1e3)
print(A_regen/A_harv)

gp['ms']=2.5

cl=np.array([[0.17,0.35,0.7],[0.3,0.6,0.8],[0.5,0.9,1],[0,0.75,0],[0.6,1,0]])

plt.close('all');
fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,7.5)); ms=1.25; lw=0.5;
ax.plot(tv,dNFD['Area']/1e3,'-bo',mfc=cl[0,:],mec=cl[0,:],color=cl[0,:],lw=lw,ms=ms,label='Harvest area (NFD)')
ax.plot(tv,dCC['Area Harvested']/1e3,'-bo',mfc=cl[1,:],mec=cl[1,:],color=cl[1,:],lw=lw,ms=ms,label='Harvest area (consolidated cutblocks database)')
ax.plot(tv,(metaOP['ts']['OPENING_GROSS_AREA']-metaOP['ts']['Area Reserves'])/1e3,'-gs',mfc=cl[2,:],mec=cl[2,:],color=cl[2,:],lw=lw,ms=ms,label='Gross area minus reserves (RESULTS)')
ax.plot(tv,metaAT['ts']['Area Regen Total']/1e3,'-r^',mfc=cl[3,:],mec=cl[3,:],color=cl[3,:],lw=lw,ms=ms,label='Area planted (RESULTS)')
ax.set(position=[0.085,0.125,0.88,0.84],xticks=np.arange(1800,2120,5),yticks=np.arange(0,400,25),ylabel='Area harvested (Thousand ha yr$^-$$^1$)',xlabel='Time, years',ylim=[0,300],xlim=[1949.5,2021.5])
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=1.5)
ax.legend(loc='upper left',facecolor=[1,1,1],frameon=False)
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Harvest\Harvest Area\AreaHarvestedBC','png',900)

#%% NSR analysis

# Maximum of harvest area estimates
H_max=np.nanmax(np.column_stack((dNFD['Area'],dCC['Area Harvested'],metaOP['ts']['OPENING_GROSS_AREA']-metaOP['ts']['Area Reserves'])),axis=1)

# Import NOSE
nose=gu.ipickle(r'D:\Data\FCI_Projects\BCFCS_NOSE_Sparse\Inputs\AIL.pkl')
A_regen_wo_nose=np.zeros(tv.size)
iT=np.where( (tv>=nose['Year'][0]) & (tv<=nose['Year'][-1]) )[0]
A_regen_wo_nose[iT]=nose['Area All']

plt.close('all')
#plt.plot(tv,,'og-')
plt.plot(tv,H_max-metaAT['ts']['Area Regen Total'],'-ko')

NSR=H_max-metaAT['ts']['Area Regen Total']
NSR_wo_nose=H_max-metaAT['ts']['Area Regen Total']+A_regen_wo_nose

plt.close('all');
fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(15,7));
ax[0].plot(tv,np.zeros(tv.size),'-k',lw=0.5)
ax[0].fill_between(tv,NSR/1e3,NSR_wo_nose/1e3,color=[0.9,0.9,0.9],alpha=1,lw=0,label='NOSE programs')
ax[0].plot(tv,NSR_wo_nose/1e3,'-go',lw=0.5,ms=1.25,label='Without NOSE')
ax[0].plot(tv,NSR/1e3,'-bo',mfc=cl[0,:],mec=cl[0,:],color=cl[0,:],lw=0.5,ms=1.25,label='Actual')
ax[0].set(xticks=np.arange(1800,2120,10),yticks=np.arange(-125,275,25),ylabel='Change in harvest NSR area (ha x 1000 yr$^-$$^1$)',xlabel='Time, years',ylim=[-125,250],xlim=[1949.5,2021.5])
ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].tick_params(length=1.5)
ax[0].legend(loc='lower left',facecolor=[1,1,1],frameon=False)

ax[1].fill_between(tv,np.cumsum(NSR/1e6),np.cumsum(NSR_wo_nose/1e6),color=[0.9,0.9,0.9],alpha=1,lw=0)
ax[1].plot(tv,np.cumsum(NSR_wo_nose)/1e6,'-g',lw=0.5,ms=2,label='Without NOSE')
ax[1].plot(tv,np.cumsum((H_max-metaAT['ts']['Area Regen Total'])/1e6),'-b',mfc=cl[0,:],mec=cl[0,:],color=cl[0,:],ms=2,label='Harvest area (NFD)')
ax[1].set(xticks=np.arange(1800,2120,10),yticks=np.arange(0,10000,1),ylabel='Cumulative harvest NSR area (Mha)',xlabel='Time, years',ylim=[0,7],xlim=[1949.5,2021.5])
ax[1].yaxis.set_ticks_position('both'); ax[1].xaxis.set_ticks_position('both'); ax[1].tick_params(length=1.5)
gu.axletters(ax,plt,0.04,0.92,FontColor=gp['cla'],LetterStyle='Caps',FontWeight='Bold')
plt.tight_layout()
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Harvest\Harvest Area\NSR_Harvest_Summary','png',900)


