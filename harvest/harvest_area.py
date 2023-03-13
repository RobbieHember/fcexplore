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

#%% Import data

gp=gu.SetGraphics('Manuscript')

# For CRS
gdf_bm=gpd.read_file(r'C:\Users\rhember\Documents\Data\Basemaps\Basemaps.gdb',layer='NRC_POLITICAL_BOUNDARIES_1M_SP')

# Initialize dictionary
d={}
d['tv']=np.arange(1950,2022,1)

#%% National Forest Database

d0=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\Harvest\Harvest Area\NFD - Area harvested by ownership and harvesting method - EN FR.xlsx')
d['Area Harvested NFD']=np.nan*np.ones(d['tv'].size)
for iT in range(d['tv'].size):
    ind=np.where( (d0['Jurisdiction']=='British Columbia') & (d0['Year']==d['tv'][iT]) )[0]
    if ind.size==0:
        continue
    d['Area Harvested NFD'][iT]=np.nansum(d0['Area (hectares)'][ind])

#%% Annual harvest area from consolidated cutblock DB

fin=r'C:\Users\rhember\Documents\Data\ForestInventory\Disturbances\20210401\Disturbances.gdb'

d['Area Harvested CC']=np.zeros(d['tv'].size)
#d['Area Harvested CC From Geom']=np.zeros(d['tv'].size)
with fiona.open(fin,layer='VEG_CONSOLIDATED_CUT_BLOCKS_SP') as source:
    for feat in source:
        if feat['geometry']==None:
            continue
        iT=np.where(d['tv']==feat['properties']['HARVEST_YEAR'])[0]
        if iT.size==0:
            continue
        d['Area Harvested CC'][iT]=d['Area Harvested CC'][iT]+feat['properties']['AREA_HA']

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
metaOP['ts']['OPENING_GROSS_AREA']=np.zeros(d['tv'].size)
metaOP['ts']['GEOMETRY_Area']=np.zeros(d['tv'].size)
metaOP['ts']['Area Reserves']=np.zeros(d['tv'].size)
for iT in range(d['tv'].size):
    ind=np.where( (metaOP['gdf']['DENUDATION_1_COMPLETION_Year']==d['tv'][iT]) & (metaOP['gdf']['DENUDATION_1_DISTURBANCE_CODE']=='L') | \
                 (metaOP['gdf']['DENUDATION_1_COMPLETION_Year']==d['tv'][iT]) & (metaOP['gdf']['DENUDATION_1_DISTURBANCE_CODE']=='S') | \
                 (metaOP['gdf']['DENUDATION_1_COMPLETION_Year']==d['tv'][iT]) & (metaOP['gdf']['DENUDATION_1_DISTURBANCE_CODE']=='R') | \
                 (metaOP['gdf']['DENUDATION_2_COMPLETION_Year']==d['tv'][iT]) & (metaOP['gdf']['DENUDATION_2_DISTURBANCE_CODE']=='L') | \
                 (metaOP['gdf']['DENUDATION_2_COMPLETION_Year']==d['tv'][iT]) & (metaOP['gdf']['DENUDATION_2_DISTURBANCE_CODE']=='S') | \
                 (metaOP['gdf']['DENUDATION_2_COMPLETION_Year']==d['tv'][iT]) & (metaOP['gdf']['DENUDATION_2_DISTURBANCE_CODE']=='R') )[0]
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

# Get a complete list of FSC
#metaAT['gdf'].keys()
#list(np.unique(metaAT['gdf']['SILV_FUND_SOURCE_CODE'][metaAT['gdf']['SILV_FUND_SOURCE_CODE']!=None]))

metaAT['gdf']['Year']=np.zeros(metaAT['gdf']['OPENING_ID'].size)
for i in range(metaAT['gdf']['OPENING_ID'].size):
    try:
        metaAT['gdf']['Year'][i]=np.array(metaAT['gdf']['ATU_COMPLETION_DATE'][i][0:4]).astype(float)
    except:
        pass

# There is veritually no direct seeding so I stopped doing planting seperately
metaAT['ts']={}
metaAT['ts']['Area Regen Total']=np.zeros(d['tv'].size)
metaAT['ts']['% Missing']=np.zeros(d['tv'].size)
#metaAT['ts']['Area PL']=np.zeros(d['tv'].size)
#metaAT['ts']['Area DS']=np.zeros(d['tv'].size)
for iT in range(d['tv'].size):
    #ind=np.where( (metaAT['gdf']['Year']==d['tv'][iT]) & (metaAT['gdf']['SILV_BASE_CODE']=='PL') & (metaAT['gdf']['SILV_METHOD_CODE']!='LAYOT') )[0]
    #metaAT['ts']['Area PL'][iT]=np.sum(metaAT['gdf']['ACTUAL_TREATMENT_AREA'][ind])
    #ind=np.where( (metaAT['gdf']['Year']==d['tv'][iT]) & (metaAT['gdf']['SILV_BASE_CODE']=='DS') & (metaAT['gdf']['SILV_METHOD_CODE']!='GS') )[0]
    #metaAT['ts']['Area DS'][iT]=np.sum(metaAT['gdf']['ACTUAL_TREATMENT_AREA'][ind])
    ind=np.where( (metaAT['gdf']['Year']==d['tv'][iT]) & (metaAT['gdf']['RESULTS_IND']=='Y') & (metaAT['gdf']['SILV_BASE_CODE']=='PL') & (metaAT['gdf']['SILV_METHOD_CODE']!='LAYOT') | \
        (metaAT['gdf']['Year']==d['tv'][iT]) & (metaAT['gdf']['RESULTS_IND']=='Y') & (metaAT['gdf']['SILV_BASE_CODE']=='DS') & (metaAT['gdf']['SILV_TECHNIQUE_CODE']!='GS'))[0]
    if ind.size==0:
        continue
    metaAT['ts']['Area Regen Total'][iT]=np.nansum(metaAT['gdf']['ACTUAL_TREATMENT_AREA'][ind])
    ind2=np.where(np.isnan(metaAT['gdf']['ACTUAL_TREATMENT_AREA'][ind])==True)[0]
    metaAT['ts']['% Missing'][iT]=ind2.size/ind.size*100

# Check nans - very little!
#plt.plot(d['tv'],metaAT['ts']['% Missing'],'-ko')

# Populate
d['Area Harvested RESULTS']=metaOP['ts']['OPENING_GROSS_AREA']-metaOP['ts']['Area Reserves']
d['Area Planted RESULTS']=metaAT['ts']['Area Regen Total']

#%% Save time series of harvest area

# Maximum of harvest area estimates
d['Area Harvested Max']=np.nanmax(np.column_stack( (d['Area Harvested NFD'],d['Area Harvested CC'],d['Area Harvested RESULTS']) ),axis=1)

gu.opickle(r'C:\Users\rhember\Documents\Data\Harvest\Harvest Area\HarvestAreaBC.pkl',d)



