#%% Import modules

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
#import fcgadgets.macgyver.util_query_gdb as qgdb
import fcgadgets.bc1ha.bc1ha_util as u1ha
import fcgadgets.macgyver.util_fcs_graphs as ufcs

#%% Import data

meta=u1ha.Init()
gp=gu.SetGraphics('Manuscript')

# Initialize dictionary
d={}
d['tv']=np.arange(1950,2023,1)

#%% National Forest Database

d0=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\Harvest\Harvest Area\NFD - Area Harvested by ownership and harvesting method - EN FR.xlsx')
d['Area Harv NFD']=np.nan*np.ones(d['tv'].size)
for iT in range(d['tv'].size):
    ind=np.where( (d0['Jurisdiction']=='British Columbia') & (d0['Year']==d['tv'][iT]) )[0]
    if ind.size==0:
        continue
    d['Area Harv NFD'][iT]=np.nansum(d0['Area (hectares)'][ind])

#%% Annual harvest area from consolidated cutblock DB

lNam='VEG_CONSOLIDATED_CUT_BLOCKS_SP'
vNam='HARVEST_YEAR'
nPack=6
tv,N=u1ha.TimeSeriesOccurrenceFromPackedEvents(meta,lNam,vNam,nPack)
d['Area Harv CC']=np.nan*np.ones(d['tv'].size)
for iT in range(d['tv'].size):
    ind=np.where(tv==d['tv'][iT])[0]
    if N[ind]>0:
        d['Area Harv CC'][iT]=N[ind]


#%% Harvest from opening layer

# # Import opening layer
# metaOP={}
# metaOP['Path']=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20220422\Results.gdb'
# metaOP['Layer']='RSLT_OPENING_SVW'; # fiona.listlayers(op['Path'])
# metaOP['crs']=meta['Geos']['crs'].crs
# metaOP['Keep Geom']='Off'
# metaOP['Select Openings']=np.array([])
# metaOP['SBC']=np.array([])
# metaOP['FSC']=np.array([])
# metaOP['ROI']=[]
# metaOP['gdf']=qgdb.Query_Openings(metaOP,[])

# # Import forest Cover Reserves
# metaFCR={}
# metaFCR['Path']=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20220422\Results.gdb'
# metaFCR['Layer']='RSLT_FOREST_COVER_RESERVE_SVW'; # fiona.listlayers(metaFCR['Path'])
# metaFCR['crs']=meta['Geos']['crs'].crs
# metaFCR['Keep Geom']='Off'
# metaFCR['Select Openings']=np.array([])
# metaFCR['SBC']=np.array([])
# metaFCR['FSC']=np.array([])
# metaFCR['ROI']=[]
# metaFCR['gdf']=qgdb.Query_Openings(metaFCR,[])

# # Add reserve areas to opening layer
# metaOP['gdf']['Area Reserves']=np.zeros(metaOP['gdf']['OPENING_ID'].size)
# for i in range(metaFCR['gdf']['OPENING_ID'].size):
#     ind=np.where(metaOP['gdf']['OPENING_ID']==metaFCR['gdf']['OPENING_ID'][i])[0]
#     if ind.size==0:
#         print('Could not find a match')
#         continue
#     metaOP['gdf']['Area Reserves'][ind]=metaOP['gdf']['Area Reserves'][ind]+metaFCR['gdf']['GEOMETRY_Area'][i]

# # Extract year from date
# metaOP['gdf']['DENUDATION_1_COMPLETION_Year']=np.zeros(metaOP['gdf']['OPENING_ID'].size)
# metaOP['gdf']['DENUDATION_2_COMPLETION_Year']=np.zeros(metaOP['gdf']['OPENING_ID'].size)
# for i in range(metaOP['gdf']['OPENING_ID'].size):
#     try:
#         metaOP['gdf']['DENUDATION_1_COMPLETION_Year'][i]=np.array(metaOP['gdf']['DENUDATION_1_COMPLETION_DATE'][i][0:4]).astype(float)
#     except:
#         pass
#     try:
#         metaOP['gdf']['DENUDATION_2_COMPLETION_Year'][i]=np.array(metaOP['gdf']['DENUDATION_2_COMPLETION_DATE'][i][0:4]).astype(float)
#     except:
#         pass
# #u=np.unique(metaOP['gdf']['DENUDATION_1_DISTURBANCE_CODE'])

# metaOP['ts']={}
# metaOP['ts']['OPENING_GROSS_AREA']=np.zeros(d['tv'].size)
# metaOP['ts']['GEOMETRY_Area']=np.zeros(d['tv'].size)
# metaOP['ts']['Area Reserves']=np.zeros(d['tv'].size)
# for iT in range(d['tv'].size):
#     ind=np.where( (metaOP['gdf']['DENUDATION_1_COMPLETION_Year']==d['tv'][iT]) & (metaOP['gdf']['DENUDATION_1_DISTURBANCE_CODE']=='L') | \
#                  (metaOP['gdf']['DENUDATION_1_COMPLETION_Year']==d['tv'][iT]) & (metaOP['gdf']['DENUDATION_1_DISTURBANCE_CODE']=='S') | \
#                  (metaOP['gdf']['DENUDATION_1_COMPLETION_Year']==d['tv'][iT]) & (metaOP['gdf']['DENUDATION_1_DISTURBANCE_CODE']=='R') | \
#                  (metaOP['gdf']['DENUDATION_2_COMPLETION_Year']==d['tv'][iT]) & (metaOP['gdf']['DENUDATION_2_DISTURBANCE_CODE']=='L') | \
#                  (metaOP['gdf']['DENUDATION_2_COMPLETION_Year']==d['tv'][iT]) & (metaOP['gdf']['DENUDATION_2_DISTURBANCE_CODE']=='S') | \
#                  (metaOP['gdf']['DENUDATION_2_COMPLETION_Year']==d['tv'][iT]) & (metaOP['gdf']['DENUDATION_2_DISTURBANCE_CODE']=='R') )[0]
#     metaOP['ts']['OPENING_GROSS_AREA'][iT]=np.nansum(metaOP['gdf']['OPENING_GROSS_AREA'][ind])
#     metaOP['ts']['GEOMETRY_Area'][iT]=np.nansum(metaOP['gdf']['GEOMETRY_Area'][ind])
#     metaOP['ts']['Area Reserves'][iT]=np.nansum(metaOP['gdf']['Area Reserves'][ind]/10000)

#%% Annual harvest area from RESULTS

# metaAT={}
# metaAT['Path']=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20220422\Results.gdb'
# metaAT['Layer']='RSLT_ACTIVITY_TREATMENT_SVW'; # fiona.listlayers(op['Path'])
# metaAT['crs']=meta['Geos']['crs'].crs
# metaAT['Keep Geom']='Off'
# metaAT['Select Openings']=np.array([])
# metaAT['SBC']=np.array([])
# metaAT['FSC']=np.array([])
# metaAT['ROI']=[]
# metaAT['gdf']=qgdb.Query_Openings(metaAT,[])

# # Get a complete list of FSC
# #metaAT['gdf'].keys()
# #list(np.unique(metaAT['gdf']['SILV_FUND_SOURCE_CODE'][metaAT['gdf']['SILV_FUND_SOURCE_CODE']!=None]))

# metaAT['gdf']['Year']=np.zeros(metaAT['gdf']['OPENING_ID'].size)
# for i in range(metaAT['gdf']['OPENING_ID'].size):
#     try:
#         metaAT['gdf']['Year'][i]=np.array(metaAT['gdf']['ATU_COMPLETION_DATE'][i][0:4]).astype(float)
#     except:
#         pass

# # There is veritually no direct seeding so I stopped doing planting seperately
# metaAT['ts']={}
# metaAT['ts']['Area Regen Total']=np.zeros(d['tv'].size)
# metaAT['ts']['% Missing']=np.zeros(d['tv'].size)
# #metaAT['ts']['Area PL']=np.zeros(d['tv'].size)
# #metaAT['ts']['Area DS']=np.zeros(d['tv'].size)
# for iT in range(d['tv'].size):
#     #ind=np.where( (metaAT['gdf']['Year']==d['tv'][iT]) & (metaAT['gdf']['SILV_BASE_CODE']=='PL') & (metaAT['gdf']['SILV_METHOD_CODE']!='LAYOT') )[0]
#     #metaAT['ts']['Area PL'][iT]=np.sum(metaAT['gdf']['ACTUAL_TREATMENT_AREA'][ind])
#     #ind=np.where( (metaAT['gdf']['Year']==d['tv'][iT]) & (metaAT['gdf']['SILV_BASE_CODE']=='DS') & (metaAT['gdf']['SILV_METHOD_CODE']!='GS') )[0]
#     #metaAT['ts']['Area DS'][iT]=np.sum(metaAT['gdf']['ACTUAL_TREATMENT_AREA'][ind])
#     ind=np.where( (metaAT['gdf']['Year']==d['tv'][iT]) & (metaAT['gdf']['RESULTS_IND']=='Y') & (metaAT['gdf']['SILV_BASE_CODE']=='PL') & (metaAT['gdf']['SILV_METHOD_CODE']!='LAYOT') | \
#         (metaAT['gdf']['Year']==d['tv'][iT]) & (metaAT['gdf']['RESULTS_IND']=='Y') & (metaAT['gdf']['SILV_BASE_CODE']=='DS') & (metaAT['gdf']['SILV_TECHNIQUE_CODE']!='GS'))[0]
#     if ind.size==0:
#         continue
#     metaAT['ts']['Area Regen Total'][iT]=np.nansum(metaAT['gdf']['ACTUAL_TREATMENT_AREA'][ind])
#     ind2=np.where(np.isnan(metaAT['gdf']['ACTUAL_TREATMENT_AREA'][ind])==True)[0]
#     metaAT['ts']['% Missing'][iT]=ind2.size/ind.size*100

# # Check nans - very little!
# #plt.plot(d['tv'],metaAT['ts']['% Missing'],'-ko')

# # Populate
# d['Area Harv RESULTS']=metaOP['ts']['OPENING_GROSS_AREA']-metaOP['ts']['Area Reserves']
# d['Area Planted RESULTS']=metaAT['ts']['Area Regen Total']

#%% Harvest area from cruise comp

dCru=gu.ipickle(r'C:\Users\rhember\Documents\Data\ECAS\Received 2023-04-04\rh_COMP_DATA_CleanCombined.pkl')

d['Area Harv Cruise']=np.nan*np.ones(d['tv'].size)
for iT in range(d['tv'].size):
    ind=np.where( (dCru['Year']==d['tv'][iT]) )[0]
    if d['tv'][iT]<2015:
        continue
    d['Area Harv Cruise'][iT]=np.nansum(dCru['NET_AREA'][ind])

#%% Harvest area from NTEM

zH=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Disturbances\\Harvest_NTEM_Year.tif')
zH=zH['Data'][zH['Data']>0]
d['Area Harv NTEM']=np.nan*np.ones(d['tv'].size)
for iT in range(d['tv'].size):
    ind=np.where(zH==d['tv'][iT])
    if ind[0].size!=0:
        d['Area Harv NTEM'][iT]=ind[0].size

#%% Save time series of harvest area

# Maximum of harvest area estimates
#d['Area Harv Max']=np.nanmax(np.column_stack( (d['Area Harv NFD'],d['Area Harv CC'],d['Area Harv RESULTS']) ),axis=1)

gu.opickle(r'C:\Users\rhember\Documents\Data\Harvest\Harvest Area\HarvestAreaBC.pkl',d)

#%% Plot time series

meta['Graphics']['Print Figures']='On'
meta['Graphics'][ 'Print Figure Path']=r'C:\Users\rhember\OneDrive - Government of BC\Figures\Harvest\Harvest Area'
#Plot_HarvestAreaTimeSeries(meta,[1990,2023])
ufcs.Plot_HarvestAreaTimeSeries(meta,[1990,2023])

