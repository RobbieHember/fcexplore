'''

QUERY NUTRIENT APPLICATION

'''

#%% Import modules

import numpy as np
import pandas as pd
import geopandas as gpd
import fiona
import matplotlib.pyplot as plt
import numpy.matlib as ml
from shapely.geometry import Polygon,Point
import time
import gc as garc
from fcgadgets.macgyver import utilities_general as gu
from fcgadgets.macgyver import utilities_gis as gis
from fcgadgets.macgyver import utilities_inventory as invu
from fcgadgets.cbrunner import cbrun_utilities
from fcgadgets.bc1ha import bc1ha_utilities as bc1hau
from fcgadgets.macgyver import utilities_query_gdb as qgdb

# Set figure properties
gp=gu.SetGraphics('Manuscript')

#%% Path

pth=r'D:\Data\FCI_Projects\BCFCS_NM_1k_Sparse'
pth_fig=r'C:\Users\rhember\OneDrive - Government of BC\Figures\BCFCS_NM_1k_Sparse'

#%% Query geodatabase of openings with non-obligation stand establishment

# Takes 30 min to import geodatabases

pth_RESULTS=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20220422\Results.gdb'

# Get spatial reference system
bm=gpd.read_file(r'C:\Users\rhember\Documents\Data\Basemaps\bc_land.shp')
crs=bm.crs; del bm

# Get unique set of openings with NOSE
atu={}
atu['Path']=pth_RESULTS
atu['Layer']='RSLT_ACTIVITY_TREATMENT_SVW'; # fiona.listlayers(atu['Path'])
atu['crs']=crs
atu['Keep Geom']='Off'
atu['Select Openings']=np.array([])
atu['SBC']=np.array(['FE'])
atu['FSC']=np.array([])
atu['ROI']=[]
atu['gdf']=qgdb.Query_Openings(atu,[])

# Unique openings
uOID=np.unique(atu['gdf']['OPENING_ID'])

# Get ATU planting with spatial
atuWS={}
atuWS['Path']=pth_RESULTS
atuWS['Layer']='RSLT_ACTIVITY_TREATMENT_SVW'; # fiona.listlayers(atuWS['Path'])
atuWS['crs']=crs
atuWS['Keep Geom']='On'
atuWS['Select Openings']=np.array([])
atuWS['SBC']=np.array(['FE'])
atuWS['FSC']=np.array([])
atuWS['ROI']=[]
atuWS['gdf']=qgdb.Query_Openings(atuWS,[])

# Get opening layer data for unique openings with NOSE
op={}
op['Path']=pth_RESULTS
op['Layer']='RSLT_OPENING_SVW'; # fiona.listlayers(op['Path'])
op['crs']=crs
op['Keep Geom']='On'
op['Select Openings']=uOID
op['SBC']=np.array([])
op['FSC']=np.array([])
op['ROI']=[]
op['gdf']=qgdb.Query_Openings(op,[])

# Get forest cover layer data for unique openings with NOSE
fcsilv={}
fcsilv['Path']=pth_RESULTS
fcsilv['Layer']='RSLT_FOREST_COVER_SILV_SVW'; # fiona.listlayers(fcsilv['Path'])
fcsilv['crs']=crs
fcsilv['Keep Geom']='On'
fcsilv['Select Openings']=uOID
fcsilv['SBC']=np.array([])
fcsilv['FSC']=np.array([])
fcsilv['ROI']=[]
fcsilv['gdf']=qgdb.Query_Openings(fcsilv,[])

#%% Query geodatabase of openings with fert

# Where possible replace opening layer geometry with planting geometry or
# forest cover with artificial

# This is a mask, but polygon areas should be bigger than area treated
# due to reliance on opening spatial where planting spatial is missing

List=[None]*int(1e5)
cnt=0
for i in range(uOID.size):

    # Try to use ATU planting spatial
    ind=np.where(atuWS['gdf']['OPENING_ID']==uOID[i])[0]
    if ind.size>0:
        for j in range(ind.size):
            feat={}
            feat['properties']={}
            feat['geometry']=atuWS['gdf'].loc[ind[j],'geometry']
            List[cnt]=feat
            cnt=cnt+1
            #break
    else:

        ind_atu=np.where(atu['gdf']['OPENING_ID']==uOID[i])[0]

        # Default to artifical forest cover layer spatial
        ind=np.where(fcsilv['gdf']['OPENING_ID']==uOID[i])[0]
        A_big0=np.sum(fcsilv['gdf']['GEOMETRY_Area'][ind].values)/10000
        if ind.size>0:
            for j in range(ind.size):
                if (fcsilv['gdf'].loc[ind[j],'STOCKING_TYPE_CODE']=='ART'):

                    A_atu=atu['gdf']['ACTUAL_TREATMENT_AREA'][ind_atu[0]]
                    # A_atu/A_big0

                    feat={}
                    feat['properties']={}
                    feat['geometry']=fcsilv['gdf'].loc[ind[j],'geometry']

                    # Shrink to match treatment area
                    a1=gpd.GeoDataFrame.from_features([feat],crs=op['crs'])
                    A_old=a1.area[0]/10000
                    rat0=A_old/A_big0

                    bin=np.arange(-500,10,10)
                    A_bin=np.zeros(bin.size)
                    for iBin in range(bin.size):
                        a2=a1.buffer(bin[iBin])
                        A_bin[iBin]=a2.area[0]/10000
                    rat_bin=A_bin/A_atu

                    AD=np.abs(rat_bin-rat0)
                    iMinD=np.where( AD==np.min(AD) )[0]
                    if iMinD.size>1:
                        iMinD=iMinD[-1]
                    a2=gpd.GeoDataFrame(geometry=a1.buffer(bin[iMinD]))
                    A_new=a2.area[0]/10000
                    rat_new=A_atu/A_new
                    #A_new/A_old

                    flg=0
                    if flg==1:
                        fig,ax=plt.subplots(1)
                        a1.plot(ax=ax)
                        #a3.plot(ax=ax,lw=1,color='r',facecolor=None)

                    feat['geometry']=a2.loc[0,'geometry']

                    List[cnt]=feat
                    cnt=cnt+1
                    #ind_atu=np.where(atu['gdf']['OPENING_ID']==uOID[i])[0]
                    #atu['gdf']['ACTUAL_TREATMENT_AREA'][ind_atu]
        else:

            # Default to opening
            ind=np.where(op['gdf']['OPENING_ID']==uOID[i])[0]
            A_big0=np.sum(op['gdf']['GEOMETRY_Area'][ind].values)/10000
            for j in range(ind.size):

                A_atu=atu['gdf']['ACTUAL_TREATMENT_AREA'][ind_atu[0]]

                feat={}
                feat['properties']={}
                feat['geometry']=op['gdf'].loc[ind[j],'geometry']

                # Shrink to match treatment area
                a1=gpd.GeoDataFrame.from_features([feat],crs=op['crs'])
                A_old=a1.area[0]/10000
                rat0=A_old/A_big0

                bin=np.arange(-500,10,10)
                A_bin=np.zeros(bin.size)
                for iBin in range(bin.size):
                    a2=a1.buffer(bin[iBin])
                    A_bin[iBin]=a2.area[0]/10000
                rat_bin=A_bin/A_atu

                AD=np.abs(rat_bin-rat0)
                iMinD=np.where( AD==np.min(AD) )[0]
                if iMinD.size>1:
                    iMinD=iMinD[-1]
                a2=gpd.GeoDataFrame(geometry=a1.buffer(bin[iMinD]))
                A_new=a2.area[0]/10000
                rat_new=A_atu/A_new

                flg=0
                if flg==1:
                    fig,ax=plt.subplots(1)
                    a1.plot(ax=ax)
                    a2.plot(ax=ax,lw=1,color='r',facecolor=None)
                    #a3.plot(ax=ax,lw=1,color='r',facecolor=None)

                feat['geometry']=a2.loc[0,'geometry']

                List[cnt]=feat
                cnt=cnt+1

# Truncate
List=List[0:cnt]

# Convert to geodatabase
gdf_NA=gpd.GeoDataFrame.from_features(List,crs=op['crs'])

ind=np.where(gdf_NA['geometry'].is_valid==True)[0]
gdf_NA=gdf_NA.loc[ind]
gdf_NA=gdf_NA.reset_index()

gdf_NA.plot(linewidth=1,color='r')

# Save
gdf_NA.to_file(pth + '\\Geospatial\\NA_polygons.geojson',driver='GeoJSON')

#%% Annual treatment area treated by funding source

# Get unique funding source codes
uFSC=np.unique(atu['gdf']['SILV_FUND_SOURCE_CODE'][np.where(atu['gdf']['SILV_FUND_SOURCE_CODE']!=None)])

ail={}
ail['Year']=np.arange(1970,2021,1)
ail['Area Total']=np.zeros(ail['Year'].size)
ail['Area By FSC']=np.zeros((ail['Year'].size,uFSC.size))

# Calculate area by funding source code
for iT in range(ail['Year'].size):

    ind=np.where( (atu['gdf']['Year']==ail['Year'][iT]) & (atu['gdf']['SILV_BASE_CODE']=='FE') & (atu['gdf']['SILV_TECHNIQUE_CODE']=='CA') & (atu['gdf']['SILV_METHOD_CODE']=='HELI') & \
        (atu['gdf']['RESULTS_IND']=='Y') & (atu['gdf']['ACTUAL_TREATMENT_AREA']!=None) & (atu['gdf']['ATU_COMPLETION_DATE']!=None) & (atu['gdf']['SILV_FUND_SOURCE_CODE']!=None) )[0]

    ail['Area Total'][iT]=ail['Area Total'][iT]+np.sum(atu['gdf']['ACTUAL_TREATMENT_AREA'][ind])
    for iFSC in range(uFSC.size):
        ind=np.where( (atu['gdf']['Year']==ail['Year'][iT]) & (atu['gdf']['SILV_FUND_SOURCE_CODE']==uFSC[iFSC]) )[0]
        ail['Area By FSC'][iT,iFSC]=ail['Area By FSC'][iT,iFSC]+np.sum(atu['gdf']['ACTUAL_TREATMENT_AREA'][ind])

# Save
gu.opickle(pth + '\\Inputs\\AIL.pkl',ail)
#plt.plot(ail['Year'],ail['Area Total'],'-o')

#%%

A=0.0
for i in range(len(gdf_NA)):
    try:
        A=A+gdf_NA.loc[i,'geometry'].area/10000
    except:
        pass
print(A)
print(np.sum(ail['Area Total']))
A/np.sum(ail['Area Total'])

#%% Plot time series of area treated by funding source

cl=np.random.random((uFSC.size,3))

plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(15.5,6.5));
A_cumu=np.zeros(ail['Year'].size)
for iFSC in range(uFSC.size):
    plt.bar(ail['Year'],ail['Area By FSC'][:,iFSC]/1000,0.8,bottom=A_cumu,facecolor=cl[iFSC,:],label=uFSC[iFSC])
    A_cumu=A_cumu+ail['Area By FSC'][:,iFSC]/1000

ax.set(position=[0.06,0.12,0.92,0.86],xticks=np.arange(1950,2020+1,5),ylabel='Treatment area (hectares x 1000)',xlabel='Time, years',
       xlim=[ail['Year'][0]-0.75,ail['Year'][-1]+0+.75],ylim=[0,45])
plt.legend(frameon=False,loc='upper left',facecolor=[1,1,1],labelspacing=0.25,ncol=3)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=1.5)
gu.PrintFig(pth_fig + '\\AIL_ByFundingSource_All','png',900)