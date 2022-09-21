'''

Query forest cover reserves

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
from fcgadgets.macgyver import utilities_query_gdbs as qgdb
from fcgadgets.cbrunner import cbrun_utilities as cbu

#%% Import a GDB with SRS for BC

gdf_bm=gpd.read_file(r'C:\Users\rhember\Documents\Data\Basemaps\bc_land.shp')

#%% Project parameters

meta={}
meta['Path']=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20220422\Results.gdb'
meta['Layer']='RSLT_FOREST_COVER_RESERVE_SVW'; # fiona.listlayers(meta['Path'])
meta['crs']=gdf_bm.crs
meta['keep_geom']='Off'

#%% Have a quick look at layer to see what the fields are

with fiona.open(meta['Path'],layer=meta['Layer']) as source:
    for feat in source:
        break
feat['properties']

#%% Import forest cover reserve

d=qgdb.Query_FC_Reserves(meta)

#%% Import openings 





#%% Project parameters

meta={}
meta['Path']=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20220422\Results.gdb'
meta['Layer']='RSLT_OPENING_SVW'; # fiona.listlayers(meta['Path'])
meta['crs']=gdf_bm.crs
meta['Keep Geom']='Off'
meta['Select Openings']=np.array([])

d=qgdb.Query_Openings(meta)





with fiona.open(meta2['Path'],layer=meta2['Layer']) as source:
    for feat in source:
        break
feat['properties']



#%% Analysis

tv=np.arange(1950,2023,1)