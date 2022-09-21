'''

Get planting and surveys that overlay with wildfires

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
from fcgadgets.cbrunner import cbrun_utilities as cbu

#%% Get planting and surveys that overlay with wildfire
    
gdf_su=GetSurveySpatial()

gdf_pl=GetPlanting()

gdf_wf=GetWildfire()

gdf_pl_wf=gpd.overlay(gdf_pl,gdf_wf,how='intersection')

gdf_su_wf=gpd.overlay(gdf_su,gdf_wf,how='intersection')

gdf_su_wf_notpl=gpd.overlay(gdf_su_wf,gdf_pl,how='difference')

gdf_su_wf_notpl.plot()
gdf_su_wf.plot()
gdf_pl.plot()

gdf_pl_wf.to_file(filename=r'C:\Users\rhember\Documents\Data\FCI_Projects\FCI_RollupFCI_Inv\Geospatial\pl_wf.geojson',driver='GeoJSON')
gdf_su_wf_notpl.to_file(filename=r'C:\Users\rhember\Documents\Data\FCI_Projects\FCI_RollupFCI_Inv\Geospatial\su_wf_notpl.geojson',driver='GeoJSON')
gdf_su_wf.to_file(filename=r'C:\Users\rhember\Documents\Data\FCI_Projects\FCI_RollupFCI_Inv\Geospatial\su_wf.geojson',driver='GeoJSON')
gdf_pl.to_file(filename=r'C:\Users\rhember\Documents\Data\FCI_Projects\FCI_RollupFCI_Inv\Geospatial\pl.geojson',driver='GeoJSON')