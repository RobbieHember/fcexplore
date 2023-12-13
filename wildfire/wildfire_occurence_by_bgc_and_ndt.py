'''

WILDFIRE OCCURRENCE STATS AND SCENARIOS BY BGC ZONE and NDT ZONE

'''

#%% Import Python modules

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import geopandas as gpd
import time
import gc as garc
import scipy.stats as stats
import statsmodels.api as sm
import statsmodels.formula.api as smf
import fcgadgets.macgyver.utilities_general as gu
import fcgadgets.macgyver.utilities_gis as gis

#%% Project assumptions

tv_obs=np.arange(1920,2022,1)
tv_scn=np.arange(-2000,2201,1)
tv_long=np.arange(1000,2022,1)

#%% Import zones

#lut=u1ha.Import_BC1ha_LUTs()
d=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\lut_bgcz_ndt_combo.xlsx')

zZone=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\bgcz_ndt_combo.tif')

#%% Calculate observed historical probability of occurrence from wildfire perimiter

# Get indices to each zone
iZone=[]
for iZ in range(d['ID'].size):
    iZone.append(np.where(zZone['Data']==d['ID'][iZ]))

# Import wildfire perimiters by year and sum area burned
d['Area Burned']=np.zeros(d['ID'].size)
for iT in range(tv_obs.size):
    print(tv_obs[iT])
    zFire=gis.OpenGeoTiff(r'C:\Users\rhember\Documents\Data\BC1ha\Disturbances\PROT_HISTORICAL_FIRE_POLYS_SP_' + str(tv_obs[iT]) + '.tif')
    for iZ in range(d['ID'].size):
        zFire0=zFire['Data'][iZone[iZ]]
        ind=np.where( (zFire0>0) )[0]
        if ind.size==0:
            continue
        d['Area Burned'][iZ]=d['Area Burned'][iZ]+ind.size

# Calculate annual probability of wildfire occurrence
d['Po (%/yr)']=d['Area Burned']/(tv_obs.size*d['Area'])*100

# Save
df=pd.DataFrame(d)
df.to_excel(r'C:\Users\rhember\Documents\Data\BC1ha\VRI\lut_bgcz_ndt_combo.xlsx')

