'''

HARVEST STATS AND SCENARIOS

See documentation.

'''

#%% Import Python modules

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pickle
import geopandas as gpd
import openpyxl
import time
import gc as garc
import scipy.stats as stats
import statsmodels.api as sm
from numpy import matlib as mb
from shapely.geometry import Point
from fcgadgets.utilities import utilities_general as gu
from fcgadgets.utilities import utilities_gis as gis
from fcgadgets.utilities import utilities_inventory as invu
from fcgadgets.taz import wildfire_stat_models as wfsm

#%% Path of file to store stats and scenarios

PathData=r'G:\My Drive\Data\Wildfire\Wildfire_Stats_Scenarios_By_BGCZ\Wildfire_Stats_Scenarios_By_BGCZ.pkl'
PathFigures=r'G:\My Drive\Figures\Wildfire\Wildfire_Stats_Sceanrios_By_BGCZ'

#%% Set figure properties

params=gu.Import_GraphicsParameters('spyder_fs7')
plt.rcParams.update(params)

#%% Import data from BC state of forest report

df=pd.read_csv(r'G:\My Drive\Data\Harvest\SummaryDataChangeInTimberHarvest\bctimberharvest.csv')

plt.plot(df['Year'],df['Total_harvest_millions_m3'],'-o')

#%% Historical time series of harvest rates

tv=np.arange(1800,2018,1)

f1=0.0004*10**((tv-1841)/100)
f2=(1/(1+np.exp(0.075*(tv-1955))))

plt.close('all')
plt.plot(tv,f1*f2,'-')
plt.grid()

f1=0.0015*27**((tv-1900)/100)
f2=(1/(1+np.exp(0.12*(tv-1945))))
plt.plot(tv,f1*f2,'--')

#%%


# Correct it to blend into observations
Adjuster=np.ones(tv.size)
it=np.where(tv>=1955)[0]
for i in range(it.size):
    Adjuster[it[i]]=0.93*Adjuster[it[i]-1]
yhat_adj=yhat*Adjuster

plt.plot(tv,yhat_adj,'--')

# Save
df=pd.DataFrame(data=np.column_stack((tv,yhat_adj)),columns=['Time','Probability of Harvest'])
df.to_excel(r'G:\My Drive\Data\Harvest\Historical BC Harvesting\HarvestHistoricalProbabilitySimple.xlsx')






















