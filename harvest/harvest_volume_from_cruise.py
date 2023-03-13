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

gp=gu.SetGraphics('Manuscript')

#%% Prepare

tv=np.arange(2010,2022,1)

#%%

d={}
for iT in range(tv.size):
    df=pd.read_csv(r'C:\Users\rhember\Documents\Data\ForestInventory\Timber Cruises\InteriorStoneQuery ' + str(tv[iT]) + '.csv')
    d0=df.to_dict(orient='list')
    for k in d0.keys():
        idx=np.char.find(k,')')
        k2=k[idx+2:]
        if k2 not in d.keys():
            d[k2]=d0[k]
        else:
            d[k2]=np.append(d[k2],d0[k])

for k in d.keys():
    if d[k].dtype=='float64':
        d[k]=np.nan_to_num(d[k])

#%% Derived variables

# volume per hectare for the cutting permit (warning:  this may be a mixture of total VPH and Coniferous VPH)
# If you want to be sure of getting con. Vph then calculate it from NCV and net merch area
d['V Net Tot (m3/ha)']=d['CP VPH']
d['V Net Conif (m3/ha)']=d['NCV']/d['Tot Merch Area']
d['V Net Decid (m3/ha)']=d['Decid. Vol.']/d['Tot Merch Area']
d['V Net From Types (m3/ha)']=d['GSCC Vol/Ha']+d['GSPC Vol/Ha']+d['CHGCC Vol/Ha']+d['CHGPC Vol/Ha']+d['SCC Vol/Ha']

ind=np.where( (d['V Net Tot (m3/ha)']>0) & (d['V Net Conif (m3/ha)']>=0) & (d['V Net Decid (m3/ha)']>=0) & (d['V Net From Types (m3/ha)']>0) )[0]

print(np.mean(d['V Net Tot (m3/ha)'][ind]))
print(np.mean(d['V Net Conif (m3/ha)'][ind]))
print(np.mean(d['V Net Decid (m3/ha)'][ind]))
print(np.mean(d['V Net From Types (m3/ha)'][ind]))

#%% How many
d['Mark']=d['Mark'].astype('object')

u=np.unique(d['Mark'])

N=np.zeros(u.size)
for i in range(u.size):
    ind=np.where(d['Mark']==u[i])[0]
    d['CP'][ind]
    N[i]=ind.size

    plt.plot(d['Albers X Coordinate'][ind],d['Albers Y Coordinate'][ind],'b.')
np.mean(N)

#%%
np.unique(d['Point of App.'])



plt.hist(d['GSCC Vol/Ha'])
plt.hist(d['GSPC Vol/Ha'])
plt.hist(d['CHGCC Vol/Ha'])
plt.hist(d['CHGPC Vol/Ha'])
plt.hist(d['SCC Vol/Ha'])

plt.hist(d['V Tot (m3/ha)'])



(32) Decid. Vol.
(33) Salvage



