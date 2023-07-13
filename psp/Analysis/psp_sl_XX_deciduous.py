
"""
PSP - DESCRIPTIVE STATISTICS
"""

#%% Import modules

import numpy as np
import gc as garc
import time
import matplotlib.pyplot as plt
import geopandas as gpd
import pandas as pd
import fiona
import statsmodels.formula.api as smf
import fcgadgets.macgyver.utilities_general as gu
import fcgadgets.macgyver.utilities_gis as gis
import fcexplore.psp.Processing.psp_utilities as ugp
from fcgadgets.bc1ha import bc1ha_utilities as u1ha

#%% Set figure properties

gp=gu.SetGraphics('Manuscript')

#%% Import data

metaGP,gplt=ugp.ImportPSPs(type='Stand')

#%%

ind=np.where( (gplt['PTF CNY']>=0) & (gplt['Deciduous L %N t0']>=0) & (gplt['Ecozone BC L1']==metaGP['LUT']['Ecozone BC L1']['IDF']) )[0]
bw=10; bin=np.arange(5,105,bw)
N,mu,med,sig,se=gu.discres(gplt['Deciduous L %BA t0'][ind],gplt['Ctot L t0'][ind],bw,bin)
plt.plot(bin,mu,'bo')


#%%

bw=20; bin=np.arange(20,300,bw)
plt.close('all')
ind=np.where( (gplt['PTF CNY']>=0) & (gplt['Ecozone BC L1']!=metaGP['LUT']['Ecozone BC L1']['CWH']) & (gplt['Deciduous L %N t0']<50) )[0]
N,mu,med,sig,se=gu.discres(gplt['Age Mean t0'][ind],gplt['Ctot L t0'][ind],bw,bin)
plt.plot(bin,mu,'-bo')
ind=np.where( (gplt['PTF CNY']>=0) & (gplt['Ecozone BC L1']!=metaGP['LUT']['Ecozone BC L1']['CWH']) & (gplt['Deciduous L %N t0']>50) )[0]
N,mu,med,sig,se=gu.discres(gplt['Age Mean t0'][ind],gplt['Ctot L t0'][ind],bw,bin)
plt.plot(bin,mu,'-gs')



bw=20; bin=np.arange(20,300,bw)
plt.close('all')
ind=np.where( (gplt['PTF CNY']>=0) & (gplt['Ecozone BC L1']==metaGP['LUT']['Ecozone BC L1']['IDF']) & (gplt['Deciduous L %N t0']<50) )[0]
N,mu,med,sig,se=gu.discres(gplt['Age Mean t0'][ind],gplt['Ctot L t0'][ind],bw,bin)
plt.plot(bin,mu,'-bo')
ind=np.where( (gplt['PTF CNY']>=0) & (gplt['Ecozone BC L1']==metaGP['LUT']['Ecozone BC L1']['IDF']) & (gplt['Deciduous L %N t0']>50) )[0]
N,mu,med,sig,se=gu.discres(gplt['Age Mean t0'][ind],gplt['Ctot L t0'][ind],bw,bin)
plt.plot(bin,mu,'-gs')

ind=np.where( (gplt['PTF CNY']==1) & (gplt['Ecozone BC L1']==metaGP['LUT']['Ecozone BC L1']['IDF']) & (gplt['Deciduous L %N t0']<50) )[0]
print(np.nanmean(gplt['Ctot Mort+Lost'][ind]))
ind=np.where( (gplt['PTF CNY']==1) & (gplt['Ecozone BC L1']==metaGP['LUT']['Ecozone BC L1']['IDF']) & (gplt['Deciduous L %N t0']>50) )[0]
print(np.nanmean(gplt['Ctot Mort+Lost'][ind]))


plt.plot(bin,mu,'-bo')
ind=np.where( (gplt['PTF CNY']>=0) & (gplt['Ecozone BC L1']==metaGP['LUT']['Ecozone BC L1']['IDF']) & (gplt['Deciduous L %N t0']>50) )[0]
N,mu,med,sig,se=gu.discres(gplt['Age Med t0'][ind],gplt['Ctot L t0'][ind],bw,bin)
plt.plot(bin,mu,'-gs')


bw=10; bin=np.arange(0,200,bw)
plt.close('all')
ind=np.where( (gplt['PTF CNY']>=0) & (gplt['Ecozone BC L1']==metaGP['LUT']['Ecozone BC L1']['IDF']) & (gplt['Deciduous L %N t0']<50) )[0]
N,mu,med,sig,se=gu.discres(gplt['cwd_gs_n'][ind],gplt['Ctot L t0'][ind],bw,bin)
plt.plot(bin,mu,'-bo')
ind=np.where( (gplt['PTF CNY']>=0) & (gplt['Ecozone BC L1']==metaGP['LUT']['Ecozone BC L1']['IDF']) & (gplt['Deciduous L %N t0']>50) )[0]
N,mu,med,sig,se=gu.discres(gplt['cwd_gs_n'][ind],gplt['Ctot L t0'][ind],bw,bin)
plt.plot(bin,mu,'-gs')

#%% Export summary by Plot Type

for k1 in metaGP['LUT']['Ecozone BC L1'].keys():
    d={}
    ind=np.where( (gplt['PTF CNY']==1) & (gplt['Ecozone BC L1']==metaGP['LUT']['Ecozone BC L1'][k1]) & (gplt['Deciduous L %N t0']<50) )[0]
    for k in gplt.keys():
        d[k]=np.round(np.nanmean(gplt[k][ind]),decimals=2)
    ind=np.where( (gplt['PTF CNY']==1) & (gplt['Ecozone BC L1']==metaGP['LUT']['Ecozone BC L1'][k1]) & (gplt['Deciduous L %N t0']>50) )[0]
    for k in gplt.keys():
        d[k]=np.append(d[k],np.round(np.nanmean(gplt[k][ind]),decimals=2))
    df=pd.DataFrame(d,index=[0,1])
    df.to_excel(r'C:\Users\rhember\Documents\Data\GroundPlots\PSP-NADB2\Processed\SummarySL_ByDecidFrac50_' + k1 + '.xlsx')