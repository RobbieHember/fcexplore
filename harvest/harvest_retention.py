#%% Import modules

from os import listdir
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.io
import copy
import geopandas as gpd
import fiona
import scipy.io as spio
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats, linalg

import fcgadgets.macgyver.utilities_general as gu
import fcgadgets.macgyver.utilities_gis as gis
import fcgadgets.bc1ha.bc1ha_utilities as u1ha

#%% Import data

meta=u1ha.Init()
meta=u1ha.ImportLUTs(meta)
gp=gu.SetGraphics('Manuscript')
zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])

vList=['harv_yearlast_cc','biomass_glob','bgcz']
z=u1ha.Import_Raster(meta,[],vList)

#%% Buffer around harvest

Mask=np.zeros(zRef['Data'].shape,dtype='int8')
Mask[z['harv_yearlast_cc']['Data']>0]=1
MaskB=gis.BufferRasterMask(Mask,1)

#%%

ikp=np.where( (z['bgcz']['Data']==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['CWH']) & (z['harv_yearlast_cc']['Data']>0) & (z['harv_yearlast_cc']['Data']<2010) )
ikp=np.where( (z['bgcz']['Data']==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['CWH']) & (z['harv_yearlast_cc']['Data']>0) & (z['harv_yearlast_cc']['Data']<2010) )

bw=10; bin=np.arange(-30,250+bw,bw); x=2010-z['harv_yearlast_cc']['Data'][ikp]
y=z['biomass_glob']['Data'][ikp]

N,mu,med,sig,se=gu.discres(x,y,bw,bin)
plt.close('all')
plt.plot(bin,mu,'bo')

#%%

ikp=np.where( (z['bgcz']['Data']==meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE']['CWH']) & (z['harv_yearlast_cc']['Data']==1990) )

plt.hist(z['biomass_glob']['Data'][ikp])



