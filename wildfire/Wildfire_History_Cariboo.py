
#%% Import modules

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
from scipy import stats
from fcgadgets.macgyver import utilities_general as gu

#%%

d=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\Wildfire\Cariboo Dendro Wildfire Occurrence\Harvey et al 2017 fig 2.xlsx')

d['Year']=np.round(d['Year'])
d['Site']=np.round(d['Site'])

tv=np.arange(1750,2010)

#%%

u=np.unique(d['Site'])

v=np.zeros((tv.size,u.size))
for i in range(u.size):
    ind=np.where(d['Site']==u[i])[0]
    ic,ia,ib=np.intersect1d(tv,d['Year'][ind],return_indices=True)
    v[ia,i]=1

#%%

#plt.close('all')
fig,ax=plt.subplots(1)
plt.plot(tv,gu.movingave(np.sum(v,axis=1)/v.shape[1],30,'Historical'),lw=2)
