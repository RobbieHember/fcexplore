
'''
HARVEST INFO BY TIMBER MARK - ANALYSIS
'''

#%% Import modules

from os import listdir
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from fcgadgets.macgyver import utilities_general as gu
import fiona
from fcgadgets.macgyver import utilities_inventory as invu

#%% Parameters

gp=gu.SetGraphics('Manuscript')
gradeL=['1','2','3','4','5','6','7','8','B','C','D','E','F','G','H','L','M','U','W','X','Y','Z']

#%% Import data

dTM=gu.ipickle(r'C:\Users\rhember\Documents\Data\Harvest\HBS\HavestSummary_ByTM.pkl')

#%% Derived variables

# Volume felled (Volume delivered + Waste)
#dTM['V Felled m3/ha']=dTM['V Logs Delivered Abs m3/ha']+dTM['V NonLog Delivered Abs m3/ha']+dTM['Waste Total m3/ha']
dTM['V Felled m3/ha']=dTM['V Logs Delivered Abs m3/ha']+dTM['Waste Total m3/ha']

# Waste ratio
dTM['WR']=dTM['Waste Total m3/ha']/dTM['V Felled m3/ha']*100

# Sawlog recovery ratio
dTM['Sawlog Ratio']=(dTM['V Logs Delivered Grade 1 Abs m3/ha']+dTM['V Logs Delivered Grade 2 Abs m3/ha'])/dTM['V Logs Delivered m3/ha']*100
#dTM['Sawlog Ratio']=(dTM['V Logs Delivered Grade 1 Abs m3/ha']+dTM['V Logs Delivered Grade 2 Abs m3/ha'])/dTM['V Felled m3/ha']*100

# Region
dTM['Reg']=np.array(['Interior' for _ in range(dTM['TM'].size)],dtype=object)
ind=np.where( (dTM['District']=='DCR') | (dTM['District']=='DNI') | (dTM['District']=='DSI') | \
             (dTM['District']=='DSC') | (dTM['District']=='DSQ') | (dTM['District']=='DQC') | \
             (dTM['District']=='DCK') )[0]
dTM['Reg'][ind]='Coast'

#%% Relationship between net volumes from different sources

list(dTM.keys())

ind=np.where( (dTM['Cruise_V Net (m3/ha)']>0) & (dTM['V Logs Delivered Abs m3/ha']<1500) )[0]

x=dTM['Cruise_V Net (m3/ha)'][ind]
y=dTM['V Logs Delivered Abs m3/ha'][ind]
rs,txt=gu.GetRegStats(x,y)

plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(8,8))
ax.plot(x,y,'b.',mew=0.25,ms=3,mfc='b',mec='w')
ax.plot([0,2000],[0,2000],'k-',lw=2,color=[0.8,0.8,0.8])
ax.text(1300,1300,'1:1',fontsize=12,ha='center',va='center')
ax.plot(rs['xhat'],rs['yhat'],'r--',lw=1.5,color=[0,0,0])
ax.text(1400,200,txt,fontsize=8,ha='right',va='center')
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=2)
ax.set(xlabel='Net vol from cruise (m3/ha)',ylabel='Net volume delivered (m3/ha)',xlim=[0,1500],ylim=[0,1500])
plt.tight_layout()
#gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Harvest\Timber Cruise Info\TimeSeries','png',900)

#%%

list(dTM.keys())

ind=np.where( (dTM['Reg']=='Interior') & (dTM['Cruise_Year']>=2014) & (dTM['Cruise_V Net (m3/ha)']>0) & (dTM['V Logs Delivered Abs m3/ha']<1500) )[0]

x=dTM['Cruise_Pct Dead Net'][ind]

#y=dTM['Cruise_V Gross (m3/ha)'][ind]
#y=dTM['Cruise_PCT Net'][ind]
#y=dTM['Sawlog Ratio'][ind]
#y=dTM['WR'][ind]
#y=dTM['Waste Total m3/ha'][ind]
#y=dTM['Waste Sawlog m3/ha'][ind]
#y=dTM['Waste GradeY/4 m3/ha'][ind]
#y=dTM['Waste Accumulation %'][ind]
#y=dTM['Waste Dispersed %'][ind]
#y=dTM['Waste Standing %'][ind]
#y=dTM['Cruise_Age'][ind]

plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(8,8))

y=dTM['Cruise_V Gross (m3/ha)'][ind]
bw=20; bin=np.arange(0,120,bw)
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
ax.plot(bin,mu,'go',ms=5,mfc='w',mec='g')

y=dTM['Cruise_V Net (m3/ha)'][ind]
bw=20; bin=np.arange(0,120,bw)
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
ax.plot(bin,mu,'rs',ms=5,mfc='w',mec='r')

y=dTM['V Logs Delivered Abs m3/ha'][ind]
bw=20; bin=np.arange(0,120,bw)
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
ax.plot(bin,mu,'b^',ms=4,mfc='w',mec='b')

ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=2)
ax.set(xlabel='Dead volume component (%)',ylabel='Net volume delivered (m3/ha)')
plt.tight_layout()
#gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Harvest\Timber Cruise Info\TimeSeries','png',900)