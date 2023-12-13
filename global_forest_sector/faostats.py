
'''
Global Forest sector - Mine FAOSTATS for info

'''

#%% Prepare session

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd
import openpyxl
import copy
import gc as garc
import time
from scipy import stats
from fcgadgets.macgyver import utilities_general as gu

#%% Set figure properties

fs=6
params={'font.sans-serif':'Arial','font.size':fs,'axes.edgecolor':'black','axes.labelsize':fs,'axes.labelcolor':'black','axes.titlesize':fs,'axes.linewidth':0.5,'lines.linewidth':0.5,
        'text.color':'black','xtick.color':'black','xtick.labelsize':fs,'xtick.major.width':0.5,'xtick.major.size':3,'xtick.direction':'in','ytick.color':'black','ytick.labelsize':fs,
        'ytick.major.width':0.5,'ytick.major.size':3,'ytick.direction':'in','legend.fontsize':fs,'savefig.dpi':300,'savefig.transparent':True}
plt.rcParams.update(params)

#%% Import data

# Crosswalk between countries and regions
dfCW=pd.read_excel(r'C:\Users\rhember\Documents\Data\FAOSTAT\CountryRegionCrosswalk.xlsx')

# Forest area
dFA=pd.read_excel(r'C:\Users\rhember\Documents\Data\FAOSTAT\FAOSTAT Summary.xlsx',sheet_name='Area Forest Raw')

# Forest area primary forest
dFAP=pd.read_excel(r'C:\Users\rhember\Documents\Data\FAOSTAT\FAOSTAT Summary.xlsx',sheet_name='Area Primary Forest')

# Forest conversion area
dFCA=pd.read_excel(r'C:\Users\rhember\Documents\Data\FAOSTAT\FAOSTAT Summary.xlsx',sheet_name='Forest Conversion Raw')

# Forest products
dFP=pd.read_excel(r'C:\Users\rhember\Documents\Data\FAOSTAT\FAOSTAT Summary.xlsx',sheet_name='Products Raw')


#%% Compile regions

rL=['Europe','Northern America','South America','Asia','Oceania','Africa']
d={}
for r in rL:
    d[r]={'Area Forest':0.0,'Area Primary Forest':0.0,'Area Deforestation':0.0,'Roundwood Prod':0.0,'Industrial Roundwood Prod':0.0,'Sawnwood Prod':0.0}

# Add forest area (ha)
for i in range(len(dFA)):
    if dFA.loc[i,'Year']!=2020:
        continue
    ind=np.where(dfCW['name']==dFA.loc[i,'Area'])[0]
    if ind.size==0:
        continue
    if np.isnan(dFA.loc[i,'Value'])==True:
        continue
    reg=dfCW.loc[ind,'region'].values[0]
    sreg=dfCW.loc[ind,'sub-region'].values[0]
    ireg=dfCW.loc[ind,'intermediate-region'].values[0]
    if reg not in rL:
        if sreg not in rL:
            if ireg not in rL:
                continue
    if (sreg=='Northern America'):
        reg=sreg
    elif (ireg=='South America'):
        reg=ireg
    d[reg]['Area Forest']=d[reg]['Area Forest']+dFA.loc[i,'Value']*1000
    
# Add primary forest area (2017 is the last year with estimates)
for i in range(len(dFAP)):
    if dFAP.loc[i,'Year']!=2017:
        continue
    ind=np.where(dfCW['name']==dFAP.loc[i,'Area'])[0]
    if ind.size==0:
        continue
    if np.isnan(dFAP.loc[i,'Value'])==True:
        continue
    reg=dfCW.loc[ind,'region'].values[0]
    sreg=dfCW.loc[ind,'sub-region'].values[0]
    ireg=dfCW.loc[ind,'intermediate-region'].values[0]
    if reg not in rL:
        if sreg not in rL:
            if ireg not in rL:
                continue
    if (sreg=='Northern America'):
        reg=sreg
    elif (ireg=='South America'):
        reg=ireg
    d[reg]['Area Primary Forest']=d[reg]['Area Primary Forest']+dFAP.loc[i,'Value']*1000
 
# Add area deforestation
for i in range(len(dFCA)):
    if (dFCA.loc[i,'Year']!=2020) | (dFCA.loc[i,'Element']!='Area'):
        continue
    ind=np.where(dfCW['name']==dFCA.loc[i,'Area'])[0]
    if ind.size==0:
        continue
    if np.isnan(dFCA.loc[i,'Value'])==True:
        continue
    reg=dfCW.loc[ind,'region'].values[0]
    sreg=dfCW.loc[ind,'sub-region'].values[0]
    ireg=dfCW.loc[ind,'intermediate-region'].values[0]
    if reg not in rL:
        if sreg not in rL:
            if ireg not in rL:
                continue
    if (sreg=='Northern America'):
        reg=sreg
    elif (ireg=='South America'):
        reg=ireg
    d[reg]['Area Deforestation']=d[reg]['Area Deforestation']+dFCA.loc[i,'Value']*1000
      
# Add roundwood production (m3/yr)
for i in range(len(dFP)):
    if (dFP.loc[i,'Item']!='Roundwood') | (dFP.loc[i,'Element']!='Production'):
        continue
    ind=np.where(dfCW['name']==dFP.loc[i,'Area'])[0]
    if ind.size==0:
        continue
    reg=dfCW.loc[ind,'region'].values[0]
    sreg=dfCW.loc[ind,'sub-region'].values[0]
    ireg=dfCW.loc[ind,'intermediate-region'].values[0]
    if reg not in rL:
        if sreg not in rL:
            if ireg not in rL:
                continue
    if (sreg=='Northern America'):
        reg=sreg
    elif (ireg=='South America'):
        reg=ireg
    d[reg]['Roundwood Prod']=d[reg]['Roundwood Prod']+dFP.loc[i,'Y2020']

# Add industrial roundwood production (m3/yr)
for i in range(len(dFP)):
    if (dFP.loc[i,'Item']!='Industrial roundwood') | (dFP.loc[i,'Element']!='Production'):
        continue
    ind=np.where(dfCW['name']==dFP.loc[i,'Area'])[0]
    if ind.size==0:
        continue
    reg=dfCW.loc[ind,'region'].values[0]
    sreg=dfCW.loc[ind,'sub-region'].values[0]
    ireg=dfCW.loc[ind,'intermediate-region'].values[0]
    if reg not in rL:
        if sreg not in rL:
            if ireg not in rL:
                continue
    if (sreg=='Northern America'):
        reg=sreg
    elif (ireg=='South America'):
        reg=ireg
    d[reg]['Industrial Roundwood Prod']=d[reg]['Industrial Roundwood Prod']+dFP.loc[i,'Y2020']
       
# Add sawnwood production (m3/yr)
for i in range(len(dFP)):
    if (dFP.loc[i,'Item']!='Sawnwood') | (dFP.loc[i,'Element']!='Production'):
        continue
    ind=np.where(dfCW['name']==dFP.loc[i,'Area'])[0]
    if ind.size==0:
        continue
    reg=dfCW.loc[ind,'region'].values[0]
    sreg=dfCW.loc[ind,'sub-region'].values[0]
    ireg=dfCW.loc[ind,'intermediate-region'].values[0]
    if reg not in rL:
        if sreg not in rL:
            if ireg not in rL:
                continue
    if (sreg=='Northern America'):
        reg=sreg
    elif (ireg=='South America'):
        reg=ireg
    d[reg]['Sawnwood Prod']=d[reg]['Sawnwood Prod']+dFP.loc[i,'Y2020']

# Add fuelwood
for reg in d.keys():
    d[reg]['Fuelwood']=d[reg]['Roundwood Prod']-d[reg]['Industrial Roundwood Prod']

# Save to spreadsheet
df=pd.DataFrame(d)
df.to_excel(r'C:\Users\rhember\Documents\Data\FAOSTAT\SummaryStats.xlsx')



#%% Plot

fw=np.array([])
irw=np.array([])
for reg in d.keys():
    fw=np.append(fw,d[reg]['Fuelwood'])
    irw=np.append(irw,d[reg]['Industrial Roundwood Prod'])

plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,6))
ax.bar(np.arange(fw.size),fw/1e6,fc=[0.7,0.8,1])
ax.bar(np.arange(fw.size),irw/1e6,bottom=fw/1e6,fc=[0.8,1,0.7])
ax.set(position=[0.09,0.12,0.88,0.82],xticks=np.arange(fw.size),xticklabels=np.array(list(d.keys())),ylabel='Million m3/year')
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=1.5)
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Global Forest Sector\RegionalSummaryHarvest','png',500)
    