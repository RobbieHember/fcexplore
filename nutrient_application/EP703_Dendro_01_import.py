'''
EP703 RETROSPECTIVE MONITORING STUDY
'''

#%% Import modules

from os import listdir
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from fcgadgets.macgyver import utilities_general as gu

from scipy.optimize import curve_fit
import scipy.io as spio
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats, linalg

#%% Graphics parameters

fs=6
params={'font.sans-serif':'Arial','font.size':fs,'axes.edgecolor':'black','axes.labelsize':fs,'axes.labelcolor':'black','axes.titlesize':fs,'axes.linewidth':0.5,'lines.linewidth':0.5,
        'text.color':'black','xtick.color':'black','xtick.labelsize':fs,'xtick.major.width':0.5,'xtick.major.size':3,'xtick.direction':'in','ytick.color':'black','ytick.labelsize':fs,
        'ytick.major.width':0.5,'ytick.major.size':3,'ytick.direction':'in','legend.fontsize':fs,'savefig.dpi':300,'savefig.transparent':True}
plt.rcParams.update(params)

#%% Preapre project

meta={}
meta['Paths']={}
meta['Paths']['Inputs']=r'C:\Users\rhember\Documents\Data\TreeRings\EP703\Processed\Cores'

#%%

# Import tree-level observations
dTL=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\TreeRings\EP703\Processed\Cores\EP703_RM_Master_stacked.xlsx')

# Import tree summary file 
dTS=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\TreeRings\EP703\Raw\Field Data received 20210125\RM Tree Summary.xlsx')

# Import plot summary file
dPS=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\TreeRings\EP703\Raw\Field Data received 20210125\RM Plot Summary.xlsx')


#%% Add derived variables

# Add unique tree ID
dTL['ID_Tree Unique']=np.zeros(dTL['TRW'].size,dtype=int)
cnt=0
u=np.unique(np.column_stack((dTL['ID_Inst'],dTL['ID_Plot'],dTL['ID_Tree'])),axis=0)
for iU in range(u.shape[0]):
    ind=np.where( (dTL['ID_Inst']==u[iU,0]) & (dTL['ID_Plot']==u[iU,1]) & (dTL['ID_Tree']==u[iU,2]) )[0]
    dTL['ID_Tree Unique'][ind]=cnt
    cnt=cnt+1

# Add info from site and tree summary files and add some derived variables
dTL['BGC SS']=np.zeros(dTL['TRW'].size)
dTL['Dob 2020']=np.zeros(dTL['TRW'].size)
dTL['H 2020']=np.zeros(dTL['TRW'].size)
dTL['Year First']=np.zeros(dTL['TRW'].size)
dTL['Year Last']=np.zeros(dTL['TRW'].size)
dTL['Time since first ring']=np.zeros(dTL['TRW'].size)
dTL['Dib']=np.zeros(dTL['TRW'].size)
dTL['Dib Last']=np.zeros(dTL['TRW'].size)
dTL['Bsw']=np.zeros(dTL['TRW'].size)
dTL['Gsw']=np.zeros(dTL['TRW'].size)
dTL['RGR']=np.zeros(dTL['TRW'].size)
dTL['TRW Standardized RR']=np.zeros(dTL['TRW'].size)
dTL['BAI Standardized RR']=np.zeros(dTL['TRW'].size)
dTL['Bsw Standardized RR']=np.zeros(dTL['TRW'].size)
dTL['Gsw Standardized AD']=np.zeros(dTL['TRW'].size)
dTL['Gsw Standardized RR']=np.zeros(dTL['TRW'].size)
dTL['RGR Standardized RR']=np.zeros(dTL['TRW'].size)

uTID=np.unique(dTL['ID_Tree Unique'])
for iU in range(uTID.size):
    
    # Index to unique tree ID
    indUT=np.where(dTL['ID_Tree Unique']==uTID[iU])[0]
    
    # Define first and last year of measurement for each tree
    iWithData=np.where(dTL['Dib'][indUT]>-1.0)[0]
    dTL['Year First'][indUT]=dTL['Year'][indUT[iWithData[0]]]
    dTL['Year Last'][indUT]=dTL['Year'][indUT[iWithData[-1]]]
    dTL['Time since first ring'][indUT[iWithData]]=np.arange(1,iWithData.size+1,1)
    
    # Index to site summary info
    indPS=np.where( (dPS['ID_Inst']==dTL['ID_Inst'][indUT[0]]) & (dPS['ID_Plot']==dTL['ID_Plot'][indUT[0]]) )[0]
    
    # Site series
    dTL['BGC SS'][indUT]=dPS['BGC site series'][indPS]
    
    # Index to tree summary file
    indTS=np.where( (dTS['ID_Inst']==dTL['ID_Inst'][indUT[0]]) & (dTS['ID_Plot']==dTL['ID_Plot'][indUT[0]]) & (dTS['ID_Tree']==dTL['ID_Tree'][indUT[0]]) )[0]
        
    if indTS.size>0:
        try:
            dTL['Dob 2020'][indUT]=dTS['DBH_2020 (cm)'][indTS].astype(np.float)   
        except:
            print('Something odd about DBH 2020')
        dTL['H 2020'][indUT]=dTS['Ht_2020 (m)'][indTS][0]
    else:
        print('No crosswalk found with tree summary')
        continue

    # Inside-bark diameter (cm)
    dTL['Dib'][indUT]=0.1*np.cumsum(2*np.nan_to_num(dTL['TRW'][indUT]))
    
    # Inside-bark diameter for last year (cm) 
    # *** Only populated if last year >= 2019 (for comparison with 2020 tape measurement)
    if dTL['Year'][indUT[iWithData[-1]]]>=2019:
        dTL['Dib Last'][indUT[iWithData[-1]]]=dTL['Dib'][indUT[iWithData[-1]]]
    
    # Stemwood biomass from allometric function (kgDM) - DBH only (Ung et al. 2008)
    dTL['Bsw'][indUT]=0.0204*dTL['Dib'][indUT]**2.6974
    
    # Stemwood biomass growth (kgDM/yr)
    dTL['Gsw'][indUT[1:]]=np.diff(dTL['Bsw'][indUT])
    
    # Relative growth rate
    dTL['RGR'][indUT[1:]]=np.log(dTL['Bsw'][indUT[1:]])-np.log(dTL['Bsw'][indUT[0:-1]])

    # Standardized growth relative to mean growth during a period leading up to N application
    ind0=np.where( (dTL['ID_Tree Unique']==uTID[iU]) & (dTL['Year']>=1971-5) & (dTL['Year']<=1970) )[0]
    dTL['TRW Standardized RR'][indUT]=dTL['TRW'][indUT]/np.mean(dTL['TRW'][ind0])
    dTL['BAI Standardized RR'][indUT]=dTL['BAI'][indUT]/np.mean(dTL['BAI'][ind0])
    dTL['Bsw Standardized RR'][indUT]=dTL['Bsw'][indUT]/np.mean(dTL['Bsw'][ind0])
    dTL['Gsw Standardized AD'][indUT]=dTL['Gsw'][indUT]-np.mean(dTL['Gsw'][ind0])
    dTL['Gsw Standardized RR'][indUT]=dTL['Gsw'][indUT]/np.maximum(0.1,np.mean(dTL['Gsw'][ind0]))
    dTL['RGR Standardized RR'][indUT]=dTL['RGR'][indUT]/np.mean(dTL['RGR'][ind0])
    
# Change name of treatment
dTL['Treatment']=dTL['Treatment.EP703']

# Delete fields
del dTL['Unnamed: 0'],dTL['Treatment.EP703']

# Define unique installations
uInst=np.unique(dTL['ID_Inst'])

# Define unique site series
uSS=np.unique(dTL['BGC SS'])

#%% Comparison between Dob and Dib from cores

# Isolate data
ikp=np.where( (dTL['Dob 2020']>0) & (dTL['Dib Last']>0) )[0]

# Fit linear best fit relationship
y=dTL['Dib Last'][ikp]
x=dTL['Dob 2020'][ikp]
x1=sm.tools.tools.add_constant(x)
md=sm.OLS(y,x1).fit()
md.summary()
xhat=np.linspace(np.min(x1[:,1]),np.max(x1[:,1]),10)
yhat=md.predict(np.c_[np.ones(xhat.size),xhat])


# Plot
plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(8.5,8.5))
ax.plot([0,90],[0,90],'k-',lw=2,color=[0.8,0.8,0.8])
ax.plot(x,y,'o',ms=4,mec='w',mfc='k')
ax.plot(xhat,yhat,'r-',lw=1.25)
ax.set(xlim=[0,80],ylim=[0,80],xlabel='Tape-based Dob in 2020 (cm)',ylabel='Core-based Dib in 2019 (cm)')
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')
plt.tight_layout()
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\EP703\TreeRings\DiameterComparison','png',150)

