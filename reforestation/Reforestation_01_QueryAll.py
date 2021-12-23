'''

REFORESTATION (ALL)

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
from fcgadgets.cbrunner import cbrun_utilities

#%% Set figure properties

fs=6
params={'font.sans-serif':'Arial','font.size':fs,'axes.edgecolor':'black','axes.labelsize':fs,'axes.labelcolor':'black','axes.titlesize':fs,'axes.linewidth':0.5,'lines.linewidth':0.5,
        'text.color':'black','xtick.color':'black','xtick.labelsize':fs,'xtick.major.width':0.5,'xtick.major.size':3,'xtick.direction':'in','ytick.color':'black','ytick.labelsize':fs,
        'ytick.major.width':0.5,'ytick.major.size':3,'ytick.direction':'in','legend.fontsize':fs,'savefig.dpi':300,'savefig.transparent':True}
plt.rcParams.update(params)

#%% Define path

path=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20210930\Results.gdb'

#%% Query AT layer for all planting entries

# Initialize 
dAT={}
dAT['Year']=np.zeros(1000000)
dAT['Area']=np.zeros(1000000)
dAT['FSC']=np.array(['empty' for _ in range(dAT['Year'].size)],dtype=object)  
dAT['ACTIVITY_TREATMENT_UNIT_ID']=np.zeros(1000000)
dAT['HasSpatial']=np.zeros(1000000)

# Import data
cnt=0
with fiona.open(path,layer='RSLT_ACTIVITY_TREATMENT_SVW') as source:        
    for feat in source:        
        
        prp=feat['properties']            
        
        if (prp['SILV_BASE_CODE']=='PL') & (prp['SILV_TECHNIQUE_CODE']!='SE') & (prp['SILV_TECHNIQUE_CODE']!='CG') & (prp['SILV_METHOD_CODE']!='LAYOT') & \
            (prp['RESULTS_IND']=='Y') & (prp['ACTUAL_TREATMENT_AREA']!=None) & (prp['ATU_COMPLETION_DATE']!=None):    
            
            Year=int(prp['ATU_COMPLETION_DATE'][0:4]) 
            dAT['Year'][cnt]=Year
            dAT['Area'][cnt]=prp['ACTUAL_TREATMENT_AREA']
            
            if (prp['SILV_FUND_SOURCE_CODE']!=None):
                dAT['FSC'][cnt]=prp['SILV_FUND_SOURCE_CODE']
            
            dAT['ACTIVITY_TREATMENT_UNIT_ID'][cnt]=prp['ACTIVITY_TREATMENT_UNIT_ID']
            
            if feat['geometry']!=None:
                dAT['HasSpatial'][cnt]=1
            
            cnt=cnt+1

# Remove excess zeros
for k in dAT.keys():
    dAT[k]=dAT[k][0:cnt]

#%% Annual treatment area (by funding source)

# Get unique funding source codes
uFSC=np.unique(dAT['FSC'])

# Last year of completed record
YearLastComp=2021

ail={}
ail['Year']=np.arange(1970,YearLastComp+1,1)
ail['AreaTotal']=np.zeros(ail['Year'].size)
ail['AreaByFSC']=np.zeros((ail['Year'].size,uFSC.size))

# Calculate area by funding source code
for iAT in range(dAT['Year'].size):
    #if dAT['HasSpatial'][i]==0:
    #    continue
    it=np.where(ail['Year']==dAT['Year'][iAT])[0]
    ail['AreaTotal'][it]=ail['AreaTotal'][it]+dAT['Area'][iAT]
    for iFSC in range(uFSC.size):
        if dAT['FSC'][iAT]==uFSC[iFSC]:
            ail['AreaByFSC'][it,iFSC]=ail['AreaByFSC'][it,iFSC]+dAT['Area'][iAT]

#plt.plot(ail['Year'],ail['AreaTotal'],'o')

#%% Group funding sources

# Define which funding source codes are licensee vs. non-ob
ListOfNonObFSC=['FTL','FTM','RBM','RBL','FR','VG','FIL','FID','FIM','S','FRP','XXX','O','GFS','IR','FES','FCE','FCM']
ListOfLicenseeFSC=['BCT','LFP''IA','IR','VOI','SBF']

# Categorize by funding program
y=[None]*12
c=-1

itLast10=np.where( (ail['Year']>=YearLastComp-9) & (ail['Year']<=YearLastComp) )[0]
itLast4=np.where( (ail['Year']>=YearLastComp-3) & (ail['Year']<=YearLastComp) )[0]
itLast2=np.where( (ail['Year']>=YearLastComp-1) & (ail['Year']<=YearLastComp) )[0]


ind=np.where( (uFSC!='BCT') & (uFSC!='FTL') & (uFSC!='FTM') & (uFSC!='FCE') & (uFSC!='FCM') & (uFSC!='RBM') & (uFSC!='RBL') & (uFSC!='FR') & (uFSC!='VG') & \
             (uFSC!='FIL') & (uFSC!='FID') & (uFSC!='FIM') & (uFSC!='S') & (uFSC!='FES') & (uFSC!='XXX') & (uFSC!='LFP') & (uFSC!='GFS') & \
             (uFSC!='IR') & (uFSC!='O') & (uFSC!='FRP') )[0]
c=c+1; y[c]={}; y[c]['Name']='Other'; y[c]['Area']=np.sum(ail['AreaByFSC'][:,ind],axis=1); y[c]['Color']=[0.85,0.85,0.85]
uFSC[ind]
#np.sum(ail['AreaByFSC'][:,ind],axis=0,dtype=int)
#ail_oth_Last10=np.mean(y[c]['Area'][itLast10])
#ail_oth_Last4=np.mean(y[c]['Area'][itLast4])

ind=np.where( (uFSC=='FTL') | (uFSC=='FTM') )[0]
c=c+1; y[c]={}; y[c]['Name']='FFT'; y[c]['Area']=np.sum(ail['AreaByFSC'][:,ind],axis=1); y[c]['Color']=[0,0.6,1]
#ail_fft_Last10=np.mean(y[c]['Area'][itLast10])
#ail_fft_Last4=np.mean(y[c]['Area'][itLast4])

ind=np.where( (uFSC=='BCT') )[0]
c=c+1; y[c]={}; y[c]['Name']='BCTS'; y[c]['Area']=np.sum(ail['AreaByFSC'][:,ind],axis=1); y[c]['Color']=[0,0.25,0]

ind=np.where( (uFSC=='RBM') | (uFSC=='RBL') )[0]
c=c+1; y[c]={}; y[c]['Name']='FRBC'; y[c]['Area']=np.sum(ail['AreaByFSC'][:,ind],axis=1); y[c]['Color']=[0,0.4,0]

ind=np.where( (uFSC=='FR') )[0]
c=c+1; y[c]={}; y[c]['Name']='FRDA'; y[c]['Area']=np.sum(ail['AreaByFSC'][:,ind],axis=1); y[c]['Color']=[0.5,0.5,0.7]

ind=np.where( (uFSC=='VG') )[0]
c=c+1; y[c]={}; y[c]['Name']='Ministry Outstanding'; y[c]['Area']=np.sum(ail['AreaByFSC'][:,ind],axis=1); y[c]['Color']=[0.9,0.9,0.75]

ind=np.where( (uFSC=='FIL') | (uFSC=='FID') | (uFSC=='FIM') )[0]
c=c+1; y[c]={}; y[c]['Name']='FIA'; y[c]['Area']=np.sum(ail['AreaByFSC'][:,ind],axis=1); y[c]['Color']=[0.5,1,1]

ind=np.where( (uFSC=='S') )[0]
c=c+1; y[c]={}; y[c]['Name']='Section 88'; y[c]['Area']=np.sum(ail['AreaByFSC'][:,ind],axis=1); y[c]['Color']=[0.95,0.8,0.4]

ind=np.where( (uFSC=='FRP') )[0]
c=c+1; y[c]={}; y[c]['Name']='FRP Section 108'; y[c]['Area']=np.sum(ail['AreaByFSC'][:,ind],axis=1); y[c]['Color']=[0,0,0.75]

ind=np.where( (uFSC=='XXX') | (uFSC=='O') | (uFSC=='GFS') | (uFSC=='IR') )[0]
c=c+1; y[c]={}; y[c]['Name']='Royalties (XXX, O, GFS, IR)'; y[c]['Area']=np.sum(ail['AreaByFSC'][:,ind],axis=1); y[c]['Color']=[1,0.4,0.85]

ind=np.where( (uFSC=='FES') )[0]
c=c+1; y[c]={}; y[c]['Name']='FES'; y[c]['Area']=np.sum(ail['AreaByFSC'][:,ind],axis=1); y[c]['Color']=[1,1,0]

ind=np.where( (uFSC=='FCE') | (uFSC=='FCM') )[0]
c=c+1; y[c]={}; y[c]['Name']='FCI'; y[c]['Area']=np.sum(ail['AreaByFSC'][:,ind],axis=1); y[c]['Color']=[0,0.9,0]
#ail_fci_Last4=np.mean(y[c]['Area'][itLast4])


#%% Plot time series of area treated by funding source

plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(15.5,6.5)); 
Tot=0
for i in range(len(y)):    
    plt.bar(ail['Year'],y[i]['Area']/1000,0.8,bottom=Tot,facecolor=y[i]['Color'],label=y[i]['Name'])
    Tot=Tot+y[i]['Area']/1000

ax.set(position=[0.06,0.12,0.92,0.86],ylabel='Treatment area (hectares x 1000)',xlabel='Time, years',xlim=[ail['Year'][0]-0.75,ail['Year'][-1]+2+.75],xticks=np.arange(1950,2035,5))
plt.legend(frameon=False,facecolor=[1,1,1],labelspacing=0.25)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')
#gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Reforestation\ReforestationNonOb\RefNonOb_AIL_ByFS_IncludingOblig','png',900)

#%% Random sample 

# Subsampling rate
ail['samp_rate_mp']=0.2

# Get IDs for random sample that will be included
ail['id_atu_subsample']=np.zeros(100000,dtype=np.int32)
ail['AreaSubSample']=np.zeros(ail['Year'].size)
n=np.zeros(ail['Year'].size)
cnt=0
for iY in range(ail['Year'].size):
    indY=np.where(dAT['Year']==ail['Year'][iY])[0]
    r=np.random.permutation(indY.size)
    nKeep=int(np.ceil(ail['samp_rate_mp']*indY.size))
    for iMP in range(nKeep):
        ail['id_atu_subsample'][cnt]=dAT['ACTIVITY_TREATMENT_UNIT_ID'][indY[r[iMP]]]
        ail['AreaSubSample'][iY]=ail['AreaSubSample'][iY]+dAT['Area'][indY[r[iMP]]]
        cnt=cnt+1
    n[i]=nKeep

ail['id_atu_subsample']=ail['id_atu_subsample'][0:cnt]

# Annual area expansion factor
ail['AEF_MP']=ail['AreaTotal']/ail['AreaSubSample']
plt.plot(ail['Year'],ail['AEF_MP'],'-o')

# Save total area time series for area expansion factor
gu.opickle(r'D:\Data\FCI_Projects\SummaryReforestation\Inputs\AnnualImplementationLevel.pkl',ail)


# ail0=gu.ipickle(r'D:\Data\FCI_Projects\SummaryReforestation\Inputs\AnnualImplementationLevel.pkl')
#plt.plot(ail0['Year'],ail0['AEF_MP'],'-o')