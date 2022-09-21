
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
from fcgadgets.cbrunner import cbrun_utilities as cbu

#%% Start year

YrStart=2007

#%% Graphics parameters

fs=7
params={'font.sans-serif':'Arial','font.size':fs,'axes.edgecolor':'black','axes.labelsize':fs,'axes.labelcolor':'black','axes.titlesize':fs,'axes.linewidth':0.5,'lines.linewidth':0.5,
        'text.color':'black','xtick.color':'black','xtick.labelsize':fs,'xtick.major.width':0.5,'xtick.major.size':3,'xtick.direction':'in','ytick.color':'black','ytick.labelsize':fs,
        'ytick.major.width':0.5,'ytick.major.size':3,'ytick.direction':'in','legend.fontsize':fs,'savefig.dpi':300,'savefig.transparent':True}
plt.rcParams.update(params)

#%% Import forest tenure layer

# List layers
pthin=r'C:\Users\rhember\Documents\Data\ForestInventory\LandUse\20220422\LandUse.gdb'
# fiona.listlayers(pthin)
    
N=5000000
dFT={}
dFT['TIMBER_MARK']=np.array(['' for _ in range(N)],dtype=object)
dFT['OPENING_ID']=np.zeros(N)
dFT['PLANNED_GROSS_BLOCK_AREA']=np.zeros(N)
dFT['PLANNED_NET_BLOCK_AREA']=np.zeros(N)
dFT['DISTURBANCE_GROSS_AREA']=np.zeros(N)
dFT['DISTURBANCE_START_YEAR']=np.zeros(N)
dFT['DISTURBANCE_END_YEAR']=np.zeros(N)
#dFT['FEATURE_AREA']=np.zeros(N)
#dFT['CLIENT_NUMBER']=np.array(['' for _ in range(N)],dtype=object)
#dFT['ADMIN_DISTRICT_CODE']=np.array(['' for _ in range(N)],dtype=object)
cnt=0
with fiona.open(pthin,layer='FTEN_CUT_BLOCK_OPEN_ADMIN') as source:
    for feat in source:
        p=feat['properties']

        if p['TIMBER_MARK']==None:
            continue
        
        if p['DISTURBANCE_START_DATE']==None:
            Year=0
        else:
            Year=int(p['DISTURBANCE_START_DATE'][0:4])
        if p['DISTURBANCE_END_DATE']==None:
            Year2=0
        else:
            Year2=int(p['DISTURBANCE_END_DATE'][0:4])    
        #if Year<YrStart:
        #    continue
        dFT['TIMBER_MARK'][cnt]=p['TIMBER_MARK']
        #dFT['CLIENT_NUMBER'][cnt]=p['CLIENT_NUMBER']
        #dFT['ADMIN_DISTRICT_CODE'][cnt]=p['ADMIN_DISTRICT_CODE']
        #dFT['FEATURE_AREA'][cnt]=p['FEATURE_AREA']/10000
        dFT['DISTURBANCE_GROSS_AREA'][cnt]=p['DISTURBANCE_GROSS_AREA']
        dFT['PLANNED_NET_BLOCK_AREA'][cnt]=p['PLANNED_NET_BLOCK_AREA']
        if p['OPENING_ID']!=None:
            dFT['OPENING_ID'][cnt]=p['OPENING_ID']
        dFT['DISTURBANCE_START_YEAR'][cnt]=Year
        dFT['DISTURBANCE_END_YEAR'][cnt]=Year2
        cnt=cnt+1

# Truncate
for k in dFT.keys():
    dFT[k]=dFT[k][0:cnt]
 
# QA
flg=0
if flg==1:
    ind=np.where(dFT['TIMBER_MARK']=='61/206')[0]
    dFT['OPENING_ID'][ind]
    dFT['PLANNED_NET_BLOCK_AREA'][ind]
    dFT['DISTURBANCE_GROSS_AREA'][ind]
    np.sum(dFT['PLANNED_NET_BLOCK_AREA'][ind])

#%% Database by unique timber marks
    
# Unique timber marks
uTM=np.unique(dFT['TIMBER_MARK'])

dTM={}
dTM['TM']=uTM
dTM['N Openings']=np.zeros(uTM.size)
dTM['DISTURBANCE_GROSS_AREA']=np.zeros(uTM.size)
dTM['PLANNED_NET_BLOCK_AREA']=np.zeros(uTM.size)
dTM['Year First']=np.zeros(uTM.size)
dTM['Year Last']=np.zeros(uTM.size)
dTM['District']=np.array(['' for _ in range(dTM['TM'].size)],dtype=object)
#dTM['ADMIN_DISTRICT_CODE']=np.array(['' for _ in range(uTM.size)],dtype=object)
#dTM['FEATURE_AREA']=np.zeros(uTM.size)

for iTM in range(uTM.size):
    ind=np.where(dFT['TIMBER_MARK']==uTM[iTM])[0]
    #dTM['CLIENT_NUMBER'][iTM]=dFT['CLIENT_NUMBER'][ind]
    dTM['N Openings'][iTM]=np.unique(dFT['OPENING_ID'][ind]).size
    #dTM['ADMIN_DISTRICT_CODE'][iTM]=dFT['ADMIN_DISTRICT_CODE'][ind[0]]
    #dTM['FEATURE_AREA'][iTM]=np.round(np.sum(dFT['FEATURE_AREA'][ind]))
    dTM['DISTURBANCE_GROSS_AREA'][iTM]=np.round(np.nansum(dFT['DISTURBANCE_GROSS_AREA'][ind]))
    dTM['PLANNED_NET_BLOCK_AREA'][iTM]=np.round(np.sum(dFT['PLANNED_NET_BLOCK_AREA'][ind]))
    dTM['Year First'][iTM]=np.min(dFT['DISTURBANCE_START_YEAR'][ind])
    dTM['Year Last'][iTM]=np.max(dFT['DISTURBANCE_END_YEAR'][ind])
 
# Plot
flg=0
if flg==1:
    tv=np.arange(1990,2022,1)
    y=np.zeros(tv.size)
    for iT in range(tv.size):
        ind=np.where(dFT['Year First']==tv[iT])[0]
        y[iT]=np.sum(dFT['FEATURE_AREA'][ind])
    plt.plot(tv,y,'-bo')
    
#%% Import HBS 
# The system appears to be down often
# -> Date of Invoice, "Billing processed" only 
# Can only do 12months at a time

# Definition of volume types (from Ross Pavan)
#   1.	Normal:  the volume delivered to the sawmill (the volume calculation is based on weight to volume ratios)
#   2.	Cruise Based:  the volume is derived from the cruise compilation (the billed volume never exceeds the cruise volume for a cutting permit, for a timber sale license the billed volume will never exceed the cruise volume plus any additional volume added after the fact (example: external R/W volume)
#   3.	Waste: the volume is derived from the waste assessment of a scale based cutting authority
#   4.	Beachcomb: volume scaled because of the beachcombing efforts of an individual or salvage logger 


##------------------------------------------------------------------------------
## Populate from annual summary files
#binY=np.arange(YrStart,2021,1)
#N=2000000
#cnt=0
#tm_hbs=np.array(['' for _ in range(N)],dtype=object)
#for iY in range(binY.size):    
#    print(binY[iY])    
#    df=pd.read_csv(r'C:\Users\rhember\Documents\Data\Harvest\HBS\HBS ' + str(binY[iY]) + '.csv',header=None)    
#    tm=df.iloc[:,0].to_numpy()
#    tm=tm.astype(str)
#    tm=np.unique(tm)
#    tm_hbs[cnt:cnt+tm.size]=tm
#    cnt=cnt+tm.size
#
#tm_hbs=tm_hbs[0:cnt]
#tm_flg=np.zeros(tm_hbs.size)
#
#for i in range(tm_hbs.size):
#    ind=np.where(dTM['TM']==tm_hbs[i])[0]
#    if ind.size>0:
#        tm_flg[i]=1
##------------------------------------------------------------------------------

# Grades
gradeL=['1','2','3','4','5','6','7','8','B','C','D','E','F','G','H','L','M','U','W','X','Y','Z']    
    
# Initialize
dTM['Tenure Type']=np.array(['' for _ in range(dTM['TM'].size)],dtype=object)
dTM['V Logs Delivered m3']=np.zeros(dTM['TM'].size)
dTM['V Logs Delivered Abs m3']=np.zeros(dTM['TM'].size)
dTM['V Logs Delivered m3/ha']=np.zeros(dTM['TM'].size)
dTM['V Logs Delivered Abs m3/ha']=np.zeros(dTM['TM'].size)
dTM['V Logs Cruise m3/ha']=np.zeros(dTM['TM'].size)
#dTM['V Logs Delivered m3 By Grade']=np.zeros((dTM['TM'].size,len(gradeL)))
dTM['V Logs Waste m3/ha']=np.zeros(dTM['TM'].size)
for Grade in gradeL:
    dTM['V Logs Delivered Grade ' + Grade + ' m3/ha']=np.zeros(dTM['TM'].size)
    dTM['V Logs Delivered Grade ' + Grade + ' Abs m3/ha']=np.zeros(dTM['TM'].size)
    dTM['V Logs Delivered Grade ' + Grade + ' Abs m3']=np.zeros(dTM['TM'].size)
dTM['Stump Logs Delivered $']=np.zeros(dTM['TM'].size)
dTM['Stump Logs Delivered Abs $']=np.zeros(dTM['TM'].size)
dTM['V NonLog Delivered m3/ha']=np.zeros(dTM['TM'].size)
dTM['V NonLog Delivered Abs m3/ha']=np.zeros(dTM['TM'].size)
dTM['V NonLog Delivered Abs m3']=np.zeros(dTM['TM'].size)

# Populate from annual summary files
binY=np.arange(YrStart,2021,1)
for iY in range(binY.size):
    
    print(binY[iY])
    
    df=pd.read_csv(r'C:\Users\rhember\Documents\Data\Harvest\HBS\HBS ' + str(binY[iY]) + '.csv',header=None)
    
    dHB={}
    for i in range(len(df.columns)):
        dHB[i]=df.iloc[:,i].to_numpy()
    
    # Remove unnecessary data to speed this up
    #ind=np.where( (dHB[6]=='Logs') )[0]
    #for k in dHB.keys():
    #    dHB[k]=dHB[k][ind]
    
    for iTM in range(uTM.size):
        #print(iTM)
        
        ind=np.where( (dHB[0]==uTM[iTM]) )[0]
        
        if ind.size==0:
            continue
        
        Mat=dHB[6][ind]
        Type=dHB[9][ind]
        V=dHB[11][ind]
        Stump=dHB[12][ind]
        Vabs=np.abs(dHB[11][ind])
        Stump_abs=np.abs(dHB[12][ind])
        Grade=dHB[7][ind]
    
        # logs delivered
        ind2=np.where( (Mat=='Logs') )[0]
        if ind2.size>0:
            dTM['Tenure Type'][iTM]=dTM['Tenure Type'][iTM]+dHB[19][ind[ind2][0]]
            dTM['V Logs Delivered m3'][iTM]=np.round(dTM['V Logs Delivered m3'][iTM]+np.sum(V[ind2]),decimals=0)
            dTM['V Logs Delivered Abs m3'][iTM]=np.round(dTM['V Logs Delivered Abs m3'][iTM]+np.sum(Vabs[ind2]),decimals=0)
            dTM['V Logs Delivered m3/ha'][iTM]=np.round(dTM['V Logs Delivered m3/ha'][iTM]+np.sum(V[ind2]),decimals=0)
            dTM['V Logs Delivered Abs m3/ha'][iTM]=np.round(dTM['V Logs Delivered Abs m3/ha'][iTM]+np.sum(Vabs[ind2]),decimals=0)            
            dTM['Stump Logs Delivered $'][iTM]=np.round(dTM['Stump Logs Delivered $'][iTM]+np.sum(Stump[ind2]),decimals=0)
            dTM['Stump Logs Delivered Abs $'][iTM]=np.round(dTM['Stump Logs Delivered Abs $'][iTM]+np.sum(Stump_abs[ind2]),decimals=0)
                
        # Logs delivered by grade
        for iGrade in range(len(gradeL)):
            ind2=np.where( (Mat=='Logs') & (Grade==gradeL[iGrade])  )[0]
            if ind2.size>0:
                #dTM['V Logs Delivered m3 By Grade'][iTM,iGrade]=np.round(dTM['V Logs Delivered m3 By Grade'][iTM,iGrade]+np.sum(V[ind2]),decimals=0)        
                dTM['V Logs Delivered Grade ' + gradeL[iGrade] + ' m3/ha'][iTM]=np.round(dTM['V Logs Delivered Grade ' + gradeL[iGrade] + ' m3/ha'][iTM]+np.sum(V[ind2]))
                dTM['V Logs Delivered Grade ' + gradeL[iGrade] + ' Abs m3/ha'][iTM]=np.round(dTM['V Logs Delivered Grade ' + gradeL[iGrade] + ' Abs m3/ha'][iTM]+np.sum(Vabs[ind2]))
                dTM['V Logs Delivered Grade ' + gradeL[iGrade] + ' Abs m3'][iTM]=np.round(dTM['V Logs Delivered Grade ' + gradeL[iGrade] + ' Abs m3'][iTM]+np.sum(Vabs[ind2]))
        
        # Logs Waste
        ind2=np.where( (Type=='WA') & (Mat=='Logs') | (Type=='WU') & (Mat=='Logs') )[0]
        if ind2.size>0:
            dTM['V Logs Waste m3/ha'][iTM]=np.round(dTM['V Logs Waste m3/ha'][iTM]+np.sum(V[ind2]),decimals=0)

        # Non-logs delivered
        ind2=np.where( (Mat!='Logs') )[0]
        if ind2.size>0:
            dTM['V NonLog Delivered m3/ha'][iTM]=np.round(dTM['V NonLog Delivered m3/ha'][iTM]+np.sum(V[ind2]),decimals=0)
            dTM['V NonLog Delivered Abs m3/ha'][iTM]=np.round(dTM['V NonLog Delivered Abs m3/ha'][iTM]+np.sum(Vabs[ind2]),decimals=0)
            dTM['V NonLog Delivered Abs m3'][iTM]=np.round(dTM['V NonLog Delivered Abs m3'][iTM]+np.sum(Vabs[ind2]),decimals=0)
        
dTM['V Logs Delivered m3/ha']=dTM['V Logs Delivered m3/ha']/dTM['DISTURBANCE_GROSS_AREA']
dTM['V Logs Delivered Abs m3/ha']=dTM['V Logs Delivered Abs m3/ha']/dTM['DISTURBANCE_GROSS_AREA']
dTM['V Logs Waste m3/ha']=dTM['V Logs Waste m3/ha']/dTM['DISTURBANCE_GROSS_AREA']
dTM['V Logs Cruise m3/ha']=dTM['V Logs Cruise m3/ha']/dTM['DISTURBANCE_GROSS_AREA']
dTM['V NonLog Delivered m3/ha']=dTM['V NonLog Delivered m3/ha']/dTM['DISTURBANCE_GROSS_AREA']
dTM['V NonLog Delivered Abs m3/ha']=dTM['V NonLog Delivered Abs m3/ha']/dTM['DISTURBANCE_GROSS_AREA']

for Grade in gradeL:
    dTM['V Logs Delivered Grade ' + Grade + ' m3/ha']=dTM['V Logs Delivered Grade ' + Grade + ' m3/ha']/dTM['DISTURBANCE_GROSS_AREA']
    dTM['V Logs Delivered Grade ' + Grade + ' Abs m3/ha']=dTM['V Logs Delivered Grade ' + Grade + ' Abs m3/ha']/dTM['DISTURBANCE_GROSS_AREA']

#%% Import Waste summaries
# Don't download more than 3 years at a time
# Don't download the regions - tempting but too big

pthinW=r'C:\Users\rhember\Documents\Data\Waste Wood\FromWasteSystem'
fL=listdir(pthinW)

N=2000000
dW={}
dW['TM']=np.array(['' for _ in range(N)],dtype=object)
dW['Net Area Ha']=np.zeros(N)
dW['Total A. Sawlog  Volume m3']=np.zeros(N)
dW['A. Grade Y/4 Volume m3']=np.zeros(N)
dW['Av. + Unav. All Grades Volume m3']=np.zeros(N)
dW['Disp_A']=np.zeros(N)
dW['Disp_V']=np.zeros(N)
dW['Accu_A']=np.zeros(N)
dW['Accu_V']=np.zeros(N)
dW['Stan_A']=np.zeros(N)
dW['Stan_V']=np.zeros(N)
dW['Waste Billing']=np.zeros(N)

cnt=0

for iF in range(len(fL)):
    
    df=pd.read_excel(pthinW + '\\' + fL[iF],skiprows=list(range(0,14)))
        
    dW['TM'][cnt:cnt+len(df)]=df['TM'].to_numpy()
    dW['Net Area Ha'][cnt:cnt+len(df)]=df['Net Area Ha'].to_numpy()
    dW['Total A. Sawlog  Volume m3'][cnt:cnt+len(df)]=df['Total A. Sawlog  Volume m3'].to_numpy()
    dW['A. Grade Y/4 Volume m3'][cnt:cnt+len(df)]=df['A. Grade Y/4 Volume m3'].to_numpy()
    dW['Av. + Unav. All Grades Volume m3'][cnt:cnt+len(df)]=df['Av. + Unav. All Grades Volume m3'].to_numpy()
    dW['Disp_A'][cnt:cnt+len(df)]=df['Area\nHa'].to_numpy()
    dW['Disp_V'][cnt:cnt+len(df)]=df[' Volume m3'].to_numpy()
    dW['Accu_A'][cnt:cnt+len(df)]=df['Area\nHa.1'].to_numpy()
    dW['Accu_V'][cnt:cnt+len(df)]=df[' Volume m3.1'].to_numpy()
    dW['Stan_A'][cnt:cnt+len(df)]=df['Area\nHa.2'].to_numpy()
    dW['Stan_V'][cnt:cnt+len(df)]=df[' Volume m3.2'].to_numpy()
    dW['Waste Billing'][cnt:cnt+len(df)]=df['Waste $Billing'].to_numpy()    
    cnt=cnt+len(df)

# Truncate
for k in dW.keys():
    dW[k]=dW[k][0:cnt]

#%% Add waste data to TM database

# Initialize
dTM['Waste N Entries']=np.zeros(dTM['TM'].size)

dTM['Waste Net Area Tot']=np.zeros(dTM['TM'].size)
dTM['Waste Sawlog m3/ha']=np.zeros(dTM['TM'].size)
dTM['Waste GradeY/4 m3/ha']=np.zeros(dTM['TM'].size)
dTM['Waste Total m3/ha']=np.zeros(dTM['TM'].size)

dTM['Waste Dispersed m3/ha']=np.zeros(dTM['TM'].size)
dTM['Waste Accumulation m3/ha']=np.zeros(dTM['TM'].size)
dTM['Waste Standing m3/ha']=np.zeros(dTM['TM'].size)
dTM['Waste Bill $']=np.zeros(dTM['TM'].size)

for iTM in range(uTM.size):
    
    ind=np.where( (dW['TM']==uTM[iTM]) )[0]
    
    dTM['Waste N Entries'][iTM]=ind.size
    if ind.size>0:
        
        A_tot=np.sum(dW['Net Area Ha'][ind])
        
        dTM['Waste Net Area Tot'][iTM]=A_tot
        
        dTM['Waste Sawlog m3/ha'][iTM]=np.round(np.sum(dW['Total A. Sawlog  Volume m3'][ind])/A_tot,decimals=0)
        dTM['Waste GradeY/4 m3/ha'][iTM]=np.round(np.sum(dW['A. Grade Y/4 Volume m3'][ind])/A_tot,decimals=0)
        dTM['Waste Total m3/ha'][iTM]=np.round(np.sum(dW['Av. + Unav. All Grades Volume m3'][ind])/A_tot,decimals=0)
        
        dTM['Waste Dispersed m3/ha'][iTM]=np.round(np.sum(dW['Disp_V'][ind])/A_tot,decimals=0)
        dTM['Waste Accumulation m3/ha'][iTM]=np.round(np.sum(dW['Accu_V'][ind])/A_tot,decimals=0)
        dTM['Waste Standing m3/ha'][iTM]=np.round(np.sum(dW['Stan_V'][ind])/A_tot,decimals=0)
        dTM['Waste Bill $'][iTM]=np.sum(dW['Waste Billing'][ind])

# Percent fate of waste
tot=dTM['Waste Dispersed m3/ha']+dTM['Waste Accumulation m3/ha']+dTM['Waste Standing m3/ha']
dTM['Waste Dispersed %']=np.round(dTM['Waste Dispersed m3/ha']/tot*100)
dTM['Waste Accumulation %']=np.round(dTM['Waste Accumulation m3/ha']/tot*100)
dTM['Waste Standing %']=np.round(dTM['Waste Standing m3/ha']/tot*100)

#dTM['Waste Sawlog m3/ha']=dTM['Waste Sawlog m3/ha']/dTM['Waste Net Area Tot']
#dTM['Waste GradeY/4 m3/ha']=dTM['Waste GradeY/4 m3/ha']/dTM['Waste Net Area Tot']
#dTM['Waste Total m3/ha']=dTM['Waste Total m3/ha']/dTM['Waste Net Area Tot']

# Calculate volume using waste net area
#dTM['V Delivered m3/ha (A from WS)']=np.round(dTM['V Delivered m3']/dTM['Waste Net Area Tot'],decimals=0)
#ind=np.where( (dTM['V Delivered m3/ha (A from WS)']<-10000) | (dTM['V Delivered m3/ha (A from WS)']>10000) )[0]
#dTM['V Delivered m3/ha (A from WS)'][ind]=0
 
#%% Keep a complete list of unique openings for crosswalk to VRI

# No need to keep earlier timber marks
ind=np.where(dTM['Year First']>=2006)[0]
uTM_Recent=np.unique(dTM['TM'][ind])

N=2000000
dOp={}
dOp['TM']=np.array(['' for _ in range(N)],dtype=object)
dOp['OPENING_ID']=np.zeros(N)
dOp['PLANNED_NET_BLOCK_AREA']=np.zeros(N)
cnt=0
for iTM in range(uTM_Recent.size):
    ind=np.where(dFT['TIMBER_MARK']==uTM_Recent[iTM])[0]
    dOp['TM'][cnt:cnt+ind.size]=uTM_Recent[iTM]
    dOp['OPENING_ID'][cnt:cnt+ind.size]=dFT['OPENING_ID'][ind]
    dOp['PLANNED_NET_BLOCK_AREA'][cnt:cnt+ind.size]=dFT['PLANNED_NET_BLOCK_AREA'][ind]
    cnt=cnt+ind.size

# Truncate
for k in dOp.keys():
    dOp[k]=dOp[k][0:cnt]
       
#%% Import RESULTS OP layer to get stand age

# List layers
pthin=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20210930\Results.gdb'
# fiona.listlayers(pthin)

dOp['PREV_AGE_CLASS_CODE']=np.zeros(dOp['OPENING_ID'].size)
dOp['DISTRICT_CODE']=np.array(['' for _ in range(dOp['OPENING_ID'].size)],dtype=object)

cnt=0
with fiona.open(pthin,layer='RSLT_OPENING_SVW') as source:
    for feat in source:
        p=feat['properties']
        #break
        
        ind=np.where(dOp['OPENING_ID']==p['OPENING_ID'])[0]
        if ind.size==0:
            continue
        dOp['PREV_AGE_CLASS_CODE'][ind]=p['PREV_AGE_CLASS_CODE']
        dOp['DISTRICT_CODE'][ind]=p['DISTRICT_CODE']        
        print(cnt)
        cnt=cnt+1

# Convert age class to age (years)
AgeClassAge=[10,30,50,70,90,110,130,195,250]  
dOp['Age RSLTS']=np.nan*np.ones(dOp['OPENING_ID'].size)
for i in range(len(AgeClassAge)):
    ind=np.where(dOp['PREV_AGE_CLASS_CODE']==np.array(i-1))[0]
    dOp['Age RSLTS'][ind]=AgeClassAge[i]

# Add to TM database
dTM['Age RSLTS WA']=np.zeros(dTM['TM'].size)
dTM['District']=np.array(['' for _ in range(dTM['TM'].size)],dtype=object)
for iTM in range(uTM.size):
    ind=np.where( (dOp['TM']==uTM[iTM]) )[0]
    if ind.size==0:
        continue
    dTM['Age RSLTS WA'][iTM]=np.round(np.sum(dOp['Age RSLTS'][ind]*dOp['PLANNED_NET_BLOCK_AREA'][ind])/np.sum(dOp['PLANNED_NET_BLOCK_AREA'][ind]),decimals=0)
    dTM['District'][iTM]=dOp['DISTRICT_CODE'][ind[0]]

#%% Save
    
gu.opickle(r'C:\Users\rhember\Documents\Data\Waste Wood\WasteSummary.pkl',dTM)

gu.opickle(r'C:\Users\rhember\Documents\Data\Waste Wood\WasteSummary_UniqueOpenings.pkl',dOp)
      
#%% Export to spreadsheet for review

#del dTM['V Logs Delivered m3 By Grade']
df=pd.DataFrame(dTM)
df.to_excel(r'C:\Users\rhember\Documents\Data\Waste Wood\WasteSummary.xlsx',index=False)

#%% Add insect and fire damage

dTM=gu.ipickle(r'C:\Users\rhember\Documents\Data\Waste Wood\WasteSummary.pkl')

geos=gu.ipickle(r'D:\Data\FCI_Projects\SummaryWaste\Geospatial\geos.pkl')
dOpSG=gu.ipickle(r'D:\Data\FCI_Projects\SummaryWaste\Geospatial\RSLT_OPENING_SVW.pkl')
dPestSG=gu.ipickle(r'D:\Data\FCI_Projects\SummaryWaste\Geospatial\PEST_INFESTATION_POLY.pkl')
dFireSG=gu.ipickle(r'D:\Data\FCI_Projects\SummaryWaste\Geospatial\PROT_HISTORICAL_FIRE_POLYS_SP.pkl')
dVRISG=gu.ipickle(r'D:\Data\FCI_Projects\SummaryWaste\Geospatial\VEG_COMP_LYR_R1_POLY.pkl')

meta={'Paths':{}}
meta['Paths']['Model Code']=r'C:\Users\rhember\Documents\Code_Python\fcgadgets\cbrunner'
lut=invu.Load_LUTs(meta)

#%% Get beetle damage 

# LUT for severity to mortality
Mort_LUT=np.array([0,0,20,40,70,5,90])
TimeWindow=12.0

dTM['D_Beetle']=np.zeros(dTM['TM'].size)
for iTM in range(uTM.size):
    
    if dTM['Year First'][iTM]<2001:
        continue
    
    try:
        id=lut['LUT']['OP']['TIMBER_MARK'][ dTM['TM'][iTM] ][0]
        #break
    except:
        continue
   
    ind1=np.where(dTM['TM']==uTM[iTM])[0]
    ind2=np.where(dOpSG['TIMBER_MARK']==id)[0]
    if ind2.size==0:
        continue
    
    #if ind2.size>0:
    #    break
    idx2sxy=dOpSG['IdxToSXY'][ind2]
    
    # Index to pest
    indP=np.nonzero(np.in1d(dPestSG['IdxToSXY'],idx2sxy))[0]
    id_spc=dPestSG['PEST_SPECIES_CODE'][indP]
    id_sev=dPestSG['PEST_SEVERITY_CODE'][indP]
    yr=dPestSG['CAPTURE_YEAR'][indP]
    
    # Conver severity to mortality
    #lut['LUT']['Pest']['PEST_SEVERITY_CODE']
    Mort=Mort_LUT[id_sev]
    
    # Is it beetles?
    flg_B=np.zeros(indP.size)
    ind=np.where( (id_spc>=32) & (id_spc<=42) & (dTM['Year First'][iTM]-yr<TimeWindow) )[0]
    flg_B[ind]=1
    
    # Sum the average per-hectare mortality due to beetles
    dTM['D_Beetle'][iTM]=np.sum(flg_B*Mort)/indP.size


#%%
    
# plt.hist(dTM['D_Beetle'])

gu.opickle(r'C:\Users\rhember\Documents\Data\Waste Wood\WasteSummary.pkl',dTM)


