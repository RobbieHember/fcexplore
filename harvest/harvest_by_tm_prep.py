
'''
SUMMARIZE HARVEST VOLUME AND WASTE FROM HBS AND WASTE SYSTEM
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

YrStart=2007

# Grades
gradeL=['1','2','3','4','5','6','7','8','B','C','D','E','F','G','H','L','M','U','W','X','Y','Z']

# Graphics parameters
gp=gu.SetGraphics('Manuscript')

# Paths
Paths={}
Paths['Dist']=r'C:\Users\rhember\Documents\Data\ForestInventory\Disturbances\20230501\Disturbances.gdb'
Paths['Results']=r'C:\Users\rhember\Documents\Data\ForestInventory\Results\20230430\Results.gdb'
Paths['HBS']=r'C:\Users\rhember\Documents\Data\Harvest\HBS'
Paths['Waste']=r'C:\Users\rhember\Documents\Data\Waste Wood\FromWasteSystem'
# List layers # fiona.listlayers(Paths['Dist'])

#%% Import forest tenure layer

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
with fiona.open(Paths['Dist'],layer='FTEN_CUT_BLOCK_OPEN_ADMIN') as source:
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

#%% Create database by unique timber marks

# Indices to unique timber marks
iTM=gu.IndicesFromUniqueArrayValues(dFT['TIMBER_MARK'])
n=len(iTM)

dTM={}
dTM['TM']=np.array(['' for _ in range(n)],dtype=object)
dTM['N Openings']=np.zeros(n)
dTM['DISTURBANCE_GROSS_AREA']=np.zeros(n)
dTM['PLANNED_NET_BLOCK_AREA']=np.zeros(n)
dTM['Year First']=np.zeros(n)
dTM['Year Last']=np.zeros(n)
dTM['District']=np.array(['' for _ in range(n)],dtype=object)
#dTM['ADMIN_DISTRICT_CODE']=np.array(['' for _ in range(uTM.size)],dtype=object)
#dTM['FEATURE_AREA']=np.zeros(uTM.size)
cnt=0
for k in iTM.keys():
    #ind=np.where(dFT['TIMBER_MARK']==uTM[iTM])[0]
    dTM['TM'][cnt]=k
    dTM['N Openings'][cnt]=np.unique(dFT['OPENING_ID'][iTM[k]]).size
    dTM['DISTURBANCE_GROSS_AREA'][cnt]=np.round(np.nansum(dFT['DISTURBANCE_GROSS_AREA'][iTM[k]]))
    dTM['PLANNED_NET_BLOCK_AREA'][cnt]=np.round(np.sum(dFT['PLANNED_NET_BLOCK_AREA'][iTM[k]]))
    dTM['Year First'][cnt]=np.min(dFT['DISTURBANCE_START_YEAR'][iTM[k]])
    dTM['Year Last'][cnt]=np.max(dFT['DISTURBANCE_END_YEAR'][iTM[k]])
    #dTM['CLIENT_NUMBER'][iTM]=dFT['CLIENT_NUMBER'][ind]
    #dTM['ADMIN_DISTRICT_CODE'][iTM]=dFT['ADMIN_DISTRICT_CODE'][ind[0]]
    #dTM['FEATURE_AREA'][iTM]=np.round(np.sum(dFT['FEATURE_AREA'][ind]))
    cnt=cnt+1

# Unique timber marks
uTM=np.unique(dTM['TM'])

# Plot
flg=0
if flg==1:
    tv=np.arange(1990,2023,1)
    y=np.zeros(tv.size)
    for iT in range(tv.size):
        ind=np.where(dFT['Year First']==tv[iT])[0]
        y[iT]=np.sum(dFT['FEATURE_AREA'][ind])
    plt.plot(tv,y,'-bo')

#%% Import HBS
# The system appears to be down often
# -> Date of Invoice, "Billing processed" only
# Can only do 12 months at a time

# Definition of volume types (from Ross Pavan)
#   1.	Normal:  the volume delivered to the sawmill (the volume calculation is based on weight to volume ratios)
#   2.	Cruise Based:  the volume is derived from the cruise compilation (the billed volume never exceeds the cruise volume for a cutting permit, for a timber sale license the billed volume will never exceed the cruise volume plus any additional volume added after the fact (example: external R/W volume)
#   3.	Waste: the volume is derived from the waste assessment of a scale based cutting authority
#   4.	Beachcomb: volume scaled because of the beachcombing efforts of an individual or salvage logger

# Initialize
dTM['Tenure Type']=np.array(['' for _ in range(dTM['TM'].size)],dtype=object)
dTM['Tenure Holder']=np.array(['' for _ in range(dTM['TM'].size)],dtype=object)
dTM['Year Min']=np.zeros(dTM['TM'].size)
dTM['Year Max']=np.zeros(dTM['TM'].size)
dTM['V Logs Delivered m3']=np.zeros(dTM['TM'].size)
dTM['V Logs Delivered Abs m3']=np.zeros(dTM['TM'].size)
dTM['V Logs Delivered m3/ha']=np.zeros(dTM['TM'].size)
dTM['V Logs Delivered Abs m3/ha']=np.zeros(dTM['TM'].size)
dTM['V Logs Cruise m3/ha']=np.zeros(dTM['TM'].size)
for Grade in gradeL:
    dTM['V Logs Delivered Grade ' + Grade + ' m3/ha']=np.zeros(dTM['TM'].size)
    dTM['V Logs Delivered Grade ' + Grade + ' Abs m3/ha']=np.zeros(dTM['TM'].size)
    dTM['V Logs Delivered Grade ' + Grade + ' Abs m3']=np.zeros(dTM['TM'].size)
dTM['V Logs Waste m3/ha']=np.zeros(dTM['TM'].size)
dTM['Stump Logs Delivered $']=np.zeros(dTM['TM'].size)
dTM['Stump Logs Delivered Abs $']=np.zeros(dTM['TM'].size)
dTM['V NonLog Delivered m3/ha']=np.zeros(dTM['TM'].size)
dTM['V NonLog Delivered Abs m3/ha']=np.zeros(dTM['TM'].size)
dTM['V NonLog Delivered Abs m3']=np.zeros(dTM['TM'].size)
dTM['V NonLog Waste m3/ha']=np.zeros(dTM['TM'].size)

# Populate from annual summary files
binY=np.arange(YrStart,2023,1)
for iY in range(binY.size):

    print(binY[iY])

    df=pd.read_csv(Paths['HBS'] + '\\HBS ' + str(binY[iY]) + '.csv',header=None)

    dHB={}
    for i in range(len(df.columns)):
        dHB[i]=df.iloc[:,i].to_numpy()

    # Create indices for each TM
    dHB[0]=dHB[0].astype('U')
    tm=gu.IndicesFromUniqueArrayValues(dHB[0])

    for iTM in range(uTM.size):

        if uTM[iTM] not in tm.keys():
            # TM not active this year, continue
            continue

        ind=tm[uTM[iTM]]

        Mat=dHB[6][ind]
        Type=dHB[9][ind]
        V=dHB[11][ind]
        V_abs=np.abs(dHB[11][ind])
        Stump=dHB[12][ind]
        Stump_abs=np.abs(dHB[12][ind])
        Grade=dHB[7][ind]

        # logs delivered
        ind2=np.where( (Mat=='Logs') )[0]
        if ind2.size>0:

            if dTM['Year Min'][iTM]==0:
                dTM['Year Min'][iTM]=binY[iY]
            elif binY[iY]<=dTM['Year Min'][iTM]:
                dTM['Year Min'][iTM]=binY[iY]

            if dTM['Year Max'][iTM]==0:
                dTM['Year Max'][iTM]=binY[iY]
            elif binY[iY]>=dTM['Year Max'][iTM]:
                dTM['Year Max'][iTM]=binY[iY]

            dTM['Tenure Type'][iTM]=dHB[19][ind[ind2][0]]
            dTM['Tenure Holder'][iTM]=dHB[26][ind[ind2][0]]
            dTM['V Logs Delivered m3'][iTM]=np.round(dTM['V Logs Delivered m3'][iTM]+np.sum(V[ind2]),decimals=0)
            dTM['V Logs Delivered Abs m3'][iTM]=np.round(dTM['V Logs Delivered Abs m3'][iTM]+np.sum(V_abs[ind2]),decimals=0)
            dTM['V Logs Delivered m3/ha'][iTM]=np.round(dTM['V Logs Delivered m3/ha'][iTM]+np.sum(V[ind2]),decimals=0)
            dTM['V Logs Delivered Abs m3/ha'][iTM]=np.round(dTM['V Logs Delivered Abs m3/ha'][iTM]+np.sum(V_abs[ind2]),decimals=0)
            dTM['Stump Logs Delivered $'][iTM]=np.round(dTM['Stump Logs Delivered $'][iTM]+np.sum(Stump[ind2]),decimals=0)
            dTM['Stump Logs Delivered Abs $'][iTM]=np.round(dTM['Stump Logs Delivered Abs $'][iTM]+np.sum(Stump_abs[ind2]),decimals=0)

        # Logs delivered by grade
        for iGrade in range(len(gradeL)):
            ind2=np.where( (Mat=='Logs') & (Grade==gradeL[iGrade])  )[0]
            if ind2.size>0:
                #dTM['V Logs Delivered m3 By Grade'][iTM,iGrade]=np.round(dTM['V Logs Delivered m3 By Grade'][iTM,iGrade]+np.sum(V[ind2]),decimals=0)
                dTM['V Logs Delivered Grade ' + gradeL[iGrade] + ' m3/ha'][iTM]=np.round(dTM['V Logs Delivered Grade ' + gradeL[iGrade] + ' m3/ha'][iTM]+np.sum(V[ind2]))
                dTM['V Logs Delivered Grade ' + gradeL[iGrade] + ' Abs m3/ha'][iTM]=np.round(dTM['V Logs Delivered Grade ' + gradeL[iGrade] + ' Abs m3/ha'][iTM]+np.sum(V_abs[ind2]))
                dTM['V Logs Delivered Grade ' + gradeL[iGrade] + ' Abs m3'][iTM]=np.round(dTM['V Logs Delivered Grade ' + gradeL[iGrade] + ' Abs m3'][iTM]+np.sum(V_abs[ind2]))

        # Logs Waste
        ind2=np.where( (Type=='WA') & (Mat=='Logs') | (Type=='WU') & (Mat=='Logs') )[0]
        if ind2.size>0:
            dTM['V Logs Waste m3/ha'][iTM]=np.round(dTM['V Logs Waste m3/ha'][iTM]+np.sum(V[ind2]),decimals=0)

        # Non-logs delivered
        ind2=np.where( (Mat!='Logs') )[0]
        if ind2.size>0:
            dTM['V NonLog Delivered m3/ha'][iTM]=np.round(dTM['V NonLog Delivered m3/ha'][iTM]+np.sum(V[ind2]),decimals=0)
            dTM['V NonLog Delivered Abs m3/ha'][iTM]=np.round(dTM['V NonLog Delivered Abs m3/ha'][iTM]+np.sum(V_abs[ind2]),decimals=0)
            dTM['V NonLog Delivered Abs m3'][iTM]=np.round(dTM['V NonLog Delivered Abs m3'][iTM]+np.sum(V_abs[ind2]),decimals=0)

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
# Just inlcude status = "Billing Issued"
# Don't download more than 3 years at a time
# Don't download the regions - tempting but too big

fL=listdir(Paths['Waste'])

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

    df=pd.read_excel(Paths['Waste'] + '\\' + fL[iF],skiprows=list(range(0,14)))
    #a=gu.ReadExcel(Paths['Waste'] + '\\' + fL[iF],skiprows=list(range(0,14)))

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

# for iTM in range(uTM.size):
#     ind=np.where( (dW['TM']==uTM[iTM]) )[0]

# Create indices for each TM
dW['TM']=dW['TM'].astype('U')
tm=gu.IndicesFromUniqueArrayValues(dW['TM'])

for iTM in range(uTM.size):

    if uTM[iTM] not in tm.keys():
        # TM not active this year, continue
        continue

    ind=tm[uTM[iTM]]

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

#%% Add timber cruise data

dTC=gu.ipickle(r'C:\Users\rhember\Documents\Data\ECAS\Received 2023-04-04\rh_COMP_DATA_CleanCombined.pkl')

vL=['NET_AREA','Year','Age','V Gross (m3/ha)','V Net (m3/ha)','Pct Dead Net','PCT Net','PCT_DIST','PCT_DECAY','PCT_WASTE','PCT_WASTE_BILL','PCT_BREAK','PCT_DWB','PCT_NET_2NDGROWTH','PCT_INSECT_VOL','PCT_NO_INSECT_M3','PCT_GREEN_INSECT_M3', \
    'PCT_RED_INSECT_M3','PCT_GREY_INSECT_M3','PCT_OTHER_INSECT_M3','X_DEFOLIATOR_LIVE_CAMB_PCT','Y_DEFOLIATOR_DEAD_CAMB_PCT']
for v in vL:
    dTM['Cruise_' + v]=np.zeros(dTM['TM'].size)

uTM=np.unique(dTC['PRIMARY_MARK'])
for iTM in range(uTM.size):
    ind1=np.where(dTC['PRIMARY_MARK']==uTM[iTM])[0]
    ind2=np.where(dTM['TM']==uTM[iTM])[0]
    for v in vL:
        dTM['Cruise_' + v][ind2]=dTC[v][ind1]

#%% Save

gu.opickle(r'C:\Users\rhember\Documents\Data\Harvest\HBS\HavestSummary_ByTM.pkl',dTM)

#gu.opickle(r'C:\Users\rhember\Documents\Data\Waste Wood\WasteSummary_UniqueOpenings.pkl',dOp)

#%%

list(dTM.keys())

ind=np.where( (dTM['Cruise_V Net (m3/ha)']>0) )[0]

plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(8,8))
ax.plot(dTM['Cruise_V Net (m3/ha)'][ind],dTM['V Logs Delivered Abs m3/ha'][ind],'b.',mew=0.25,ms=3,mfc='b',mec='w')
ax.plot([0,2000],[0,2000],'k-',lw=2,color=[0.8,0.8,0.8])
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=2)
ax.set(xlabel='Net vol from cruise (m3/ha)',ylabel='Net volume delivered (m3/ha)',xlim=[0,1500],ylim=[0,1500])
plt.tight_layout()
#gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Harvest\Timber Cruise Info\TimeSeries','png',900)

#%% Export to spreadsheet for review

#del dTM['V Logs Delivered m3 By Grade']
df=pd.DataFrame(dTM)
df.to_excel(r'C:\Users\rhember\Documents\Data\Waste Wood\WasteSummary.xlsx',index=False)

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
# fiona.listlayers(Paths['Results'])

dOp['PREV_AGE_CLASS_CODE']=np.zeros(dOp['OPENING_ID'].size)
dOp['DISTRICT_CODE']=np.array(['' for _ in range(dOp['OPENING_ID'].size)],dtype=object)

cnt=0
with fiona.open(Paths['Results'],layer='RSLT_OPENING_SVW') as source:
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



