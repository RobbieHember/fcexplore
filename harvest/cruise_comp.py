'''
TIMBER CRUISE COMPILIATION
'''
#%% Import modules
import numpy as np
import matplotlib.pyplot as plt
import fcgadgets.macgyver.util_general as gu

#%% Import data
gp=gu.SetGraphics('Manuscript')
dI=gu.ReadCSV(r'C:\Data\ECAS\Received 2023-04-04\rh_COMP_DATA_INTERIOR.csv')
dC=gu.ReadCSV(r'C:\Data\ECAS\Received 2023-04-04\rh_COMP_DATA_COAST.csv')
dAge=gu.ReadCSV(r'C:\Data\ECAS\Received 2023-04-04\rh_DAT_PLOT_AGES.csv')

#%% Combine coast and interior
d={}
for k in dI.keys():
    try:
        d[k]=np.append(dI[k],dC[k])
    except:
        d[k]=np.append(dI[k],np.nan*np.ones(dC['PCT_DIST'].size))

for k in dC.keys():
    if k in d.keys():
        continue
    try:
        d[k]=np.append(np.nan*np.ones(dI['PCT_DIST'].size),dC[k])
    except:
        pass

list(d.keys())

# QA:
#for k in d.keys():
#    print(d[k].size)

#%% Clean variables with mixed strings
vL=['PCT_DIST','PCT_DECAY','PCT_WASTE','PCT_WASTE_BILL','PCT_BREAK','PCT_DWB','PCT_NET_2NDGROWTH','PCT_INSECT_VOL','PCT_NO_INSECT_M3','PCT_GREEN_INSECT_M3', \
    'PCT_RED_INSECT_M3','PCT_GREY_INSECT_M3','PCT_OTHER_INSECT_M3','X_DEFOLIATOR_LIVE_CAMB_PCT','Y_DEFOLIATOR_DEAD_CAMB_PCT']
for v in vL:
    for i in range(d[v].size):
        try:
            d[v][i].astype('float')
        except:
            d[v][i]='0'
    d[v]=d[v].astype('float')
    ind=np.where( (d[v]<0) | (d[v]>100) )[0]
    d[v][ind]=np.nan

#%% Add basic derived variables (more below)

d['V Gross (m3/ha)']=d['GROSS_M3']/d['NET_AREA']
d['V Net (m3/ha)']=d['NET_M3']/d['NET_AREA']
d['Pct Dead Net']=d['NET_DEAD_M3']/d['NET_M3']*100
d['PCT Net']=d['V Net (m3/ha)']/d['V Gross (m3/ha)']*100

# Percent conifer
u=np.unique(d['PRIMARY_MARK'])
d['PCT CON']=np.nan*np.ones(d['NET_AREA'].size)
for i in range(u.size):
    ind1=np.where( (dC['PRIMARY_MARK']==u[i]) & (dC['SPECIES']=='ALL') )[0]
    ind2=np.where( (dC['PRIMARY_MARK']==u[i]) & (dC['SPECIES']=='CON') )[0]
    d['PCT CON'][ind1]=np.nanmean(d['V Gross (m3/ha)'][ind2])/np.nanmean(d['V Gross (m3/ha)'][ind1])*100
d['PCT CON']=gu.Clamp(d['PCT CON'],0,100)

#%% Collapse down to ALL species and exclude reductions

ind=np.where( (d['SPECIES']=='ALL') & (d['REDN']!='Y') )[0]
for k in d.keys():
    d[k]=d[k][ind]

#%% Some TMs have duplication due to differences in PARENT_ID -> remove

flg=np.ones(d['SPECIES'].size)
uTM=np.unique(d['PRIMARY_MARK'])
for iTM in range(uTM.size):
    ind=np.where(d['PRIMARY_MARK']==uTM[iTM])[0]
    if ind.size>1:
        flg[ind[0:-1]]=0

ind=np.where(flg==1)[0]
for k in d.keys():
    d[k]=d[k][ind]

#%% Derived variables

d['Age']=np.nan*np.ones(d['NET_AREA'].size)
for i in range(d['NET_AREA'].size):
    ind=np.where( (dAge['ECAS_ID']==d['ECAS_ID'][i]) & (dAge['AGE_TENS']>0) )[0]
    d['Age'][i]=np.mean(dAge['AGE_TENS'][ind]*10)

d['Year']=np.nan*np.ones(d['NET_AREA'].size)
for i in range(d['NET_AREA'].size):
    try:
        d['Year'][i]=d['COMP_DATE'][i][-4:]
    except:
        d['Year'][i]=2000+np.array(d['COMP_DATE'][i][-2::],dtype='float')
        #print(d['Year'][i]) # Make sure it's working

#%% Save cleaned version

gu.opickle(r'C:\Users\rhember\Documents\Data\ECAS\Received 2023-04-04\rh_COMP_DATA_CleanCombined.pkl',d)

# d=gu.ipickle(r'C:\Users\rhember\Documents\Data\ECAS\Received 2023-04-04\rh_COMP_DATA_CleanCombined.pkl')

y=np.zeros(tv.size)
for iT in range(tv.size):
    ind=np.where(d['Year']==tv[iT])[0]
    y[iT]=ind.size
plt.plot(tv,y,'ko')

#%% Explore

flg=0
if flg==1:
    # TIMBER MARK LIST APPEARS TO BE THE SAME
    u=np.unique(d['PRIMARY_MARK'])
    print(u.size)

    ind=np.where( (d['PRIMARY_MARK']==u[10]) )[0]

    list(d['TIMBER_MARK_LIST'][ind])
    d['MARK_COUNT'][ind]
    d['LOCATION'][ind]
    d['CP'][ind] # CPs don't change within TM
    list(d['SPECIES'][ind])
    d['Age'][ind]
    d['PCT_BLOWDOWN'][ind]
    d['PCT_DIST'][ind]
    d['AVG_SNAG_DBH'][ind]

    ind=np.where( (d['SPECIES']=='ALL') )[0]

    plt.close('all')
    plt.hist(d['V Net (m3/ha)'][ind],np.arange(0,1000,10))

#%% Histograms

plt.close('all')
fig,ax=plt.subplots(2,3,figsize=gu.cm2inch(15,8))
ax[0,0].hist(d['V Gross (m3/ha)'],np.arange(0,1220,20))
ax[0,0].set(xlabel='Gross volume (m3/ha)',ylabel='Frequency')
ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=2)
ax[0,1].hist(d['V Net (m3/ha)']/d['V Gross (m3/ha)']*100,np.arange(0,105,5))
ax[0,1].set(xlabel='Net volume fraction (%)',ylabel='Frequency')
ax[0,2].hist(d['Pct Dead Net'],np.arange(0,105,5))
ax[0,2].set(xlabel='Dead volume fraction (%)',ylabel='Frequency')
ax[1,0].hist(d['PCT CON'],np.arange(0,105,5))
ax[1,0].set(xlabel='Conifer fraction (%)',ylabel='Frequency')
ax[1,1].hist(d['PCT_NET_2NDGROWTH'],np.arange(0,105,5)) #
ax[1,1].set(xlabel='Second growth fraction (%)',ylabel='Frequency')
ax[1,2].hist(d['PCT_NET_IMM'],np.arange(0,105,5))
ax[1,2].set(xlabel='Immature fraction (%)',ylabel='Frequency')
plt.tight_layout()
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Harvest\Timber Cruise Info\Histograms','png',900)

#%% Time series

tv=np.arange(2000,2022,1)
dt={'Mean':{},'Sum':{}}
for k in d.keys():
    dt['Mean'][k]=np.zeros(tv.size)
    dt['Sum'][k]=np.zeros(tv.size)
    for iT in range(tv.size):
        ind=np.where( (d['Year']==tv[iT]) )[0]
        try:
            dt['Mean'][k][iT]=np.nanmean(d[k][ind])
            dt['Sum'][k][iT]=np.nansum(d[k][ind])
        except:
            pass

iT=np.where(tv>=2014)[0]

plt.close('all')
fig,ax=plt.subplots(2,3,figsize=gu.cm2inch(15,8))
ax[0,0].plot(tv[iT],dt['Mean']['V Gross (m3/ha)'][iT],'-bo')
ax[0,0].set(xlabel='Time, years',ylabel='Gross volume (m3/ha)')
ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=2)
ax[0,1].plot(tv[iT],(dt['Mean']['V Net (m3/ha)'][iT]/dt['Mean']['V Gross (m3/ha)'][iT])*100,'-bo')
ax[0,1].set(xlabel='Time, years',ylabel='Net volume fraction (%)')
ax[0,2].plot(tv[iT],dt['Mean']['Pct Dead Net'][iT],'-bo')
ax[0,2].set(xlabel='Time, years',ylabel='Dead volume fraction (%)')
ax[1,0].plot(tv[iT],dt['Sum']['NET_AREA'][iT]/1e3,'-bo')
ax[1,0].set(xlabel='Time, years',ylabel='Net area (Kha)')
ax[1,1].plot(tv[iT],dt['Mean']['PCT_NET_2NDGROWTH'][iT],'-bo') #
ax[1,1].set(xlabel='Time, years',ylabel='Second growth fraction (%)')
ax[1,2].plot(tv[iT],dt['Mean']['PCT_NET_IMM'][iT],'-bo')
ax[1,2].set(xlabel='Time, years',ylabel='Immature fraction (%)')
plt.tight_layout()
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Harvest\Timber Cruise Info\TimeSeries','png',900)

#%%

plt.close('all')
plt.plot(tv,dt['Age'],'-bo')

#%%

plt.close('all')
plt.plot(tv,dt['Pct Dead Net'],'-bo')

#%%

plt.close('all')
plt.plot(tv,dt['PCT_DWB'],'-bo')

#%%

#ind=np.where( (dI['SPECIES']=='ALL') )[0]

u=np.unique(d['Age I'])
yG=np.zeros(u.size)
yN=np.zeros(u.size)
for i in range(u.size):
    ind=np.where( (d['Age I']==u[i]) & (d['SPECIES']=='ALL') )[0]
    yG[i]=np.nanmean(d['V Gross (m3/ha)'][ind])
    yN[i]=np.nanmean(d['V Net (m3/ha)'][ind])


plt.close('all')
plt.plot(u,yG,'bo')
plt.plot(u,yN,'gs')

plt.close('all')
plt.plot(u,yN/yG,'gs')

