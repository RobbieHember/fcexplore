
#%% Import data

dS=gu.ReadExcel(r'D:\Data\FCI_Projects\SummaryReforestation\Inputs\SummaryAttributes_BySXY.xlsx')


#%% Add planted species and planting density

uSXY=np.unique(dS['ID_SXY'])
dS['Spc']=np.array(['' for _ in range(dS['PL_S1_CD'].size)],dtype=object)
dS['PD']=np.zeros(dS['ID_SXY'].size)
for iSXY in range(uSXY.size):    
    ind=np.where( dS['ID_SXY']==uSXY[iSXY] )[0]
    iSpc=np.where( (dS['ID_SXY']==uSXY[iSXY]) & (dS['PL_S1_CD']!=' ') )[0]
    if iSpc.size>0:
        dS['Spc'][ind]=str(dS['PL_S1_CD'][iSpc[0]])
    iPD=np.where( (dS['ID_SXY']==uSXY[iSXY]) & (dS['PL SPH']>300) )[0]
    if iPD.size>0:
        dS['PD'][ind]=np.mean(dS['PL SPH'][iPD])

dS['Surv FG']=dS['FC SPH FG']/dS['PD']
dS['Surv TWS']=dS['FC SPH TWS']/dS['PD']

#%% Define sample

iBGC=np.where(dS['BGCz']=='SBS')[0]
dS2=dS.copy()
for k in dS.keys():
    dS2[k]=dS[k][iBGC]

ind=np.where(dS['Spc']=='FDI')[0]
dS2=dS.copy()
for k in dS.keys():
    dS2[k]=dS[k][ind]

#%% Get time profiles
    
rs={}
rs['Time']=np.arange(1,51,1)
nS=int(20000)
rs['FC SPH TOT']=np.nan*np.ones((rs['Time'].size,nS))
rs['FC SPH TWS']=np.nan*np.ones((rs['Time'].size,nS))
rs['FC SPH FG']=np.nan*np.ones((rs['Time'].size,nS))
rs['FC CC']=np.nan*np.ones((rs['Time'].size,nS))
rs['Surv TWS']=np.nan*np.ones((rs['Time'].size,nS))
rs['Surv FG']=np.nan*np.ones((rs['Time'].size,nS))

uSXY=np.unique(dS2['ID_SXY'])
cnt=0
for iSXY in range(uSXY.size):    
    #ind=np.where(dS2['ID_SXY']==uSXY[iSXY])[0]    
    iH=np.where( (dS2['ID_SXY']==uSXY[iSXY]) & (dS2['Dist Type']=='Harvest') )[0]    
    if (iH.size==0):
        continue
    
    iFC=np.where( (dS2['ID_SXY']==uSXY[iSXY]) & (dS2['FC Status']!='') )[0]
    #iFC=np.where( (dS2['ID_SXY']==uSXY[iSXY]) & (dS2['FC SPH FG']>0) )[0]
    if (iFC.size==0):
        continue
    
    print(iSXY)
    for j in range(iFC.size):
        dt=np.floor(dS2['Year'][iFC[j]]-dS2['Year'][iH[0]])
        it=np.where(rs['Time']==dt)[0]
        rs['FC CC'][it,cnt]=dS2['FC CC'][iFC[j]]
        rs['FC SPH TOT'][it,cnt]=dS2['FC SPH TOT'][iFC[j]]
        rs['FC SPH TWS'][it,cnt]=dS2['FC SPH TWS'][iFC[j]]
        rs['FC SPH FG'][it,cnt]=dS2['FC SPH FG'][iFC[j]]  
        rs['Surv TWS'][it,cnt]=dS2['Surv TWS'][iFC[j]]  
        rs['Surv FG'][it,cnt]=dS2['Surv FG'][iFC[j]]  
    cnt=cnt+1
    
for k in rs.keys():
    if k=='Time':
        continue
    ind=np.where(rs[k]<0)
    rs[k][ind]=np.nan
N_Sites=cnt    

#%% Plot

plt.close('all')
plt.figure(1)
plt.plot(rs['Time'],np.sum(~np.isnan(rs['FC SPH FG']),axis=1)/N_Sites,'-o')

plt.close('all')
plt.figure(1)
plt.plot(rs['Time'],np.nanmedian(rs['Surv TWS'],axis=1),'o')
plt.plot(rs['Time'],np.nanmedian(rs['Surv FG'],axis=1),'s')

plt.close('all')
plt.figure(1)
plt.plot(rs['Time'],np.nanpercentile(rs['Surv TWS'],90,axis=1),'s')
plt.plot(rs['Time'],np.nanpercentile(rs['Surv TWS'],10,axis=1),'v')

plt.close('all')
plt.figure(1)
plt.plot(rs['Time'],np.nanpercentile(rs['Surv FG'],90,axis=1),'s')
plt.plot(rs['Time'],np.nanpercentile(rs['Surv FG'],50,axis=1),'v')

plt.figure(2)
plt.plot(rs['Time'],np.nanmean(rs['FC SPH TWS'],axis=1),'o')
plt.plot(rs['Time'],np.nanpercentile(rs['FC SPH TWS'],90,axis=1),'s')
plt.plot(rs['Time'],np.nanpercentile(rs['FC SPH TWS'],10,axis=1),'v')

plt.figure(3)
plt.plot(rs['Time'],np.nanmean(rs['FC SPH FG'],axis=1),'o')
plt.plot(rs['Time'],np.nanpercentile(rs['FC SPH FG'],90,axis=1),'s')
plt.plot(rs['Time'],np.nanpercentile(rs['FC SPH FG'],10,axis=1),'v')


plt.figure(4)
plt.plot(rs['Time'],np.nanmean(rs['FC CC'],axis=1),'o')
plt.plot(rs['Time'],np.nanpercentile(rs['FC CC'],90,axis=1),'s')
plt.plot(rs['Time'],np.nanpercentile(rs['FC CC'],10,axis=1),'v')


plt.figure(5)
plt.plot(rs['Time'],np.nanmedian(rs['FC SPH TOT'],axis=1),'o')
plt.plot(rs['Time'],np.nanmean(rs['FC SPH TWS'],axis=1),'s')
plt.plot(rs['Time'],np.nanmean(rs['FC SPH FG'],axis=1),'^')


