# -*- coding: utf-8 -*-
"""
Created on Fri Oct 15 09:28:32 2021

@author: RHEMBER
"""

#%% Calculate responses

vd='RGR'

dC=np.zeros(uInst.size)
dF=np.zeros(uInst.size)
for iI in range(uInst.size):
    indC0=np.where( (dTL['ID_Inst']==uInst[iI]) & (dTL['Year']>=1971-5) & (dTL['Year']<=1970) & (dTL['Treatment']=='C') & (dTL[vd]>-1) & (dTL[vd]<100) )[0]
    indC1=np.where( (dTL['ID_Inst']==uInst[iI]) & (dTL['Year']>=1972) & (dTL['Year']<=1971+8) & (dTL['Treatment']=='C') & (dTL[vd]>-1) & (dTL[vd]<100) )[0]
    indF0=np.where( (dTL['ID_Inst']==uInst[iI]) & (dTL['Year']>=1971-5) & (dTL['Year']<=1970) & (dTL['Treatment']=='F1') & (dTL[vd]>-1) & (dTL[vd]<100) )[0]
    indF1=np.where( (dTL['ID_Inst']==uInst[iI]) & (dTL['Year']>=1972) & (dTL['Year']<=1971+8) & (dTL['Treatment']=='F1') & (dTL[vd]>-1) & (dTL[vd]<100) )[0]
    yC0=np.nanmean(dTL[vd][indC0])
    yC1=np.nanmean(dTL[vd][indC1])
    yF0=np.nanmean(dTL[vd][indF0])
    yF1=np.nanmean(dTL[vd][indF1])
    dC[iI]=yC1/yC0
    dF[iI]=yF1/yF0

plt.close('all')
plt.bar(np.arange(dC.size)-0.25,dC,0.4)
plt.bar(np.arange(dC.size)+0.25,dF,0.4)
#plt.bar(np.arange(dC.size)+0.25,dF-dC,0.4)

print(np.mean(dC))
print(np.mean(dF))
print(np.mean(dF)/np.mean(dC))

  
#%% Time response

vd='TRW Standardized RR'
vd='BAI Standardized RR'

bin=np.arange(1965,2022,1)
yC=np.zeros(bin.size)
yCP=np.zeros(bin.size)
yF=np.zeros(bin.size)
for i in range(bin.size):
    
    ind=np.where( (dTL['Year']==bin[i]) & (dTL['Treatment']=='C') & (dTL[vd]>0) & (dTL[vd]<2000) )[0]
    yC[i]=np.nanmedian(dTL[vd][ind])
    
    #ind=np.where( (dTL['ID_Inst']==id) & (dTL['Year']==bin[i]) & (dTL['Treatment']!='C') & (dTL['Treatment']!='F1') )[0]
    ind=np.where( (dTL['Year']==bin[i]) & (dTL['BGC SS']==1) & (dTL[vd]>0) & (dTL[vd]<2000) )[0]
    yCP[i]=np.nanmedian(dTL[vd][ind])
    
    ind=np.where( (dTL['Year']==bin[i]) & (dTL['Treatment']=='F1') & (dTL[vd]>0) & (dTL[vd]<2000) )[0]
    yF[i]=np.nanmedian(dTL[vd][ind])

plt.close('all')
plt.plot(bin,yC,'-bo')
#plt.plot(bin,yCP,'-cd')
plt.plot(bin,yF,'-gs')  

plt.figure()
plt.plot(bin,yF/yC,'-bo')
plt.plot(bin,yF/yCP,'-gs')

#%% Time response

vd='TRW Standardized RR'
vd='RGR Standardized RR'

id=uInst[7]

bin=np.arange(1950,2022,1)
yC=np.zeros(bin.size)
yCP=np.zeros(bin.size)
yF=np.zeros(bin.size)
for i in range(bin.size):
    ind=np.where( (dTL['ID_Inst']==id) & (dTL['Year']==bin[i]) & (dTL['Treatment']=='C') & (dTL[vd]>0) & (dTL[vd]<2000) )[0]
    yC[i]=np.nanmedian(dTL[vd][ind])
    #ind=np.where( (dTL['ID_Inst']==id) & (dTL['Year']==bin[i]) & (dTL['Treatment']!='C') & (dTL['Treatment']!='F1') )[0]
    ind=np.where( (dTL['ID_Inst']==id) & (dTL['Year']==bin[i]) & (dTL['BGC SS']==6) & (dTL[vd]>0) & (dTL[vd]<2000) )[0]
    yCP[i]=np.nanmedian(dTL[vd][ind])
    ind=np.where( (dTL['ID_Inst']==id) & (dTL['Year']==bin[i]) & (dTL['Treatment']=='F1') & (dTL[vd]>0) & (dTL[vd]<2000) )[0]
    yF[i]=np.nanmedian(dTL[vd][ind])

plt.close('all')
plt.plot(bin,yC,'-bo')
#plt.plot(bin,yCP,'-cd')
plt.plot(bin,yF,'-gs')

plt.figure()
plt.plot(bin,yF/yC,'-bo')

#%%

plt.close('all')
plt.plot(bin,yF-yCP,'-bo')
plt.plot(bin,yF-yC,'-cd')

#%% Age response

vd='TRW Standardized RR'
vd='Gsw Standardized RR'

bin=np.arange(1,80,1)
yC=np.zeros(bin.size)
yCP=np.zeros(bin.size)
yF=np.zeros(bin.size)
for i in range(bin.size):
    ind=np.where( (dTL['Time since first ring']==bin[i]) & (dTL['Treatment']=='C') )[0]
    yC[i]=np.nanmedian(dTL[vd][ind])
    #ind=np.where( (dTL['Age']==bin[i]) & (dTL['Treatment']!='C') & (dTL['Treatment']!='F1') )[0]
    ind=np.where( (dTL['Time since first ring']==bin[i]) & (dTL['BGC SS']==3) )[0]
    yCP[i]=np.nanmedian(dTL[vd][ind])
    ind=np.where( (dTL['Time since first ring']==bin[i]) & (dTL['Treatment']=='F1') )[0]
    yF[i]=np.nanmedian(dTL[vd][ind])

plt.close('all')
plt.plot(bin,yC,'-bo')
#plt.plot(bin,yCP,'-cd')
plt.plot(bin,yF,'-gs')

#%% Age response

id=uInst[5]

bin=np.arange(1,80,1)
yC=np.zeros(bin.size)
yCP=np.zeros(bin.size)
yF=np.zeros(bin.size)
for i in range(bin.size):
    ind=np.where( (dTL['ID_Inst']==id) & (dTL['Age']==bin[i]) & (dTL['Treatment']=='C') )[0]
    yC[i]=np.nanmean(dTL['Gsw Standardized RR'][ind])
    #ind=np.where( (dTL['ID_Inst']==id) & (dTL['Age']==bin[i]) & (dTL['Treatment']!='C') & (dTL['Treatment']!='F1') )[0]
    ind=np.where( (dTL['ID_Inst']==id) & (dTL['Age']==bin[i]) & (dTL['BGC SS']==3) )[0]
    yCP[i]=np.nanmean(dTL['Gsw Standardized RR'][ind])
    ind=np.where( (dTL['ID_Inst']==id) & (dTL['Age']==bin[i]) & (dTL['Treatment']=='F1') )[0]
    yF[i]=np.nanmean(dTL['Gsw Standardized RR'][ind])

plt.close('all')
plt.plot(bin,yC,'-bo')
#plt.plot(bin,yCP,'-cd')
plt.plot(bin,yF,'-gs')

#%% Biomass response

id=uInst[0]

bin=np.arange(5,150,5)
yC=np.zeros(bin.size)
yCP=np.zeros(bin.size)
yF=np.zeros(bin.size)
for i in range(bin.size):
    ind=np.where( (dTL['ID_Inst']==id) & (np.abs(dTL['Bsw']-bin[i])<=2) & (dTL['Treatment']=='C') )[0]
    yC[i]=np.nanmean(dTL['Gsw'][ind])
    ind=np.where( (dTL['ID_Inst']==id) & (np.abs(dTL['Bsw']-bin[i])<=2) & (dTL['Treatment']!='C') & (dTL['Treatment']!='F1') )[0]
    yCP[i]=np.nanmean(dTL['Gsw'][ind])
    ind=np.where( (dTL['ID_Inst']==id) & (np.abs(dTL['Bsw']-bin[i])<=2) & (dTL['Treatment']=='F1') )[0]
    yF[i]=np.nanmean(dTL['Gsw'][ind])

plt.close('all')
plt.plot(bin,yC,'-bo')
#plt.plot(bin,yCP,'-cd')
plt.plot(bin,yF,'-gs')

#%% Comparison of growth trends across site series

vd='Gsw'
tv=np.arange(1950,2022,1)
y=np.zeros((tv.size,uSS.size))
for iSS in range(uSS.size):
    for iT in range(tv.size):
        ind=np.where( (dTL['Year First']<1966) & (dTL['ID_Inst']==uInst[0]) & (dTL['Year']==tv[iT]) & (dTL['BGC SS']==uSS[iSS]) & (dTL['Treatment']!='F1') & (dTL[vd]>0) & (dTL[vd]<2000) )[0]
        y[iT,iSS]=np.nanmedian(dTL[vd][ind])
    
plt.close('all')
plt.plot(tv,y[:,0],'-ko')
plt.plot(tv,y[:,1],'-ro')
plt.plot(tv,y[:,2],'-bo')
plt.plot(tv,y[:,3],'-gs')

#%%

vd='Gsw'
tv=np.arange(1950,2022,1)
y=np.zeros((tv.size,uSS.size))
for iSS in range(uSS.size):
    for iT in range(tv.size):
        ind=np.where( (dTL['Year First']<1966) & (dTL['Year']==tv[iT]) & (dTL['BGC SS']==uSS[iSS]) & (dTL['Treatment']!='F1') & (dTL[vd]>0) & (dTL[vd]<2000) )[0]
        y[iT,iSS]=np.nanmedian(dTL[vd][ind])
    #y[:,iSS]=y[:,iSS]-gu.movingave(y[:,iSS],5,'Centre')
    
plt.close('all')
plt.plot(tv,y[:,0],'-ko')
plt.plot(tv,y[:,1],'-ro')
plt.plot(tv,y[:,2],'-bo')
plt.plot(tv,y[:,3],'-gs')

#%%

tv=np.arange(1,80,1)
y=np.zeros((tv.size,uSS.size))
for iSS in range(uSS.size):
    for iT in range(tv.size):
        ind=np.where( (dTL['Year First']<1966) & (dTL['Time since first ring']==tv[iT]) & (dTL['BGC SS']==uSS[iSS]) & (dTL['Treatment']!='F1') )[0]
        y[iT,iSS]=np.nanmedian(dTL['Gsw Standardized RR'][ind])
    
plt.close('all')
plt.plot(tv,y[:,0],'-ko')
plt.plot(tv,y[:,1],'-ro')
plt.plot(tv,y[:,2],'-bo')
plt.plot(tv,y[:,3],'-gs')


#%% Comparison between Dob and H in 2020

# Isolate data
iC=np.where( (dTL['Dob 2020']>0) & (dTL['H 2020']>0) & (dTL['Treatment']=='C') )[0]
iF1=np.where( (dTL['Dob 2020']>0) & (dTL['H 2020']>0) & (dTL['Treatment']=='F1') )[0]

# Fit linear best fit relationship
yC=dTL['H 2020'][iC]
xC=dTL['Dob 2020'][iC]
x1C=sm.tools.tools.add_constant(xC)
mdC=sm.OLS(yC,x1C).fit()
mdC.summary()
xhat=np.linspace(np.min(x1C[:,1]),np.max(x1C[:,1]),10)
yhatC=mdC.predict(np.c_[np.ones(xhat.size),xhat])

yF1=dTL['H 2020'][iF1]
xF1=dTL['Dob 2020'][iF1]
x1F1=sm.tools.tools.add_constant(xF1)
mdF1=sm.OLS(yF1,x1F1).fit()
mdF1.summary()
yhatF1=mdF1.predict(np.c_[np.ones(xhat.size),xhat])


# Plot
plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(8.5,8.5))
#ax.plot([0,90],[0,90],'k-',lw=2,color=[0.8,0.8,0.8])
ax.plot(xC,yC,'o',ms=4,mec='w',mfc='b')
ax.plot(xF1,yF1,'o',ms=4,mec='w',mfc='r')

ax.plot(xhat,yhatC,'b-',lw=1.25)
ax.plot(xhat,yhatF1,'r--',lw=1.25)
ax.set(xlim=[0,80],ylim=[0,80],xlabel='Inside-bark diameter from cores 2019 (cm)',ylabel='Outside-bark diameter from 2020 (cm)')
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')
plt.tight_layout()