'''

EP703 RETROSPECTIVE MONITORING STUDY - ANALYSIS - EVALUATING THE RM METHOD

'''

#%% Import modules

from os import listdir
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from fcgadgets.macgyver import utilities_general as gu
import scipy.io
import copy
from scipy.optimize import curve_fit
import scipy.io as spio
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats, linalg
from matplotlib.patches import Rectangle
from scipy.io import loadmat

#%% Save figures?

flag_savefigs='On'

#%% Import data

# Run import script

# Import EP703 tree-level data
dEPTL=gu.ipickle(r'C:\Users\rhember\Documents\Data\EP703\Processed\EP703_TL.pkl')
ID_FD_EPTL=1

#%% Define reference periods before and after N application

t_ref0=[1970-4,1970] # 5 years
t_ref1=[1972,1972+9] # 10 years

#%% Custom function for calculating the C.I.s for response ratios

def GetCIsForRR(lo0,hi0,lo1,hi1):
    tmp=np.array([lo1/lo0,lo1/hi0,hi1/lo0,hi1/hi0]).T
    lo=np.min(tmp,axis=0)
    hi=np.max(tmp,axis=0)    
    return lo,hi

# Custom function for calculating the C.I.s for % differences
def GetCIsForDR(lo0,hi0,lo1,hi1):
    a=(lo1-lo0)/lo0*100
    b=(lo1-hi0)/hi0*100 
    c=(hi1-lo0)/lo0*100 
    d=(hi1-hi0)/hi0*100
    tmp=np.array([ a,b,c,d ]).T
    lo=np.min(tmp,axis=1)
    hi=np.max(tmp,axis=1)    
    return lo,hi

#%% Calculate error variance in ratio

ind=np.where( (dTL['Treatment']!='F') )[0]   
uT=np.unique(dTL['ID_Tree_Unique'][ind])

var='TRW'

dS={}
dS['All']=np.zeros(uT.size)
dS['C']=np.zeros(uT.size)
dS['SS3']=np.zeros(uT.size)
dS['SS5']=np.zeros(uT.size)
for iT in range(uT.size):
    ind0=np.where( (dTL['ID_Tree_Unique']==uT[iT]) & (dTL['Year']>=t_ref0[0]) & (dTL['Year']<=t_ref0[1]) & (dTL[var]>-1) & (dTL[var]<5000) & (dTL['QA Summary']==1) )[0]
    ind1=np.where( (dTL['ID_Tree_Unique']==uT[iT]) & (dTL['Year']>=t_ref1[0]) & (dTL['Year']<=t_ref1[1]) & (dTL[var]>-1) & (dTL[var]<5000) & (dTL['QA Summary']==1) )[0]
    dS['All'][iT]=np.nanmean(dTL[var][ind1])/np.nanmean(dTL[var][ind0])
    ind0=np.where( (dTL['ID_Tree_Unique']==uT[iT]) & (dTL['Treatment']=='C') & (dTL['Year']>=t_ref0[0]) & (dTL['Year']<=t_ref0[1]) & (dTL[var]>-1) & (dTL[var]<5000) & (dTL['QA Summary']==1) )[0]
    ind1=np.where( (dTL['ID_Tree_Unique']==uT[iT]) & (dTL['Treatment']=='C') & (dTL['Year']>=t_ref1[0]) & (dTL['Year']<=t_ref1[1]) & (dTL[var]>-1) & (dTL[var]<5000) & (dTL['QA Summary']==1) )[0]
    dS['C'][iT]=np.nanmean(dTL[var][ind1])/np.nanmean(dTL[var][ind0])
    ind0=np.where( (dTL['ID_Tree_Unique']==uT[iT]) & (dTL['BGC SS']==3) & (dTL['Year']>=t_ref0[0]) & (dTL['Year']<=t_ref0[1]) & (dTL[var]>-1) & (dTL[var]<5000) & (dTL['QA Summary']==1) )[0]
    ind1=np.where( (dTL['ID_Tree_Unique']==uT[iT]) & (dTL['BGC SS']==3) & (dTL['Year']>=t_ref1[0]) & (dTL['Year']<=t_ref1[1]) & (dTL[var]>-1) & (dTL[var]<5000) & (dTL['QA Summary']==1) )[0]
    dS['SS3'][iT]=np.nanmean(dTL[var][ind1])/np.nanmean(dTL[var][ind0])
    ind0=np.where( (dTL['ID_Tree_Unique']==uT[iT]) & (dTL['BGC SS']==5) & (dTL['Year']>=t_ref0[0]) & (dTL['Year']<=t_ref0[1]) & (dTL[var]>-1) & (dTL[var]<5000) & (dTL['QA Summary']==1) )[0]
    ind1=np.where( (dTL['ID_Tree_Unique']==uT[iT]) & (dTL['BGC SS']==5) & (dTL['Year']>=t_ref1[0]) & (dTL['Year']<=t_ref1[1]) & (dTL[var]>-1) & (dTL[var]<5000) & (dTL['QA Summary']==1) )[0]
    dS['SS5'][iT]=np.nanmean(dTL[var][ind1])/np.nanmean(dTL[var][ind0])

bw=0.1; bin=np.arange(0,2+bw,bw)
y=np.zeros(bin.size)
for i in range(bin.size):
    ind=np.where(np.abs(dS['All']-bin[i])<=bw/2)[0]
    y[i]=ind.size

mu=np.nanmean(dS['All'])
sd=np.nanstd(dS['All'])
se=sd/np.sqrt(dS['All'].size)
cv=sd/mu*100

mu=np.nanmean(dS['SS3'])
sd=np.nanstd(dS['SS3'])
se=sd/np.sqrt(dS['SS3'].size)
cv=sd/mu*100

plt.close('all'); 
fig,ax=plt.subplots(1,figsize=gu.cm2inch(12,7)); 
#h,bins=np.histogram(dS['D'],bins=np.arange(0,2,0.1))
#ax.bar(bins[1:],h,1,fc=cl1,label='Control plot');
ax.bar(bin,y,1,fc=[0.29,0.47,0.74],label='Control plot');

#%% Calculate the variation in ring width after application (for unfertilized trees)

dS={}
dS['D']=np.array([])
for iI in range(uInst.size):
    ind=np.where( (dTL['ID_Inst']==uInst[iI]) & (dTL['Treatment']!='F') & (dTL['Year']>=t_ref1[0]) & (dTL['Year']<=t_ref1[1]) & (dTL[var]>-1) & (dTL[var]<5000) & (dTL['QA Summary']==1) )[0]
    uT=np.unique(dTL['ID_Tree_Unique'][ind])
    for iR in range(10):
        ids=np.random.choice(uT.size,2,replace=False)
    
        ind1=np.where( (dTL['ID_Inst']==uInst[iI]) & (dTL['ID_Tree_Unique']==uT[ids[0]]) & (dTL['Treatment']!='F') & (dTL['Year']>=t_ref1[0]) & (dTL['Year']<=t_ref1[1]) & (dTL[var]>-1) & (dTL[var]<5000) & (dTL['QA Summary']==1) )[0]
        ind2=np.where( (dTL['ID_Inst']==uInst[iI]) & (dTL['ID_Tree_Unique']==uT[ids[1]]) & (dTL['Treatment']!='F') & (dTL['Year']>=t_ref1[0]) & (dTL['Year']<=t_ref1[1]) & (dTL[var]>-1) & (dTL[var]<5000) & (dTL['QA Summary']==1) )[0]
        dS['D']=np.append(dS['D'],(np.nanmean(dTL[var][ind2])-np.nanmean(dTL[var][ind1]))/np.nanmean(dTL[var][ind1])*100)

np.mean(np.abs(dS['D']))

#%% Response (following Ballard and Majid 1984)
    
# Initialize
varL=['TRW','BAI','Gsw','RGR','TRW S NE DPLR','TRW RR','BAI RR','TRW S NE DPLR RR']
trL=['D','R']

dBM={}
for var in varL:
    dBM[var]={}
    for tr in trL:        
        dBM[var][tr]={}
        dBM[var][tr]['Mean']=np.zeros(uInst.size)
        dBM[var][tr]['Median']=np.zeros(uInst.size)

# Populate
for iI in range(uInst.size):
    for var in varL:
        
        # Indices
        indC0=np.where( (dTL['ID_Inst']==uInst[iI]) & (dTL['Year']>=t_ref0[0]) & (dTL['Year']<=t_ref0[1]) & (dTL['Treatment']=='C') & (dTL[var]>-1) & (dTL[var]<5000) & (dTL['QA Summary']==1) )[0]
        indC1=np.where( (dTL['ID_Inst']==uInst[iI]) & (dTL['Year']>=t_ref1[0]) & (dTL['Year']<=t_ref1[1]) & (dTL['Treatment']=='C') & (dTL[var]>-1) & (dTL[var]<5000) & (dTL['QA Summary']==1) )[0]
        indF0=np.where( (dTL['ID_Inst']==uInst[iI]) & (dTL['Year']>=t_ref0[0]) & (dTL['Year']<=t_ref0[1]) & (dTL['Treatment']=='F1') & (dTL[var]>-1) & (dTL[var]<5000) & (dTL['QA Summary']==1) )[0]
        indF1=np.where( (dTL['ID_Inst']==uInst[iI]) & (dTL['Year']>=t_ref1[0]) & (dTL['Year']<=t_ref1[1]) & (dTL['Treatment']=='F1') & (dTL[var]>-1) & (dTL[var]<5000) & (dTL['QA Summary']==1) )[0]
        indSS30=np.where( (dTL['ID_Inst']==uInst[iI]) & (dTL['Year']>=t_ref0[0]) & (dTL['Year']<=t_ref0[1]) & (dTL['BGC SS']==3) & (dTL[var]>-1) & (dTL[var]<5000) & (dTL['Treatment']!='C') & (dTL['Treatment']!='F1') & (dTL['QA Summary']==1) )[0]
        indSS31=np.where( (dTL['ID_Inst']==uInst[iI]) & (dTL['Year']>=t_ref1[0]) & (dTL['Year']<=t_ref1[1]) & (dTL['BGC SS']==3) & (dTL[var]>-1) & (dTL[var]<5000) & (dTL['Treatment']!='C') & (dTL['Treatment']!='F1') & (dTL['QA Summary']==1) )[0]
        indSS50=np.where( (dTL['ID_Inst']==uInst[iI]) & (dTL['Year']>=t_ref0[0]) & (dTL['Year']<=t_ref0[1]) & (dTL['BGC SS']==5) & (dTL[var]>-1) & (dTL[var]<5000) & (dTL['Treatment']!='C') & (dTL['Treatment']!='F1') & (dTL['QA Summary']==1) )[0]
        indSS51=np.where( (dTL['ID_Inst']==uInst[iI]) & (dTL['Year']>=t_ref1[0]) & (dTL['Year']<=t_ref1[1]) & (dTL['BGC SS']==5) & (dTL[var]>-1) & (dTL[var]<5000) & (dTL['Treatment']!='C') & (dTL['Treatment']!='F1') & (dTL['QA Summary']==1) )[0]
        indSS60=np.where( (dTL['ID_Inst']==uInst[iI]) & (dTL['Year']>=t_ref0[0]) & (dTL['Year']<=t_ref0[1]) & (dTL['BGC SS']==6) & (dTL[var]>-1) & (dTL[var]<5000) & (dTL['Treatment']!='C') & (dTL['Treatment']!='F1') & (dTL['QA Summary']==1) )[0]
        indSS61=np.where( (dTL['ID_Inst']==uInst[iI]) & (dTL['Year']>=t_ref1[0]) & (dTL['Year']<=t_ref1[1]) & (dTL['BGC SS']==6) & (dTL[var]>-1) & (dTL[var]<5000) & (dTL['Treatment']!='C') & (dTL['Treatment']!='F1') & (dTL['QA Summary']==1) )[0]
    
        # Median 
        yC0=np.nanmedian(dTL[var][indC0]); yC1=np.nanmedian(dTL[var][indC1]); yF0=np.nanmedian(dTL[var][indF0]); yF1=np.nanmedian(dTL[var][indF1]); ySS30=np.nanmedian(dTL[var][indSS30]); ySS31=np.nanmedian(dTL[var][indSS31]); ySS50=np.nanmedian(dTL[var][indSS50]); ySS51=np.nanmedian(dTL[var][indSS51]); ySS60=np.nanmedian(dTL[var][indSS60]); ySS61=np.nanmedian(dTL[var][indSS61]);
        dBM[var]['D']['Median'][iI]=yF1-yF0
        dBM[var]['R']['Median'][iI]=yF1-yF0*(yC1/yC0)
        
        # Mean
        muC0=np.nanmean(dTL[var][indC0]); muC1=np.nanmean(dTL[var][indC1]); muF0=np.nanmean(dTL[var][indF0]); muF1=np.nanmean(dTL[var][indF1]); muSS30=np.nanmean(dTL[var][indSS30]); muSS31=np.nanmean(dTL[var][indSS31]); muSS50=np.nanmean(dTL[var][indSS50]); muSS51=np.nanmean(dTL[var][indSS51]); muSS60=np.nanmean(dTL[var][indSS60]); muSS61=np.nanmean(dTL[var][indSS61]);
        dBM[var]['D']['Mean'][iI]=yF1-yF0
        dBM[var]['R']['Mean'][iI]=yF1-yF0*(yC1/yC0)

# Plot
stat='Mean'        
var='TRW'
plt.close('all'); 
fig,ax=plt.subplots(1,figsize=gu.cm2inch(8,5)); 
ax.plot([-1,11],[0,0],'k-')
ax.bar(np.arange(uInst.size),dBM[var]['R'][stat],0.6,fc=[0.7,0.7,0.7])
ax.bar(uInst.size,np.mean(dBM[var]['R'][stat]),0.6,fc=[0.35,0.35,0.35])
ax.set(position=[0.15,0.15,0.8,0.8],xlim=[-0.5,uInst.size+0.5],xticks=np.arange(uInst.size+1),xticklabels=np.append(uInst,'All'),ylim=[-0.2,0.5], \
       yticks=np.arange(-0.3,0.6,0.1),ylabel='$\itR$$_{t1-10}$',xlabel='Installation ID')
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=1.5)
if flag_savefigs=='On':
    gu.PrintFig(meta['Paths']['Figures'] + '\\Response_BM84' + var,'png',900)
    fig.savefig(meta['Paths']['Figures'] + '\\Response_BM84_' + var + '.svg',format='svg',dpi=1200)


#%% Calculate response ratio

# Initialize
varL=['TRW','BAI','Gsw','RGR','TRW S NE DPLR','TRW RR','BAI RR','TRW S NE DPLR RR']
treatL=['C','F1','SS3','SS5','SS6']
dBM={}
for var in varL:
    dBM[var]={}
    for tn in treatL:
        dBM[var][tn]={}
        dBM[var][tn]['Mean']=np.zeros(uInst.size)
        dBM[var][tn]['Median']=np.zeros(uInst.size)
        #dBM[var][tn]['SD']=np.zeros(uInst.size)
        dBM[var][tn]['CI Lo']=np.zeros(uInst.size)
        dBM[var][tn]['CI Hi']=np.zeros(uInst.size)

# Populate
for iI in range(uInst.size):
    for var in varL:
        
        # Indices
        indC0=np.where( (dTL['ID_Inst']==uInst[iI]) & (dTL['Year']>=t_ref0[0]) & (dTL['Year']<=t_ref0[1]) & (dTL['Treatment']=='C') & (dTL[var]>-1) & (dTL[var]<5000) & (dTL['QA Summary']==1) )[0]
        indC1=np.where( (dTL['ID_Inst']==uInst[iI]) & (dTL['Year']>=t_ref1[0]) & (dTL['Year']<=t_ref1[1]) & (dTL['Treatment']=='C') & (dTL[var]>-1) & (dTL[var]<5000) & (dTL['QA Summary']==1) )[0]
        indF0=np.where( (dTL['ID_Inst']==uInst[iI]) & (dTL['Year']>=t_ref0[0]) & (dTL['Year']<=t_ref0[1]) & (dTL['Treatment']=='F1') & (dTL[var]>-1) & (dTL[var]<5000) & (dTL['QA Summary']==1) )[0]
        indF1=np.where( (dTL['ID_Inst']==uInst[iI]) & (dTL['Year']>=t_ref1[0]) & (dTL['Year']<=t_ref1[1]) & (dTL['Treatment']=='F1') & (dTL[var]>-1) & (dTL[var]<5000) & (dTL['QA Summary']==1) )[0]
        indSS30=np.where( (dTL['ID_Inst']==uInst[iI]) & (dTL['Year']>=t_ref0[0]) & (dTL['Year']<=t_ref0[1]) & (dTL['BGC SS']==3) & (dTL[var]>-1) & (dTL[var]<5000) & (dTL['Treatment']!='C') & (dTL['Treatment']!='F1') & (dTL['QA Summary']==1) )[0]
        indSS31=np.where( (dTL['ID_Inst']==uInst[iI]) & (dTL['Year']>=t_ref1[0]) & (dTL['Year']<=t_ref1[1]) & (dTL['BGC SS']==3) & (dTL[var]>-1) & (dTL[var]<5000) & (dTL['Treatment']!='C') & (dTL['Treatment']!='F1') & (dTL['QA Summary']==1) )[0]
        indSS50=np.where( (dTL['ID_Inst']==uInst[iI]) & (dTL['Year']>=t_ref0[0]) & (dTL['Year']<=t_ref0[1]) & (dTL['BGC SS']==5) & (dTL[var]>-1) & (dTL[var]<5000) & (dTL['Treatment']!='C') & (dTL['Treatment']!='F1') & (dTL['QA Summary']==1) )[0]
        indSS51=np.where( (dTL['ID_Inst']==uInst[iI]) & (dTL['Year']>=t_ref1[0]) & (dTL['Year']<=t_ref1[1]) & (dTL['BGC SS']==5) & (dTL[var]>-1) & (dTL[var]<5000) & (dTL['Treatment']!='C') & (dTL['Treatment']!='F1') & (dTL['QA Summary']==1) )[0]
        indSS60=np.where( (dTL['ID_Inst']==uInst[iI]) & (dTL['Year']>=t_ref0[0]) & (dTL['Year']<=t_ref0[1]) & (dTL['BGC SS']==6) & (dTL[var]>-1) & (dTL[var]<5000) & (dTL['Treatment']!='C') & (dTL['Treatment']!='F1') & (dTL['QA Summary']==1) )[0]
        indSS61=np.where( (dTL['ID_Inst']==uInst[iI]) & (dTL['Year']>=t_ref1[0]) & (dTL['Year']<=t_ref1[1]) & (dTL['BGC SS']==6) & (dTL[var]>-1) & (dTL[var]<5000) & (dTL['Treatment']!='C') & (dTL['Treatment']!='F1') & (dTL['QA Summary']==1) )[0]
    
        # Median 
        yC0=np.nanmedian(dTL[var][indC0]); yC1=np.nanmedian(dTL[var][indC1]); yF0=np.nanmedian(dTL[var][indF0]); yF1=np.nanmedian(dTL[var][indF1]); ySS30=np.nanmedian(dTL[var][indSS30]); ySS31=np.nanmedian(dTL[var][indSS31]); ySS50=np.nanmedian(dTL[var][indSS50]); ySS51=np.nanmedian(dTL[var][indSS51]); ySS60=np.nanmedian(dTL[var][indSS60]); ySS61=np.nanmedian(dTL[var][indSS61]);
        dBM[var]['C']['Median'][iI]=yC1/yC0
        dBM[var]['F1']['Median'][iI]=yF1/yF0
        dBM[var]['SS3']['Median'][iI]=ySS31/ySS30
        dBM[var]['SS5']['Median'][iI]=ySS51/ySS50
        dBM[var]['SS6']['Median'][iI]=ySS61/ySS60
        
        # Mean
        muC0=np.nanmean(dTL[var][indC0]); muC1=np.nanmean(dTL[var][indC1]); muF0=np.nanmean(dTL[var][indF0]); muF1=np.nanmean(dTL[var][indF1]); muSS30=np.nanmean(dTL[var][indSS30]); muSS31=np.nanmean(dTL[var][indSS31]); muSS50=np.nanmean(dTL[var][indSS50]); muSS51=np.nanmean(dTL[var][indSS51]); muSS60=np.nanmean(dTL[var][indSS60]); muSS61=np.nanmean(dTL[var][indSS61]);
        dBM[var]['C']['Mean'][iI]=muC1/muC0
        dBM[var]['F1']['Mean'][iI]=muF1/muF0
        dBM[var]['SS3']['Mean'][iI]=muSS31/muSS30
        dBM[var]['SS5']['Mean'][iI]=muSS51/muSS50
        dBM[var]['SS6']['Mean'][iI]=muSS61/muSS60
        
        # Error (1 x S.E.)
        multip=1
        eC0=multip*np.nanstd(dTL[var][indC0])/np.sqrt(indC0.size); 
        eC1=multip*np.nanstd(dTL[var][indC1])/np.sqrt(indC1.size);         
        loC0=muC0-eC0
        hiC0=muC0+eC0
        loC1=muC1-eC1
        hiC1=muC1+eC1
        dBM[var]['C']['CI Lo'][iI],dBM[var]['C']['CI Hi'][iI]=GetCIsForRR(loC0,hiC0,loC1,hiC1)
        
        eF0=multip*np.nanstd(dTL[var][indF0])/np.sqrt(indF0.size); 
        eF1=multip*np.nanstd(dTL[var][indF1])/np.sqrt(indF1.size); 
        loF0=muF0-eF0
        hiF0=muF0+eF0
        loF1=muF1-eF1
        hiF1=muF1+eF1
        dBM[var]['F1']['CI Lo'][iI],dBM[var]['F1']['CI Hi'][iI]=GetCIsForRR(loF0,hiF0,loF1,hiF1)
        
        eSS30=multip*np.nanstd(dTL[var][indSS30])/np.sqrt(indSS30.size); 
        eSS31=multip*np.nanstd(dTL[var][indSS31])/np.sqrt(indSS31.size); 
        loSS30=muSS30-eSS30
        hiSS30=muSS30+eSS30
        loSS31=muSS31-eSS31
        hiSS31=muSS31+eSS31
        dBM[var]['SS3']['CI Lo'][iI],dBM[var]['SS3']['CI Hi'][iI]=GetCIsForRR(loSS30,hiSS30,loSS31,hiSS31)
        
        eSS50=multip*np.nanstd(dTL[var][indSS50])/np.sqrt(indSS50.size); 
        eSS51=multip*np.nanstd(dTL[var][indSS51])/np.sqrt(indSS51.size); 
        loSS50=muSS50-eSS50
        hiSS50=muSS50+eSS50
        loSS51=muSS51-eSS51
        hiSS51=muSS51+eSS51
        dBM[var]['SS5']['CI Lo'][iI],dBM[var]['SS5']['CI Hi'][iI]=GetCIsForRR(loSS50,hiSS50,loSS51,hiSS51)
        
        eSS60=multip*np.nanstd(dTL[var][indSS60])/np.sqrt(indSS60.size); 
        eSS61=multip*np.nanstd(dTL[var][indSS61])/np.sqrt(indSS61.size); 
        loSS60=muSS60-eSS60
        hiSS60=muSS60+eSS60
        loSS61=muSS61-eSS61
        hiSS61=muSS61+eSS61
        dBM[var]['SS6']['CI Lo'][iI],dBM[var]['SS6']['CI Hi'][iI]=GetCIsForRR(loSS60,hiSS60,loSS61,hiSS61)

#%% Calculate time series of response ratio

varL=['TRW RR','TRW S NE DPLR RR','BAI RR','Gsw RR','RGR RR']
treatL=['C','F','SS3','SS5']
statL=['Mean','Median','SE']

# Time period of analysis
binT=np.arange(1950,2022,1)

# Initialize dictionary of results
dTS={}
for v in varL:
    dTS[v]={}
    for treat in treatL:
        dTS[v][treat]={}
        for s in statL:
            dTS[v][treat][s]=np.zeros(binT.size)

# Populate dictionary of results
for vd in varL:
    for iT in range(binT.size):
        ind=np.where( (dTL['Year']==binT[iT]) & (dTL['Treatment']=='C') & (dTL[vd]>0) & (dTL[vd]<2000) & (dTL['QA Summary']==1) )[0]
        dTS[vd]['C']['Mean'][iT]=np.nanmean(dTL[vd][ind])
        dTS[vd]['C']['Median'][iT]=np.nanmedian(dTL[vd][ind])
        dTS[vd]['C']['SE'][iT]=np.nanstd(dTL[vd][ind])/np.sqrt(ind.size)
        ind=np.where( (dTL['Year']==binT[iT]) & (dTL['Treatment']=='F1') & (dTL[vd]>0) & (dTL[vd]<2000) & (dTL['QA Summary']==1) )[0]
        dTS[vd]['F']['Mean'][iT]=np.nanmean(dTL[vd][ind])
        dTS[vd]['F']['Median'][iT]=np.nanmedian(dTL[vd][ind])
        dTS[vd]['F']['SE'][iT]=np.nanstd(dTL[vd][ind])/np.sqrt(ind.size)
        ind=np.where( (dTL['Year']==binT[iT]) & (dTL['BGC SS']==3) & (dTL[vd]>0) & (dTL[vd]<2000) & (dTL['Treatment']!='C') & (dTL['Treatment']!='F1') & (dTL['QA Summary']==1) )[0]
        dTS[vd]['SS3']['Mean'][iT]=np.nanmean(dTL[vd][ind])
        dTS[vd]['SS3']['Median'][iT]=np.nanmean(dTL[vd][ind])
        dTS[vd]['SS3']['SE'][iT]=np.nanstd(dTL[vd][ind])/np.sqrt(ind.size)
        ind=np.where( (dTL['Year']==binT[iT]) & (dTL['BGC SS']==5) & (dTL[vd]>0) & (dTL[vd]<2000) & (dTL['Treatment']!='C') & (dTL['Treatment']!='F1') & (dTL['QA Summary']==1) )[0]
        #ind=np.where( (dTL['Year']==binT[iT]) & (dTL['BGC SS Comb']==99) & (dTL[vd]>0) & (dTL[vd]<2000) & (dTL['Treatment']!='C') & (dTL['Treatment']!='F1') & (dTL['QA Summary']==1) )[0]
        dTS[vd]['SS5']['Mean'][iT]=np.nanmean(dTL[vd][ind])
        dTS[vd]['SS5']['Median'][iT]=np.nanmedian(dTL[vd][ind])
        dTS[vd]['SS5']['SE'][iT]=np.nanstd(dTL[vd][ind])/np.sqrt(ind.size)
    
#%% Bar charts of change in C and F1 plots (between pre- and post-application reference periods)

var='TRW'
stat='Mean'

cl1=[0.8,0.8,0.8]; cl2=[0.55,0.55,0.55]

plt.close('all'); 
fig,ax=plt.subplots(3,1,figsize=gu.cm2inch(8,11)); 

# Both C and F
ax[0].bar(np.arange(uInst.size)-0.15,dBM[var]['C'][stat],0.3,fc=cl1,label='Unfertilized');
ax[0].bar(np.arange(uInst.size)+0.15,dBM[var]['F1'][stat],0.3,fc=cl2,label='Fertilized');
ax[0].errorbar(np.arange(uInst.size)-0.15,dBM[var]['C'][stat],yerr=[ dBM[var]['C'][stat]-dBM[var]['C']['CI Lo'], dBM[var]['C']['CI Hi']-dBM[var]['C'][stat] ],color='k',ls='',capsize=2)
ax[0].errorbar(np.arange(uInst.size)+0.15,dBM[var]['F1'][stat],yerr=[ dBM[var]['F1'][stat]-dBM[var]['F1']['CI Lo'], dBM[var]['F1']['CI Hi']-dBM[var]['C'][stat] ],color='k',ls='',capsize=2)
#ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both')
ax[0].set(position=[0.16,0.7,0.8,0.29],xlim=[-0.5,uInst.size-0.5],xticks=np.arange(uInst.size),xticklabels='',ylim=[0.5,1.35], \
  yticks=np.arange(0.4,2,0.2),ylabel='$\itr_f$ or $\itr_u$') # '$\itw$ (mm yr$^-$$^1$)')
ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].tick_params(length=1.5)
ax[0].legend(loc='upper left',frameon=False,fontsize=6,bbox_to_anchor=(0.55,0.45,0.5,0.5))

# Absolute differences
da_be=dBM[var]['F1'][stat]-dBM[var]['C'][stat]
da_lo,da_hi=gu.GetCIsFromDifference(dBM[var]['C']['CI Lo'],dBM[var]['C']['CI Hi'],dBM[var]['F1']['CI Lo'],dBM[var]['F1']['CI Hi'])

da_be_mu=np.mean(da_be)
da_be_se=1*np.nanstd(da_be)/np.sqrt(da_be_mu.size)

da_be=np.append(da_be,da_be_mu)
da_lo=np.append(da_lo,da_be_mu-da_be_se)
da_hi=np.append(da_hi,da_be_mu+da_be_se)

#plt.close('all'); 
#fig,ax=plt.subplots(1,figsize=gu.cm2inch(14.5,5)); 
ax[1].plot([-1,20],[0,0],'k-')
ax[1].bar(np.arange(uInst.size+1),da_be,0.5,fc=cl2)
ax[1].errorbar(np.arange(uInst.size+1),da_be,yerr=[ da_be-da_lo, da_hi-da_be ],color='k',ls='',capsize=2)
#ax[1].yaxis.set_ticks_position('both'); ax[1].xaxis.set_ticks_position('both')
ax[1].set(position=[0.16,0.39,0.8,0.29],xlim=[-0.5,uInst.size+0.5],xticks=np.arange(uInst.size+1),xticklabels='', \
  yticks=np.arange(-0.2,0.4,0.1),ylabel='$\itI$$_{a,t1-10}$')
ax[1].yaxis.set_ticks_position('both'); ax[1].xaxis.set_ticks_position('both'); ax[1].tick_params(length=1.5)

# Relative differences
dr_be=(dBM[var]['F1'][stat]-dBM[var]['C'][stat])/dBM[var]['C'][stat]*100
dr_lo,dr_hi=GetCIsForDR(dBM[var]['C']['CI Lo'],dBM[var]['C']['CI Hi'],dBM[var]['F1']['CI Lo'],dBM[var]['F1']['CI Hi'])
dr_be_mu=np.mean(dr_be)
dr_be_se=1*np.nanstd(dr_be)/np.sqrt(dr_be_mu.size)
dr_be=np.append(dr_be,dr_be_mu)
dr_lo=np.append(dr_lo,dr_be_mu-dr_be_se)
dr_hi=np.append(dr_hi,dr_be_mu+dr_be_se)

#plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(14.5,5)); 
ax[2].plot([-1,20],[0,0],'k-')
ax[2].bar(np.arange(uInst.size+1),dr_be,0.5,fc=cl2)
ax[2].errorbar(np.arange(uInst.size+1),dr_be,yerr=[ dr_be-dr_lo, dr_hi-dr_be ],color='k',ls='',capsize=2)
#ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')
ax[2].set(position=[0.16,0.08,0.8,0.29],xlim=[-0.5,uInst.size+0.5],xticks=np.arange(uInst.size+1),xticklabels=np.append(uInst,'All'), \
  ylim=[-25,50],yticks=np.arange(-20,50,10),ylabel='$\itI$$_{r,t1-10}$ (%)',xlabel='Installation ID')
ax[2].yaxis.set_ticks_position('both'); ax[2].xaxis.set_ticks_position('both'); ax[2].tick_params(length=1.5)
gu.axletters(ax,plt,0.0185,0.875,FontWeight='Bold') # ,LetterStyle='Caps'
if flag_savefigs=='On':
    gu.PrintFig(meta['Paths']['Figures'] + '\\Response_' + var,'png',900)
    fig.savefig(meta['Paths']['Figures'] + '\\Response_' + var + '.svg',format='svg',dpi=1200)

#%% Plot time series

# Central tendancy
stat='Mean'
#stat='Median'

# Variable
vd='TRW RR'
#vd='BAI RR'
#vd='TRW S NE DPLR RR'
#vd='RGR RR'
#vd='Gsw RR'

# Relative response (%)
be=(dTS[vd]['F'][stat]-dTS[vd]['C'][stat])/dTS[vd]['C'][stat]*100
lo,hi=GetCIsForDR(dTS[vd]['C'][stat]-dTS[vd]['C']['SE'],dTS[vd]['C'][stat]+dTS[vd]['C']['SE'],dTS[vd]['F'][stat]-dTS[vd]['F']['SE'],dTS[vd]['F'][stat]+dTS[vd]['F']['SE'])

# Absolute respones (mm)
#be=dTS[vd]['F'][stat]-dTS[vd]['C'][stat]

# Mean ten-year response
ind=np.where( (binT>=t_ref1[0]) & (binT<=t_ref1[1]) )[0]
mu=np.round(np.mean(be[ind]),decimals=1)

# Time period
it=np.where( (binT>=1965) & (binT<=2020) )[0]

# Plot
plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(12,8)); 
ax.add_patch(Rectangle([1970.5-5,-100],5,250,fc=[0.93,0.95,1],ec="none"))
ax.add_patch(Rectangle([1971.5,-100],10,250,fc=[1,0.94,0.94],ec="none"))
#ax.annotate('Application',xy=(1971.9,-19),xytext=(1971.9,41.5),
#    arrowprops={'color':'red','arrowstyle':'->'},ha='center',color='r');
ax.plot(binT[it],np.zeros(it.size),'k-',lw=1,color='k')
ax.plot(binT[it],be[it],'-ko',ms=3.5,mfc='w',mec='k',lw=1,color='k')
ax.errorbar(binT[it],be[it],yerr=[ be[it]-lo[it], hi[it]-be[it] ],color=[0.8,0.8,0.8],ls='',capsize=0.5,elinewidth=2)
ax.text(1977,32,'$\itI$$_{r,t1-10}$ = ' + str(mu) + '%',fontsize=10,color=[0.7,0.3,0.3])
ax.text(1965.25,-14,'Before',fontsize=9,style='italic',weight='bold',color=[0.05,0.1,0.75]) # $\itt$$_{ref}$
ax.text(1974.5,-14,'After',fontsize=9,style='italic',weight='bold',color=[0.5,0,0]) # $\itt$$_{1-10}$
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=1.5)
ax.set(position=[0.14,0.11,0.83,0.86],xlim=[binT[it][0]-0.5,binT[it][-1]+0.5],xticks=np.arange(1970,2025,10),xlabel='Time, years', \
       ylim=[-20,45],yticks=np.arange(-20,60,5),ylabel='Relative response (%)')

#if flag_savefigs=='On':
#    gu.PrintFig(meta['Paths']['Figures'] + '\\' + vd + '_vs_Time','png',900)
#    fig.savefig(meta['Paths']['Figures'] + '\\' + vd + '_vs_Time.svg',format='svg',dpi=1200)

def FittingCurves():
    
    # Fitting
    xhat=np.arange(1960,2020,1)
    iFit=np.where( (binT>=1972) & (binT<=2020) )[0]
    
    def Func(x,a,b,c):
        y= ( a+(b-a)*np.exp(-c*(x-1971)) )
        return y
    
    flg=0
    if flg==1:
        b=[0.1,3,0.2]
        plt.close('all')
        plt.plot(xhat,Func(xhat,b[0],b[1],b[2]),'b-',lw=1)
    
        poptG,pcovG=curve_fit(Func,binT[iFit],be[iFit],[0.1,3,0.15])
        yhat=Func(xhat,poptG[0],poptG[1],poptG[2]) # ,poptG[3]
        ax.plot(xhat,yhat,'b--',lw=1)
    
    def Func(x,a,b,c):
        y=np.maximum(0,a*(b/(c*(b-c))*(np.exp(-c*(x-1971))-np.exp(-b*(x-1971)))))
        return y
    
    #b=[3,1.05,0.2]
    #plt.close('all')
    #plt.plot(xhat,Func(xhat,b[0],b[1],b[2]),'b-',lw=1)
    poptG,pcovG=curve_fit(Func,binT[iFit],be[iFit],[3,1.05,0.2])
    yhat=Func(xhat,poptG[0],poptG[1],poptG[2]) # ,poptG[3]
    ax.plot(xhat,yhat,'b--',lw=1)
    
    yhat=Func(xhat,0.4,poptG[1],poptG[2]) # ,poptG[3]
    ax.plot(xhat,yhat,'r--',lw=1)
    
    def Func(x,a,b,c):
        y=1-( ( (1-a) * x**b) / ( c + (x)**b ) )
        return y



#%% Plot time series of detrended responses (long time series)

# Central tendancy
stat='Mean'

# Time period
it=np.where( (binT>=1965) & (binT<=2010) )[0]

plt.close('all'); 
fig,ax=plt.subplots(2,1,figsize=gu.cm2inch(12,10)); 

v='TRW RR'
#v='RGR RR'

ax[0].plot(binT[it],np.zeros(it.size),'k-',lw=1,color='k')
be=(dTS[v]['F'][stat][it]-dTS[v]['C'][stat][it])/dTS[v]['C'][stat][it]*100
ax[0].plot(binT[it],be,'-ko',ms=3.5,mfc='w',mec='k',lw=1,color='k',label='Original control')
be=(dTS[v]['F'][stat][it]-dTS[v]['SS3'][stat][it])/dTS[v]['SS3'][stat][it]*100
ax[0].plot(binT[it],be,'-rs',ms=3.5,mfc='w',mec='r',lw=1,color='r',label='Dry control')
be=(dTS[v]['F'][stat][it]-dTS[v]['SS5'][stat][it])/dTS[v]['SS5'][stat][it]*100
ax[0].plot(binT[it],be,'-bd',ms=3.5,mfc='w',mec='b',lw=1,color='b',label='Wet control')
ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].tick_params(length=1.5)
ax[0].set(position=[0.11,0.54,0.87,0.44],xlim=[binT[it][0]-0.5,binT[it][-1]+0.5],xticks=np.arange(1965,2025,5),ylim=[-45,45], \
  yticks=np.arange(-50,100,10),xticklabels=[''],ylabel='$\itI_{r,t}$ (%)')
ax[0].legend(loc='lower left',frameon=False,ncol=3,fontsize=6)

v='TRW S NE DPLR RR'

be=(dTS[v]['F'][stat][it]-dTS[v]['C'][stat][it])/dTS[v]['C'][stat][it]*100
ax[1].add_patch(Rectangle([binT[it[0]]-10,-1*np.std(be)],100,2*np.std(be),fc=[0.875,0.875,0.875],ec="none"))
ax[1].plot(binT[it],np.zeros(it.size),'k-',lw=1,color='k')

ax[1].plot(binT[it],be,'-ko',ms=3.5,mfc='w',mec='k',lw=1,color='k')
be=(dTS[v]['F'][stat][it]-dTS[v]['SS3'][stat][it])/dTS[v]['SS3'][stat][it]*100
ax[1].plot(binT[it],be,'-rs',ms=3.5,mfc='w',mec='r',lw=1,color='r')
be=(dTS[v]['F'][stat][it]-dTS[v]['SS5'][stat][it])/dTS[v]['SS5'][stat][it]*100
ax[1].plot(binT[it],be,'-bd',ms=3.5,mfc='w',mec='b',lw=1,color='b')
ax[1].yaxis.set_ticks_position('both'); ax[1].xaxis.set_ticks_position('both'); ax[1].tick_params(length=1.5)
ax[1].set(position=[0.11,0.08,0.87,0.44],xlim=[binT[it][0]-0.5,binT[it][-1]+0.5],xticks=np.arange(1965,2025,5),ylim=[-45,45], \
  yticks=np.arange(-50,100,10),ylabel='$\itI_{r,t}$ (%)',xlabel='Time, years')

gu.axletters(ax,plt,0.024,0.9,FontWeight='Bold',Labels=['Without detrending','With detrending'],LabelSpacer=0.04) # ,LetterStyle='Caps'

if flag_savefigs=='On':
    gu.PrintFig(meta['Paths']['Figures'] + '\\GrowthIndex_vs_Time2','png',900)

#%% Plot time series (short detrending period)

# Central tendancy
stat='Mean'

# Time period
it=np.where( (binT>=1965) & (binT<=1981) )[0]
itR=np.where( (binT[it]>=1972) & (binT[it]<=1981) )[0]

plt.close('all'); 
fig,ax=plt.subplots(3,1,figsize=gu.cm2inch(8.5,12)); 

v='TRW RR'#v='RGR RR'
ax[0].plot(binT[it],np.zeros(it.size),'k-',lw=1,color='k')
be=(dTS[v]['F'][stat][it]-dTS[v]['C'][stat][it])/dTS[v]['C'][stat][it]*100
ax[0].plot(binT[it],be,'-ko',ms=3.5,mfc='w',mec='k',lw=1,color='k',label='Original control')
be=(dTS[v]['F'][stat][it]-dTS[v]['SS3'][stat][it])/dTS[v]['SS3'][stat][it]*100
ax[0].plot(binT[it],be,'-rs',ms=3.5,mfc='w',mec='r',lw=1,color='r',label='Dry control')
be=(dTS[v]['F'][stat][it]-dTS[v]['SS5'][stat][it])/dTS[v]['SS5'][stat][it]*100
ax[0].plot(binT[it],be,'-bd',ms=3.5,mfc='w',mec='b',lw=1,color='b',label='Wet control')
ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].tick_params(length=1.5)
ax[0].set(position=[0.15,0.7,0.82,0.29],xlim=[binT[it][0]-0.5,binT[it][-1]+0.5],xticks=np.arange(1965,2025,2),ylim=[-30,45], \
  yticks=np.arange(-50,100,10),xticklabels=[''],ylabel='$\itI_{r,t}$ (%)')
ax[0].legend(loc='lower left',frameon=False,ncol=3,fontsize=6)

#------------------------------------------------------------------------------
# Short fit (all 17 years)
rsShort=[]
ax[1].plot(binT[it],np.zeros(it.size),'k-',lw=1,color='k')

y=(dTS[v]['F'][stat][it]-dTS[v]['C'][stat][it])/dTS[v]['C'][stat][it]*100
x=binT[it]
x1=sm.tools.tools.add_constant(x)
md=sm.OLS(y,x1).fit()
md.summary()
xhat=np.linspace(np.min(x1[:,1]),np.max(x1[:,1]),y.shape[0])
yhat=md.predict(np.c_[np.ones(xhat.size),xhat])
yd=y-yhat
rsShort.append(np.mean(yd[itR]))
ax[1].plot(binT[it],yd,'-ko',ms=3.5,mfc='w',mec='k',lw=1,color='k')

y=(dTS[v]['F'][stat][it]-dTS[v]['SS3'][stat][it])/dTS[v]['SS3'][stat][it]*100
x=binT[it]
x1=sm.tools.tools.add_constant(x)
md=sm.OLS(y,x1).fit()
md.summary()
xhat=np.linspace(np.min(x1[:,1]),np.max(x1[:,1]),y.shape[0])
yhat=md.predict(np.c_[np.ones(xhat.size),xhat])
yd=y-yhat
rsShort.append(np.mean(yd[itR]))
ax[1].plot(binT[it],yd,'-rs',ms=3.5,mfc='w',mec='r',lw=1,color='r')

y=(dTS[v]['F'][stat][it]-dTS[v]['SS5'][stat][it])/dTS[v]['SS5'][stat][it]*100
x=binT[it]
x1=sm.tools.tools.add_constant(x)
md=sm.OLS(y,x1).fit()
md.summary()
xhat=np.linspace(np.min(x1[:,1]),np.max(x1[:,1]),y.shape[0])
yhat=md.predict(np.c_[np.ones(xhat.size),xhat])
yd=y-yhat
rsShort.append(np.mean(yd[itR]))
ax[1].plot(binT[it],yd,'-bd',ms=3.5,mfc='w',mec='b',lw=1,color='b')

ax[1].yaxis.set_ticks_position('both'); ax[1].xaxis.set_ticks_position('both'); ax[1].tick_params(length=1.5)
ax[1].set(position=[0.15,0.39,0.82,0.29],xlim=[binT[it][0]-0.5,binT[it][-1]+0.5],xticks=np.arange(1965,2025,2),ylim=[-30,45], \
  yticks=np.arange(-50,100,10),xticklabels=[''],ylabel='$\itI_{r,t}$ (%)')

#------------------------------------------------------------------------------
# Short fit (excluding first 5 years after application)

rsShortE=[]
ax[2].add_patch(Rectangle([1971.5,-100],5,250,fc=[0.85,0.85,0.85],ec="none"))
ax[2].plot(binT[it],np.zeros(it.size),'k-',lw=1,color='k')
it2=np.where( (binT[it]>=1960) & (binT[it]<=1970) | (binT[it]>=1976) & (binT[it]<=1981) )[0]

y=(dTS[v]['F'][stat][it]-dTS[v]['C'][stat][it])/dTS[v]['C'][stat][it]*100
y1=y[it2]
x=binT[it2]
x1=sm.tools.tools.add_constant(x)
md=sm.OLS(y1,x1).fit()
md.summary()
xhat=np.linspace(np.min(x1[:,1]),np.max(x1[:,1]),y.shape[0])
yhat=md.predict(np.c_[np.ones(xhat.size),xhat])
yd=y-yhat
rsShortE.append(np.mean(yd[itR]))
ax[2].plot(binT[it],yd,'-ko',ms=3.5,mfc='w',mec='k',lw=1,color='k')

y=(dTS[v]['F'][stat][it]-dTS[v]['SS3'][stat][it])/dTS[v]['SS3'][stat][it]*100
y1=y[it2]
x=binT[it2]
x1=sm.tools.tools.add_constant(x)
md=sm.OLS(y1,x1).fit()
md.summary()
xhat=np.linspace(np.min(x1[:,1]),np.max(x1[:,1]),y.shape[0])
yhat=md.predict(np.c_[np.ones(xhat.size),xhat])
yd=y-yhat
rsShortE.append(np.mean(yd[itR]))
ax[2].plot(binT[it],yd,'-rs',ms=3.5,mfc='w',mec='r',lw=1,color='r')

y=(dTS[v]['F'][stat][it]-dTS[v]['SS5'][stat][it])/dTS[v]['SS5'][stat][it]*100
y1=y[it2]
x=binT[it2]
x1=sm.tools.tools.add_constant(x)
md=sm.OLS(y1,x1).fit()
md.summary()
xhat=np.linspace(np.min(x1[:,1]),np.max(x1[:,1]),y.shape[0])
yhat=md.predict(np.c_[np.ones(xhat.size),xhat])
yd=y-yhat
rsShortE.append(np.mean(yd[itR]))
ax[2].plot(binT[it],yd,'-bd',ms=3.5,mfc='w',mec='b',lw=1,color='b')

ax[2].yaxis.set_ticks_position('both'); ax[1].xaxis.set_ticks_position('both'); ax[2].tick_params(length=1.5)
ax[2].set(position=[0.15,0.08,0.82,0.29],xlim=[binT[it][0]-0.5,binT[it][-1]+0.5],xticks=np.arange(1965,2025,2),ylim=[-30,45], \
  yticks=np.arange(-50,100,10),ylabel='$\itI_{r,t}$ (%)',xlabel='Time, years')

gu.axletters(ax,plt,0.028,0.9,FontWeight='Bold',Labels=['Without detrending','With detrending (full 17-year training period)','With detrending (fragmented training period)'],LabelSpacer=0.06) # ,LetterStyle='Caps'

if flag_savefigs=='On':
    gu.PrintFig(meta['Paths']['Figures'] + '\\GrowthIndex_vs_Time_ShortDetrendingPeriod','png',900)


#%% Bar chart comparison of response using original and post facto controls

stat='Mean'

sig_mult=2.0

rs={}

plt.close('all'); fig,ax=plt.subplots(2,2,figsize=gu.cm2inch(15,7.5)); 

var='TRW RR'
rs[var]={}
rs[var]['Mean']=np.zeros(3)
rs[var]['SE']=np.zeros(3)

rs[var]['Mean'][0]=np.mean((dBM[var]['F1'][stat]-dBM[var]['C'][stat])/dBM[var]['C'][stat]*100)
rs[var]['SE'][0]=sig_mult*np.std((dBM[var]['F1'][stat]-dBM[var]['C'][stat])/dBM[var]['C'][stat]*100)/np.sqrt(dBM[var]['F1'][stat].size)

rs[var]['Mean'][1]=np.mean((dBM[var]['F1'][stat]-dBM[var]['SS3'][stat])/dBM[var]['SS3'][stat]*100)
rs[var]['SE'][1]=sig_mult*np.std((dBM[var]['F1'][stat]-dBM[var]['SS3'][stat])/dBM[var]['SS3'][stat]*100)/np.sqrt(dBM[var]['F1'][stat].size)

rs[var]['Mean'][2]=np.nanmean((dBM[var]['F1'][stat]-dBM[var]['SS5'][stat])/dBM[var]['SS5'][stat]*100)
rs[var]['SE'][2]=sig_mult*np.nanstd((dBM[var]['F1'][stat]-dBM[var]['SS5'][stat])/dBM[var]['SS5'][stat]*100)/np.sqrt(dBM[var]['F1'][stat].size)

ax[0,0].plot([0,5],[0,0],'k-')
ax[0,0].bar([1,2,3],rs[var]['Mean'],0.75,fc=[0.7,0.7,0.7])
ax[0,0].errorbar([1,2,3],rs[var]['Mean'],yerr=[ rs[var]['SE'], rs[var]['SE'] ],color=[0,0,0],ls='',capsize=3,elinewidth=1)
ax[0,0].set(position=[0.08,0.56,0.4,0.4],ylim=[-5,25],yticks=np.arange(-5,40,5),ylabel='$\itI_{r,t1-10}$ (%)',xlim=[0.5,3.5], \
  xticks=np.arange(1,4,1),xticklabels=[''])
ax[0,0].yaxis.set_ticks_position('both'); ax[0,0].xaxis.set_ticks_position('both'); ax[0,0].tick_params(length=1.5)

var='TRW S NE DPLR RR'
rs[var]={}
rs[var]['Mean']=np.zeros(3)
rs[var]['SE']=np.zeros(3)

rs[var]['Mean'][0]=np.mean((dBM[var]['F1'][stat]-dBM[var]['C'][stat])/dBM[var]['C'][stat]*100)
rs[var]['SE'][0]=sig_mult*np.std((dBM[var]['F1'][stat]-dBM[var]['C'][stat])/dBM[var]['C'][stat]*100)/np.sqrt(dBM[var]['F1'][stat].size)

rs[var]['Mean'][1]=np.mean((dBM[var]['F1'][stat]-dBM[var]['SS3'][stat])/dBM[var]['SS3'][stat]*100)
rs[var]['SE'][1]=sig_mult*np.std((dBM[var]['F1'][stat]-dBM[var]['SS3'][stat])/dBM[var]['SS3'][stat]*100)/np.sqrt(dBM[var]['F1'][stat].size)

rs[var]['Mean'][2]=np.nanmean((dBM[var]['F1'][stat]-dBM[var]['SS5'][stat])/dBM[var]['SS5'][stat]*100)
rs[var]['SE'][2]=sig_mult*np.nanstd((dBM[var]['F1'][stat]-dBM[var]['SS5'][stat])/dBM[var]['SS5'][stat]*100)/np.sqrt(dBM[var]['F1'][stat].size)

ax[0,1].plot([0,5],[0,0],'k-')
ax[0,1].bar([1,2,3],rs[var]['Mean'],0.75,fc=[0.7,0.7,0.7])
ax[0,1].errorbar([1,2,3],rs[var]['Mean'],yerr=[ rs[var]['SE'], rs[var]['SE'] ],color=[0,0,0],ls='',capsize=3,elinewidth=1)
ax[0,1].set(position=[0.58,0.58,0.4,0.4],ylim=[-5,25],yticks=np.arange(-5,40,5),ylabel='$\itI_{r,t1-10}$ (%)', \
  xlim=[0.5,3.5],xticks=np.arange(1,4,1),xticklabels=[''])
ax[0,1].yaxis.set_ticks_position('both'); ax[0,1].xaxis.set_ticks_position('both'); ax[0,1].tick_params(length=1.5)

#------------------------------------------------------------------------------

ax[1,0].plot([0,5],[0,0],'k-')
ax[1,0].bar([1,2,3],rsShort,0.75,fc=[0.7,0.7,0.7])
#ax[1,0].errorbar([1,2,3],rs[var]['Mean'],yerr=[ rs[var]['SE'], rs[var]['SE'] ],color=[0,0,0],ls='',capsize=3,elinewidth=1)
ax[1,0].set(position=[0.08,0.09,0.4,0.4],ylim=[-5,25],yticks=np.arange(-5,40,5),ylabel='$\itI_{r,t1-10}$ (%)', \
  xlim=[0.5,3.5],xticks=np.arange(1,4,1),xticklabels=['Original\ncontrol','Dry\ncontrol','Wet\ncontrol'])
ax[1,0].yaxis.set_ticks_position('both'); ax[1,0].xaxis.set_ticks_position('both'); ax[1,0].tick_params(length=1.5)

ax[1,1].plot([0,5],[0,0],'k-')
ax[1,1].bar([1,2,3],rsShortE,0.75,fc=[0.7,0.7,0.7])
#ax[1,0].errorbar([1,2,3],rs[var]['Mean'],yerr=[ rs[var]['SE'], rs[var]['SE'] ],color=[0,0,0],ls='',capsize=3,elinewidth=1)
ax[1,1].set(position=[0.58,0.09,0.4,0.4],ylim=[-5,25],yticks=np.arange(-5,40,5),ylabel='$\itI_{r,t1-10}$ (%)', \
  xlim=[0.5,3.5],xticks=np.arange(1,4,1),xticklabels=['Original\ncontrol','Dry\ncontrol','Wet\ncontrol'])
ax[1,1].yaxis.set_ticks_position('both'); ax[1,1].xaxis.set_ticks_position('both'); ax[1,1].tick_params(length=1.5)


gu.axletters(ax,plt,0.025,0.87,FontWeight='Bold',Labels=['Without detrending','Detrended (full time series)','Detrended (short)','Detrend (short excluding first 5 years)'],LabelSpacer=0.05) # ,LetterStyle='Caps'


if flag_savefigs=='On':
    gu.PrintFig(meta['Paths']['Figures'] + '\\Evaluation of Detrending Bar Chart','png',900)
    #fig.savefig(meta['Paths']['Figures'] + '\\Pseudocontrols test.svg',format='svg',dpi=1200)




#%% Compare core responses with responses from the original trial

dEP={}
dEP['G_BA_DA']=np.zeros(uInst.size)
dEP['G_BA_DR']=np.zeros(uInst.size)
for iInst in range(uInst.size):
    
    # Index to plots of ith installation
    indInst0=np.where( (dTL['ID_Inst']==uInst[iInst]) & (dTL['ID_Plot']<100) )[0]
    
    #uPlot=np.unique(dTL['ID_Plot'][indInst0])
    
    iC_t0=np.where((dEPTL['ID_Spc_t0']==ID_FD_EPTL) & (dEPTL['ID_Site']==uInst[iInst]) & (dEPTL['TSF_t0']==0) & (dEPTL['CrownClass_t0']<=2) & (dEPTL['N_Dose']==0) & (dEPTL['Mort']==0) & (dEPTL['Rec']==0) )[0]
    iC_t1=np.where((dEPTL['ID_Spc_t0']==ID_FD_EPTL) & (dEPTL['ID_Site']==uInst[iInst]) & (dEPTL['TSF_t0']==9) & (dEPTL['CrownClass_t0']<=2) & (dEPTL['N_Dose']==0) & (dEPTL['Mort']==0) & (dEPTL['Rec']==0) )[0]
    G_BA_C=np.sum(dEPTL['BA_t0'][iC_t1])-np.sum(dEPTL['BA_t0'][iC_t0])/9
    
    iF_t0=np.where((dEPTL['ID_Spc_t0']==ID_FD_EPTL) & (dEPTL['ID_Site']==uInst[iInst]) & (dEPTL['TSF_t0']==0) & (dEPTL['CrownClass_t0']<=2) & (dEPTL['N_Dose']==225) & (dEPTL['Mort']==0) & (dEPTL['Rec']==0) )[0]
    iF_t1=np.where((dEPTL['ID_Spc_t0']==ID_FD_EPTL) & (dEPTL['ID_Site']==uInst[iInst]) & (dEPTL['TSF_t0']==9) & (dEPTL['CrownClass_t0']<=2) & (dEPTL['N_Dose']==225) & (dEPTL['Mort']==0) & (dEPTL['Rec']==0) )[0]
    G_BA_F=np.sum(dEPTL['BA_t0'][iF_t1])-np.sum(dEPTL['BA_t0'][iF_t0])/9
    
    dEP['G_BA_DA'][iInst]=G_BA_F-G_BA_C
    dEP['G_BA_DR'][iInst]=(G_BA_F-G_BA_C)/G_BA_C*100

#var='TRW'
var='TRW RR'
#var='BAI'
stat='Mean'

#da_be=dBM[var]['F1'][stat]-dBM[var]['C'][stat]
#da_lo,da_hi=gu.GetCIsFromDifference(dBM[var]['C']['CI Lo'],dBM[var]['C']['CI Hi'],dBM[var]['F1']['CI Lo'],dBM[var]['F1']['CI Hi'])

dr_be=(dBM[var]['F1'][stat]-dBM[var]['C'][stat])/dBM[var]['C'][stat]*100
dr_lo,dr_hi=GetCIsForDR(dBM[var]['C']['CI Lo'],dBM[var]['C']['CI Hi'],dBM[var]['F1']['CI Lo'],dBM[var]['F1']['CI Hi'])

plt.close('all')
plt.plot(dEP['G_BA_DR'],dr_be,'o')

#%% Compare variation in EP703 data
 
dI={}; N=5000; cnt=0
dI['G_DA_T']=np.nan*np.empty(N)
dI['G_DR_T']=np.nan*np.empty(N)
dI['G_DA_C']=np.nan*np.empty(N)
dI['G_DR_C']=np.nan*np.empty(N)
for iInst in range(uInst.size):
    
    # Index to plots of ith installation
    ind0=np.where( (dTL['ID_Inst']==uInst[iInst]) & (dTL['Treatment']=='C') | (dTL['ID_Inst']==uInst[iInst]) & (dTL['Treatment']=='F1') )[0]
    
    uPlot=np.unique(dTL['ID_Plot'][ind0])
    
    for iPlot in range(uPlot.size):
        
        ind1=np.where( (dTL['ID_Inst']==uInst[iInst]) & (dTL['ID_Plot']==uPlot[iPlot]) & (dTL['QA Summary']==1) )[0]
        
        uTree=np.unique(dTL['ID_Tree'][ind1])
        
        for iTree in range(uTree.size):
            
            # Tape
            tsf2=9
            indT0=np.where( (dEPTL['ID_Site']==uInst[iInst]) & (dEPTL['ID_Plot']==uPlot[iPlot]) & (dEPTL['ID_Tree']==uTree[iTree])& (dEPTL['TSF_t0']==0) )[0]
            indT1=np.where( (dEPTL['ID_Site']==uInst[iInst]) & (dEPTL['ID_Plot']==uPlot[iPlot]) & (dEPTL['ID_Tree']==uTree[iTree])& (dEPTL['TSF_t0']==tsf2) )[0]
            yr0=dEPTL['Year_t0'][indT0]
            
            indC0=np.where( (dTL['ID_Inst']==uInst[iInst]) & (dTL['ID_Plot']==uPlot[iPlot]) & (dTL['ID_Tree']==uTree[iTree]) & (dTL['Year']==yr0) )[0]
            indC1=np.where( (dTL['ID_Inst']==uInst[iInst]) & (dTL['ID_Plot']==uPlot[iPlot]) & (dTL['ID_Tree']==uTree[iTree]) & (dTL['Year']==yr0+tsf2-1) )[0] 
            #indC0=np.where( (dGF['ID_Inst']==uInst[iInst]) & (dGF['ID_Plot']==uPlot[iPlot]) & (dGF['ID_Tree']==uTree[iTree]) & (dGF['Year']==yr0) )[0]
            #indC1=np.where( (dGF['ID_Inst']==uInst[iInst]) & (dGF['ID_Plot']==uPlot[iPlot]) & (dGF['ID_Tree']==uTree[iTree]) & (dGF['Year']==yr0+tsf2-1) )[0]
            
            if (indC0.size==0) | (indT0.size==0) | (indC1.size==0) | (indT1.size==0):
                print('Missing')
                continue
            
            dI['G_DA_T'][cnt]=(dEPTL['D_t0'][indT1]-dEPTL['D_t0'][indT0])/tsf2
            dI['G_DR_T'][cnt]=(dEPTL['D_t0'][indT1]-dEPTL['D_t0'][indT0])/dEPTL['D_t0'][indT0]*100/tsf2
            dI['G_DA_C'][cnt]=(dTL['Dib'][indC1]-dTL['Dib'][indC0])/tsf2
            dI['G_DR_C'][cnt]=(dTL['Dib'][indC1]-dTL['Dib'][indC0])/dTL['Dib'][indC0]/tsf2
            #dI['G_DA_C'][cnt]=(dGF['Dib'][indC1]-dGF['Dib'][indC0])/tsf2
            #dI['G_DR_C'][cnt]=(dGF['Dib'][indC1]-dGF['Dib'][indC0])/dGF['Dib'][indC0]/tsf2
            cnt=cnt+1

for k in dI.keys():
    dI[k]=dI[k][0:cnt]

#%%

ikp=np.where( (dI['G_DA_T']!=0) & (dI['G_DA_C']!=0) & (np.isnan(dI['G_DA_C']+dI['G_DA_T'])==False) )[0] 

y=dI['G_DA_C'][ikp]
x=dI['G_DA_T'][ikp]

x1=sm.tools.tools.add_constant(x)
md=sm.OLS(y,x1).fit()
md.summary()
xhat=np.linspace(np.min(x1[:,1]),np.max(x1[:,1]),10)
yhat=md.predict(np.c_[np.ones(xhat.size),xhat])

yhat1=md.predict(x1)
E_mu=np.mean(yhat1-y)
E_mu_R=np.mean( (yhat1-y)/y*100 )
AE_mu=np.mean(np.abs(yhat1-y))
AE_mu_R=np.mean( np.abs((yhat1-y)/y*100) )

# Plot
plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(8.5,8.5))
ax.plot([0,1],[0,1],'k-',lw=2,color=[0.8,0.8,0.8])
ax.plot(x,y,'o',ms=4,mec='w',mfc='k')
ax.plot(xhat,yhat,'r-',lw=1.25)
ax.set(xlim=[0,1],ylim=[0,1],xlabel='Diameter increment from tape (cm/yr)',ylabel='Diameter increment from core (cm/yr)')
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=1.5)
plt.tight_layout()
#gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\EP703\TreeRings\DiameterComparison','png',150)


np.std(dI['G_DA_T'][ikp])
np.std(dI['G_DA_C'][ikp])

np.std(dI['G_DA_T'][ikp])/np.mean(dI['G_DA_T'][ikp])*100
np.std(dI['G_DA_C'][ikp])/np.mean(dI['G_DA_C'][ikp])*100


#%% Age response

def AnalyzeAgeResponse():
    # Choose variable to analyze
    vd='TRW'
    #vd='TRW RR'
    #vd='BAI'
    #vd='Bsw'
    #vd='Gsw RR'
    #vd='RGR RR'
    #vd='TRW S NE DPLR RR'
    
    # Time period of analysis
    bin=np.arange(1,101,1)
    
    # Initialize dictionary of results
    dR={}
    treatL=['C','F','SS3','SS5','SS6']
    varL=['Mean','SE']
    for treat in treatL:
        dR[treat]={}
        for var in varL:
            dR[treat][var]=np.zeros(bin.size)
    
    # Populate dictionary of results
    for i in range(bin.size):
        ind=np.where( (dTL['Age']==bin[i]) & (dTL['Treatment']=='C') & (dTL[vd]>0) & (dTL[vd]<2000) & (dTL['QA Summary']==1) )[0]
        dR['C']['Mean'][i]=np.nanmean(dTL[vd][ind])
        dR['C']['SE'][i]=np.nanstd(dTL[vd][ind])/np.sqrt(ind.size)
        ind=np.where( (dTL['Age']==bin[i]) & (dTL['Treatment']=='F1') & (dTL[vd]>0) & (dTL[vd]<2000) & (dTL['QA Summary']==1) )[0]
        dR['F']['Mean'][i]=np.nanmean(dTL[vd][ind])    
        dR['F']['SE'][i]=np.nanstd(dTL[vd][ind])/np.sqrt(ind.size)
        ind=np.where( (dTL['Age']==bin[i]) & (dTL['BGC SS']==3) & (dTL[vd]>0) & (dTL[vd]<2000) & (dTL['Treatment']!='C') & (dTL['Treatment']!='F1') & (dTL['QA Summary']==1) )[0]
        dR['SS3']['Mean'][i]=np.nanmean(dTL[vd][ind])
        ind=np.where( (dTL['Age']==bin[i]) & (dTL['BGC SS']==5) & (dTL[vd]>0) & (dTL[vd]<2000) & (dTL['Treatment']!='C') & (dTL['Treatment']!='F1') & (dTL['QA Summary']==1) )[0]
        #ind=np.where( (dTL['Year']==bin[i]) & (dTL['BGC SS Comb']==99) & (dTL[vd]>0) & (dTL[vd]<2000) & (dTL['Treatment']!='C') & (dTL['Treatment']!='F1') & (dTL['QA Summary']==1) )[0]
        dR['SS5']['Mean'][i]=np.nanmean(dTL[vd][ind])
        ind=np.where( (dTL['Age']==bin[i]) & (dTL['BGC SS']==6) & (dTL[vd]>0) & (dTL[vd]<2000) & (dTL['Treatment']!='C') & (dTL['Treatment']!='F1') & (dTL['QA Summary']==1) )[0]
        dR['SS6']['Mean'][i]=np.nanmean(dTL[vd][ind])

    
#%%

def fAge_NE(x,a,b,c):
    y=a*np.exp(-b*x)-c*x
    return y

poptG,pcovG=curve_fit(fAge_NE,bin[10:],dR['C']['Mean'][10:],[7,0.04,0.01])
dR['C']['yhat']=fAge_NE(bin,poptG[0],poptG[1],poptG[2]) 

poptG,pcovG=curve_fit(fAge_NE,bin[10:],dR['SS3']['Mean'][10:],[7,0.04,0.01])
dR['SS3']['yhat']=fAge_NE(bin,poptG[0],poptG[1],poptG[2]) 

poptG,pcovG=curve_fit(fAge_NE,bin[10:],dR['SS5']['Mean'][10:],[7,0.04,0.01])
dR['SS5']['yhat']=fAge_NE(bin,poptG[0],poptG[1],poptG[2]) 
  
plt.close('all')
#plt.plot(bin,dR['C']['Mean'],'-ko')
plt.plot(bin,dR['C']['yhat'],'b-',lw=1)
plt.plot(bin,dR['SS3']['yhat'],'r--',lw=1)
plt.plot(bin,dR['SS5']['yhat'],'g--',lw=1)

#%%

d=gu.ReadExcel(r'C:\Users\rhember\Documents\Data\TreeRings\EP703\TIPSY Runs\TipsyRunsDouglasFirWidth.xlsx')

fig,ax=plt.subplots(1,figsize=gu.cm2inch(8.5,8.5))
#plt.plot(d['Abh'],d['Gw all 25'],'.r-')
#plt.plot(d['Abh'],d['Gw all 30'],'b-')
#plt.plot(d['Abh'],d['Gw all 35'],'g-')

poptG,pcovG=curve_fit(fAge_NE,d['Abh'][6:],d['Gw all 25'][6:],[7,0.04,0.01])
yhat=fAge_NE(bin,poptG[0],poptG[1],poptG[2]) 
plt.plot(bin,yhat,'r-')

poptG,pcovG=curve_fit(fAge_NE,d['Abh'][6:],d['Gw all 30'][6:],[7,0.04,0.01])
yhat=fAge_NE(bin,poptG[0],poptG[1],poptG[2]) 
plt.plot(bin,yhat,'b-')

poptG,pcovG=curve_fit(fAge_NE,d['Abh'][6:],d['Gw all 35'][6:],[7,0.04,0.01])
yhat=fAge_NE(bin,poptG[0],poptG[1],poptG[2]) 
plt.plot(bin,yhat,'g-')


#%%

pth=r'E:\Data\ForestInventory\PSP-NADB\Data\02_ProcessingSteps\BC\Level2\FI_PSP_L2_BC-DB.mat'

d=gu.ImportMat(pth,'tobs')
for k in d.keys():
    d[k]=d[k].flatten()/100

ind=np.where( (d['ID_Spc']==25) & (d['EcoZone_BC_L2']>=34) & (d['EcoZone_BC_L2']<=35) )[0]
for k in d.keys():
    d[k]=d[k][ind]

# Fix stand age variable in BC

u=np.unique(d['ID_Plot'])
d['AgeFixed_t0']=d['Age_t0_Stand']
d['AgeFixed_t1']=d['Age_t1_Stand']
for i in range(u.size):
    ind0=np.where( (d['ID_Plot']==u[i]) )[0]
    uT=np.unique(d['ID_Tree'][ind0])
    for j in range(uT.size):
        ind1=np.where( (d['ID_Plot']==u[i]) & (d['ID_Tree']==uT[j]) )[0]
        if ind1.size==1:
            continue
        d['AgeFixed_t0'][ind1]=d['Age_t0_Stand'][ind1[0]]+np.append(0,np.cumsum(d['DT'][ind1][0:-1]))
        d['AgeFixed_t1'][ind1]=d['Age_t0_Stand'][ind1[0]]+np.cumsum(d['DT'][ind1])  

d['AgeFixed']=(d['AgeFixed_t0']+d['AgeFixed_t1'])/2

#%%
    
ind=np.where( (d['Dam_G']>0) & (d['Age_t0_Stand']>0) )[0]

x=d['AgeFixed_t0'][ind]
y=d['Dam_G'][ind]/2*10

xbin=np.arange(5,205,5)
ybin=np.zeros(xbin.size)
for i in range(xbin.size):
    ind=np.where(np.abs(x-xbin[i])<=2.5)[0]
    ybin[i]=np.nanmean(y[ind])

plt.close('all')
plt.plot(x,y,'k.')
plt.plot(xbin,ybin,'ks')


#%% Modelling Age and N application effect

def Func2(x,a,b,c,d,e,f,g):
    fA=a*np.exp(-b*x[:,0])-c*x[:,0]
    #bia=np.minimum(1,np.maximum(0,x-b[5]))
    #fN= ( b[2]+(b[3]-b[2])*np.exp(-b[4]*(x-b[5])) )
    #fN=np.maximum(0,b[2]*b[3]/(b[4]*(b[3]-b[5]))*(np.exp(-b[5]*(x[:,1]-1971))-np.exp(-b[3]*(x[:,1]-1971))))
    fN=np.maximum(0,d*e/(f*(e-g))*(np.exp(-g*(x[:,1]-1971))-np.exp(-e*(x[:,1]-1971))))
    fN=x[:,2]*fN
    y=fA+fN
    return y

plt.close('all')

xhat=np.arange(1,101,1)
b=[8,0.05,0.01]
plt.plot(xhat,fAge_NE(xhat,b[0],b[1],b[2]),'-c',lw=2.5)

xhat2=np.column_stack((xhat,np.arange(1945,1945+xhat.size,1),np.ones(xhat.size)))
b=[8,0.05,0.01, 0.15,0.5,0.2,0.15]
plt.plot(xhat,Func2(xhat2,b[0],b[1],b[2],b[3],b[4],b[5],b[6]),'-bs',lw=1.5)

ind=np.where(dTL['Year']==1971)[0]
np.mean(dTL['Age'][ind])

#%% Plot Sample size

dSS=np.zeros((binT.size,4)); 
for iT in range(binT.size):
    ind=np.where( (dTL['Year']==binT[iT]) & (dTL['Treatment']=='C') & (dTL['TRW']>0) & (dTL['TRW']<2000) & (dTL['QA Summary']==1) )[0]; dSS[iT,0]=ind.size
    ind=np.where( (dTL['Year']==binT[iT]) & (dTL['Treatment']=='F1') & (dTL['TRW']>0) & (dTL['TRW']<2000) & (dTL['QA Summary']==1) )[0]; dSS[iT,1]=ind.size
    ind=np.where( (dTL['Year']==binT[iT]) & (dTL['BGC SS']==3) & (dTL['TRW']>0) & (dTL['TRW']<2000) & (dTL['QA Summary']==1) )[0]; dSS[iT,2]=ind.size
    ind=np.where( (dTL['Year']==binT[iT]) & (dTL['BGC SS']==5) & (dTL['TRW']>0) & (dTL['TRW']<2000) & (dTL['QA Summary']==1) )[0]; dSS[iT,3]=ind.size

plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(12,8)); 
plt.plot(binT,dSS,lw=1)


#%%

ind=np.where( (dTL['BGC SS']==1) & (dTL['H 2020']>0) )[0]
np.mean(dTL['H 2020'][ind])

ind=np.where( (dTL['BGC SS']==3) & (dTL['H 2020']>0) )[0]
np.mean(dTL['H 2020'][ind])

ind=np.where( (dTL['BGC SS']==5) & (dTL['H 2020']>0) )[0]
np.mean(dTL['H 2020'][ind])



##%% Modelling
#
## Add binary indicator of application
#dTL['BIN']=np.zeros(dTL['ID_Tree'].size)
#ind=np.where( (dTL['Treatment']=='F1') & (dTL['Year']>=1972) )[0]
#dTL['BIN'][ind]=1
#
#dTL['TSF']=np.maximum(0,dTL['Year']-1971)
# 
## Prepare dataframe for R
#ikp=np.where( (dTL['Treatment']=='C') & (dTL['QA Summary']==1) & (dTL['TRW']>0) & (dTL['Age']>10) | \
#             (dTL['Treatment']=='F1') & (dTL['QA Summary']==1) & (dTL['TRW']>0) & (dTL['Age']>10) )[0]
#
#d={}
#for k in dTL.keys():
#    d[k]=dTL[k][ikp]
#
#df=pd.DataFrame(d)
#df.to_excel(r'C:\Users\rhember\Documents\Data\TreeRings\EP703\ForR_Complete.xlsx')