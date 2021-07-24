'''
PROVINGICAL FERTILIZATION PROGRAM - UTILITIES
'''

#%% Import modules

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
import gc as garc
import scipy.stats as stats
from fcgadgets.macgyver import utilities_general as gu
from fcgadgets.macgyver import utilities_gis as gis
from fcgadgets.macgyver import utilities_inventory as invu
from fcgadgets.cbrunner import cbrun_utilities as cbu

#%% Percent of areas that have been harvested

def Calc_PercentOfStandsHarvested(meta,dmec):

    y=np.zeros(meta['Project']['N Stand Full'])
    
    for iStand in range(meta['Project']['N Stand Full']):        
    
        if dmec[iStand]==None:
            continue  
    
        dmec[iStand]['Year']=np.floor(dmec[iStand]['Year'])
    
        # Index to fertilization
        iFert=np.where( (dmec[iStand]['ID_Type']==meta['LUT']['Dist']['Fertilization Aerial']) )[0]
    
        if iFert.size==0: 
            continue

        iFert=iFert[-1]
        
        iHarv=np.where( (dmec[iStand]['Year']>=dmec[iStand]['Year'][iFert]) & (dmec[iStand]['ID_Type']==meta['LUT']['Dist']['Harvest']) | (dmec[iStand]['Year']>=dmec[iStand]['Year'][iFert]) & (dmec[iStand]['ID_Type']==meta['LUT']['Dist']['Harvest Salvage']) )[0]
    
        if iHarv.size>0:
            y[iStand]=1
    
    txt='Total number of sampled hectares = ' + str(y.size)
    print(txt)
    
    txt='Number of sampled hectares that were harvested = ' + str(np.sum(y))
    print(txt)   
    
    txt='Percent of stands harvested = ' + str(np.round(np.sum(y)/y.size*100))
    print(txt)

#%% Age at time of fertilization (from VRI)

def Calc_AgeAtTimeOfFert_FromVRI(meta,dmec,vri):
    AgeAtFert=np.nan*np.ones(1000000)
    cnt=0
    cnt_full=0
    for iStand in range(meta['Project']['N Stand Full']):
            
        if dmec[iStand]==None:
            continue    
    
        dmec[iStand]['Year']=np.floor(dmec[iStand]['Year'])
    
        indFert=np.where( (dmec[iStand]['ID_Type']==meta['LUT']['Dist']['Fertilization Aerial']) )[0]
        cnt_full=cnt_full+indFert.size
    
        if indFert.size==0:
            continue    
    
        Age=vri['PROJ_AGE_1'][iStand]
    
        if Age<0:
            continue
    
        for iFert in range(indFert.size):
            
            dYear=vri['REFERENCE_YEAR'][iStand]-dmec[iStand]['Year'][indFert[iFert]]
            AgeAtFert[cnt]=Age-dYear 
            cnt=cnt+1

    #print(cnt)
    #print(cnt/cnt_full*100)

    AgeAtFert=AgeAtFert[0:cnt]
    
    txt='Mean age at time of fert (from VRI) = ' + str(np.round(np.mean(AgeAtFert)))
    print(txt)
    
    txt='Median age at time of fert (from VRI) = ' + str(np.round(np.median(AgeAtFert)))
    print(txt)
    
    return AgeAtFert

#%% Age at time of fertilization (from disturbance databases)

def Calc_AgeAtTimeOfFert_FromDisturbanceDBs(meta,dmec):
    
    AgeAtFert=np.nan*np.ones(1000000)
    cnt=0
    cnt_full=0
    for iStand in range(meta['Project']['N Stand Full']):
        
        if dmec[iStand]==None:
            continue    
    
        dmec[iStand]['Year']=np.floor(dmec[iStand]['Year'])
        
        indFert=np.where( (dmec[iStand]['ID_Type']==meta['LUT']['Dist']['Fertilization Aerial']) )[0]
    
        cnt_full=cnt_full+indFert.size
    
        if indFert.size==0:
            continue
    
        for iFert in range(indFert.size):
    
            iHarv=np.where( (dmec[iStand]['Year']<dmec[iStand]['Year'][indFert[iFert]]) & (dmec[iStand]['ID_Type']==meta['LUT']['Dist']['Harvest']) | (dmec[iStand]['Year']<dmec[iStand]['Year'][indFert[iFert]]) & (dmec[iStand]['ID_Type']==meta['LUT']['Dist']['Harvest Salvage']) )[0]
            
            iDist=np.where( (dmec[iStand]['Year']<dmec[iStand]['Year'][indFert[iFert]]) & (dmec[iStand]['ID_Type']==meta['LUT']['Dist']['Harvest']) | 
                (dmec[iStand]['Year']<dmec[iStand]['Year'][indFert[iFert]]) & (dmec[iStand]['ID_Type']==meta['LUT']['Dist']['Knockdown']) |
                (dmec[iStand]['Year']<dmec[iStand]['Year'][indFert[iFert]]) & (dmec[iStand]['ID_Type']==meta['LUT']['Dist']['Wildfire']) |
                (dmec[iStand]['Year']<dmec[iStand]['Year'][indFert[iFert]]) & (dmec[iStand]['ID_Type']==meta['LUT']['Dist']['IBM']) & (dmec[iStand]['MortalityFactor']>80) |
                (dmec[iStand]['Year']<dmec[iStand]['Year'][indFert[iFert]]) & (dmec[iStand]['ID_Type']==meta['LUT']['Dist']['IBS']) & (dmec[iStand]['MortalityFactor']>80) |
                (dmec[iStand]['Year']<dmec[iStand]['Year'][indFert[iFert]]) & (dmec[iStand]['ID_Type']==meta['LUT']['Dist']['IBB']) & (dmec[iStand]['MortalityFactor']>80) |
                (dmec[iStand]['Year']<dmec[iStand]['Year'][indFert[iFert]]) & (dmec[iStand]['ID_Type']==meta['LUT']['Dist']['IBD']) & (dmec[iStand]['MortalityFactor']>80) )[0]
        
            #iBurn=np.where( (dmec[iStand]['Year']<dmec[iStand]['Year'][indFert[iFert]]) & (dmec[iStand]['ID_Type']==meta['LUT']['Dist']['Harvest']) | (dmec[iStand]['Year']<dmec[iStand]['Year'][indFert[iFert]]) & (dmec[iStand]['ID_Type']==meta['LUT']['Dist']['Harvest Salvage']) )[0]
        
            if iDist.size==0:
                continue
        
            #AgeAtFert[cnt]=dmec[iStand]['Year'][indFert[iFert]]-dmec[iStand]['Year'][iHarv[0]]
            AgeAtFert[cnt]=dmec[iStand]['Year'][indFert[iFert]]-dmec[iStand]['Year'][iDist[0]]
            cnt=cnt+1

    AgeAtFert=AgeAtFert[0:cnt]
    #print(cnt)
    #print(cnt/cnt_full*100)
    txt='Mean age at time of fert (from Disturbance DBs) = ' + str(np.round(np.nanmean(AgeAtFert)))
    print(txt)
    
    txt='Median age at time of fert (from Disturbance DBs) = ' + str(np.round(np.nanmedian(AgeAtFert)))
    print(txt)
    
    return AgeAtFert

#%% Histogram of age at time of fertilization

def Plot_AgeAtFert(AgeAtFert):

    bw=5; 
    bin=np.arange(10,160,bw)
    bA=np.zeros(bin.size)
    for i in range(bin.size):
        ind=np.where( np.abs(AgeAtFert-bin[i])<=bw/2 )[0]
        bA[i]=ind.size
    bA=bA/np.sum(bA)*100

    plt.close('all');
    fig,ax=plt.subplots(1,figsize=gu.cm2inch(12,5));
    ax.bar(bin,bA,bw*0.8,ec='None') #'o',ms=6,mec='k',mfc='w',lw=1)    
    #ax[0].set(position=[0.13,0.13,0.82,0.82],xlim=[-2.5,67.5],xticks=np.arange(0,100,5),xlabel='Stand age, years',ylabel='Frequency (%)')
    ax.set(xlim=[-10,120],xticks=np.arange(0,170,10),xlabel='Stand age at time of application, years',ylabel='Frequency (%)')
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')
    #gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Fertilization\NutrientManagementSummary\HarvestRate','png',500)

    return

#%% Annual probability of harvest vs. time since fertilization

def Calc_AnnProbHarvest_Vs_TimeSinceFert(meta,dmec):
    
    # Time since fertilization
    tsf=np.arange(0,51,1)

    # Occurrence
    Occ=np.zeros((tsf.size,500000))
    cnt=0
    for iStand in range(meta['Project']['N Stand Full']):        
    
        if dmec[iStand]==None:
            cnt=cnt+1
            continue    
    
        dmec[iStand]['Year']=np.floor(dmec[iStand]['Year'])
    
        iFert=np.where( (dmec[iStand]['ID_Type']==meta['LUT']['Dist']['Fertilization Aerial']) )[0]
    
        if iFert.size==0:
            cnt=cnt+1
            continue    
    
        iFert=iFert[-1]
    
        #if (dmec[iStand]['Year'][iA[0]]==np.min(dmec[iStand]['Year'])):
        #    # No previous events
        #    d['No history']=d['No history']+1
    
        # Harvesting
        iHarv=np.where( (dmec[iStand]['Year']>=dmec[iStand]['Year'][iFert]) & (dmec[iStand]['ID_Type']==meta['LUT']['Dist']['Harvest']) | (dmec[iStand]['Year']>=dmec[iStand]['Year'][iFert]) & (dmec[iStand]['ID_Type']==meta['LUT']['Dist']['Harvest Salvage']) )[0]
    
        if iHarv.size==0:
            cnt=cnt+1
            continue
    
        iHarv=iHarv[-1]
    
        #indF=np.where( (dmec[iStand]['Year']>=dmec[iStand]['Year'][iFert]) & (dmec[iStand]['ID_Type']==meta['LUT']['Dist']['Harvest']) | (dmec[iStand]['Year']>=dmec[iStand]['Year'][iFert]) & (dmec[iStand]['ID_Type']==meta['LUT']['Dist']['Harvest Salvage']) )[0]
        #if indF.size>0:
    
        DeltaT=dmec[iStand]['Year'][iHarv]-dmec[iStand]['Year'][iFert]
        iTSF=np.where(tsf==np.round(DeltaT))[0]
        Occ[iTSF,cnt]=1
    
        cnt=cnt+1

    Occ=Occ[:,0:cnt]
    Pa=np.sum(Occ,axis=1)/Occ.shape[1]*100
    
    #print(cnt)
    #print(len(dmec))

    return tsf,Occ,Pa,cnt

#%%
    
def Plot_AnnProbHarvest(tsf,Pa):
    
    plt.close('all');
    fig,ax=plt.subplots(1,figsize=gu.cm2inch(12,6.5));
    ax.bar(tsf,Pa,0.8,ec='None') #'o',ms=6,mec='k',mfc='w',lw=1)

    #plt.close('all');
    #fig,ax=plt.subplots(1,figsize=gu.cm2inch(7.8,6));
    #ax[1].bar(bin,bA,0.8*bw,ec='None',fc=[0.8,0.82,0.8]) #'o',ms=6,mec='k',mfc='w',lw=1)    
    ax.set(xlim=[-0.75,50],xticks=np.arange(-5,100,5),xlabel='Time since application, years',ylabel='Probability of harvest (% yr$^-$$^1$)')
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')

    #gu.axletters(ax,plt,0.03,0.85)
    #gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Fertilization\NutrientManagementSummary\HarvestRate','png',500)
    
    return

#%% Time between harvest and fertilization

def Calc_TimeBtwnFertAndHarvest(meta,dmec):
    
    tbfh=np.nan*np.ones(len(dmec))

    for iStand in range(meta['Project']['N Stand Full']):        
        if dmec[iStand]==None:
            continue
    
        # Fertilization
        iFert=np.where( (dmec[iStand]['ID_Type']==meta['LUT']['Dist']['Fertilization Aerial']) )[0]
    
        if iFert.size==0: 
            continue
    
        iFert=iFert[-1]
    
        # Harvesting
        iHarv=np.where( (dmec[iStand]['Year']>=dmec[iStand]['Year'][iFert]) & (dmec[iStand]['ID_Type']==meta['LUT']['Dist']['Harvest']) )[0]
    
        if iHarv.size==0:
            continue
    
        iHarv=iHarv[-1]
        tbfh[iStand]=dmec[iStand]['Year'][iHarv]-dmec[iStand]['Year'][iFert]
    
    return tbfh

#%% 
    
def Plot_TimeBtwnFertAndHarvest(tbfh):
    
    bw=1; 
    bin=np.arange(0,50,bw)
    y=np.zeros(bin.size)
    for i in range(bin.size):
        ind=np.where(np.round(tbfh)==bin[i])[0]
        y[i]=ind.size
    
    plt.close('all');
    fig,ax=plt.subplots(1,figsize=gu.cm2inch(12,6.5));
    ax.bar(bin,y,0.8,ec='None') #'o',ms=6,mec='k',mfc='w',lw=1)
    ax.set(position=[0.08,0.12,0.84,0.82],xlim=[-0.85,50],xticks=np.arange(0,55,5),xlabel='Time since application, years',ylabel='Num. stands harvested') # (% yr$^-$$^1$)
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')
    
    return

#%% Species composition

def Extract_SpeciesComposition(meta,ba):
    
    cd=np.array([])
    nt=np.array([])
    for i in range(4):
        tmp=ba['VRI_LIVE_STEMS_PER_HA']*(ba['Spc_Pct' + str(i+1)].astype(float)/100)
        ind=np.where(tmp>0)[0]
        cd=np.append(cd,ba['Spc_CD' + str(i+1)][ind].astype(float))    
        nt=np.append(nt,tmp[ind])

    uSpc=np.unique(cd[ (cd>0) & (nt>0) ])
    namSpc=np.array([])
    ntSpc=np.zeros(uSpc.size)
    for i in range(uSpc.size):
        ind=np.where(cd==uSpc[i])[0]
        ntSpc[i]=np.sum(nt[ind])
        namSpc=np.append(namSpc,cbu.lut_n2s(meta['LUT']['VRI']['SPECIES_CD_1'],uSpc[i]))
    
    pSpc=ntSpc/np.sum(ntSpc)*100
    iSort=np.flip(np.argsort(pSpc))
    namSpc=namSpc[iSort]
    pSpc=pSpc[iSort]
    ntSpc=ntSpc[iSort]
    
    return namSpc,ntSpc,pSpc

#%%

def Plot_SpeciesComp(namSpc,pSpc,stop_at):
    
    plt.close('all');
    fig,ax=plt.subplots(1,figsize=gu.cm2inch(15.5,6));
    ax.bar(namSpc[:stop_at],pSpc[:stop_at],0.8,ec='None') #'o',ms=6,mec='k',mfc='w',lw=1)
    ax.set(position=[0.08,0.12,0.84,0.82],ylabel='Percent of sample') # (% yr$^-$$^1$)
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')
    
    return

#%% BGC zone distirbution

def Calc_BGCZone_Composiiton(meta,ba):
    
    ind=np.where(ba['BEC_ZONE_CODE']>0)[0]
    u=np.unique(ba['BEC_ZONE_CODE'][ind])
    
    n=np.zeros(u.size)
    nam=np.array([])
    for i in range(u.size):
        ind=np.where(ba['BEC_ZONE_CODE']==u[i])[0]
        n[i]=ind.size
        nam=np.append(nam,cbu.lut_n2s(meta['LUT']['VRI']['BEC_ZONE_CODE'],u[i]))
    
    iSort=np.flip(np.argsort(n))
    n=n[iSort]
    nam=nam[iSort]
    
    return nam,n

#%% 

def Plot_BGC_Zone_Composition(nam,n):
    
    plt.close('all');
    fig,ax=plt.subplots(1,figsize=gu.cm2inch(15.5,6));
    ax.bar(nam,n,0.8,ec='None') #'o',ms=6,mec='k',mfc='w',lw=1)
    ax.set(position=[0.08,0.12,0.84,0.82],ylabel='Percent of sample') # (% yr$^-$$^1$)
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')
    
    return

#%% Site index

def Calc_SiteIndex(meta,ba):
    
    ind=np.where( (ba['SI']>0) & (ba['FIZ']==1) )[0]
    txt='Mean site index (coast) = ' + str(np.round(np.nanmean(ba['SI'][ind])))
    print(txt)
    
    ind=np.where( (ba['SI']>0) & (ba['FIZ']==2) )[0]
    txt='Mean site index (interior) = ' + str(np.round(np.nanmean(ba['SI'][ind])))
    print(txt)
    
    ind=np.where( (ba['SI']>0) )[0]
    txt='Mean site index (all) = ' + str(np.round(np.nanmean(ba['SI'][ind])))
    print(txt)
    
    return

#%% Identify what distrubances precede fertilization 

def Calc_PrecedingDisturbance(meta,dmec):

    d={}
    d['No history']=0
    d['Harvest']=0
    d['Wildfire, stand replacing']=0
    d['Wildfire, partial']=0
    d['Insects, stand replacing']=0
    d['Insects, partial']=0

    AgeAtFert=np.nan*np.ones(1e6)
    cnt=0

    for iStand in range(meta['Project']['N Stand Full']):
        
        if dmec[iStand]==None:
            continue    
    
        dmec[iStand]['Year']=np.floor(dmec[iStand]['Year'])
        
        indFert=np.where( (dmec[iStand]['ID_Type']==meta['LUT']['Dist']['Fertilization Aerial']) )[0]
    
        if indFert.size==0: 
            continue    
    
        for iFert in range(indFert.size):
    
            if (dmec[iStand]['Year'][indFert[iFert]]==np.min(dmec[iStand]['Year'])):
                # No previous events
                d['No history']=d['No history']+1
    
        # Harvesting
        ind=np.where( (dmec[iStand]['Year']<=dmec[iStand]['Year'][iA[0]]) & (dmec[iStand]['MortalityFactor']==100) & np.isin(dmec[iStand]['ID_Type'],[3,7,12,19,20,21,22,23,25]) )[0]
        if ind.size>0:
            if (dmec[iStand]['ID_Type'][ind[-1]]==meta['LUT']['Dist']['Harvest']):
                d['Harvest']=d['Harvest']+1
                AgeAtFert[iStand]=dmec[iStand]['Year'][iA[0]]-dmec[iStand]['Year'][ind[-1]]
            elif (dmec[iStand]['ID_Type'][ind[-1]]==meta['LUT']['Dist']['Wildfire']):
                d['Wildfire, stand replacing']=d['Wildfire, stand replacing']+1  
            else:
                d['Insects, stand replacing']=d['Insects, stand replacing']+1
        else:
            ind=np.where( (dmec[iStand]['Year']<=dmec[iStand]['Year'][iA[0]]) & (dmec[iStand]['MortalityFactor']<100) & np.isin(dmec[iStand]['ID_Type'],[3,7,12,19,20,21,22,23,25]) )[0]
            if ind.size>0:
                if (dmec[iStand]['ID_Type'][ind[-1]]==meta['LUT']['Dist']['Wildfire']):
                    d['Wildfire, partial']=d['Wildfire, partial']+1  
                elif (np.isin(dmec[iStand]['ID_Type'][ind[-1]],[12,13,14,19,20,21,22,23])):
                    d['Insects, partial']=d['Insects, partial']+1

    return d