'''

FOREST NUTRIENT APPLICATION DATABASE - UTILITIES

'''



#%% Import modules

import os
import numpy as np
import pandas as pd
import geopandas as gpd
import fiona
import matplotlib.pyplot as plt
import gc
from fcgadgets.macgyver import utilities_general as gu
from fcgadgets.macgyver import utilities_gis as gis
from sklearn.utils import resample
import statsmodels.api as sm

#%% Set figure properties

#params={'font.sans-serif':'Calibri',
#        'font.size':11,
#        'axes.edgecolor':'black',
#        'axes.labelsize':11,
#        'axes.labelcolor':'black',
#        'axes.titlesize':7,
#        'axes.linewidth':0.5,        
#        'lines.linewidth':0.5,
#        'text.color':'black',
#        'xtick.color':'black',        
#        'xtick.labelsize':11,
#        'xtick.major.width':0.5,
#        'xtick.major.size':3,
#        'xtick.direction':'in',
#        'ytick.color':'black',
#        'ytick.labelsize':11,
#        'ytick.major.width':0.5,
#        'ytick.major.size':3,
#        'ytick.direction':'in',
#        'legend.fontsize':9,        
#        'savefig.dpi':300,
#        'savefig.transparent':True}
#plt.rcParams.update(params)

#%% Import database

def ImportDB(path):
    
    #path=r'C:\Users\rhember\Documents\Data\Nutrient Addition Experiments Database\Nutrient Addition Experiments Database.xlsx'
    d=gu.ReadExcel(path,'Table1')

    # Add derived variables
    d['N Total']=d['N Num Apps']*d['N Dose Per App (kgN/ha)']

    return d

#%% Plot stemwood net growth response of coastal Douglas-fir

def Plot_StemwoodNetGrowth_FDC_1(d):

    plt.close('all');
    fig,ax=plt.subplots(1,figsize=gu.cm2inch(20,11)); 
    ms=6; xmx=7; 
    ax.plot([0,xmx],[0,xmx],'k-',color=[0.82,0.82,0.82],linewidth=2)

    # Add Shawnigan Lake - 1 app
    ikp=np.where( (d['Species leading']=='Pseudotsuga menziesii var. menziesii') &
             (d['N Num Apps']==1) & 
             (d['Experiment Name']=='Shawnigan Lake') &
             (d['Stemwood Biomass Growth Net Units Combined']=='MgC/ha/yr') & 
             (d['Stemwood Biomass Growth Net Control Combined']>-999) & 
             (d['Stemwood Biomass Growth Net Control Combined']<30) & 
             (d['Duration (years)']==15) &
             (d['P Num Apps']>=0) &
             (d['N Total Mean Annual (kgN/ha/yr)']>=-999) & 
             (d['Exclusion Reason']=='None') )[0]

    x0=d['Stemwood Biomass Growth Net Control Combined'][ikp]
    y0=d['Stemwood Biomass Growth Net Fertilized Combined'][ikp]
    Da=np.round(np.mean(y0-x0),decimals=2)
    Dr=int(np.median((y0-x0)/x0*100))
    ax.plot(x0,y0,'ko',markersize=ms,markerfacecolor=[0.65,0.85,1],markeredgecolor='k',linewidth=0.25,label='Shawnigan Lake, year 15, 1 application')

    # Add Shawnigan Lake -2 app
    ikp=np.where( (d['Species leading']=='Pseudotsuga menziesii var. menziesii') &
             (d['N Num Apps']==2) & 
             (d['Experiment Name']=='Shawnigan Lake') &
             (d['Stemwood Biomass Growth Net Units Combined']=='MgC/ha/yr') & 
             (d['Stemwood Biomass Growth Net Control Combined']>-999) & 
             (d['Stemwood Biomass Growth Net Control Combined']<30) & 
             (d['Duration (years)']==15) &
             (d['P Num Apps']>=0) &
             (d['N Total Mean Annual (kgN/ha/yr)']>=-999) & 
             (d['Exclusion Reason']=='None') )[0]

    x0=d['Stemwood Biomass Growth Net Control Combined'][ikp]
    y0=d['Stemwood Biomass Growth Net Fertilized Combined'][ikp]
    Da=np.round(np.mean(y0-x0),decimals=2)
    Dr=int(np.median((y0-x0)/x0*100))
    ax.plot(x0,y0,'ks',markersize=ms,markerfacecolor=[0.65,0.85,1],markeredgecolor='k',linewidth=0.25,label='Shawnigan Lake, year 15, 2 applications')

    # Add Shawnigan Lake - 3 app
    ikp=np.where( (d['Species leading']=='Pseudotsuga menziesii var. menziesii') &
             (d['N Num Apps']==3) & 
             (d['Experiment Name']=='Shawnigan Lake') &
             (d['Stemwood Biomass Growth Net Units Combined']=='MgC/ha/yr') & 
             (d['Stemwood Biomass Growth Net Control Combined']>-999) & 
             (d['Stemwood Biomass Growth Net Control Combined']<30) & 
             (d['Duration (years)']==15) &
             (d['P Num Apps']>=0) &
             (d['N Total Mean Annual (kgN/ha/yr)']>=-999) & 
             (d['Exclusion Reason']=='None') )[0]

    x0=d['Stemwood Biomass Growth Net Control Combined'][ikp]
    y0=d['Stemwood Biomass Growth Net Fertilized Combined'][ikp]
    Da=np.round(np.mean(y0-x0),decimals=2)
    Dr=int(np.median((y0-x0)/x0*100))
    ax.plot(x0,y0,'kd',markersize=ms,markerfacecolor=[0.65,0.85,1],markeredgecolor='k',linewidth=0.25,label='Shawnigan Lake, year 15, 3 applications')

    # Add Shawnigan Lake - 1 app
    ikp=np.where( (d['Species leading']=='Pseudotsuga menziesii var. menziesii') &
             (d['N Num Apps']==1) & 
             (d['Experiment Name']=='Shawnigan Lake') &
             (d['Stemwood Biomass Growth Net Units Combined']=='MgC/ha/yr') & 
             (d['Stemwood Biomass Growth Net Control Combined']>-999) & 
             (d['Stemwood Biomass Growth Net Control Combined']<30) & 
             (d['Duration (years)']==40) &
             (d['P Num Apps']>=0) &
             (d['N Total Mean Annual (kgN/ha/yr)']>=-999) & 
             (d['Exclusion Reason']=='None') )[0]

    x0=d['Stemwood Biomass Growth Net Control Combined'][ikp]
    y0=d['Stemwood Biomass Growth Net Fertilized Combined'][ikp]
    Da=np.round(np.mean(y0-x0),decimals=2)
    Dr=int(np.median((y0-x0)/x0*100))
    ax.plot(x0,y0,'ko',markersize=ms,markerfacecolor=[0.35,0.45,1],markeredgecolor='k',linewidth=0.25,label='Shawnigan Lake, year 40, 1 application')

    # Add Shawnigan Lake -2 app
    ikp=np.where( (d['Species leading']=='Pseudotsuga menziesii var. menziesii') &
             (d['N Num Apps']==2) & 
             (d['Experiment Name']=='Shawnigan Lake') &
             (d['Stemwood Biomass Growth Net Units Combined']=='MgC/ha/yr') & 
             (d['Stemwood Biomass Growth Net Control Combined']>-999) & 
             (d['Stemwood Biomass Growth Net Control Combined']<30) & 
             (d['Duration (years)']==40) &
             (d['P Num Apps']>=0) &
             (d['N Total Mean Annual (kgN/ha/yr)']>=-999) & 
             (d['Exclusion Reason']=='None') )[0]

    x0=d['Stemwood Biomass Growth Net Control Combined'][ikp]
    y0=d['Stemwood Biomass Growth Net Fertilized Combined'][ikp]
    Da=np.round(np.mean(y0-x0),decimals=2)
    Dr=int(np.median((y0-x0)/x0*100))
    ax.plot(x0,y0,'ks',markersize=ms,markerfacecolor=[0.35,0.45,1],markeredgecolor='k',linewidth=0.25,label='Shawnigan Lake, year 40, 2 applications')

    # Add Shawnigan Lake - 3 app
    ikp=np.where( (d['Species leading']=='Pseudotsuga menziesii var. menziesii') &
             (d['N Num Apps']==3) & 
             (d['Experiment Name']=='Shawnigan Lake') &
             (d['Stemwood Biomass Growth Net Units Combined']=='MgC/ha/yr') & 
             (d['Stemwood Biomass Growth Net Control Combined']>-999) & 
             (d['Stemwood Biomass Growth Net Control Combined']<30) & 
             (d['Duration (years)']==40) &
             (d['P Num Apps']>=0) &
             (d['N Total Mean Annual (kgN/ha/yr)']>=-999) & 
             (d['Exclusion Reason']=='None') )[0]

    x0=d['Stemwood Biomass Growth Net Control Combined'][ikp]
    y0=d['Stemwood Biomass Growth Net Fertilized Combined'][ikp]
    Da=np.round(np.mean(y0-x0),decimals=2)
    Dr=int(np.median((y0-x0)/x0*100))
    ax.plot(x0,y0,'kd',markersize=ms,markerfacecolor=[0.35,0.45,1],markeredgecolor='k',linewidth=0.25,label='Shawnigan Lake, year 40, 3 applications')


    # Add EP703
    ikp=np.where( (d['Species leading']=='Pseudotsuga menziesii var. menziesii') &
             (d['Experiment Name']=='EP703') &
             (d['Stemwood Biomass Growth Net Units Combined']=='MgC/ha/yr') & 
             (d['Stemwood Biomass Growth Net Control Combined']>-999) & 
             (d['Stemwood Biomass Growth Net Control Combined']<30) & 
             (d['P Num Apps']>=0) &
             (d['Duration (years)']==3) &
             (d['N Total Mean Annual (kgN/ha/yr)']>=-999) & 
             (d['Exclusion Reason']=='None') )[0]

    x0=d['Stemwood Biomass Growth Net Control Combined'][ikp]
    y0=d['Stemwood Biomass Growth Net Fertilized Combined'][ikp]
    Da=np.round(np.mean(y0-x0),decimals=2)
    Dr=int(np.median((y0-x0)/x0*100))
    ax.plot(x0,y0,'ko',markersize=ms,markerfacecolor=[0.5,0.9,0.25],markeredgecolor='k',linewidth=0.25,label='EP703, year 3')

    ikp=np.where( (d['Species leading']=='Pseudotsuga menziesii var. menziesii') &
             (d['Experiment Name']=='EP703') &
             (d['Duration (years)']==12) &
             (d['Stemwood Biomass Growth Net Units Combined']=='MgC/ha/yr') & 
             (d['Stemwood Biomass Growth Net Control Combined']>-999) & 
             (d['Stemwood Biomass Growth Net Control Combined']<30) & 
             (d['Duration (years)']>=0) &
             (d['Duration (years)']<=120) &
             (d['Thinned']!='NoNo') &
             (d['P Num Apps']>=0) &
             (d['N Total Mean Annual (kgN/ha/yr)']>=-999) & 
             (d['Exclusion Reason']=='None') )[0]

    x0=d['Stemwood Biomass Growth Net Control Combined'][ikp]
    y0=d['Stemwood Biomass Growth Net Fertilized Combined'][ikp]
    Da=np.round(np.mean(y0-x0),decimals=2)
    Dr=int(np.median((y0-x0)/x0*100))
    ax.plot(x0,y0,'k^',markersize=ms,markerfacecolor=[0.5,0.9,0.25],markeredgecolor='k',linewidth=0.25,label='EP703, year 12')

    # Add DF49
    ikp=np.where( (d['Species leading']=='Pseudotsuga menziesii var. menziesii') &
             (d['Experiment Name']=='DF49') &
             (d['Stemwood Biomass Growth Net Units Combined']=='MgC/ha/yr') & 
             (d['Stemwood Biomass Growth Net Control Combined']>-999) & 
             (d['Stemwood Biomass Growth Net Control Combined']<30) & 
             (d['Duration (years)']>=0) &
             (d['Duration (years)']<=120) &
             (d['Thinned']!='NoNo') &
             (d['P Num Apps']>=0) &
             (d['N Total Mean Annual (kgN/ha/yr)']>=-999) & 
             (d['Exclusion Reason']=='None') )[0]

    x0=d['Stemwood Biomass Growth Net Control Combined'][ikp]
    y0=d['Stemwood Biomass Growth Net Fertilized Combined'][ikp]
    Da=np.round(np.mean(y0-x0),decimals=2)
    Dr=int(np.median((y0-x0)/x0*100))
    ax.plot(x0,y0,'ko',markersize=ms,markerfacecolor=[1,0.65,0.85],markeredgecolor='k',linewidth=0.25,label='DF49, year 5')

    # Add Cedar River
    ikp=np.where( (d['Species leading']=='Pseudotsuga menziesii var. menziesii') &
             (d['Site Name']=='Cedar River') &
             (d['Stemwood Biomass Growth Net Units Combined']=='MgC/ha/yr') & 
             (d['Stemwood Biomass Growth Net Control Combined']>-999) & 
             (d['Stemwood Biomass Growth Net Control Combined']<30) & 
             (d['Duration (years)']>=0) &
             (d['Duration (years)']<=120) &
             (d['Thinned']!='NoNo') &
             (d['P Num Apps']>=0) &
             (d['N Total Mean Annual (kgN/ha/yr)']>=-999) & 
             (d['Exclusion Reason']=='None') )[0]

    x0=d['Stemwood Biomass Growth Net Control Combined'][ikp]
    y0=d['Stemwood Biomass Growth Net Fertilized Combined'][ikp]
    Da=np.round(np.mean(y0-x0),decimals=2)
    Dr=int(np.median((y0-x0)/x0*100))
    ax.plot(x0,y0,'kd',markersize=ms,markerfacecolor=[1,0.75,0.25],markeredgecolor='k',linewidth=0.25,label='Cedar River, year 24')

    # Add Wind River
    ikp=np.where( (d['Species leading']=='Pseudotsuga menziesii var. menziesii') &
             (d['Site Name']=='Wind River') &
             (d['Stemwood Biomass Growth Net Units Combined']=='MgC/ha/yr') & 
             (d['Stemwood Biomass Growth Net Control Combined']>-999) & 
             (d['Stemwood Biomass Growth Net Control Combined']<30) & 
             (d['Duration (years)']>=0) &
             (d['Duration (years)']<=120) &
             (d['Thinned']!='NoNo') &
             (d['P Num Apps']>=0) &
             (d['N Total Mean Annual (kgN/ha/yr)']>=-999) & 
             (d['Exclusion Reason']=='None') )[0]

    x0=d['Stemwood Biomass Growth Net Control Combined'][ikp]
    y0=d['Stemwood Biomass Growth Net Fertilized Combined'][ikp]
    Da=np.round(np.mean(y0-x0),decimals=2)
    Dr=int(np.median((y0-x0)/x0*100))
    ax.plot(x0,y0,'ko',markersize=ms,markerfacecolor=[0.8,0,0],markeredgecolor='k',linewidth=0.25,label='Wind River, year 15')

    ax.set(position=[0.08,0.11,0.55,0.87],xlim=[0,xmx],ylim=[0,xmx],xlabel='Net growth of stemwood biomass, control (MgC ha$^-$$^1$ yr$^-$$^1$)',ylabel='Net growth of stemwood biomass, fertilized (MgC ha$^-$$^1$ yr$^-$$^1$)')
    ax.legend(title='Legend',bbox_to_anchor=(1.05, 1),loc='upper left',frameon=False)
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')

    # All trials for global regression
    ikp=np.where( (d['Species leading']=='Pseudotsuga menziesii var. menziesii') &
             (d['N Num Apps']>=0) &
             (d['Stemwood Biomass Growth Net Units Combined']=='MgC/ha/yr') & 
             (d['Stemwood Biomass Growth Net Control Combined']>-999) & 
             (d['Stemwood Biomass Growth Net Control Combined']<30) & 
             (d['Duration (years)']>=0) &
             (d['Duration (years)']<=120) &
             (d['Thinned']!='NoNo') &
             (d['P Num Apps']>=0) &
             (d['N Total Mean Annual (kgN/ha/yr)']>=-999) & 
             (d['Exclusion Reason']=='None') )[0]

    x0=d['Stemwood Biomass Growth Net Control Combined'][ikp]
    y0=d['Stemwood Biomass Growth Net Fertilized Combined'][ikp]
    Da=np.round(np.mean(y0-x0),decimals=2)
    Dr=int(np.median((y0-x0)/x0*100))

    # With intercept
    x1=sm.tools.tools.add_constant(x0)
    md=sm.OLS(y0,x1).fit()
    xhat=np.arange(0,7.1,0.1)
    yhat=md.predict(np.c_[np.ones(xhat.size),xhat])
    plt.plot(xhat,yhat,'k--',linewidth=1)
    md.summary()

    ax.text(0.95*xmx,0.1*xmx,'Global relationship:\n' + 
        'y = ' + str(np.round(md.params[1],decimals=2)) + 'x + ' + str(np.round(md.params[0],decimals=2)) + '\n' +
        'Mean difference = ' + str(Da) + ' (MgC ha$^-$$^1$ yr$^-$$^1$)\n' +
        'Median relative difference = ' + str(Dr) + ' (%)',ha='right')

    return

#%% Plot second summary for coastal Douglas-fir

def Plot_StemwoodNetGrowth_FDC_2(d,units):

    plt.close('all');
    fig,ax=plt.subplots(1,figsize=gu.cm2inch(16,14)); 
    
    fs2=11
    dy=25
    
    # EP703
    ikp=np.where( (d['Species leading']=='Pseudotsuga menziesii var. menziesii') &
             (d['Experiment Name']=='EP703') &
             (d['Stemwood Biomass Growth Net Units Combined']=='MgC/ha/yr') & 
             (d['Stemwood Biomass Growth Net Control Combined']>-999) & 
             (d['Thinned']!='NoNo') &
             (d['Exclusion Reason']=='None') )[0]

    for i in range(ikp.size):
        if d['N Num Apps'][ikp[i]]==1:
            smb='o'
        elif d['N Num Apps'][ikp[i]]==2:
            smb='s'
        elif d['N Num Apps'][ikp[i]]==3:
            smb='d'    
        if i==0:
            ax.plot(d['Duration (years)'][ikp[i]],d['N Total'][ikp[i]],smb,color=[0.4,0.7,0],markeredgecolor='k',linewidth=0.25,label='EP703, 1 application')
        else:
            ax.plot(d['Duration (years)'][ikp[i]],d['N Total'][ikp[i]],smb,color=[0.4,0.7,0],markeredgecolor='k',linewidth=0.25)
        
        if units=='Actual':
            y=str(np.round(d['Stemwood Biomass Growth Net DA Combined'][ikp[i]],decimals=2))
        else:
            y=str(int(d['Stemwood Biomass Growth Net DR Combined (%)'][ikp[i]]))
        ax.text(d['Duration (years)'][ikp[i]],d['N Total'][ikp[i]]+dy,y,ha='center',va='center',color=[0.4,0.7,0],fontsize=fs2)

    # Shawnigan Lake
    ikp=np.where( (d['Species leading']=='Pseudotsuga menziesii var. menziesii') &
             (d['Experiment Name']=='Shawnigan Lake') &
             (d['Stemwood Biomass Growth Net Units Combined']=='MgC/ha/yr') & 
             (d['Stemwood Biomass Growth Net Control Combined']>-999) & 
             (d['Thinned']!='NoNo') &
             (d['Exclusion Reason']=='None') )[0]

    flg1=0; flg2=0; flg3=0;
    for i in range(ikp.size): 
        if d['N Num Apps'][ikp[i]]==1:
            smb='o'
        elif d['N Num Apps'][ikp[i]]==2:
            smb='s'
        elif d['N Num Apps'][ikp[i]]==3:
            smb='d' 
        
        if (smb=='o') & (flg1==0):
            ax.plot(d['Duration (years)'][ikp[i]],d['N Total'][ikp[i]],smb,color=[0.4,0.6,0.85],markeredgecolor='k',linewidth=0.25,label='Shawnigan Lake, 1 application')
            flg1=1
        elif (smb=='s') & (flg2==0):
            ax.plot(d['Duration (years)'][ikp[i]],d['N Total'][ikp[i]],smb,color=[0.4,0.6,0.85],markeredgecolor='k',linewidth=0.25,label='Shawnigan Lake, 2 applications')
            flg2=1
        elif (smb=='d') & (flg3==0):
            ax.plot(d['Duration (years)'][ikp[i]],d['N Total'][ikp[i]],smb,color=[0.4,0.6,0.85],markeredgecolor='k',linewidth=0.25,label='Shawnigan Lake, 3 applications')    
            flg3=1
        else:
            ax.plot(d['Duration (years)'][ikp[i]],d['N Total'][ikp[i]],smb,color=[0.4,0.6,0.85],markeredgecolor='k',linewidth=0.25)
        
        if units=='Actual':
            y=str(np.round(d['Stemwood Biomass Growth Net DA Combined'][ikp[i]],decimals=2))
        else:
            y=str(int(d['Stemwood Biomass Growth Net DR Combined (%)'][ikp[i]]))
        ax.text(d['Duration (years)'][ikp[i]],d['N Total'][ikp[i]]-dy,y,ha='center',va='center',color=[0.4,0.6,0.85],fontsize=fs2)

    # DF49
    ikp=np.where( (d['Species leading']=='Pseudotsuga menziesii var. menziesii') &
             (d['Experiment Name']=='DF49') &
             (d['Stemwood Biomass Growth Net Units Combined']=='MgC/ha/yr') & 
             (d['Stemwood Biomass Growth Net Control Combined']>-999) & 
             (d['Thinned']!='NoNo') &
             (d['Exclusion Reason']=='None') )[0]

    for i in range(ikp.size): 
        if d['N Num Apps'][ikp[i]]==1:
            smb='o'
        elif d['N Num Apps'][ikp[i]]==2:
            smb='s'
        elif d['N Num Apps'][ikp[i]]==3:
            smb='d' 
        ax.plot(d['Duration (years)'][ikp[i]],d['N Total'][ikp[i]],smb,color=[0.7,0.2,1],markeredgecolor='k',linewidth=0.25,label='DF49, 1 application')
        
        if units=='Actual':
            y=str(np.round(d['Stemwood Biomass Growth Net DA Combined'][ikp[i]],decimals=2))
        else:
            y=str(int(d['Stemwood Biomass Growth Net DR Combined (%)'][ikp[i]]))
        ax.text(d['Duration (years)'][ikp[i]],d['N Total'][ikp[i]]-dy,y,ha='center',va='center',color=[0.7,0.2,1],fontsize=fs2)

    # Cedar River
    ikp=np.where( (d['Species leading']=='Pseudotsuga menziesii var. menziesii') &
             (d['Site Name']=='Cedar River') &
             (d['Stemwood Biomass Growth Net Units Combined']=='MgC/ha/yr') & 
             (d['Stemwood Biomass Growth Net Control Combined']>-999) & 
             (d['Thinned']!='NoNo') &
             (d['Exclusion Reason']=='None') )[0]

    for i in range(ikp.size): 
        if d['N Num Apps'][ikp[i]]==1:
            smb='o'
        elif d['N Num Apps'][ikp[i]]==2:
            smb='s'
        elif d['N Num Apps'][ikp[i]]==3:
            smb='d' 
        ax.plot(d['Duration (years)'][ikp[i]],d['N Total'][ikp[i]],smb,color=[1,0.6,0],markeredgecolor='k',linewidth=0.25,label='Cedar River, 3 applications')
        
        if units=='Actual':
            y=str(np.round(d['Stemwood Biomass Growth Net DA Combined'][ikp[i]],decimals=2))
        else:
            y=str(int(d['Stemwood Biomass Growth Net DR Combined (%)'][ikp[i]]))
        ax.text(d['Duration (years)'][ikp[i]],d['N Total'][ikp[i]]+dy,y,ha='center',va='center',color=[1,0.6,0],fontsize=fs2)

    # Wind River
    ikp=np.where( (d['Species leading']=='Pseudotsuga menziesii var. menziesii') &
             (d['Site Name']=='Wind River') &
             (d['Stemwood Biomass Growth Net Units Combined']=='MgC/ha/yr') & 
             (d['Stemwood Biomass Growth Net Control Combined']>-999) & 
             (d['Thinned']!='NoNo') &
             (d['Exclusion Reason']=='None') )[0]

    for i in range(ikp.size): 
        if d['N Num Apps'][ikp[i]]==1:
            smb='o'
        elif d['N Num Apps'][ikp[i]]==2:
            smb='s'
        elif d['N Num Apps'][ikp[i]]==3:
            smb='d' 
        if i==0:
            ax.plot(d['Duration (years)'][ikp[i]],d['N Total'][ikp[i]],smb,color=[0.9,0,0],markeredgecolor='k',linewidth=0.25,label='Wind River, 1 application')
        else:
            ax.plot(d['Duration (years)'][ikp[i]],d['N Total'][ikp[i]],smb,color=[0.9,0,0],markeredgecolor='k',linewidth=0.25)
        
        if units=='Actual':
            y=str(np.round(d['Stemwood Biomass Growth Net DA Combined'][ikp[i]],decimals=2))
        else:
            y=str(int(d['Stemwood Biomass Growth Net DR Combined (%)'][ikp[i]]))
        ax.text(d['Duration (years)'][ikp[i]],d['N Total'][ikp[i]]+dy,y,ha='center',va='center',color=[0.9,0,0],fontsize=fs2)

    ax.set(position=[0.1,0.08,0.88,0.9],xlim=[0,50],ylim=[0,1000],xlabel='Duration, years',ylabel='Total amount added (kgN ha $^-$$^1$)')
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')
    ax.legend(title='Legend',bbox_to_anchor=(1.05, 1),loc='upper left',frameon=False)
    
    return

#%% Pinus contorta

def Plot_StemwoodNetGrowth_PLI_1(d):
    
    ikp=np.where( (d['Species leading']=='Pinus contorta') &
             (d['Stemwood Biomass Growth Net Units Combined']=='MgC/ha/yr') & 
             (d['Stemwood Biomass Growth Net Control Combined']>-999) & 
             (d['Stemwood Biomass Growth Net Control Combined']<10) & 
             (d['Duration (years)']>=0) &
             (d['Duration (years)']<=120) &
             (d['Thinned']!='NoNo') &
             (d['P Num Apps']>=0) &
             (d['N Total Mean Annual (kgN/ha/yr)']>-999) & 
             (d['Exclusion Reason']=='None') )[0]

    x0=d['Stemwood Biomass Growth Net Control Combined'][ikp]
    y0=d['Stemwood Biomass Growth Net Fertilized Combined'][ikp]
    Da=np.round(np.mean(y0-x0),decimals=2)
    Dr=int(np.median((y0-x0)/x0*100))

    # Regression with intercept
    x1=sm.tools.tools.add_constant(x0)
    md=sm.OLS(y0,x1).fit()
    xhat=np.arange(0,7.1,0.1)
    yhat=md.predict(np.c_[np.ones(xhat.size),xhat])
    md.summary()

    plt.close('all');
    fig,ax=plt.subplots(1,figsize=gu.cm2inch(12,12)); 
    ms=5; xmx=1.3*np.max([np.max(y0),np.max(x0)])
    
    ax.plot([0,xmx],[0,xmx],'k-',color=[0.82,0.82,0.82],linewidth=2)
    ax.plot(x0,y0,'ko',markersize=ms,markerfacecolor='w',markeredgecolor='k')
    ax.plot(xhat,yhat,'k--',linewidth=1)
    ax.set(position=[0.08,0.08,0.9,0.9],xlim=[0,xmx],ylim=[0,xmx],xlabel='Net growth of stemwood biomass, control (MgC ha$^-$$^1$ yr$^-$$^1$)',ylabel='Net growth of stemwood biomass, fertilized (MgC ha$^-$$^1$ yr$^-$$^1$)')
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')
    
    ax.text(0.95*xmx,0.1*xmx,'Global relationship:\n' + 
        'y = ' + str(np.round(md.params[1],decimals=2)) + 'x + ' + str(np.round(md.params[0],decimals=2)) + '\n' +
        'Mean difference = ' + str(Da) + ' (MgC ha$^-$$^1$ yr$^-$$^1$)\n' +
        'Median relative difference = ' + str(Dr) + ' (%)',ha='right')
    
    return

#%% Plot stemwood net growth response of Picea glauca

def Plot_StemwoodNetGrowth_SW_1(d):

    ikp=np.where( (d['Species leading']=='Picea glauca') &
             (d['Stemwood Biomass Growth Net Units Combined']=='MgC/ha/yr') & 
             (d['Stemwood Biomass Growth Net Control Combined']>-999) & 
             (d['Stemwood Biomass Growth Net Control Combined']<10) & 
             (d['Duration (years)']>=0) &
             (d['Duration (years)']<=120) &    
             (d['P Num Apps']>=0) &
             (d['N Total Mean Annual (kgN/ha/yr)']>-999) & 
             (d['Exclusion Reason']=='None') )[0]

    x0=d['Stemwood Biomass Growth Net Control Combined'][ikp]
    y0=d['Stemwood Biomass Growth Net Fertilized Combined'][ikp]
    Da=np.round(np.mean(y0-x0),decimals=2)
    Dr=int(np.median((y0-x0)/x0*100))

    try:
        # With intercept
        x1=sm.tools.tools.add_constant(x0)
        md=sm.OLS(y0,x1).fit()
        xhat=np.arange(0,7.1,0.1)
        yhat=md.predict(np.c_[np.ones(xhat.size),xhat])
        plt.plot(xhat,yhat,'r--',linewidth=1)
        md.summary()
    except:
        pass

    plt.close('all');
    fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,12)); 
    ms=5; xmx=1.3*np.max([np.max(y0),np.max(x0)])
    ax.plot([0,xmx],[0,xmx],'k-',color=[0.82,0.82,0.82],linewidth=2)
    
    for i in range(ikp.size):
        tn=d['Treatment Note 1'][ikp[i]]
        if tn=='ON1':
            ax.plot(x0[i],y0[i],'kd',markersize=ms,markerfacecolor=[0,0.5,0],markeredgecolor='k',linewidth=0.25,label='Optimal Nutrition 1')
        elif tn=='ON2':
            ax.plot(x0[i],y0[i],'kd',markersize=ms,markerfacecolor=[0.45,1,0],markeredgecolor='k',linewidth=0.25,label='Optimial Nutrition 2')
        elif tn=='NB':
            ax.plot(x0[i],y0[i],'ko',markersize=ms,markerfacecolor=[0,1,1],markeredgecolor='k',linewidth=0.25,label='N + Boron')
        elif tn=='Complete':
            ax.plot(x0[i],y0[i],'ks',markersize=ms,markerfacecolor=[0,0,1],markeredgecolor='k',linewidth=0.25,label='Complete Mix')    
        

    ax.set(position=[0.08,0.08,0.9,0.9],xlim=[0,xmx],ylim=[0,xmx],xlabel='Net growth of stemwood biomass, control (MgC ha$^-$$^1$ yr$^-$$^1$)',ylabel='Net growth of stemwood biomass, fertilized (MgC ha$^-$$^1$ yr$^-$$^1$)')
    ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')
    ax.legend(title='Legend',bbox_to_anchor=(1.05, 1),loc='upper left',frameon=False)
    
    ax.text(0.95*xmx,0.1*xmx,'Global relationship:\n' + 
        'Mean difference = ' + str(Da) + ' (MgC ha$^-$$^1$ yr$^-$$^1$)\n' +
        'Median relative difference = ' + str(Dr) + ' (%)',ha='right')

    return
