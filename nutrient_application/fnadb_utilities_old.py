'''

NUTC DATABASE ANALYSIS

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

#%% Set paths

PathFigures=r'G:\My Drive\Figures\NutC Database'

#%% Set figure properties

params={'font.sans-serif':'Calibri',
        'font.size':11,
        'axes.edgecolor':'black',
        'axes.labelsize':11,
        'axes.labelcolor':'black',
        'axes.titlesize':7,
        'axes.linewidth':0.5,        
        'lines.linewidth':0.5,
        'text.color':'black',
        'xtick.color':'black',        
        'xtick.labelsize':11,
        'xtick.major.width':0.5,
        'xtick.major.size':3,
        'xtick.direction':'in',
        'ytick.color':'black',
        'ytick.labelsize':11,
        'ytick.major.width':0.5,
        'ytick.major.size':3,
        'ytick.direction':'in',
        'legend.fontsize':7,        
        'savefig.dpi':300,
        'savefig.transparent':True}
plt.rcParams.update(params)

#%% Import NutC database

path=r'C:\Users\rhember\Documents\Data\Nutrient Addition Experiments Database\Nutrient Addition Experiments Database.xlsx'
d=gu.ReadExcel(path,'Table1')

d['N Total']=d['N Num Apps']*d['N Dose Per App (kgN/ha)']


#%% Global relative net stemwood biomass response

ikp=np.where( (d['Stemwood Volume Growth Net DR (%)']>-100) & 
             (d['Stemwood Volume Growth Net DR (%)']<=3000) & 
             (np.isnan(d['Stemwood Volume Growth Net DR (%)'])==False) &
             (d['Duration (years)']>=5) &
             (d['Thinned']!='NoNo') &
             (d['P Num Apps']==0) )[0]
#(d['Biome']=='Cold Temperate')

y=d['Stemwood Volume Growth Net DR (%)'][ikp]

ss=ikp.size
med=np.median(y)
mu=np.mean(y)
sd=np.std(y)
se=np.std(y)/np.sqrt(ikp.size)

NR=10000
medR=np.nan*np.ones(NR)
for i in range(NR):
    medR[i]=np.median(resample(y,n_samples=y.size))
mu=np.median(medR)
p=np.percentile(medR,[2.5,97.5])

print(np.median(d['N Dose Per App (kgN/ha)'][ikp]))



n,bins,patches=plt.hist(y,np.arange(-100,3000,10))

plt.close('all'); 
fig,ax=plt.subplots(1,figsize=gu.cm2inch(15.5,6))
ax.bar(bins[0:-1],n,10,facecolor=[0.78,0.78,0.78])
ax.set(position=[0.07,0.16,0.89,0.78],
       xlim=[-100,400],yticks=np.arange(0,60,10),
       xscale='linear',xticks=[0,100,200,300,400],       
       ylabel='Number of instalations',xlabel='% response');
ax.yaxis.grid()
plt.text(100,44,r'Median = 25%, Mean = 75%, S.D. = 227%')
gu.PrintFig(PathFigures + '\\NutC_GlobalRelativeStemwoodResponse','png',900)


#%% Coastal Douglas-fir net stemwood biomass production

plt.close('all'); 
fig,ax=plt.subplots(1)

ikp=np.where( (d['Species leading']=='Pseudotsuga menziesii var. menziesii') & \
            (d['Stemwood Biomass Growth Net Units Combined']=='MgC/ha/yr') & \
             (np.isnan(d['Stemwood Biomass Growth Net DA Combined'])==0) & \
             (d['Reference']!='Mitchell et al 1996') & \
             (d['N Num Apps']==1) )[0]
y=d['Stemwood Biomass Growth Net DA Combined']
#y=d['Stemwood Biomass Growth Net DR Combined (%)']
x=d['Duration (years)']
ax.plot(x[ikp],y[ikp],'g^',color=[0.8,0.8,0.8])

ikp=np.where( (d['Species leading']=='Pseudotsuga menziesii var. menziesii') & \
             (d['Stemwood Biomass Growth Net Units Combined']=='MgC/ha/yr') & \
             (np.isnan(d['Stemwood Biomass Growth Net DA Combined'])==0) & \
             (d['Experiment Name']=='Shawnigan Lake') & \
             (d['Reference']!='Mitchell et al 1996') & \
             (d['N Num Apps']==1))[0]
y=d['Stemwood Biomass Growth Net DA Combined']
#y=d['Stemwood Biomass Growth Net DR Combined (%)']
x=d['Duration (years)']
ax.plot(x[ikp],y[ikp],'bo')

ikp=np.where( (d['Species leading']=='Pseudotsuga menziesii var. menziesii') & \
             (d['Stemwood Biomass Growth Net Units Combined']=='MgC/ha/yr') & \
             (np.isnan(d['Stemwood Biomass Growth Net DA Combined'])==0) & \
             (d['Experiment Name']=='EP703') & \
             (d['N Num Apps']==1))[0]
y=d['Stemwood Biomass Growth Net DA Combined']
y#=d['Stemwood Biomass Growth Net DR Combined (%)']
x=d['Duration (years)']
ax.plot(x[ikp],y[ikp],'co')

ikp=np.where( (d['Species leading']=='Pseudotsuga menziesii var. menziesii') & \
             (d['Stemwood Biomass Growth Net Units Combined']=='MgC/ha/yr') & \
             (np.isnan(d['Stemwood Biomass Growth Net DA Combined'])==0) & \
             (d['N Num Apps']>1))[0]

y=d['Stemwood Biomass Growth Net DA Combined']
y#=d['Stemwood Biomass Growth Net DR Combined (%)']
x=d['Duration (years)']
ax.plot(x[ikp],y[ikp],'rs')
ax.set(xlabel='Time since application',ylabel='Stemwood biomass growth (Mg C ha-1 yr-1)')

#%% Coastal Douglas-fir net stemwood biomass production

plt.close('all'); 
fig,ax=plt.subplots(1)

ikp=np.where( (d['Species leading']=='Pseudotsuga menziesii var. menziesii') & \
             (d['Stemwood Biomass Growth Net Units Combined']=='MgC/ha/yr') & \
             (np.isnan(d['Stemwood Biomass Growth Net DA Combined'])==0) & \
             (d['Reference']!='Mitchell et al 1996') & \
             (d['N Num Apps']>0))[0]

y=d['Stemwood Biomass Growth Net DA Combined'][ikp]
#x=d['Duration (years)']
x=d['N Num Apps'][ikp]*d['N Dose Per App (kgN/ha)'][ikp]
ax.plot(x,y,'bo')
ax.set(xlabel='Time since application',ylabel='Stemwood biomass growth (Mg C ha-1 yr-1)')

#%% Interior Douglas-fir net stemwood biomass growth

#nam='Tsuga heterophylla'
nam='Pseudotsuga menziesii var. menziesii'
#nam='Pseudotsuga menziesii var. glauca'
#nam='Pinus contorta'
#nam='Picea glauca'
#nam='Picea mariana'
#nam='Pinus sylvestris'
#nam='Pinus taeda'
#nam='Thuja plicata'
#nam='Picea sitchensis'
nam='Pinus ponderosa'
#nam='Picea amabilis'
#nam='Picea abies'
#nam='Pinus banksiana'

ikp=np.where( (d['Species leading']==nam) & \
             (d['Stemwood Biomass Growth Net Units Combined']=='MgC/ha/yr') & \
             (np.isnan(d['Stemwood Biomass Growth Net DA Combined'])==0) & \
             (d['Reference']!='Mitchell et al 1996') & \
             (d['N Num Apps']>0) & \
             (d['P Num Apps']>0) & \
             (d['N Dose Per App (kgN/ha)']>0) )[0]
#(d['N Type']=='Urea') & \

da=d['Stemwood Biomass Growth Net DA Combined'][ikp]
dr=d['Stemwood Biomass Growth Net DR Combined (%)'][ikp]
N=d['Num Installations'][ikp]
L=d['Duration (years)'][ikp]
Nadd=d['N Num Apps'][ikp]*d['N Dose Per App (kgN/ha)'][ikp]
da_w=np.nansum(da*N)/np.nansum(N)
dr_w=np.nansum(dr*N)/np.nansum(N)
L_w=np.nansum(L*N)/np.nansum(N)
Nadd_w=np.nansum(Nadd*N)/np.nansum(N)

print(np.nansum(N))
print(L_w)
print(Nadd_w)
print(da_w)
print(np.nanstd(da))





#%% Lodgepole pine net stemwood biomass production

plt.close('all'); 
fig,ax=plt.subplots(1)

ikp=np.where( (d['Species leading']=='Pinus contorta') & \
             (d['Stemwood Biomass Growth Net Units Combined']=='MgC/ha/yr') & \
             (np.isnan(d['Stemwood Biomass Growth Net DA Combined'])==0) & \
             (d['N Num Apps']>0))[0]

y=d['Stemwood Biomass Growth Net DA Combined'][ikp]
x=d['Duration (years)'][ikp]
#x=d['N Num Apps'][ikp]*d['N Dose Per App (kgN/ha)'][ikp]
#x=d['N Added Mean Annual (kgN/ha/yr)'][ikp]
ax.plot(x,y,'bo')
ax.set(xlabel='Time since application',ylabel='Stemwood biomass growth (Mg C ha-1 yr-1)')


#%% Interior species net stemwood biomass production

ikp=np.where( (d['Species leading']=='Pinus contorta') & \
             (d['Stemwood Biomass Growth Net Units Combined']=='MgC/ha/yr') & \
             (np.isnan(d['Stemwood Biomass Growth Net DA Combined'])==0) & \
             (d['N Num Apps']>0) |
             (d['Species leading']=='Picea glauca') & \
             (d['Stemwood Biomass Growth Net Units Combined']=='MgC/ha/yr') & \
             (np.isnan(d['Stemwood Biomass Growth Net DA Combined'])==0) & \
             (d['N Num Apps']>0) |
             (d['Species leading']=='Pseudotsuga menziesii var. glauca') & \
             (d['Stemwood Biomass Growth Net Units Combined']=='MgC/ha/yr') & \
             (np.isnan(d['Stemwood Biomass Growth Net DA Combined'])==0) & \
             (d['N Num Apps']>0) )[0]

plt.close('all'); 
fig,ax=plt.subplots(1)
y#=d['Stemwood Biomass Growth Net DA Combined'][ikp]
y=d['Stemwood Biomass Growth Net DR Combined (%)'][ikp]
#x=d['Duration (years)'][ikp]
#x=d['N Num Apps'][ikp]*d['N Dose Per App (kgN/ha)'][ikp]
x=d['N Added Mean Annual (kgN/ha/yr)'][ikp]
ax.plot(x,y,'bo')
ax.set(xlabel='Time since application',ylabel='Stemwood biomass growth (Mg C ha-1 yr-1)')




#%% Net stemwood biomass response by Biome

biome=['Warm Temperate','Cold Temperate','Boreal']
SS=np.zeros((len(biome),1))
Duration=np.zeros((len(biome),2))
Na=np.zeros((len(biome),2))
DA=np.zeros((len(biome),3))
DR=np.zeros((len(biome),3))
NUE=np.zeros((len(biome),3))
for i in range(len(biome)):
    ikp=np.where( (d['Biome']==biome[i]) & (d['Stemwood Biomass Growth Net Units Combined']=='MgC/ha/yr') & (d['Stemwood Biomass Growth Net DA Combined']>-1) )[0]
    SS[i]=ikp.size
    DR[i,0]=np.nanmean(d['Stemwood Biomass Growth Net DR Combined (%)'][ikp])
    DR[i,1]=np.nanmedian(d['Stemwood Biomass Growth Net DR Combined (%)'][ikp])
    DR[i,2]=np.nanstd(d['Stemwood Biomass Growth Net DR Combined (%)'][ikp])/np.sqrt(ikp.size)
    DA[i,0]=np.nanmean(d['Stemwood Biomass Growth Net DA Combined'][ikp])
    DA[i,1]=np.nanmedian(d['Stemwood Biomass Growth Net DA Combined'][ikp])
    DA[i,2]=np.nanstd(d['Stemwood Biomass Growth Net DA Combined'][ikp])/np.sqrt(ikp.size)
    NUE[i,0]=np.nanmean(d['Stemwood Biomass Growth Net DA Combined'][ikp]*1000/d['N Added Mean Annual (kgN/ha/yr)'][ikp])
    NUE[i,1]=np.nanmedian(d['Stemwood Biomass Growth Net DA Combined'][ikp]*1000/d['N Added Mean Annual (kgN/ha/yr)'][ikp])
    NUE[i,2]=np.nanstd(d['Stemwood Biomass Growth Net DA Combined'][ikp]*1000/d['N Added Mean Annual (kgN/ha/yr)'][ikp])/np.sqrt(ikp.size)
    Duration[i,0]=np.nanmean(d['Duration (years)'][ikp])
    Duration[i,1]=np.nanmedian(d['Duration (years)'][ikp])
    Na[i,0]=np.nanmean(d['N Added Mean Annual (kgN/ha/yr)'][ikp])
    Na[i,1]=np.nanmedian(d['N Added Mean Annual (kgN/ha/yr)'][ikp])


plt.close('all'); 
fig,ax=plt.subplots(1,figsize=gu.cm2inch(15.5,10))
ax.bar([1,2,3],DA[:,1],0.8,facecolor=[0,0,0.5],yerr=DA[:,2],error_kw=dict(lw=1.5,capsize=4,capthick=1),ecolor=[0,0,0])
ax.set(position=[0.10,0.07,0.88,0.88],xlim=[0.5,3.5],yticks=np.arange(0,1.3,0.1),xticks=[1,2,3,4],ylabel='Net stemwood biomass growth (Mg C ha$^-$$^1$ yr$^-$$^1$)');
ax.set_xticklabels(['Warm temperate','Cold temperate','Boreal']);
gu.PrintFig(PathFigures + '\\NutC_StemwoodNetBiomassGrowthByBiome_MedianDA','png',900)

plt.close('all'); 
fig,ax=plt.subplots(1,figsize=gu.cm2inch(15.5,10))
ax.bar([1,2,3],DR[:,1],0.8,facecolor=[0,0,0.5],yerr=DR[:,2],error_kw=dict(lw=1.5,capsize=4,capthick=1),ecolor=[0,0,0])
ax.set(position=[0.10,0.07,0.88,0.88],xlim=[0.5,3.5],yticks=np.arange(0,100,10),xticks=[1,2,3,4],ylabel='Net stemwood biomass growth (%)');
ax.set_xticklabels(['Warm temperate','Cold temperate','Boreal']);
gu.PrintFig(PathFigures + '\\NutC_StemwoodNetBiomassGrowthByBiome_MedianDR','png',900)

plt.close('all'); 
fig,ax=plt.subplots(1,figsize=gu.cm2inch(15.5,10))
ax.bar([1,2,3],NUE[:,1],0.8,facecolor=[0,0,0.5],yerr=NUE[:,2],error_kw=dict(lw=1.5,capsize=4,capthick=1),ecolor=[0,0,0])
ax.set(position=[0.10,0.07,0.88,0.88],xlim=[0.5,3.5],yticks=np.arange(0,60,10),xticks=[1,2,3,4],ylabel='Nitrogen use efficiency (kg C kg N$^-$$^1$)');
ax.set_xticklabels(['Warm temperate','Cold temperate','Boreal']);
gu.PrintFig(PathFigures + '\\NutC_StemwoodNetBiomassGrowthByBiome_MedianNUE','png',900)



#%% Gross growth vs. Mortality response

plt.close('all'); 
fig,ax=plt.subplots(1)

ikp=np.where( (d['Stemwood Biomass Growth Gross DA Combined']>0) & (d['Stemwood Biomass Mortality DA Combined']<2) & \
              (d['Stemwood Biomass Mortality DA Combined']>0) )[0]

y=d['Stemwood Biomass Growth Gross DR Combined (%)'][ikp]
x=d['Stemwood Biomass Mortality DR Combined (%)'][ikp]
ax.plot(x,y,'bo')
ax.set(xlabel='Time since application',ylabel='Stemwood biomass growth (Mg C ha-1 yr-1)')



#%% Ratio of non-stemwood biomass to stemwood biomass response

drStemwood=d['Stemwood Biomass Growth Net DR (%)']
drBranch=d['Branch Biomass Growth Net DR (%)']
drFoliage=d['Foliage Biomass Growth Net DR (%)']
drRootC=d['Root Coarse Biomass Net Growth DR (%)']
drRootF=d['Root Fine Biomass Net Growth DR (%)']
dr=np.column_stack([drBranch,drFoliage,drRootC,drRootF])

ikp=np.where( (d['Stemwood Biomass Growth Net DR (%)']>-10) & \
    (d['N Num Apps']>0) )[0]

dr_med=np.nanmedian(dr[ikp,:],axis=0)
dr_ss=np.sum(~np.isnan(dr[ikp,:]),axis=0)
dr_se=np.nanstd(dr[ikp,:],axis=0)/np.sqrt(dr_ss)

flg=0
if flg==1:
    plt.close('all'); 
    fig,ax=plt.subplots(1,figsize=gu.cm2inch(12,8))
    ax.bar([1,2,3,4],dr_med,0.8,facecolor=[0,0,0.5],yerr=dr_se,error_kw=dict(lw=1.5,capsize=4,capthick=1),ecolor=[0,0,0])
    ax.set(position=[0.11,0.07,0.88,0.88],ylim=[0,100],xticks=[1,2,3,4],xlim=[0.5,4.5],ylabel='Median response (%)')
    ax.set_xticklabels(['Branches','Foliage','Coarse roots','Fine roots']);

rBranch=d['Branch Biomass Growth Net DR (%)']/d['Stemwood Biomass Growth Net DR Combined (%)']
#rBark=d['Bark Biomass Growth Net DR (%)']/d['Stemwood Biomass Growth Net DR Combined (%)']
rFoliage=d['Foliage Biomass Growth Net DR (%)']/d['Stemwood Biomass Growth Net DR Combined (%)']
rRootC=d['Root Coarse Biomass Net Growth DR (%)']/d['Stemwood Biomass Growth Net DR Combined (%)']
rRootF=d['Root Fine Biomass Net Growth DR (%)']/d['Stemwood Biomass Growth Net DR Combined (%)']
r=np.column_stack([rBranch,rFoliage,rRootC,rRootF])

print(np.sum(~np.isnan(r),axis=0))
ikp=np.where( (d['Stemwood Biomass Growth Net DR (%)']>-10) & \
    (d['N Num Apps']>0) & \
    (rFoliage<1000) )[0]
r_med=np.nanmedian(r[ikp,:],axis=0)
r_ss=np.sum(~np.isnan(r[ikp,:]),axis=0)
r_se=np.nanstd(r[ikp,:],axis=0)/np.sqrt(r_ss)

plt.close('all'); 
fig,ax=plt.subplots(1,figsize=gu.cm2inch(12,8))
ax.plot([0,5],[1,1],'k--',color=[0.6,0.6,0.6],linewidth=2)
ax.bar([1,2,3,4],r_med,0.8,facecolor=[0,0,0.5],yerr=r_se,error_kw=dict(lw=1.5,capsize=4,capthick=1),ecolor=[0,0,0])
ax.set(position=[0.11,0.07,0.88,0.88],xlim=[0.5,4.5],xticks=[1,2,3,4],ylim=[0,2],ylabel='Response ratio')
ax.set_xticklabels(['Branches','Foliage','Coarse roots','Fine roots']);
gu.PrintFig(PathFigures + '\\NutC_NonStemwoodToStemwoodBiomassRatio','png',900);



#%% Summary comparison

#------------------------------------------------------------------------------
cNPP=[]; 

# From demo without any explicit response of mortality or decomposition
d={}
d['Source']='This study'
d['BE']=5.3
d['CI']=np.array([np.nan,np.nan])
cNPP.append(d)

# LeBaur and Treseder (2008)
d={}
d['Source']='LeBaur and Treseder (2008)'
d['BE']=29.0
d['CI']=np.array([np.nan,np.nan])
cNPP.append(d)

# Yue et al. (2016)
d={}
d['Source']='Yue et al. (2016)'
d['BE']=52.4
d['CI']=np.array([40.6,65.2])
cNPP.append(d)

plt.close('all'); 
fig,ax=plt.subplots(1,figsize=gu.cm2inch(7,5))
y=cNPP; cm=plt.cm.get_cmap('viridis',len(y)+1)
for i in range(len(y)):
    ax.bar(i+1,y[i]['BE'],0.8,label=y[i]['Source'],facecolor=cm.colors[i,0:3]) # ,yerr=y[i]['CI'],error_kw=dict(lw=1.5,capsize=4,capthick=1),ecolor=[0,0,0]
    lab.append(y[i]['Source'])
ax.set(position=[0.15,0.03,0.82,0.94],xlim=[0.5,len(y)+.5],xticks=np.arange(1,len(y)+1,1),xticklabels='',ylabel='Response (%)')
#ax.set_xticklabels(lab);
leg=ax.legend(frameon=False,facecolor=None);
# leg.get_frame().set_linewidth(0.0)
gu.PrintFig(PathFigures + '\\NutC_CompareNPP','png',900);

#------------------------------------------------------------------------------
cM=[]; 

# From demo without any explicit response of mortality or decomposition
d={}
d['Source']='This study'
d['BE']=19.6
d['CI']=np.array([np.nan,np.nan])
cM.append(d)

# LeBaur and Treseder (2008)
d={}
d['Source']='LeBaur and Treseder (2008)'
d['BE']=29.0
d['CI']=np.array([np.nan,np.nan])
cM.append(d)

# Yue et al. (2016)
d={}
d['Source']='Yue et al. (2016)'
d['BE']=52.4
d['CI']=np.array([40.6,65.2])
cNPP.append(d)

plt.close('all'); 
fig,ax=plt.subplots(1,figsize=gu.cm2inch(7,5))
y=cNPP; cm=plt.cm.get_cmap('viridis',len(y)+1)
for i in range(len(y)):
    ax.bar(i+1,y[i]['BE'],0.8,label=y[i]['Source'],facecolor=cm.colors[i,0:3]) # ,yerr=y[i]['CI'],error_kw=dict(lw=1.5,capsize=4,capthick=1),ecolor=[0,0,0]
    lab.append(y[i]['Source'])
ax.set(position=[0.15,0.03,0.82,0.94],xlim=[0.5,len(y)+.5],xticks=np.arange(1,len(y)+1,1),xticklabels='',ylabel='Response (%)')
#ax.set_xticklabels(lab);
leg=ax.legend(frameon=False,facecolor=None);
# leg.get_frame().set_linewidth(0.0)
gu.PrintFig(PathFigures + '\\NutC_CompareNPP','png',900);


#%%

# From demo without any explicit response of mortality or decomposition
sim={}
sim['NPP']=5.3
sim['Stemwood mortality']=19.6
sim['Litterfall']=3.3
sim['Decomposition']=2.0
sim['Dead wood']=3.6
sim['Litter']=1.0
sim['Soil organic matter']=0.1

# From NutC database
ex={}

# Stemwood mortality
ikp=np.where( (d['Stemwood Biomass Mortality DA Combined']>0) & (d['Stemwood Biomass Mortality DA Combined']<5) & \
             (d['Stemwood Biomass Mortality DR Combined (%)']<1000) )[0]
da_mu=np.nanmean(d['Stemwood Biomass Mortality DA Combined'][ikp])
dr_mu=np.nanmean(d['Stemwood Biomass Mortality DA Combined'][ikp])
dr_md=np.nanmedian(d['Stemwood Biomass Mortality DR Combined (%)'][ikp])
dr_se=np.nanstd(d['Stemwood Biomass Mortality DR Combined (%)'][ikp])/np.sqrt(ikp.size)
ex['Stemwood mortality']=np.array([dr_md,dr_se])

# Litterfall
y=d['Litterfall DR (%)']
ikp=np.where( (np.isnan(y)==False) & (y>-100) & (y<1000) )[0]
dr_mu=np.nanmean(y[ikp])
dr_md=np.nanmedian(y[ikp])
dr_se=np.nanstd(y[ikp])/np.sqrt(ikp.size)
ex['Litterfall']=np.array([dr_md,dr_se])

# Fine root mortality
y=d['Root Fine Biomass Mortality DR (%)']
ikp=np.where( (np.isnan(y)==False) & (y>-100) & (y<1000) )[0]
dr_mu=np.nanmean(y[ikp])
dr_md=np.nanmedian(y[ikp])
dr_se=np.nanstd(y[ikp])/np.sqrt(ikp.size)
ex['Fine root mortality']=np.array([dr_md,dr_se])

# Decomposition
sim['Decomposition']=2.0
sim['Dead wood']=3.6
sim['Litter']=1.0
sim['Soil organic matter']=0.1


#%% N response of foliage

y=d['Foliage Biomass Growth Net DR (%)']
#y=d['Foliage Biomass Growth Gross DR (%)']

x=d['N Added Mean Annual (kgN/ha)']*d['Interval (years)']

ind=np.where( (np.isnan(y)==False) & (y>-100) )[0]

plt.close('all'); 
fig,ax=plt.subplots(1)
plt.plot(x[ind],y[ind],'o')
#ax.set(xscale='log',yscale='log')
print(ind.size)
print(np.nanmean(y[ind]))
print(np.nanmedian(y[ind]))

#%% Duration effect

plt.close('all'); 
fig,ax=plt.subplots(1)
u=np.unique(d['Species leading'])
da=np.nan*np.ones(u.size)
dur=np.nan*np.ones(u.size)
for i in range(u.size):
    ikp=np.where( (d['Species leading']==u[i]) & 
                 (d['Stemwood Biomass Growth Net DA Combined']<5) & 
                 (d['N Num Apps']==1) )[0] 
    da[i]=np.nanmean(d['Stemwood Biomass Growth Net DA Combined'][ikp])
    dur[i]=np.nanmean(d['Duration (years)'][ikp])
    plt.plot(dur[i],da[i],'o')

ind=np.where(np.isnan(da)!=True)[0]
da[ind]