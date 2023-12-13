
"""
PSP - SUMMARY
"""

#%% Import modules
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import copy
from scipy import stats
import statsmodels.formula.api as smf
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.bc1ha.bc1ha_util as u1ha
import fcexplore.field_plots.Processing.psp_util as ugp
gp=gu.SetGraphics('Manuscript')

#%% Import data
meta=u1ha.Init()
meta,gpt=ugp.ImportGroundPlotData(meta,type='Tree')

gpt['PTF CNY']=np.zeros(gpt['ID Plot'].size)
ind=np.where( (gpt['Plot Type']==meta['LUT']['GP']['Plot Type']['CMI']) |
             (gpt['Plot Type']==meta['LUT']['GP']['Plot Type']['NFI']) |
             (gpt['Plot Type']==meta['LUT']['GP']['Plot Type']['YSM']) )[0]
gpt['PTF CNY'][ind]=1

#%%
list(gpt.keys())

uS=np.unique(gpt['ID Species'])
vL=['Stand N L t0','Cag t0','Cag G','DBH t1','H t1']
d={}
d['Name']=np.array(['' for _ in range(uS.size)],dtype=object)
d['N']=np.zeros(uS.size)
for v in vL:
    d[v + '_mu']=np.zeros(uS.size)
    d[v + '_p99']=np.zeros(uS.size)
for iS in range(uS.size):    
    ind=np.where( (gpt['PTF CNY']==1) & (gpt['ID Species']==uS[iS]) & (gpt['Vital Status t0']==meta['LUT']['GP']['Vital Status']['Live']) )[0]
    if ind.size==0:
        continue
    d['Name'][iS]=u1ha.lut_n2s(meta['LUT']['GP']['Species'],uS[iS])[0]
    d['N'][iS]=ind.size
    for v in vL:
        d[v + '_mu'][iS]=np.nanmean(gpt[v][ind])
        d[v + '_p99'][iS]=np.nanpercentile(gpt[v][ind],99)

#%% Probability of mortality
d['Pm']=np.zeros(uS.size)
for iS in range(uS.size):    
    ind=np.where( (gpt['PTF CNY']==1) & (gpt['ID Species']==uS[iS]) & (gpt['Delta t']>0) & (gpt['Stand N L t0']<5000) & \
                 (gpt['DA t0']!=meta['LUT']['GP']['Damage Agent']['Fire']) & \
                 (gpt['DA t1']!=meta['LUT']['GP']['Damage Agent']['Fire']) )[0]
    if ind.size<100:
        continue
    x=np.flip(np.array([gpt['Mortality'][ind],gpt['Delta t'][ind]]).T,axis=0)
    df=pd.DataFrame(x,columns=['Mortality','DT'])
    
    # Use the formula method to perform a regression
    md=smf.logit("Mortality~DT",data=df)
    rs=md.fit(maxiter=100)    
    #print(rs.summary())
    df2=df.iloc[1].drop('Mortality')
    df2['DT']=1
    yhat=rs.predict(df2).values*100
    d['Pm'][iS]=yhat

#%%
plt.close('all')
plt.plot(d['Cag G_mu'],d['Pm'],'bo')

#%% Relationship between growth and density

ikp=np.where(d['N']>100)[0]
Space=np.sqrt(10000/d['Stand N L t0_mu'][ikp])
plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(10,10))
ax.plot(Space,d['Cag G_mu'][ikp],'bo')
for i in range(ikp.size):
    ax.text(Space[i]+0,d['Cag G_mu'][ikp[i]]+0.15,d['Name'][ikp[i]],ha='center',fontsize=7)
ax.set(xticks=np.arange(1,3,0.25),yticks=np.arange(0,10,1),xlabel='Space between trees (m)',ylabel='Aboveground growth (kgC/yr)',xlim=[1.25,2.75],ylim=[0,5])
#ax.legend(loc='upper left',frameon=False,facecolor=[1,1,1],labelspacing=0.25)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])
plt.tight_layout()
#gu.PrintFig(meta['Paths']['GP']['Figures'] + '\\GP_CompareDecidConifNetBiomassHarvestFoot_CN','png',900)


#%% Save
df=pd.DataFrame.from_dict(d)   
df.to_excel(r'C:\Data\tree_c_stats.xlsx')


#%%

s='BL'
dn=copy.deepcopy(meta['LUT']['GP']['Damage Agent'])
#ind0=np.where( (gpt['PTF CNY']==1) & (gpt['ID Species']==meta['LUT']['GP']['Species'][s]) & (gpt['Year t0'] )[0]
for k in meta['LUT']['GP']['Damage Agent'].keys():
    ind=np.where( (gpt['PTF CNY']==1) & (gpt['ID Species']==meta['LUT']['GP']['Species'][s]) & (gpt['DA t0']==meta['LUT']['GP']['Damage Agent'][k]) |
                  (gpt['PTF CNY']==1) & (gpt['ID Species']==meta['LUT']['GP']['Species'][s]) & (gpt['DA t1']==meta['LUT']['GP']['Damage Agent'][k]) )[0]
    dn[k]=ind.size#/ind0.size*100
    #cnt=cnt+1
dn

dn2=copy.deepcopy(meta['LUT']['GP']['Damage Agent'])
ind0=np.where( (gpt['PTF CNY']==1) & (gpt['ID Species']==meta['LUT']['GP']['Species'][s]) & (gpt['Mortality']==1) )[0]
for k in meta['LUT']['GP']['Damage Agent'].keys():
    ind=np.where( (gpt['PTF CNY']==1) & (gpt['ID Species']==meta['LUT']['GP']['Species'][s]) & (gpt['DA t0']==meta['LUT']['GP']['Damage Agent'][k]) & (gpt['Mortality']==1) )[0]
    dn2[k]=ind.size/ind0.size*100
    #cnt=cnt+1
dn2  

#%%
ind0=np.where( (gpt['ID Species']==meta['LUT']['GP']['Species']['SX']) )[0]
ikp=np.where((gpt['ID Species']==meta['LUT']['GP']['Species']['AT']))[0]
plt.close('all');
plt.plot(gpt['Lon'][ikp],gpt['Lat'][ikp],'b.')





