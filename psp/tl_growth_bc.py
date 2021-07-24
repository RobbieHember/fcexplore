
#%% Import modules

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy.io as spio
import statsmodels.api as sm
import statsmodels.formula.api as smf
import fcgadgets.macgyver.utilities_general as gu
from fcgadgets.psp import psp_utilities as utl

#%% Import data

q=utl.Read_L3_TL_BySpc(['PL'])

# Add derived variables
q['Ln_Csw_G']=np.log(q['Csw_G'])
q['Ln_Csw_t0']=np.log(q['Csw_t0'])
q['Ln_N_L_t0_Stand']=np.log(q['N_L_t0_Stand'])
q['Ln_Csw_L_t0_Stand']=np.log(q['Csw_L_t0_Stand'])
q['Ln_Cag_Larger_t0_Stand']=np.log(q['Cag_Larger_t0_Stand'])
q['H_t02']=q['H_t0']**2
q['ws_mjjas_n2']=q['ws_mjjas_n']**2
q['ws_mjjas_a2']=q['ws_mjjas_a']**2;

#np.percentile(q['Age_t0_Stand'],99.5)

#%% Exclusion rules

ikp=np.where( (q['ID_DB']==1) & \
             (q['EcoZone_CA_L1']!=13) & \
             (q['Mort']==0) & \
             (q['H_t0']>0) & (q['H_t0']<140) & \
             (q['Csw_t0']>0) & (q['Csw_t0']<3000) & \
             (q['Csw_G']>0) & (q['Csw_G']<100) & \
             (q['Age_t0_Stand']>0) & (q['Age_t0_Stand']<500) & \
             (q['N_L_t0_Stand']>0) & (q['N_L_t0_Stand']<10000) & \
             (q['Csw_L_t0_Stand']>0) & (q['Csw_L_t0_Stand']<2000) & \
             (q['Cag_Larger_t0_Stand']>=0) & (q['Cag_Larger_t0_Stand']<1000) & \
             (q['Ln_Cag_Larger_t0_Stand']>-1000) & (q['Ln_Cag_Larger_t0_Stand']<1000) & \
             (q['ws_mjjas_n']>=0) & (q['ws_mjjas_n']<=200) & \
             (q['ws_mjjas_a']>=-200) & (q['ws_mjjas_n']<=200) \
             )[0]

#%% Remove unnecessary variables

vs=['ID_Plot','H_t0','H_t02','Csw_G','Ln_Csw_G','Csw_t0','Ln_Csw_t0','Ln_N_L_t0_Stand','Ln_Csw_L_t0_Stand','Ln_Cag_Larger_t0_Stand','Age_t0_Stand','ws_mjjas_n','ws_mjjas_n2','ws_mjjas_a','ws_mjjas_a2']

q1={}
for key in vs:
    q1[key]=q[key][ikp]

#%% Define dpeendent variable, direct x variables

# Dependent variable string
y_s='Ln_Csw_G'

# x-factor string
xfactor_s='ws_mjjas_a'

# Independent variable direct list
xdir=['H_t0','H_t02','Age_t0_Stand','Ln_Cag_Larger_t0_Stand','Ln_Csw_L_t0_Stand','Ln_N_L_t0_Stand','ws_mjjas_n','ws_mjjas_n2','ws_mjjas_a','ws_mjjas_a2']

# Independent variable direct string
xdir_s=''
for k in xdir:
    xdir_s=xdir_s + k + ' + '

# Independent variable interactions list
xint=['H_t0','H_t02','Ln_N_L_t0_Stand','Ln_Csw_L_t0_Stand','Ln_Cag_Larger_t0_Stand','Age_t0_Stand','ws_mjjas_n']
xint_s=''
for k in xint:
    xint_s=xint_s + xfactor_s + '_x_' + k + ' + '
# Remove the last Plus symbol
xint_s=xint_s[0:-3]

#%% Calculate z-scores

q1z={}

stats={}
stats['mu']={} 
stats['sig']={}
stats['min']={}
stats['max']={}
stats['min_z']={}
stats['max_z']={}

for key in vs:    
    if (key=='ID_Plot'):
        continue
    stats['mu'][key]=np.nanmean(q1[key])
    stats['sig'][key]=np.nanstd(q1[key])
    stats['sig'][key]=np.nanstd(q1[key])
    stats['min'][key]=np.min(q1[key])
    stats['max'][key]=np.max(q1[key])
    
    q1z[key]=(q1[key]-stats['mu'][key])/stats['sig'][key]
    
    stats['min_z'][key]=np.min(q1z[key])
    stats['max_z'][key]=np.max(q1z[key])

# Check NaNs
N_nan=np.zeros(len(q1z)); cnt=0
for k in q1z.keys():
    N_nan[cnt]=np.sum(np.isnan(q1z[k])); cnt=cnt+1

#%% Add interactions

stats_int={}
stats_int['mu']={}
stats_int['sig']={}
for i in range(len(xint)):
    nam=xfactor_s + '_x_' + xint[i]
    
    tmp=q1z[xfactor_s]*q1z[xint[i]]
    
    stats_int['mu'][nam]=np.mean(tmp)
    stats_int['sig'][nam]=np.std(tmp)
    
    q1z[nam]=(tmp-stats_int['mu'][nam])/stats_int['sig'][nam]

#%% Dataframe

df=pd.DataFrame(data=q1z)
df.shape

df=df.dropna()
df.shape
df=df.reset_index(drop=True)

#%% Fit random intercept model

func_s=y_s + ' ~ ' + xdir_s + xint_s

md=smf.mixedlm(func_s,df,groups=q1['ID_Plot'])

mdf=md.fit()
print(mdf.summary())

#%% Responses to xfactor

ifactor_z=[-1.5,0,1.5]

ifactor_s='H_t0'
#ifactor_s='Ln_N_L_t0_Stand'
#ifactor_s='Ln_Csw_L_t0_Stand'
#ifactor='Ln_Cag_Larger_t0_Stand'
#ifactor='ws_mjjas_n'

y=[]
for i_ifactor in range(len(ifactor_z)):
    
    r0={}
    tmp=np.linspace(stats['min'][xfactor_s],stats['max'][xfactor_s],30)
    N=tmp.size

    # Add direct effects
    for k in xdir:
        if k=='ID_Plot':
            continue
        elif k==xfactor_s:
            r0[k]=tmp.copy()
        elif k==ifactor_s:
            r0[k]=(stats['mu'][k]+ifactor_z[i_ifactor]*stats['sig'][k])*np.ones(N) 
        else:
            r0[k]=stats['mu'][k]*np.ones(N)   
    r0[xfactor_s + '2']=r0[xfactor_s]**2
    
    if ifactor_s=='H_t0':
        r0[ifactor_s + '2']=r0[ifactor_s]**2

    # Standardize variables
    r1=r0.copy()
    for k in r1.keys():
        r1[k]=(r0[k]-stats['mu'][k])/stats['sig'][k]

    # Add interactive effects
    for i in range(len(xint)):
        nam=xfactor_s + '_x_' + xint[i]
        tmp=r1[xfactor_s]*r1[xint[i]]
        r1[nam]=(tmp-stats_int['mu'][nam])/stats_int['sig'][nam]

    y0=mdf.params['Intercept']
    for k in r1.keys():
        y0=y0+mdf.params[k]*r1[k]
    y.append(y0)

#%%
    
plt.close('all')
plt.plot(r0[xfactor_s],y[0],'-bo')
plt.plot(r0[xfactor_s],y[1],'-cs')
plt.plot(r0[xfactor_s],y[2],'-gd')





#%%
    
    
    
    
#%% Fit random intercept model

funcs='Ln_Csw_G~H_t0+H_t02+Age_t0_Stand+Ln_Cag_Larger_t0_Stand+Ln_Csw_L_t0_Stand+Ln_N_L_t0_Stand+ws_mjjas_n+ws_mjjas_n2+ws_mjjas_a+ws_mjjas_a2'
md=smf.mixedlm(funcs,df,groups=q1['ID_Plot'])
mdf=md.fit()
print(mdf.summary())

#%% Partial correlation

#cp=gu.PartialCorrelation(df.to_numpy())*100
#print(cp.astype(int))