#%% Import modules
import numpy as np
import pandas as pd
from sklearn import linear_model
import statsmodels.formula.api as smf
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.bc1ha.bc1ha_util as u1ha
gp=gu.SetGraphics('Manuscript')

#%% Import data

meta=u1ha.Init()

z0=u1ha.Import_Raster(meta,[],['refg','lc_comp1_2019','spc1_vri02','age_vri02','tdc_02','slope','tmean_ann_n','prcp_ann_n','fire_yr'],'Extract Grid')
z0['tmean_ann_n']=z0['tmean_ann_n']/10

for k in z0.keys():
    z0[k]=z0[k][0::10,0::10]

# Derived variables
fire_yr=np.zeros(z0['fire_yr'].shape)
ind=np.where( (z0['fire_yr']>=2000) ); fire_yr[ind]=1

spcg=np.ones(z0['fire_yr'].shape)
ind=np.where( (z0['spc1_vri02']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['AT']) |
             (z0['spc1_vri02']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['EP']) |
             (z0['spc1_vri02']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['ACT']) )
spcg[ind]=2

ikp=np.where( (z0['lc_comp1_2019']==1) & (z0['spc1_vri02']>0) & (z0['age_vri02']>=0) & (z0['tmean_ann_n']>=-100) & (z0['prcp_ann_n']>0) & (z0['tdc_02']>0) )

x=np.flip(np.array([fire_yr[ikp],spcg[ikp],z0['tdc_02'][ikp],z0['age_vri02'][ikp],z0['slope'][ikp],z0['tmean_ann_n'][ikp],z0['prcp_ann_n'][ikp]]).T,axis=0)
df=pd.DataFrame(x,columns=['fire_yr','spcg','tdc_02','age_vri02','slope','tmean_ann_n','prcp_ann_n'])
print(df.describe())

# Use the formula method to perform a regression
md=smf.logit("fire_yr~C(spcg)+C(tdc_02)+age_vri02+slope+tmean_ann_n+prcp_ann_n",data=df)
rs=md.fit(maxiter=100)

print(rs.summary())

df2=df.iloc[1].drop('fire_yr')
df2['spcg']=1
df2['tdc_02']=3
df2['age_vri02']=100
df2['slope']=10
df2['tmean_ann_n']=5
df2['prcp_ann_n']=1000
print(rs.predict(df2))

#%%












logr=linear_model.LogisticRegression()
logr.fit(X,y)

#predict if tumor is cancerous where the size is 3.46mm:
predicted=logr.predict(numpy.array([3.46]).reshape(-1,1))
print(predicted)

