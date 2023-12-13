'''
THIS SCRIPT IS A MESS / RETIRED
*** Comp 1 moved to a function in util ***
'''

#%% Import modules

import numpy as np
import matplotlib.pyplot as plt
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.bc1ha.bc1ha_util as u1ha
gp=gu.SetGraphics('Manuscript')

#%% Import parameters
meta=u1ha.Init()

#%% 


#%% Summarize areas

labL=list(meta['LUT']['Derived']['lcc1'].keys())
A=np.zeros(len(labL))
cnt=0
for k in meta['LUT']['Derived']['lcc1'].keys():
    ind=np.where( (zRef['Data']==1) & (zLCC1['Data']==meta['LUT']['Derived']['lcc1'][k]) )
    A[cnt]=ind[0].size/1e6
    cnt=cnt+1

# Plot bar chart of area
x=np.arange(0,A.size,1)
plt.close('all');
fig,ax=plt.subplots(1,figsize=gu.cm2inch(10,6.5));
ax.bar(x,A)
ax.set(xticks=x,xticklabels=labL,ylabel='Area (Mha)')
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=gp['tickl'])



#%% Disagreement in forest class

ind=np.where( (z['lc2']['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_2']['T']) & \
             (np.isin(z['lcc_ntem_2019']['Data'],[meta['LUT']['Derived']['lcc_ntem']['Coniferous'],meta['LUT']['Derived']['lcc_ntem']['Broadleaf'],meta['LUT']['Derived']['lcc_ntem']['Mixedwood'],meta['LUT']['Derived']['lcc_ntem']['Wetland-treed']])==True) )
print(ind[0].size/1e6)

ind=np.where( (z['lc2']['Data']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_2']['T']) & \
             (np.isin(z['lcc_ntem_2019']['Data'],[meta['LUT']['Derived']['lcc_ntem']['Coniferous'],meta['LUT']['Derived']['lcc_ntem']['Broadleaf'],meta['LUT']['Derived']['lcc_ntem']['Mixedwood'],meta['LUT']['Derived']['lcc_ntem']['Wetland-treed']])==False) )
print(ind[0].size/1e6)

ind=np.where( (z['lc2']['Data']!=meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_2']['T']) & \
             (np.isin(z['lcc_ntem_2019']['Data'],[meta['LUT']['Derived']['lcc_ntem']['Coniferous'],meta['LUT']['Derived']['lcc_ntem']['Broadleaf'],meta['LUT']['Derived']['lcc_ntem']['Mixedwood'],meta['LUT']['Derived']['lcc_ntem']['Wetland-treed']])==True) )
print(ind[0].size/1e6)

#%%

a=np.zeros(zLCC1['Data'].shape,dtype='int8');
#ind=np.where( (zLCC1['Data']==1) & (z['lc2']['Data']!=meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_2']['T']) )
#a[ind]=1

ind=np.where( (z['lc2']['Data']!=meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_2']['T']) & \
             (z['lcc_cec_2020']['Data']==meta['LUT']['Derived']['lcc_cec_c']['Forest']) & \
             (np.isin(z['lcc_ntem_2019']['Data'],[meta['LUT']['Derived']['lcc_ntem']['Coniferous'],meta['LUT']['Derived']['lcc_ntem']['Broadleaf'],meta['LUT']['Derived']['lcc_ntem']['Mixedwood'],meta['LUT']['Derived']['lcc_ntem']['Wetland-treed']])==True) )
a[ind]=1

plt.close('all')
plt.matshow(a)

labL=list(meta['LUT']['Derived']['lcc_cec_c'])
d={}
for k,i in meta['LUT']['Derived']['lcc_cec_c'].items():
    ind=np.where( (z['lcc_cec_2020']['Data']==i) & (z['lc2']['Data']!=meta['LUT']['VEG_COMP_LYR_R1_POLY']['BCLCS_LEVEL_2']['T']) & (np.isin(z['lcc_ntem_2019']['Data'],[meta['LUT']['Derived']['lcc_ntem']['Coniferous'],meta['LUT']['Derived']['lcc_ntem']['Broadleaf'],meta['LUT']['Derived']['lcc_ntem']['Mixedwood'],meta['LUT']['Derived']['lcc_ntem']['Wetland-treed']])==True) )
    d[k]=ind[0].size/1e6






#%% Change analysis

zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])

vList=['rd','lcc_cec_2010','lcc_cec_2020','harv_yr_con1','fire_yr','ibm_yr'] #'lcc1_c','gfcly','gfcly_filt',
z=u1ha.Import_Raster(meta,[],vList)

#%% Land use change analysis - BC wide

id=np.array(list(meta['LUT']['Derived']['lcc_cec_c'].values()))
# Areas at start and end
v=np.zeros( (2,id.size) )
for i in range(id.size):
    ind=np.where( (z['lcc_cec_2010']['Data']==id[i]) & (z['harv_yr_con1']['Data']==0) & (z['fire_yr']['Data']==0) & (z['ibm_yr']['Data']==0) )
    v[0,i]=ind[0].size
    ind=np.where( (z['lcc_cec_2020']['Data']==id[i]) & (z['harv_yr_con1']['Data']==0) & (z['fire_yr']['Data']==0) & (z['ibm_yr']['Data']==0) )
    v[1,i]=ind[0].size
df=pd.DataFrame(v)
df.to_excel(r'C:\Users\rhember\OneDrive - Government of BC\Manuscripts\Land Use Change Assessment\CEC Change Analysis\BC Wide\cec_c0_nodist.xlsx')

# Transitions
v=np.zeros( (id.size,id.size) )
for i in range(id.size):
    for j in range(id.size):
        ind=np.where( (z['lcc_cec_2010']['Data']==id[i]) & (z['lcc_cec_2020']['Data']==id[j]) & (z['harv_yr_con1']['Data']==0) & (z['fire_yr']['Data']==0) & (z['ibm_yr']['Data']==0) )
        v[i,j]=ind[0].size
df=pd.DataFrame(v)
df.to_excel(r'C:\Users\rhember\OneDrive - Government of BC\Manuscripts\Land Use Change Assessment\CEC Change Analysis\BC Wide\cec_c1_nodist.xlsx')

#%% Filter by regional district

ind=np.where(z['rd']['Data']==meta['LUT']['TA_REGIONAL_DISTRICTS_SVW']['REGIONAL_DISTRICT_NAME']['CAPITAL'])
for k in z.keys():
    z[k]['Data']=z[k]['Data'][ind]

#%% Land use change analysis - by regional district (based on compressed categories)

id=np.array(list(meta['LUT']['Derived']['lcc_cec_c'].values()))

rs={}
for k in meta['LUT']['TA_REGIONAL_DISTRICTS_SVW']['REGIONAL_DISTRICT_NAME'].keys():
    # Areas at start and end
    vI=np.zeros( (2,id.size) )
    for i in range(id.size):
        ind=np.where( (z['lcc_cec_2010']['Data']==id[i]) & (z['harv_yr_con1']['Data']==0) & (z['fire_yr']['Data']==0) & (z['ibm_yr']['Data']==0) & (z['rd']['Data']==meta['LUT']['TA_REGIONAL_DISTRICTS_SVW']['REGIONAL_DISTRICT_NAME'][k]) )
        vI[0,i]=ind[0].size
        ind=np.where( (z['lcc_cec_2020']['Data']==id[i]) & (z['harv_yr_con1']['Data']==0) & (z['fire_yr']['Data']==0) & (z['ibm_yr']['Data']==0) & (z['rd']['Data']==meta['LUT']['TA_REGIONAL_DISTRICTS_SVW']['REGIONAL_DISTRICT_NAME'][k]) )
        vI[1,i]=ind[0].size
    df=pd.DataFrame(vI)
    df.to_excel(r'C:\Users\rhember\OneDrive - Government of BC\Manuscripts\Land Use Change Assessment\CEC Change Analysis\cec_c0_nodist_' + k + '.xlsx')

    # Transitions
    vT=np.zeros( (id.size,id.size) )
    for i in range(id.size):
        for j in range(id.size):
            ind=np.where( (z['lcc_cec_2010']['Data']==id[i]) & (z['lcc_cec_2020']['Data']==id[j]) & (z['harv_yr_con1']['Data']==0) & (z['fire_yr']['Data']==0) & (z['ibm_yr']['Data']==0) & (z['rd']['Data']==meta['LUT']['TA_REGIONAL_DISTRICTS_SVW']['REGIONAL_DISTRICT_NAME'][k]) )
            vT[i,j]=ind[0].size
    df=pd.DataFrame(vT)
    df.to_excel(r'C:\Users\rhember\OneDrive - Government of BC\Manuscripts\Land Use Change Assessment\CEC Change Analysis\cec_c1_nodist_' + k + '.xlsx')

    # Net deforestation summary
    ind=np.array([meta['LUT']['Derived']['lcc_cec_c']['Cropland'],meta['LUT']['Derived']['lcc_cec_c']['Barren Ground'],meta['LUT']['Derived']['lcc_cec_c']['Urban']],dtype='int8')
    Ad=np.sum(vT[0,ind-1])
    Aa=np.sum(vT[ind-1,0])
    Pd=np.sum(vT[0,ind-1])/vI[0,0]*100
    Pa=np.sum(vT[ind-1,0])/vI[0,0]*100
    Pn=Pa-Pd

# ind=np.where( (z['harv_yr_con1']['Data']>0) )
# z['Data'][ind]=2
# plt.matshow(z['Data'])

# ind=np.where(z['Data']==1)
# ind[1].size/1000000

#%% Map Forest to Grassland and Forest to Shrubland
# *** Way too much to be real! ***

z=zRef.copy()
z['Data']=np.zeros(zRef['Data'].shape,dtype='int8')
ind=np.where( (z['lcc_cec_2010']['Data']==meta['LUT']['Forest']) & (z['lcc_cec_2020']['Data']==meta['LUT']['Grassland']) )
z['Data'][ind]=1
ind=np.where( (z['lcc_cec_2010']['Data']==meta['LUT']['Forest']) & (z['lcc_cec_2020']['Data']==meta['LUT']['Shrubland']) )
z['Data'][ind]=1
ind=np.where( (z['harv_yr_con1']['Data']>0) )
z['Data'][ind]=2
plt.matshow(z['Data'])




