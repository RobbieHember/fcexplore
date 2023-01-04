
import sys
import numpy as np
import gdal
import osr
import matplotlib.pyplot as plt
import geopandas as gpd
import pandas as pd
import pyproj
from pyproj import Proj, transform
import rasterio
from scipy.interpolate import griddata
from scipy.io import loadmat
import time
import statsmodels.api as sm
import pickle
import gc
from sklearn.metrics import roc_curve, auc

sys.path.append(r'I:\My Drive\Code_Python')
from CustomFunctions.basic_functions import *
from CustomFunctions.gis import *


#------------------------------------------------------------------------------
# Import aerial survey data
#------------------------------------------------------------------------------

# Import geotiff
zND1=OpenGdal(r'I:\My Drive\PEST_SPECIES_CODE_ND_year1.tif')
zH=OpenGdal(r'I:\My Drive\h1.tif')

# Projection
iPND=Proj("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +no_defs +a=6378137 +rf=298.2572221010042 +to_meter=1")
oP=Proj(init='epsg:4326')

# Sample subgridde
ivl=10
zND1.Data=zND1.Data[0::ivl,0::ivl]
zND1.X=zND1.X[0::ivl,0::ivl]
zND1.Y=zND1.Y[0::ivl,0::ivl]
zND1.m,zND1.n=zND1.Data.shape

zH.Data=zH.Data[0::ivl,0::ivl]
zH.X=zH.X[0::ivl,0::ivl]
zH.Y=zH.Y[0::ivl,0::ivl]
zH.m,zH.n=zH.Data.shape

# Sampling 
nPos=zND1.Data[zND1.Data>0].size

MaskNull=np.zeros(zND1.Data.shape)
frqNull=30
MaskNull[0::frqNull,0::frqNull]=1
nNull=MaskNull[MaskNull==1].size
print(nNull/nPos)

iMask=np.where((MaskNull==1) & (zH.Data>2) & (zH.Data<60) | (zND1.Data>0) & (zH.Data>0) & (zH.Data<60))
L=iMask[0].size

zND1.Data=zND1.Data[iMask]
zH.Data=zH.Data[iMask]

# Species mix
#u=np.unique(zSpc.Data[iMask])
#n=np.zeros(u.size)
#for i in range(u.size):
#    ind=np.where(zSpc.Data[iMask]==u[i])[0]
#    n[i]=ind.size

#------------------------------------------------------------------------------
# Import fire
#------------------------------------------------------------------------------

zFire1=OpenGdal(r'I:\My Drive\PROT_HISTORICAL_FIRE_POLYS_SP_year1.tif')
zH=OpenGdal(r'I:\My Drive\h1.tif')

tmp=zFire1.Data[zFire1.Data>0]
tv=np.arange(1950,2019,1)
A=np.zeros(tv.shape)
for i in range(tv.size):
    ind=np.where(tmp==tv[i])[0]
    A[i]=ind.size
plt.bar(tv,A)

# Sample subgridde
ivl=10
zFire1.Data=zFire1.Data[0::ivl,0::ivl]
zFire1.X=zFire1.X[0::ivl,0::ivl]
zFire1.Y=zFire1.Y[0::ivl,0::ivl]
zFire1.m,zFire1.n=zFire1.Data.shape

zH.Data=zH.Data[0::ivl,0::ivl]
zH.X=zH.X[0::ivl,0::ivl]
zH.Y=zH.Y[0::ivl,0::ivl]
zH.m,zH.n=zH.Data.shape

# Sampling 
nPos=zFire1.Data[zFire1.Data>0].size

MaskNull=np.zeros(zFire1.Data.shape)
frqNull=5
MaskNull[0::frqNull,0::frqNull]=1
nNull=MaskNull[MaskNull==1].size
print(nNull/nPos)

iMask=np.where((MaskNull==1) & (zH.Data>2) & (zH.Data<60) | (zFire1.Data>0) & (zH.Data>0) & (zH.Data<60))
L=iMask[0].size

zFire1.Data=zFire1.Data[iMask]
zH.Data=zH.Data[iMask]






#------------------------------------------------------------------------------
# Convert to matrix
#------------------------------------------------------------------------------

#tv=np.arange(1952,2019,1)
tv=np.arange(1992,2019,1)

v={}
v['ND']=np.zeros((tv.size,L))
#v['IDW']=np.zeros((tv.size,L))
#v['Fire']=np.zeros((tv.size,L))
for j in range(L):
    it=np.where(tv==zND1.Data[j])[0]
    v['ND'][it,j]=1
    #it=np.where(tv==zIDW1.Data[j])[0]
    #v['IDW'][it,j]=1
    #it=np.where(tv==zFire1.Data[j])[0]
    #v['Fire'][it,j]=1
v['h']=np.tile(zH.Data,(tv.size,1))


#------------------------------------------------------------------------------
# Create a crosswalk between BC 1ha grid and NACID grid
#------------------------------------------------------------------------------

flg=0
if flg==1:
    # Import BC 1ha grid
    ProjBC=Proj("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +no_defs +a=6378137 +rf=298.2572221010042 +to_meter=1")

    # Import NACID grid
    zNA=OpenGdal(r'E:\Data\Climate\NACID\Geotiff\NACID\NACID_ws_tmw_gs_abso_comp_hist_v1\NACID_ws_tmw_gs_abso_comp_hist_v1_2018.tif')
    ProjNA=Proj("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=40 +lon_0=-100 +x_0=0 +y_0=0 +no_defs +a=6378137 +rf=298.2572221010042 +to_meter=1")

    # Convert NACID coordinates to BC 1ha project
    xNA,yNA=pyproj.transform(ProjNA,ProjBC,zNA.X,zNA.Y)

    # Create index for NACID grid
    idxNA=np.reshape(np.arange(0,zNA.m*zNA.n,1),zNA.Data.shape,order='C')

    # Nearest-neighbour interpolation of the NACID index to the BC grid 
    # This took 6 min
    x0=xNA.flatten()
    y0=yNA.flatten()
    z0=idxNA.flatten()
    t0=time.time()
    idxBC=griddata((x0,y0),z0,(zND1.X,zND1.Y),method='nearest')
    t1=time.time(); print(t1-t0)

    # Re-projected NACID grid
    #zNA_new=zNA.Data[np.unravel_index(idxBC,zNA.Data.shape,'C')]
    #plt.close('all'); plt.imshow(zNA_new)

    # Save
    #np.save(r'I:\My Drive\idxBC',idxBC)
    np.save(r'I:\My Drive\idxBC_SubSamp10',idxBC)

# Load crosswalk
#idxBC=np.load(r'I:\My Drive\idxBC.npy')
idxBC=np.load(r'I:\My Drive\idxBC_SubSamp10.npy')


#------------------------------------------------------------------------------
# Import monthly climate data
#------------------------------------------------------------------------------

tvm=tvec.tvec('m',tv[0]-1,2019)
vm={}
vm['tmean_a']=np.nan*np.ones((tvm.shape[0],L))
vm['prcp_a']=np.nan*np.ones((tvm.shape[0],L))
vm['etp_a']=np.nan*np.ones((tvm.shape[0],L))
vm['cwd_a']=np.nan*np.ones((tvm.shape[0],L))
vm['ws_a']=np.nan*np.ones((tvm.shape[0],L))
for mo in range(0,12):
    
    z=loadmat(r'E:\Data\Climate\NACID\Grids\NACID_etp_tmw_mon_norm_1971to2000_comp_hist_v1\NACID_etp_tmw_mon_norm_1971to2000_comp_hist_v1_' + str(mo+1) + '.mat')    
    mu=z['z'][0][0][6].flatten()[idxBC]
    sf=0.01
    mu=np.reshape(mu,(zFire1.m,zFire1.n),order='C')
    etp_mu=mu[iMask].astype(float)*sf
    
    z=loadmat(r'E:\Data\Climate\NACID\Grids\NACID_eta_tmw_mon_norm_1971to2000_comp_hist_v1\NACID_eta_tmw_mon_norm_1971to2000_comp_hist_v1_' + str(mo+1) + '.mat')    
    mu=z['z'][0][0][6].flatten()[idxBC]
    sf=0.01
    mu=np.reshape(mu,(zFire1.m,zFire1.n),order='C')
    eta_mu=mu[iMask].astype(float)*sf
    
    cwd_mu=etp_mu-eta_mu
    
    z=loadmat(r'E:\Data\Climate\NACID\Grids\NACID_ws_tmw_mon_norm_1971to2000_comp_hist_v1\NACID_ws_tmw_mon_norm_1971to2000_comp_hist_v1_' + str(mo+1) + '.mat')    
    mu=z['z'][0][0][6].flatten()[idxBC]
    sf=0.01
    mu=np.reshape(mu,(zFire1.m,zFire1.n),order='C')
    ws_mu=mu[iMask].astype(float)*sf
        
    for yr in range(1951,2019,1):
        print(yr,mo+1)
        try:
            z=loadmat(r'E:\Data\Climate\NACID\Grids\NACID_etp_tmw_mon_abso_comp_hist_v1\NACID_etp_tmw_mon_abso_comp_hist_v1_' + str(yr) + '_' + str(mo+1) + '.mat')    
            data=z['z'][0][0][6].flatten()[idxBC]
            sf=0.01
            data=np.reshape(data,(zFire1.m,zFire1.n),order='C')
            etp=data[iMask].astype(float)*sf
        
            z=loadmat(r'E:\Data\Climate\NACID\Grids\NACID_eta_tmw_mon_abso_comp_hist_v1\NACID_eta_tmw_mon_abso_comp_hist_v1_' + str(yr) + '_' + str(mo+1) + '.mat')    
            data=z['z'][0][0][6].flatten()[idxBC]
            sf=0.01
            data=np.reshape(data,(zFire1.m,zFire1.n),order='C')
            eta=data[iMask].astype(float)*sf
        
            cwd=etp-eta
            
            z=loadmat(r'E:\Data\Climate\NACID\Grids\NACID_ws_tmw_mon_abso_comp_hist_v1\NACID_ws_tmw_mon_abso_comp_hist_v1_' + str(yr) + '_' + str(mo+1) + '.mat')    
            data=z['z'][0][0][6].flatten()[idxBC]
            sf=0.01
            data=np.reshape(data,(zFire1.m,zFire1.n),order='C')
            ws=data[iMask].astype(float)*sf
        
            z=loadmat(r'E:\Data\Climate\NACID\Grids\NACID_prcp_mon_anom_si_hist_v1\NACID_prcp_mon_anom_si_hist_v1_' + str(yr) + '_' + str(mo+1) + '.mat')    
            data=z['z'][0][0][0].flatten()[idxBC]
            sf=1
            data=np.reshape(data,(zFire1.m,zFire1.n),order='C')
            prcp_a=data[iMask].astype(float)*sf
        
            z=loadmat(r'E:\Data\Climate\NACID\Grids\NACID_tmean_mon_anom_si_hist_v1\NACID_tmean_mon_anom_si_hist_v1_' + str(yr) + '_' + str(mo+1) + '.mat')    
            try:
                data=z['z'][0][0][0].flatten()[idxBC]
            except:
                data=z['z'][0][0][1].flatten()[idxBC]
            sf=0.1
            data=np.reshape(data,(zFire1.m,zFire1.n),order='C')
            tmean_a=data[iMask].astype(float)*sf
        
            it=np.where((tvm[:,0]==yr) & (tvm[:,1]==mo))[0]
            vm['tmean_a'][it,:]=tmean_a
            vm['prcp_a'][it,:]=prcp_a
            vm['etp_a'][it,:]=etp-etp_mu
            vm['cwd_a'][it,:]=cwd-cwd_mu
            vm['ws_a'][it,:]=ws-ws_mu
            
        except:
            continue

vm['cmi_a']=vm['prcp_a']/30-vm['etp_a']

# Exclude outliers
ind=np.where((vm['tmean_a']<-10) | (vm['tmean_a']>10) | (vm['prcp_a']<-500))
vm['tmean_a'][ind]=np.nan
vm['prcp_a'][ind]=np.nan
vm['etp_a'][ind]=np.nan
vm['cmi_a'][ind]=np.nan
vm['cwd_a'][ind]=np.nan
vm['ws_a'][ind]=np.nan

# Save
fout=open(r'E:\Data\DisturbanceVsClimate\ND_Climate.pkl','wb')
pickle.dump(vm,fout); fout.close()

# Dependent variable
y=v['Fire'].reshape((L*tv.size,1))

# Time vectors
yr=np.concatenate((-1*np.ones((12,1)),(0*np.ones((9,1)))),axis=0)
mo=np.concatenate((np.arange(1,13,1),np.arange(1,10,1)),axis=0)



xv_str=['tmean_a','prcp_a','etp_a','cmi_a','cwd_a','ws_a']
ms=[]
for ix in range(len(xv_str)):
    d={}
    d['X Name']=xv_str[ix]
    d['prsq']=np.nan*np.ones(yr.shape)
    d['b']=np.nan*np.ones(yr.shape)
    d['pval']=np.nan*np.ones(yr.shape)
    d['aic']=np.nan*np.ones(yr.shape)
    d['auc']=np.nan*np.ones(yr.shape)
    for i in range(yr.size):
        print(ix,i)
        #it=np.where((tvm[:,0]>=1990+yr[i]) & (tvm[:,0]<=2018+yr[i]) & (tvm[:,1]==mo[i]))
        it=np.where((tvm[:,0]>=1952+yr[i]) & (tvm[:,0]<=2018+yr[i]) & (tvm[:,1]==mo[i]))
        x=vm[xv_str[ix]][it,:].reshape((L*tv.size,1))
        if np.sum(x)==0:
            x[0]=0.01
        x=(x-np.nanmean(x))/np.nanstd(x)        
        x=sm.add_constant(x)
        
        ikp=(np.isnan(x[:,1])==False) & (x[:,1]>-200) & (x[:,1]<200)
        y1=y[ikp]
        x1=x[ikp,:]
        if y1.size==0:
            continue
        sfreq=100
        logit=sm.Logit(y1[0::sfreq,0],x1[0::sfreq,:])
        mf=logit.fit()
        
        #yhat=mf.predict(x1)
        #fpr,tpr,thresholds=roc_curve(y1,yhat)
        #roc_auc=auc(fpr,tpr)
        
        d['aic'][i]=mf.aic
        d['prsq'][i]=mf.prsquared
        d['pval'][i]=mf.pvalues[1]
        d['b'][i]=mf.params[1]
        #d['auc'][i]=roc_auc
    ms.append(d.copy())




plt.close('all'); plt.plot(np.arange(-12,9,1),ms[0]['prsq'],'-o')
plt.close('all'); plt.plot(np.arange(-12,9,1),ms[0]['b'],'-o')

plt.close('all'); 
for i in range(len(xv_str)):
    plt.plot(np.arange(-12,9,1),ms[i]['b'],'-o',label=ms[i]['X Name'])
plt.legend()

plt.figure(2) 
for i in range(len(xv_str)):
    plt.plot(np.arange(-12,9,1),ms[i]['aic'],'-o',label=ms[i]['X Name'])
plt.legend()

plt.figure(3) 
for i in range(len(xv_str)):
    plt.plot(np.arange(-12,9,1),ms[i]['auc'],'-o',label=ms[i]['X Name'])
plt.legend()

plt.plot(np.arange(-12,9,1),ms[1]['b'],'-o')









zIDW1=OpenGdal(r'I:\My Drive\PEST_SPECIES_CODE_IDW_year1.tif')
tmp=zIDW1.Data[zIDW1.Data>0]
tv=np.arange(1950,2019,1)
A=np.zeros(tv.shape)
for i in range(tv.size):
    print(tv[i])
    ind=np.where(tmp==tv[i])[0]
    A[i]=ind.size
plt.bar(tv,A)








#------------------------------------------------------------------------------
# Import warm-season mean climate data
#------------------------------------------------------------------------------

tv2=np.arange(tv[0]-10,2019,1)
lg=np.arange(0,10+1,1)

v['cwd_a']=[]
v['ws_a']=[]
v['cwd_z']=[]
v['ws_z']=[]
for i in range(lg.size):
    v['cwd_a'].append(np.nan*np.ones((tv.size,L)))
    v['ws_a'].append(np.nan*np.ones((tv.size,L)))
    v['cwd_z'].append(np.nan*np.ones((tv.size,L)))
    v['ws_z'].append(np.nan*np.ones((tv.size,L)))
    

for i in range(tv2.size):  
    
    z=loadmat(r'E:\Data\Climate\NACID\Grids\NACID_cwd_tmw_gs_anom_comp_hist_v1\NACID_cwd_tmw_gs_anom_comp_hist_v1_' + str(tv2[i]) + '.mat')    
    data=z['z'][0][0][6].flatten()[idxBC]    
    data=np.reshape(data,(zND1.m,zND1.n),order='C')
    sf=0.01    
    for j in range(lg.size):
        it=np.where(tv==tv2[i]+lg[j])[0]
        if it.size!=0:
            v['cwd_a'][j][it,:]=data[iMask].astype(float)*sf

    z=loadmat(r'E:\Data\Climate\NACID\Grids\NACID_ws_tmw_gs_anom_comp_hist_v1\NACID_ws_tmw_gs_anom_comp_hist_v1_' + str(tv2[i]) + '.mat')    
    data=z['z'][0][0][6].flatten()[idxBC]    
    data=np.reshape(data,(zND1.m,zND1.n),order='C')    
    sf=1    
    for j in range(lg.size):
        it=np.where(tv==tv2[i]+lg[j])[0]
        if it.size!=0:
            v['ws_a'][j][it,:]=data[iMask].astype(float)*sf
    
    z=loadmat(r'E:\Data\Climate\NACID\Grids\NACID_cwd_tmw_gs_zsco_comp_hist_v1\NACID_cwd_tmw_gs_zsco_comp_hist_v1_' + str(tv2[i]) + '.mat')    
    data=z['z'][0][0][6].flatten()[idxBC]    
    data=np.reshape(data,(zND1.m,zND1.n),order='C')
    sf=0.01    
    for j in range(lg.size):
        it=np.where(tv==tv2[i]+lg[j])[0]
        if it.size!=0:
            v['cwd_z'][j][it,:]=data[iMask].astype(float)*sf
            
    z=loadmat(r'E:\Data\Climate\NACID\Grids\NACID_ws_tmw_gs_zsco_comp_hist_v1\NACID_ws_tmw_gs_zsco_comp_hist_v1_' + str(tv2[i]) + '.mat')    
    data=z['z'][0][0][6].flatten()[idxBC]    
    data=np.reshape(data,(zND1.m,zND1.n),order='C')
    sf=1  
    for j in range(lg.size):
        it=np.where(tv==tv2[i]+lg[j])[0]
        if it.size!=0:
            v['ws_z'][j][it,:]=data[iMask].astype(float)*sf
    
   


#------------------------------------------------------------------------------
# Satistical modelling
#------------------------------------------------------------------------------


xv_str=['cwd_a','ws_a','cwd_z','ws_z']
ms=[]
for ix in range(len(xv_str)):
    d={}
    d['X Name']=xv_str[ix]
    d['prsq']=np.nan*np.ones(lg.shape)
    d['b']=np.nan*np.ones(lg.shape)
    d['pval']=np.nan*np.ones(lg.shape)
    d['aic']=np.nan*np.ones(lg.shape)
    d['auc']=np.nan*np.ones(lg.shape)
    for i in range(lg.size):
        
        y0=v['ND'].reshape((L*tv.size,1))
        x0=v[xv_str[ix]][i].reshape((L*tv.size,1))
                       
        ikp=(x0>-200) & (x0<200)
        y1=y0[ikp]
        x1=x0[ikp]
        
        x1=(x1-np.nanmean(x1))/np.nanstd(x1)
        x1=sm.add_constant(x1)
        
        
        sfreq=1
        logit=sm.Logit(y1[0::sfreq],x1[0::sfreq,:])
        
        mf=logit.fit()
        
        #yhat=mf.predict(x1)
        #fpr,tpr,thresholds=roc_curve(y1,yhat)
        #roc_auc=auc(fpr,tpr)
        
        d['aic'][i]=mf.aic
        d['prsq'][i]=mf.prsquared
        d['pval'][i]=mf.pvalues[1]
        d['b'][i]=mf.params[1]
        #d['auc'][i]=roc_auc
    ms.append(d.copy())





plt.close('all'); 
for i in range(len(xv_str)):
    plt.plot(np.arange(-10,1,1),ms[i]['aic'],'-o',label=ms[i]['X Name'])
plt.legend()


plt.figure(2); plt.plot(np.arange(-10,1,1),ms[2]['b'],'-o')


x=v['cwd_a0'].reshape((L*tv.size,1))
y=v['ND'].reshape((L*tv.size,1))

bw0=0.1
bin0=np.arange(-2,2,bw0)
y0=np.zeros(bin0.shape)
n0=np.zeros(bin0.shape)
for i in range(bin0.size):
    ind=np.where(np.abs(x-bin0[i])<=bw0/2)[0]
    y0[i]=np.mean(y[ind])
    n0[i]=ind.size
plt.close('all');plt.plot(bin0,y0,'-o')
#plt.close('all'); plt.plot(bin0,n0,'o')




y=v['ND'].reshape((L*tv.size,1))
x=v['cwd_a0'].reshape((L*tv.size,1))
x=sm.add_constant(x)
ikp=np.where(np.abs(x)<2.2)[0]
y=y[ikp]
x=x[ikp,:]
logit=sm.Logit(y[0::10,0],x[0::10,:])
r=logit.fit()
print(r.aic)
r.summary()
r.conf_int()

# Plot the relationship
xhat=np.arange(-3,3,0.01)
xhat=sm.add_constant(xhat)
yhat=r.predict(xhat)

plt.close('all')
fig,ax=plt.subplots(figsize=(9,9))
ax.plot(xhat[:,1],yhat,'r-')
ax.set(xlim=[1100000,1600000],ylim=[450000,860000],position=[0.05,0.03,0.90,0.94],xlabel='x axis')
ax.tick_params(axis='both',which='major',labelsize=6)
plt.legend(loc='upper right')



from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression



x_train,x_test,y_train,y_test=train_test_split(x,y,test_size=0.25,random_state=0)

lr=LogisticRegression()
lr.fit(x_train,y_train)
lr.coef_

yhat_test=lr.predict(x_test)






