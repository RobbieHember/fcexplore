#%% Import modules
import os
import numpy as np
import time
import copy
import matplotlib.pyplot as plt
import pandas as pd
import pyproj
import fiona
import cv2
import geopandas as gpd
import scipy.io as spio
#from rasterio.transform import from_origin
#from shapely.geometry import Point,Polygon
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.macgyver.util_query_gdb as qgdb
import fcgadgets.bc1ha.bc1ha_util as u1ha
gp=gu.SetGraphics('Manuscript')

#%% Import paths and look-up-tables
meta=u1ha.Init()
zRef=gis.OpenGeoTiff(meta['Paths']['bc1ha Ref Grid'])

#%% Import data
vList=['refg','bsr_sc','spc1_vri23','spc1_pct_vri23']
z0=u1ha.Import_Raster(meta,[],vList,'Extract Grid')

#%% Buffer fire areas
zF=np.zeros(z0['refg'].shape,dtype='int16')
ind=np.where( (z0['bsr_sc']>=1) & (z0['bsr_sc']<=4) ); zF[ind]=1
#plt.matshow(zF,clim=[0,1])
kernel=np.ones((5,5),np.uint8)
zFB=cv2.dilate(zF,kernel,iterations=1)
zFB[zF==1]=2
#plt.matshow(zFB)

z1=copy.deepcopy(zRef)
z1['Data']=np.zeros(z0['refg'].shape,dtype='int16')
ind=np.where( (z0['spc1_vri23']==meta['LUT']['VEG_COMP_LYR_R1_POLY']['SPECIES_CD_1']['PL']) & (z0['spc1_pct_vri23']>80) & (zFB>0) )
z1['Data'][ind]=1
zFB[ind]=3

gdf=gis.DigitizeBinaryMask(z1)
gdf['Area']=gdf.geometry.area.values/10000
iKeep=np.where( (gdf['Area']>10) & (gdf['Area']<30) )[0]
xv=zRef['X'][0,:]
yv=zRef['Y'][:,0]
nt=50

#%%
tct={'h':np.zeros((iKeep.size,2)),'v':np.zeros((iKeep.size,2))}
cnt=0
for ikp in iKeep:
	g=gdf.iloc[ikp].geometry
	xc=g.representative_point().x
	yc=g.representative_point().y
	dx=np.abs(xv-xc)
	dy=np.abs(yv-yc)
	ix=np.where(dx==np.min(dx))[0][0]
	iy=np.where(dy==np.min(dy))[0][0]
	
	#z2=zFB.copy()
	#z2[iy,ix]=4
	for i in range(nt):
		if zFB[iy-i,ix]!=3:
			z2[iy-i,ix]=4
			tct['v'][cnt,0]=z0['bsr_sc'][iy-i,ix]
			break
	for i in range(nt):
		if zFB[iy+i,ix]!=3:
			z2[iy+i,ix]=4
			tct['v'][cnt,1]=z0['bsr_sc'][iy+i,ix]
			break
	for i in range(nt):
		if zFB[iy,ix-i]!=3:
			z2[iy,ix-i]=4
			tct['h'][cnt,0]=z0['bsr_sc'][iy,ix-i]
			break
	for i in range(nt):
		if zFB[iy,ix+i]!=3:
			z2[iy,ix+i]=4
			tct['h'][cnt,1]=z0['bsr_sc'][iy,ix+i]
			break
			#z0['bsr_sc'][ix-i,iy]
	cnt=cnt+1
#%%

ind1=np.where( (tct['v'][:,0]==4) & (tct['v'][:,1]==5) | (tct['v'][:,0]==4) & (tct['v'][:,1]==1) )[0]
ind2=np.where( (tct['h'][:,0]==4) & (tct['h'][:,1]==5) | (tct['h'][:,0]==4) & (tct['h'][:,1]==1) )[0]
#ind2=np.where( (tct['v'][:,0]==4) & (tct['v'][:,0]==5) | (tct['v'][:,0]==4) & (tct['v'][:,0]==1) )[0]
print(ind1.size/tct['v'].shape[0]*100)
print(ind2.size/tct['h'].shape[0]*100)


#plt.close('all'); plt.matshow(z2)
#print(tct)

#%%
w=100
z2=zFB.copy()
for i in range(0,1):
	z2[iy+i,ix-w:ix+w]=4
	z2[iy-w:iy+w,ix+i]=4
plt.close('all'); plt.matshow(z2)



