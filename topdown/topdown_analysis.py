'''
ANALYSIS OF TOP-DOWN GHG EMISSIONS FROM BYRNE ET AL 2023
'''
#%% Import modules
import numpy as np
import copy
import matplotlib.pyplot as plt
import geopandas as gpd
import pandas as pd
from shapely.geometry import Polygon,Point,box,shape
from rasterio import features
import rasterio
from rasterio.transform import from_origin
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.bc1ha.bc1ha_util as u1ha
gp=gu.SetGraphics('Manuscript')

#%% Import data

nc=gu.ReadNC(r'C:\Data\TopDown\pilot_topdown_CO2_Budget_grid_v1.nc')
list(nc.keys())

gdf_pb=gpd.read_file(r'C:\Data\Geodatabases\North America\bound_p.shp')
u_pb=gdf_pb['STATEABB'].unique()
u_pb=u_pb[u_pb!=None]
gdf_pb['ID']=np.zeros(gdf_pb['STATEABB'].size)
cnt=1
lut_pb={}
for u in u_pb:
	ind=np.where(gdf_pb['STATEABB']==u)[0]
	gdf_pb['ID'][ind]=cnt
	lut_pb[u]=cnt
	cnt=cnt+1

gdf_ez=gpd.read_file(r'C:\Data\Ecozones\nef_ca_ter_ecozone_v2_2.geojson')
uez=gdf_ez['ECOZONE_NAME_EN'].unique()
gdf_ez['ID']=np.zeros(gdf_ez['ECOZONE_NAME_EN'].size)
cnt=1
lut_ez={}
for iU in uez:
	if iU==None: continue
	ind=np.where(gdf_ez['ECOZONE_NAME_EN']==iU)[0]
	gdf_ez['ID'][ind]=cnt
	lut_ez[iU]=cnt
	cnt=cnt+1

# BC
gdf_bc=gpd.read_file(r'C:\Data\Geodatabases\North America\bound_p.shp')
gdf_bc=gdf_bc[gdf_bc['STATEABB']=='CA-BC']
gdf_bc=gdf_bc.reset_index()
gdf_bc['ID']=np.ones(gdf_bc['STATEABB'].size)

flg=1
if flg==1:
	gdf_glob=gpd.read_file(r'C:\Data\Geodatabases\Global\World_Countries_Generalized.geojson')
	u_pb=gdf_glob['COUNTRY'].unique()
	u_pb=u_pb[u_pb!=None]
	gdf_glob['ID']=np.zeros(gdf_glob['COUNTRY'].size)
	cnt=1
	lut_glob={}
	for u in u_pb:
		ind=np.where(gdf_glob['COUNTRY']==u)[0]
		gdf_glob['ID'][ind]=cnt
		lut_glob[u]=cnt
		cnt=cnt+1

#%% Rasterize

# Reference grid
#west=-180; north=90
west=-179.5; north=89.5
z=np.zeros((nc['lat'].size,nc['lon'].size),dtype='int16')
with rasterio.open(
	'C:\Data\TopDown\grid.tif',
	mode="w",
	driver="GTiff",
	height=z.shape[0],
	width=z.shape[1],
	count=1,
	dtype=z.dtype,
	crs=gdf_pb.crs,
	transform=from_origin(west,north,1,1), # Inputs (west,north,xsize,ysize)
	) as new_dataset:
	new_dataset.write(z,1)

zRef=gis.OpenGeoTiff(r'C:\Data\TopDown\grid.tif')
shapes=((geom,value) for geom, value in zip(gdf_pb['geometry'],gdf_pb['ID']))
z1=copy.deepcopy(zRef)
z1['Data']=np.zeros(zRef['Data'].shape,dtype=float)
features.rasterize(shapes=shapes,fill=0,out=z1['Data'],transform=zRef['Transform'])
#plt.close('all'); plt.matshow(pbg,clim=[0,100])
z1['Data']=z1['Data'].astype('int16')
gis.SaveGeoTiff(z1,'C:\Data\TopDown\PB2.tif')

zRef=gis.OpenGeoTiff(r'C:\Data\TopDown\grid.tif')
shapes=((geom,value) for geom, value in zip(gdf_ez['geometry'],gdf_ez['ID']))
z1=copy.deepcopy(zRef)
z1['Data']=np.zeros(zRef['Data'].shape,dtype=float)
features.rasterize(shapes=shapes,fill=0,out=z1['Data'],transform=zRef['Transform'])
#plt.close('all'); plt.matshow(z1['Data'])
z1['Data']=z1['Data'].astype('int16')
gis.SaveGeoTiff(z1,'C:\Data\TopDown\EZ.tif')

zRef=gis.OpenGeoTiff(r'C:\Data\TopDown\grid.tif')
shapes=((geom,value) for geom, value in zip(gdf_bc['geometry'].geometry.buffer(1),gdf_bc['ID']))
z1=copy.deepcopy(zRef)
z1['Data']=np.zeros(zRef['Data'].shape,dtype=float)
features.rasterize(shapes=shapes,fill=0,out=z1['Data'],transform=zRef['Transform'])
plt.close('all'); plt.matshow(z1['Data'],clim=[0,2])
z1['Data']=z1['Data'].astype('int16')
gis.SaveGeoTiff(z1,'C:\Data\TopDown\BC.tif')

zRef=gis.OpenGeoTiff(r'C:\Data\TopDown\grid.tif')
shapes=((geom,value) for geom, value in zip(gdf_glob['geometry'],gdf_glob['ID']))
z1=copy.deepcopy(zRef)
z1['Data']=np.zeros(zRef['Data'].shape,dtype=float)
features.rasterize(shapes=shapes,fill=0,out=z1['Data'],transform=zRef['Transform'])
plt.close('all'); plt.matshow(z1['Data'],clim=[0,200])
z1['Data']=z1['Data'].astype('int16')
gis.SaveGeoTiff(z1,'C:\Data\TopDown\GLOB.tif')

vL=['area','LNLGIS_NCE_meanYear','LNLGIS_dC_loss_meanYear','LNLGIS_NBE_meanYear','FF_meanYear','Crop_meanYear','Wood_meanYear','River_meanYear']
for v in vL:
	z1=copy.deepcopy(zRef)
	z1['Data']=np.flip(nc[v],axis=0).astype('float32')
	pth=r'C:\Data\TopDown' + '\\' + v + '.tif'
	gis.SaveGeoTiff(z1,pth)
	#gis.ResampleRaster(pth,sf)

#%% Open grids
vL=['area','BC','PB2','EZ','LNLGIS_NCE_meanYear','LNLGIS_dC_loss_meanYear','LNLGIS_NBE_meanYear','FF_meanYear','Crop_meanYear','Wood_meanYear','River_meanYear']
grd={}
for v in vL:
	grd[v]=gis.OpenGeoTiff(r'C:\Data\TopDown' + '\\' + v + '.tif')['Data']

#%% Global fluxes
flg=0
if flg==1:
	vL=['LNLGIS_NCE_meanYear','LNLGIS_dC_loss_meanYear','LNLGIS_NBE_meanYear','FF_meanYear','Crop_meanYear','Wood_meanYear','River_meanYear']
	for v in vL:
		print(np.sum(grd[v]*grd['area'])/1e15)

#%% Map

import matplotlib.colors
zMask=gis.OpenGeoTiff('C:\Data\TopDown\GLOB.tif')
z0=grd['Wood_meanYear']*grd['area']/1e6/1e4
bw=10; bin=np.arange(-100,110,bw)
lab=bin.astype(str)
z1=N_vis*np.ones(z0.shape,dtype='int16')
for i in range(bin.size):
	ind=np.where(np.abs(z0-bin[i])<=bw/2)
	z1[ind]=i
ind=np.where(z0<bin[0]); z1[ind]=1
ind=np.where(z0>bin[-1]); z1[ind]=i
ind=np.where(zMask['Data']==0); z1[ind]=i+1
for i in range(bin.size):
	z1[0,i]=i
lab=np.append(lab,'')
N_vis=bin.size
N_hidden=1
N_tot=N_vis+N_hidden
#cm=plt.cm.get_cmap('viridis',N_vis)
#cm=plt.cm.get_cmap('prism',N_vis)
#cm=matplotlib.colormaps['prism']
cm0=plt.cm.get_cmap("jet",N_vis)
cm=np.vstack( (cm0(0),cm0(1),cm0(2),cm0(3),cm0(4),cm0(5),cm0(6),cm0(7),cm0(8),cm0(9),cm0(10),cm0(11),cm0(12),cm0(13),cm0(14),cm0(15),cm0(16),cm0(17),cm0(18),cm0(19),cm0(20)) )
cm[10,:]=[0.85,0.85,0.85,1]
cm=np.vstack( (cm,(1,1,1,1)) )
cm=matplotlib.colors.ListedColormap(cm)

plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(24,11))
im=ax[0].matshow(z1,extent=zRef['Extent'],cmap=cm)
gdf_pb.plot(ax=ax[0],facecolor='None',linewidth=0.25,edgecolor='k')
ax[0].set(position=[0.02,0.01,0.96,0.98],xlim=[-180,180],ylim=[-75,75])
ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both');# ax[0].grid(Off); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])
ax[0].text(-172,-68,'Lateral transfers of carbon due to trade of\nwood products, 2015-2020\nUnits: tCO$_2$e/ha/year',fontsize=8,fontweight='normal')
ax[0].text(-130,5,'Positive:\nnet import',fontsize=8,fontweight='normal',va='center',fontstyle='italic')
ax[0].text(-130,-12,'Negative:\nnet export',fontsize=8,fontweight='normal',va='center',fontstyle='italic')
zmn=np.min(z1); zmx=np.max(z1); cb_ivl=(zmx-zmn)/N_tot; cb_bnd=np.arange(zmn,zmx+cb_ivl-N_hidden*cb_ivl,cb_ivl)
cb_ticks=np.arange(zmn+cb_ivl/2,N_tot-1,cb_ivl)
cb=plt.colorbar(im,cax=ax[1],cmap=cm,boundaries=cb_bnd,ticks=cb_ticks)
ax[1].set(position=[0.71,0.6,0.05,0.14])
cb.ax.set(yticklabels=lab)
cb.ax.tick_params(labelsize=8,length=0)
cb.outline.set_edgecolor('w')
for i in range(cb_bnd.size):
	ax[1].plot([0,100],[cb_bnd[i],cb_bnd[i]],'w-',linewidth=2)
ax[1].set(position=[0.1,0.23,0.0175,0.5])
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\TopDown\Global_WoodFluxCO2e_Mean2015to2020','png',900)


#%%
plt.matshow(grd['PB2'])

#%% Analysis BC
vL=['LNLGIS_NCE_meanYear','FF_meanYear','LNLGIS_dC_loss_meanYear','LNLGIS_NBE_meanYear','Crop_meanYear','Wood_meanYear','River_meanYear']
d={}
ind=np.where(grd['BC']==1)
f=grd['FF_meanYear'][ind]
A=grd['area'][ind]
print(np.sum(A)/1e6/10000)
print(np.sum(f*A)/1e12) # MtCO2e/yr
d={}
dMu={}
for v in vL:
	dMu[v]=np.nanmean(grd[v][ind])/1e6*1e4 # tCO2e/ha/yr
	d[v]=np.sum(np.sum(grd[v][ind]*grd['area'][ind]))/1e12 # MtCO2e/yr

y=np.array(list(d.values()))
xlabs=['Net carbon\nexchange\n(NBE + FF)','Fossil fuel\nemissions\n(FF)','Net change in\nbiogenic\ncarbon\n($\Delta$C = NBE\n - crops - wood\n - river)',
	   'Net biosphere\nexchange\n(NBE)\n(biomass\nburning +\ndecay - uptake)','Lateral flux\ndue to\ntrade of crops\n(crops)','Lateral flux\ndue to\ntrade of wood\n(wood)','Lateral flux\ndue to\nriver export\n(river)']
plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(14,8));
ax.plot([-2,10],[0,0],'k-')
ax.bar(np.arange(y.size),y,facecolor=[0.7,0.8,1])
ax.text(1,90,'Provincial\ninventory',ha='center',va='center',fontsize=7,fontweight='normal',color=[0,0.4,0])
ax.plot(1,64,'s',ms=4,color=[0,0.4,0])
ax.text(2,-75,'Net biomass\nproduction\nof trees\n(CMI network)',ha='center',va='center',fontsize=7,fontweight='normal',color=[0,0.4,0])
ax.plot(2,-115,'s',ms=4,color=[0,0.4,0])
ax.text(3,60,'Provincial\ninventory',ha='center',va='center',fontsize=8,fontweight='normal',color=[0,0.4,0])
ax.plot(3,35,'s',ms=4,color=[0,0.4,0])
ax.text(5,-60,'Harvest\nmortality\n(CMI network)',ha='center',va='center',fontsize=8,fontweight='normal',color=[0,0.4,0])
ax.plot(5,-25,'s',ms=4,color=[0,0.4,0])
ax.set(xticks=np.arange(y.size),xticklabels=xlabs,ylabel='GHG flux (MtCO2e/yr)',yticks=np.arange(-250,200,25),ylim=[-225,125],xlim=[-0.5,y.size-0.5])
ax.legend(loc='lower left',facecolor=[1,1,1],frameon=False);
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')
plt.tight_layout()
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\TopDown\BC_AnnualFluxCO2e_Mean2015to2020','png',900)

y=np.array(list(dMu.values()))
xlabs=['Net carbon\nexchange\n(NBE + F)','Fossil fuel\nemissions\n(FF)','Net change in\nbiogenic\ncarbon\n($\Delta$C = NBE\n - crops - wood\n - river)',
	   'Net biosphere\nexchange\n(NBE)\n(biomass\nburning +\ndecay - uptake)','Lateral flux\ndue to\ntrade of crops\n(crops)','Lateral flux\ndue to\ntrade of wood\n(wood)','Lateral flux\ndue to\nriver export\n(river)']
plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(14,8));
ax.plot([-2,10],[0,0],'k-')
ax.bar(np.arange(y.size),y,facecolor=[0.7,0.8,1])
#ax.text(1,90,'Provincial\ninventory',ha='center',va='center',fontsize=7,fontweight='normal',color=[0,0.4,0])
ax.plot(1,64,'s',ms=4,color=[0,0.4,0])
#ax.text(2,-75,'Net biomass\nproduction\nof trees\n(CMI network)',ha='center',va='center',fontsize=7,fontweight='normal',color=[0,0.4,0])
ax.plot(2,-115,'s',ms=4,color=[0,0.4,0])
#ax.text(3,60,'Provincial\ninventory',ha='center',va='center',fontsize=8,fontweight='normal',color=[0,0.4,0])
ax.plot(3,35,'s',ms=4,color=[0,0.4,0])
#ax.text(5,-60,'Harvest\nmortality\n(CMI network)',ha='center',va='center',fontsize=8,fontweight='normal',color=[0,0.4,0])
ax.plot(5,-25,'s',ms=4,color=[0,0.4,0])
ax.set(xticks=np.arange(y.size),xticklabels=xlabs,ylabel='GHG flux (tCO2e/ha/yr)',yticks=np.arange(-3,3,0.2),ylim=[-1.6,0.6],xlim=[-0.5,y.size-0.5])
plt.tight_layout()
#ax.legend(loc='lower left',facecolor=[1,1,1],frameon=False);
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\TopDown\BC_AnnualFlux_tCO2ePerHectarePerYear_Mean2015to2020','png',900)

# Time series
#vL2=['LNLGIS_dC_loss','LNLGIS_NBE','FF','Crop','Wood','River']

#%% Analysis by ecozone
vL=['LNLGIS_NCE_meanYear','LNLGIS_dC_loss_meanYear','LNLGIS_NBE_meanYear','FF_meanYear','Crop_meanYear','Wood_meanYear','River_meanYear']
#vL2=['LNLGIS_dC_loss','LNLGIS_NBE','FF','Crop','Wood','River']
d={}
dts={}
for u in uez:
	ind=np.where(grd['EZ']==lut_ez[u])
	#v='FF_meanYear'
	#f=nc[v][ind]
	#A=nc['area'][ind]
	#print(np.sum(A)/1e6/10000)
	#print(np.sum(f*A)/1e12) # MtCO2e/yr

	d[u]={}
	for v in vL:
		#d[u][v]=np.nanmean(f)/1e6*1e4 # tCO2e/ha/yr
		d[u][v]=np.sum(np.sum(grd[v][ind]*grd['area'][ind]))/1e12 # MtCO2e/yr
# 	dts[u]={}
# 	for v in vL2:
# 		dts[u][v]=np.zeros(6)
# 		for iT in range(6):
# 			#dts[u][v][iT]=np.mean(nc[v][iT,:,:][ind])/1e6*1e4
# 			dts[u][v][iT]=np.sum(nc[v][iT,:,:][ind]*nc['area'][ind])/1e12

#e='Pacific Maritime'
#e='Montane Cordillera'
#e='Boreal Shield'
e='Boreal Plains'
#e='Mixedwood Plains'
#e='Boreal Cordillera'
#e='Taiga Shield'
#e='Taiga Cordillera'
#e='Southern Arctic'
#e='Hudson Plains'
#e='Atlantic Maritime'
#e='Prairies'
a=np.array(list(d[e].values()))
plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(30,20))
plt.bar(np.arange(1,7),a)
ax.set(xticks=np.arange(1,7),xticklabels=vL)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')

#%% Analysis by PB
vL=['LNLGIS_NCE_meanYear','LNLGIS_dC_loss_meanYear','LNLGIS_NBE_meanYear','FF_meanYear','Crop_meanYear','Wood_meanYear','River_meanYear']
#vL2=['LNLGIS_dC_loss','LNLGIS_NBE','FF','Crop','Wood','River']
d={}
dts={}
for u in u_pb:
	ind=np.where(grd['PB2']==lut_pb[u])
	#v='FF_meanYear'
	f=grd['FF_meanYear'][ind]
	A=grd['area'][ind]
	#print(np.sum(A)/1e6/10000)
	#print(np.sum(f*A)/1e12) # MtCO2e/yr

	d[u]={}
	for v in vL:
		#d[u][v]=np.nanmean(grd[v][ind])/1e6*1e4 # tCO2e/ha/yr
		d[u][v]=np.sum(np.sum(grd[v][ind]*grd['area'][ind]))/1e12 # MtCO2e/yr
# 	dts[u]={}
# 	for v in vL2:
# 		dts[u][v]=np.zeros(6)
# 		for iT in range(6):
# 			#dts[u][v][iT]=np.mean(nc[v][iT,:,:][ind])/1e6*1e4
# 			dts[u][v][iT]=np.sum(nc[v][iT,:,:][ind]*nc['area'][ind])/1e12

e='India'
#e='US-OR'
a=np.array(list(d[e].values()))
plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(30,20))
plt.bar(np.arange(1,7),a)
ax.set(xticks=np.arange(1,7),xticklabels=vL)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')

#%%

plt.close('all'); plt.plot(nc['lat'],np.mean(nc['LNLGIS_NBE_meanYear'],axis=1)/1e6*1e4,'b-')





#%%

flg=0
if flg==1:
	# Confirm that it is working
	plt.close('all'); fig,ax=plt.subplots(1)
	#ax.matshow(grd['FF_meanYear'],extent=zRef['Extent'],clim=[0,1200])
	ax.matshow(grd['LNLGIS_NBE_meanYear'],extent=zRef['Extent'],clim=[-200,300])
	#ax.matshow(np.flip(nc['FF_meanYear'],axis=0),extent=zRef['Extent'],clim=[0,400])
	#ax.matshow(np.flip(nc['River_meanYear'],axis=0),extent=zRef['Extent'],clim=[-20,20])
	#ax.matshow(np.flip(nc['Wood_meanYear'],axis=0),extent=zRef['Extent'],clim=[-60,60])
	#ax.matshow(np.flip(nc['LNLG_NBE_meanYear'],axis=0),extent=zRef['Extent'],clim=[-200,300])
	gdf_pb.plot(ax=ax,facecolor='None',linewidth=1,edgecolor='w')
	gdf_ez.plot(ax=ax,facecolor='None',linewidth=0.25,edgecolor='c')

flg=0
if flg==1:
	# Confirm that it is working
	plt.close('all'); fig,ax=plt.subplots(1)
	#ax.matshow(zRef['Data'],extent=zRef['Extent'])
	ax.matshow(np.flip(nc['FF_meanYear'],axis=0),extent=zRef['Extent'],clim=[0,1200])
	#ax.matshow(np.flip(nc['LNLG_NBE_meanYear'],axis=0),extent=zRef['Extent'],clim=[-200,200])
	gdf_ez.plot(ax=ax,facecolor='None',linewidth=1,edgecolor='w')

#%% Analysis by PT
vL=['LNLGIS_dC_loss_meanYear','LNLGIS_NBE_meanYear','FF_meanYear','Crop_meanYear','Wood_meanYear','River_meanYear']
vL2=['LNLGIS_dC_loss','LNLGIS_NBE','FF','Crop','Wood','River']
d={}
dts={}
for u in u_pb:
	if u=='CA-BC':
		break
	ind=np.where(np.flip(np.flip(grd['PB2'],axis=0),axis=1)==lut_pb[u])
	v='FF_meanYear'
	f=nc[v][ind]
	A=nc['area'][ind]
	print(np.sum(A)/1e6/10000)
	print(np.sum(f*A)/1e12) # MtCO2e/yr

	d[u]={}
	for v in vL:
		#d[u][v]=np.nanmean(nc[v][ind])/1e6*1e4 # tCO2e/ha/yr
		d[u][v]=np.sum(nc[v][ind]*nc['area'][ind])/1e12 # MtCO2e/yr
	dts[u]={}
	for v in vL2:
		dts[u][v]=np.zeros(6)
		for iT in range(6):
			#dts[u][v][iT]=np.mean(nc[v][iT,:,:][ind])/1e6*1e4
			dts[u][v][iT]=np.sum(nc[v][iT,:,:][ind]*nc['area'][ind])/1e12

e='CA-BC'
y=np.array(list(d[e].values()))
plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(30,20))
plt.bar(np.arange(1,7),y)
ax.set(xticks=np.arange(1,7),xticklabels=vL)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')

fig,ax=plt.subplots(1,figsize=gu.cm2inch(30,20))
plt.plot(np.arange(2015,2021,1),dts[e]['LNLGIS_NBE'],'-bo')

fig,ax=plt.subplots(1,figsize=gu.cm2inch(30,20))
plt.plot(np.arange(2015,2021,1),dts[e]['FF'],'-bo')

#%%
ptL=['CA-NU','CA-NT','CA-YT','CA-QC','CA-NL','CA-BC','CA-AB','CA-SK','CA-MB','CA-ON','CA-NB','CA-PE','CA-NS']
a=0
for pt in ptL:
	#a=a+d[pt]['FF_meanYear']
	a=a+d[pt]['LNLGIS_NBE_meanYear']
print(a)

#%%
ptL=['US-AK','US-MN', 'US-WA','US-MT', 'US-ND', 'US-ID', 'US-WI', 'US-MI', 'US-ME','US-OR', 'US-SD', 'US-NH', 'US-VT', 'US-NY',
'US-WY', 'US-IA', 'US-NE', 'US-MA', 'US-IL', 'US-PA', 'US-CT','US-RI', 'US-CA', 'US-UT', 'US-NV', 'US-OH', 'US-IN', 'US-NJ',
'US-CO', 'US-WV', 'US-MO', 'US-KS', 'US-DE', 'US-MD', 'US-VA','US-KY', 'US-AZ', 'US-OK', 'US-NM', 'US-TN', 'US-NC', 'US-TX','US-AR', 'US-SC', 'US-AL', 'US-GA', 'US-MS', 'US-LA','US-FL']
a=0
for pt in ptL:
	a=a+d[pt]['FF_meanYear']
print(a)

#%% Analysis by ecozone
vL=['LNLGIS_dC_loss_meanYear','LNLGIS_NBE_meanYear','FF_meanYear','Crop_meanYear','Wood_meanYear','River_meanYear']
vL2=['LNLGIS_dC_loss','LNLGIS_NBE','FF','Crop','Wood','River']
d={}
dts={}
for u in uez:
	ind=np.where(np.flip(np.flip(grd['EZ'],axis=0),axis=1)==lut_ez[u])
	#v='FF_meanYear'
	#f=nc[v][ind]
	#A=nc['area'][ind]
	#print(np.sum(A)/1e6/10000)
	#print(np.sum(f*A)/1e12) # MtCO2e/yr

	d[u]={}
	for v in vL:
		#d[u][v]=np.nanmean(f)/1e6*1e4 # tCO2e/ha/yr
		d[u][v]=np.sum(np.sum(nc[v][ind]*nc['area'][ind]))/1e12 # MtCO2e/yr
	dts[u]={}
	for v in vL2:
		dts[u][v]=np.zeros(6)
		for iT in range(6):
			#dts[u][v][iT]=np.mean(nc[v][iT,:,:][ind])/1e6*1e4
			dts[u][v][iT]=np.sum(nc[v][iT,:,:][ind]*nc['area'][ind])/1e12

e='Pacific Maritime'
#e='Montane Cordillera'
#e='Boreal Shield'
#e='Boreal Plains'
#e='Mixedwood Plains'
#e='Boreal Cordillera'
#e='Taiga Shield'
#e='Taiga Cordillera'
#e='Southern Arctic'
#e='Hudson Plains'
#e='Atlantic Maritime'
#e='Prairies'
a=np.array(list(d[e].values()))
plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(30,20))
plt.bar(np.arange(1,7),a)
ax.set(xticks=np.arange(1,7),xticklabels=vL)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')

fig,ax=plt.subplots(1,figsize=gu.cm2inch(30,20))
plt.plot(np.arange(2015,2021,1),dts[e]['LNLGIS_NBE'],'-bo')


#%%

plt.close('all'); plt.matshow(grd['FF_meanYear'],clim=[0,1200])
plt.close('all'); plt.matshow(grd['Wood_meanYear'],clim=[-60,60])
#plt.close('all'); plt.matshow(np.flip(nc['LNLG_NBE_meanYear'],axis=0),clim=[-800,0])
plt.close('all'); plt.matshow(np.flip(nc['LNLG_dC_loss_meanYear'],axis=0),clim=[-800,800])
plt.close('all'); plt.matshow(np.flip(nc['Wood_meanYear'],axis=0),clim=[-100,100])
plt.close('all'); plt.matshow(np.flip(nc['River_meanYear'],axis=0),clim=[-100,100])

np.mean(nc['FF_meanYear'])*5.1e14/1e15
np.mean(nc['River_meanYear'])*5.1e14/1e15


# #%% Resample to higher resolution
# sf=5
# gis.ResampleRaster('C:\Data\TopDown\PB2.tif',sf)
# gis.ResampleRaster('C:\Data\TopDown\EZ.tif',sf)

# vL=['LNLGIS_dC_loss_meanYear','LNLGIS_NBE_meanYear','FF_meanYear','Crop_meanYear','Wood_meanYear','River_meanYear']
# for v in vL:
#  	z1=copy.deepcopy(zRef)
#  	z1['Data']=np.flip(nc[v],axis=0).astype('float32')
#  	pth=r'C:\Data\TopDown' + '\\' + v + '.tif'
#  	gis.SaveGeoTiff(z1,pth)
#  	gis.ResampleRaster(pth,sf)