'''
IMPACTS OF HARVEST RETENTION ON FOREST CARBON STOCKS IN THE BOUNDARY TIMBER SUPPLY AREA
Datasets:
    - whole stem volume from LiDAR-driven Predictive Forest Invnetory (from Geoff Quinn and Chris Butson)
    - Consolidated cutblocks database (year of harvest)
    - Harvest retention compilation 1 (see look-up table for definitions)
'''

#%% Import modules
import numpy as np
import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
from shapely.geometry import Polygon,Point
import copy
from rasterio import features
import statsmodels.formula.api as smf
import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.bc1ha.bc1ha_util as u1ha
import fcgadgets.bc1ha.bc1ha_plot as p1ha
import fcgadgets.macgyver.util_query_gdb as qgdb

#%% Import parameters
meta=u1ha.Init()
meta['Graphics']['Map']['RGSF']=1
meta['Graphics']['Map']['Fig Width']=15.5
meta['Graphics']['Map']['Side Space']=0.25
meta['Graphics']['Map']['Map Position']=[0,0,1-meta['Graphics']['Map']['Side Space']-0.01,1]
meta['Graphics']['Map']['Map Axis Vis']='off'
meta['Graphics']['Map']['Map Grid Vis']=False
meta['Graphics']['Map']['Legend X']=1-meta['Graphics']['Map']['Side Space']+meta['Graphics']['Map']['Map Position'][0]#+0.01,0.6,0.03,0.35]
meta['Graphics']['Map']['Legend Width']=0.0275
meta['Graphics']['Map']['Legend Font Size']=7
meta['Graphics']['Map']['Legend Text Space']=0.035
meta['Graphics']['Map']['Show Bound Land Mask']='On'
meta['Graphics']['Map']['Show Bound Within']='On'
meta['Graphics']['Map']['Show Lakes']='Off'
meta['Graphics']['Map']['Show Rivers']='Off'
meta['Graphics']['Map']['Show Roads']='On'
meta['Graphics']['Plot Style']='Manuscript'
meta['Graphics']['gp']=gu.SetGraphics(meta['Graphics']['Plot Style'])
meta['Graphics']['Print Figures']='On'
meta['Graphics']['Print Figure Path']=r'C:\Users\rhember\OneDrive - Government of BC\Figures\Harvest\PFI Boundary TSA'

meta['Paths']['Working']=r'C:\Data\BC20m\Boundary TSA'

# Import Growth and Yield model
dGY=gu.ReadExcel(r'C:\Users\rhember\OneDrive - Government of BC\Projects\MGEM 2023 Harvesting and Retention\TIPSY Age Response.xlsx')

#%% Define region of interest
roi={}
roi['Type']='ByTSA'
roi['Name']='Boundary TSA'
roi['List']=['Boundary TSA']
gdf=u1ha.Import_GDBs_ProvinceWide(meta)
roi=u1ha.DefineROI(meta,roi,gdf)

#%% Override the 1ha defaults with the Predictive Forest Inentory dataset from Geoff Quinn
zV=gis.OpenGeoTiff(r'C:\Data\LiDAR\Boundary TSA\From GQ 20230302\BoundaryTSA_PFI_LnTotalWSV.tif')
zV=gis.ClipRasterByXYLimits(zV,roi['grd']['xlim'],roi['grd']['ylim'])
roi['grd']=copy.deepcopy(zV)

# Convert to stemwood biomass (tC/ha)
zB=copy.deepcopy(zV)
zB['Data']=0.5*0.45*zB['Data']

#%% Rasterize TSA boundary at 20 m
df0=roi['gdf']['bound within']
df0=df0[df0.geometry!=None]
df0=df0.reset_index()
df0['ID']=1
shapes=((geom,value) for geom, value in zip(df0.geometry,df0['ID']))
z0=np.zeros(zV['Data'].shape,dtype=float)
burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=roi['grd']['Transform'])
roi['grd']['Data']=burned.astype('int8')

# Remove exterior region from volume
ind=np.where(roi['grd']['Data']==0); zV['Data'][ind]=0
ind=np.where(zV['Data']<0); zV['Data'][ind]=-99
zV['Data']=zV['Data'].astype('int16')
gis.SaveGeoTiff(zV,meta['Paths']['Working'] + '\\VolumeWholeStem.tif')

#%% Import vector data
#roi=u1ha.Import_GDB_Over_ROI(meta,roi,['fcres','cc'])

#%% Plot ROI
cm=plt.cm.get_cmap('viridis',30)
plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
im=ax[0].matshow(zV['Data'],extent=roi['grd']['Extent'],clim=[0,700],cmap=cm)
roi=p1ha.PlotVectorBaseMaps(meta,roi,ax[0])
ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])
#cc0=roi['gdf']['cc'][(roi['gdf']['cc']['HARVEST_YEAR']==2010)]
#cc0.plot(ax=ax[0],edgecolor=[0.85,0.35,1],facecolor='None',linewidth=1)
roi['gdf']['fcres']['gdf'].plot(ax=ax[0],linestyle='--',edgecolor=[0.6,1,0],facecolor='None',linewidth=1)
cb=plt.colorbar(im,cax=ax[1],cmap=cm)
ax[1].set(position=[0.76,0.5,0.03,0.35])

#%% Rasterize acquisition date
def RasterizeAcquisitionDate(meta):
	df=gpd.read_file(r'C:\Data\LiDAR\Boundary TSA\RDKB\Boundary_collection_metadata.shp')
	df['Year']=np.zeros(df['year'].size)
	for y in df['year'].unique():
		ind=np.where(df['year']==y)[0]
		df['Year'][ind]=np.array(y,dtype='float')
	shapes=((geom,value) for geom, value in zip(df.geometry,df['Year']))
	z0=np.zeros(zV['Data'].shape,dtype=float)
	burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=roi['grd']['Transform'])
	z1=copy.deepcopy(zV)
	z1['Data']=burned.astype('int16')
	gis.SaveGeoTiff(z1,meta['Paths']['Working'] + '\\Year_Acquisition.tif')
	return

#%% Rasterize harvest year from consolidated cutblocks DB
def RasterizeHarvestYear(meta):
	flg=0
	if flg==1:
		lNam='VEG_CONSOLIDATED_CUT_BLOCKS_SP'
		vNam='HARVEST_YEAR'
		ind=np.where( (meta['Geos']['Variable Info']['Layer Name']==lNam) & (meta['Geos']['Variable Info']['Variable Name']==vNam) )[0]
		pthin=meta['Paths']['GDB'][ meta['Geos']['Variable Info']['Geodatabase Name'][ind[0]] ]
		df=gpd.read_file(pthin,layer=lNam)
		df=df[df.geometry!=None]
		df=df.reset_index()
		
		# Clip to ROI
		roi['gdf']['H_CC']=df.cx[roi['grd']['xmin']:roi['grd']['xmax'],roi['grd']['ymin']:roi['grd']['ymax']]
		roi['gdf']['H_CC']=roi['gdf']['H_CC'].reset_index(drop=True)
		roi['gdf']['H_CC']=gpd.overlay(roi['gdf']['H_CC'],roi['gdf']['bound within'],how='intersection')
		roi['gdf']['H_CC'].to_file(meta['Paths']['Working'] + '\\ConsolidatedCutblocks.geojson',driver='GeoJSON')
	else:
		roi['gdf']['H_CC']=gpd.read_file(meta['Paths']['Working'] + '\\ConsolidatedCutblocks.geojson')
	
	# Rasterize
	zYearLast=copy.deepcopy(roi['grd'])
	zYearLast['Data']=np.zeros(roi['grd']['Data'].shape,dtype='int16')
	uYear=roi['gdf']['H_CC'][vNam].unique()
	tv=np.arange(np.min(uYear),np.max(uYear),1)
	for iT in range(tv.size):
	    print(tv[iT])
	    df0=roi['gdf']['H_CC'][roi['gdf']['H_CC'][vNam]==tv[iT]].copy()
	    shapes=((geom,value) for geom, value in zip(df0.geometry,df0[vNam]))
	    z0=np.zeros(roi['grd']['Data'].shape,dtype=float)
	    if len(df0)>0:
	        burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=roi['grd']['Transform'])
	    z1=copy.deepcopy(roi['grd'])
	    z1['Data']=z0.astype(meta['Geos']['Variable Info']['Precision'][ind[0]])
	    #gis.SaveGeoTiff(z1,meta['Paths']['bc1ha'] + '\\' + lNam + '\\' + vNam + '_' + str(tv[iT]) + '.tif')
	    zYearLast['Data'][burned>0]=tv[iT]
	gis.SaveGeoTiff(zYearLast,meta['Paths']['Working'] + '\\Harvest_CC_YearLast.tif')
	#plt.close('all'); plt.matshow(zYearLast['Data'],clim=[1955,2022])
	
	# Harvest year with inner buffer (to avoid edge effects)
	buf=roi['gdf']['H_CC'].geometry.buffer(-20)
	buf['ID']=np.ones(len(buf))
	shapes=((geom,value) for geom, value in zip(buf.geometry,buf['ID']))
	z0=np.zeros(roi['grd']['Data'].shape,dtype=float)
	burned=features.rasterize(shapes=shapes,fill=0,out=z0,transform=roi['grd']['Transform'])
	zH_Buf=burned.astype('int8')
	gis.SaveGeoTiff(zYearLast,meta['Paths']['Working'] + '\\Harvest_CC_YearLast_Buffer.tif')
	
	# Check that it worked
	flg=0
	if flg==1:
		Mask=np.zeros(zH['Data'].shape,dtype='int8')
		ind=np.where(zH['Data']>0); Mask[ind]=1
		ind=np.where(zH_Buf>0); Mask[ind]=2
		plt.matshow(Mask)
	return

#%% Clip and resample variables from BC1ha database
def ClipAndResampleFromBC1ha(meta):
	# Clip and resample land cover compilation 1 (2019)
	z1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\LandCoverUse\\LandCover_Comp1_2019.tif') 
	z2=gis.ClipRasterByXYLimits(z1,roi['grd']['xlim'],roi['grd']['ylim'])
	pth=meta['Paths']['Working'] + '\\LandCover_Comp1_2019.tif'
	gis.SaveGeoTiff(z2,pth)
	gis.ResampleRaster(pth,5)
	z2=gis.OpenGeoTiff(pth) 
	z2=gis.ClipToRaster(z2,roi['grd'])
	gis.SaveGeoTiff(z2,pth)
	#plt.matshow(z2['Data'])
	
	# Clip and resample harvest retention compilation 1
	z1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\Management\\HarvestRetentionComp1.tif') 
	z2=gis.ClipRasterByXYLimits(z1,roi['grd']['xlim'],roi['grd']['ylim'])
	pth=meta['Paths']['Working'] + '\\HarvestRetentionComp1.tif'
	gis.SaveGeoTiff(z2,pth)
	gis.ResampleRaster(pth,5)
	z2=gis.OpenGeoTiff(pth) 
	z2=gis.ClipToRaster(z2,roi['grd'])
	gis.SaveGeoTiff(z2,pth)
	plt.matshow(z2['Data'])
	
	# Clip and resample silviculture system
	z1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\RSLT_OPENING_SVW\\SILV_SYSTEM_CODE.tif')
	z2=gis.ClipRasterByXYLimits(z1,roi['grd']['xlim'],roi['grd']['ylim'])
	pth=meta['Paths']['Working'] + '\\SILV_SYSTEM_CODE.tif'
	gis.SaveGeoTiff(z2,pth)
	gis.ResampleRaster(pth,5)
	z2=gis.OpenGeoTiff(pth) 
	z2=gis.ClipToRaster(z2,roi['grd'])
	gis.SaveGeoTiff(z2,pth)
	plt.matshow(z2['Data'])
	
	# Clip and resample VRI volume
	z1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\LIVE_STAND_VOLUME_125.tif') 
	z2=gis.ClipRasterByXYLimits(z1,roi['grd']['xlim'],roi['grd']['ylim'])
	pth=meta['Paths']['Working'] + '\\LIVE_STAND_VOLUME_125.tif'
	gis.SaveGeoTiff(z2,pth)
	gis.ResampleRaster(pth,5)
	z2=gis.OpenGeoTiff(pth) 
	z2=gis.ClipToRaster(z2,roi['grd'])
	gis.SaveGeoTiff(z2,pth)
	#plt.matshow(z2['Data'])
	
	z1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\DEAD_STAND_VOLUME_125.tif') 
	z2=gis.ClipRasterByXYLimits(z1,roi['grd']['xlim'],roi['grd']['ylim'])
	pth=meta['Paths']['Working'] + '\\DEAD_STAND_VOLUME_125.tif'
	gis.SaveGeoTiff(z2,pth)
	gis.ResampleRaster(pth,5)
	z2=gis.OpenGeoTiff(pth) 
	z2=gis.ClipToRaster(z2,roi['grd'])
	gis.SaveGeoTiff(z2,pth)
	#plt.matshow(z2['Data'])
	
	z1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\WHOLE_STEM_BIOMASS_PER_HA.tif') 
	z2=gis.ClipRasterByXYLimits(z1,roi['grd']['xlim'],roi['grd']['ylim'])
	pth=meta['Paths']['Working'] + '\\WHOLE_STEM_BIOMASS_PER_HA.tif'
	gis.SaveGeoTiff(z2,pth)
	gis.ResampleRaster(pth,5)
	z2=gis.OpenGeoTiff(pth) 
	z2=gis.ClipToRaster(z2,roi['grd'])
	gis.SaveGeoTiff(z2,pth)
	#plt.matshow(z2['Data'])
	
	z1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\VRI 2023\\PROJ_AGE_1.tif') 
	z2=gis.ClipRasterByXYLimits(z1,roi['grd']['xlim'],roi['grd']['ylim'])
	pth=meta['Paths']['Working'] + '\\PROJ_AGE_1.tif'
	gis.SaveGeoTiff(z2,pth)
	gis.ResampleRaster(pth,5)
	z2=gis.OpenGeoTiff(pth) 
	z2=gis.ClipToRaster(z2,roi['grd'])
	gis.SaveGeoTiff(z2,pth)
	#plt.matshow(z2['Data'])
	
	z1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\OGSR_TAP_PRIORITY_DEF_AREA_SP\\PRIORITY_DEFERRAL_ID.tif') 
	z2=gis.ClipRasterByXYLimits(z1,roi['grd']['xlim'],roi['grd']['ylim'])
	pth=meta['Paths']['Working'] + '\\PRIORITY_DEFERRAL_ID.tif'
	gis.SaveGeoTiff(z2,pth)
	gis.ResampleRaster(pth,5)
	z2=gis.OpenGeoTiff(pth)
	z2=gis.ClipToRaster(z2,roi['grd'])
	gis.SaveGeoTiff(z2,pth)
	#plt.matshow(z2['Data'])
	
	# Clip and resample wildfire year
	z1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PROT_HISTORICAL_FIRE_POLYS_SP\\FIRE_YEAR_YearLast.tif')
	z2=gis.ClipRasterByXYLimits(z1,roi['grd']['xlim'],roi['grd']['ylim'])
	pth=meta['Paths']['Working'] + '\\FIRE_YEAR_YearLast.tif'
	gis.SaveGeoTiff(z2,pth)
	gis.ResampleRaster(pth,5)
	z2=gis.OpenGeoTiff(pth) 
	z2=gis.ClipToRaster(z2,roi['grd'])
	gis.SaveGeoTiff(z2,pth)
	plt.matshow(z2['Data'])
	
	# Clip and resample insect comp year
	z1=gis.OpenGeoTiff(meta['Paths']['bc1ha'] + '\\PEST_INFESTATION_POLY\\PEST_SEVERITY_CODE_IBM_YearLast.tif')
	z2=gis.ClipRasterByXYLimits(z1,roi['grd']['xlim'],roi['grd']['ylim'])
	pth=meta['Paths']['Working'] + '\\PEST_SEVERITY_CODE_IBM_YearLast.tif'
	gis.SaveGeoTiff(z2,pth)
	gis.ResampleRaster(pth,5)
	z2=gis.OpenGeoTiff(pth) 
	z2=gis.ClipToRaster(z2,roi['grd'])
	gis.SaveGeoTiff(z2,pth)
	plt.matshow(z2['Data'])
	return

#%% Import variables
zAc=gis.OpenGeoTiff(meta['Paths']['Working'] + '\\Year_Acquisition.tif')
zVL=gis.OpenGeoTiff(meta['Paths']['Working'] + '\\LIVE_STAND_VOLUME_125.tif')
zVD=gis.OpenGeoTiff(meta['Paths']['Working'] + '\\DEAD_STAND_VOLUME_125.tif')
zLCC1=gis.OpenGeoTiff(meta['Paths']['Working'] + '\\LandCover_Comp1_2019.tif')
zA_VRI=gis.OpenGeoTiff(meta['Paths']['Working'] + '\\PROJ_AGE_1.tif')
zB_VRI=gis.OpenGeoTiff(meta['Paths']['Working'] + '\\WHOLE_STEM_BIOMASS_PER_HA.tif')
zB_VRI['Data']=0.5*zB_VRI['Data']
zD=gis.OpenGeoTiff(meta['Paths']['Working'] + '\\PRIORITY_DEFERRAL_ID.tif')
zH=gis.OpenGeoTiff(meta['Paths']['Working'] + '\\Harvest_CC_YearLast.tif')
zH_Buf=gis.OpenGeoTiff(meta['Paths']['Working'] + '\\Harvest_CC_YearLast_Buffer.tif')
zR=gis.OpenGeoTiff(meta['Paths']['Working'] + '\\HarvestRetentionComp1.tif')
zSS=gis.OpenGeoTiff(meta['Paths']['Working'] + '\\SILV_SYSTEM_CODE.tif')

#%% Plot harvest
cm=plt.cm.get_cmap('viridis',30)
plt.close('all'); fig,ax=plt.subplots(1,2,figsize=gu.cm2inch(meta['Graphics']['Map']['Fig Width'],(1-meta['Graphics']['Map']['Side Space'])*meta['Graphics']['Map']['Fig Width']*roi['grd']['yxrat']))
im=ax[0].matshow(zH['Data'],extent=roi['grd']['Extent'],clim=[1900,2020],cmap=cm)
roi=p1ha.PlotVectorBaseMaps(meta,roi,ax[0])
ax[0].set(position=meta['Graphics']['Map']['Map Position'],xlim=roi['grd']['xlim'],ylim=roi['grd']['ylim'])
ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both'); ax[0].grid(meta['Graphics']['Map']['Map Grid Vis']); ax[0].axis(meta['Graphics']['Map']['Map Axis Vis'])
cc0=roi['gdf']['H_CC'][(roi['gdf']['H_CC']['HARVEST_YEAR']==2010)]
cc0.plot(ax=ax[0],edgecolor=[0.85,0.35,1],facecolor='None',linewidth=1)
#cd='G'
#cd='D'
cd='W'
#cd='R'
#cd='O'
roi['gdf']['fcres']['gdf'][ (roi['gdf']['fcres']['gdf']['SILV_RESERVE_CODE']==cd) ].plot(ax=ax[0],linestyle='--',edgecolor=[1,1,1],facecolor='None',linewidth=1)
cb=plt.colorbar(im,cax=ax[1],cmap=cm)
ax[1].set(position=[0.76,0.5,0.03,0.35])

#%%
fig,ax=plt.subplots(1)
ind=np.where( (zH['Data']==2010) & (zV['Data']>=0) )
plt.hist(zV['Data'][ind],np.arange(0,600,5))

#%%

ikp=np.where( (zV['Data']>0) & (zLCC1['Data']==1) & (zH['Data']==0) )
print(np.mean(zV['Data'][ikp]))
ikp=np.where( (zV['Data']>0) & (zLCC1['Data']==1) & (zH['Data']==0) )
print(np.mean(zB['Data'][ikp]))
ikp=np.where( (zV['Data']>0) & (zLCC1['Data']==1) & (zH['Data']>0) )
print(np.mean(zB['Data'][ikp]))

ikp=np.where( (zV['Data']>0) & (zLCC1['Data']==1) & (zH['Data']>0) )
print(np.mean(zV['Data'][ikp]))
print(np.mean(zB['Data'][ikp]))

TSH=zAc['Data']-zH['Data']
print(np.mean(TSH[ikp]))


ikp=np.where( (zLCC1['Data']==1) )
print(ikp[0].size/1e6)
ikp=np.where( (zLCC1['Data']==1) & (zH['Data']>0) )
print(ikp[0].size/1e6)
ikp=np.where( (zLCC1['Data']==1) & (zH['Data']>0) & (zSS['Data']>0) )
print(ikp[0].size/1e6)

#%% Modelling relationship between biomass and time since harvest, stratified by SSC

ikp=np.where( (zB['Data']>0) & (zLCC1['Data']==1) & (zH['Data']==0) )
mu_ref=np.mean(zB['Data'][ikp])

ikp=np.where( (zB['Data']>0) & (zH_Buf['Data']>0) & (zAc['Data']-zH['Data']>=10) & (zAc['Data']-zH['Data']<50) )
# z=np.zeros(zB['Data'].shape); z[ikp]=1; plt.matshow(z)

ivl=5
ikp[0].size
B=zB['Data'][ikp][0::ivl]
tsh=zAc['Data'][ikp][0::ivl]-zH['Data'][ikp][0::ivl]
plt.hist(tsh)
ss_id=zSS['Data'][ikp][0::ivl]
ss=np.array(['No Class' for _ in range(ss_id.size)],dtype=object)
idx=gu.IndicesFromUniqueArrayValues(ss_id)
for k in idx.keys():
	if k==0:
		continue
	ss[idx[k]]=u1ha.lut_n2s(meta['LUT']['RSLT_OPENING_SVW']['DENUDATION_1_SILV_SYSTEM_CODE'],k)[0]
u,N=gu.CountByCategories(ss,'Percent')
u,N=gu.CountByCategories(ss)

#x=np.flip(np.array([B,tsh]).T,axis=0)
df=pd.DataFrame({'B':B,'TSH':tsh,'SS':ss})
print(df.describe())
507368
# Use the formula method to perform a regression
md=smf.ols("B~TSH+C(SS)+TSH:C(SS)",data=df)
rs=md.fit(maxiter=100)
print(rs.summary())

tsh=np.arange(1,61)
lw=1.5
plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(12,10));
plt.plot([0,100],[mu_ref,mu_ref],'k-',color=[0.75,0.75,0.75],lw=2.5,label='Reference level\n(landscape mean biomass)')
y=rs.params['Intercept']+tsh*rs.params['TSH']
plt.plot(tsh,y,'k-',lw=lw,label='Clearcut with reserves')
y=rs.params['Intercept']+tsh*rs.params['TSH']+rs.params['C(SS)[T.CLEAR]']+rs.params['TSH:C(SS)[T.CLEAR]']*tsh
plt.plot(tsh,y,'b-',lw=lw,label='Clearcut')
#y=rs.params['Intercept']+tsh*rs.params['TSH']+rs.params['C(SS)[T.RETEN]']+rs.params['TSH:C(SS)[T.RETEN]']*tsh
#plt.plot(tsh,y,'g--',lw=lw,label='Retention')
#y=rs.params['Intercept']+tsh*rs.params['TSH']+rs.params['C(SS)[T.No Class]']+rs.params['TSH:C(SS)[T.No Class]']*tsh
#plt.plot(tsh,y,'r:',label='No Class')
y=rs.params['Intercept']+tsh*rs.params['TSH']+rs.params['C(SS)[T.SELEC]']+rs.params['TSH:C(SS)[T.SELEC]']*tsh
plt.plot(tsh,y,'m-.',lw=lw,label='Selection')
y=rs.params['Intercept']+tsh*rs.params['TSH']+rs.params['C(SS)[T.PATCT]']+rs.params['TSH:C(SS)[T.PATCT]']*tsh
plt.plot(tsh,y,'r-',lw=lw,label='Patch')
y=rs.params['Intercept']+tsh*rs.params['TSH']+rs.params['C(SS)[T.SEEDT]']+rs.params['TSH:C(SS)[T.SEEDT]']*tsh
plt.plot(tsh,y,'y--',lw=lw,label='Seed tree retention')
y=rs.params['Intercept']+tsh*rs.params['TSH']+rs.params['C(SS)[T.SHELT]']+rs.params['TSH:C(SS)[T.SHELT]']*tsh
plt.plot(tsh,y,'c-.',lw=lw,label='Shelterwood')
ax.set(xticks=np.arange(0,100,5),yticks=np.arange(0,400,5),ylabel='Stemwood biomass (tC/ha)',xlabel='Time since harvest, years',xlim=[0,60],ylim=[0,62])
ax.legend(loc='lower right',facecolor=[1,1,1],frameon=False);
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
plt.tight_layout()
#y=rs.params['Intercept']+tsh*rs.params['TSH']+rs.params['C(SS)[T.No Class]']+rs.params['TSH:C(SS)[T.No Class]']*tsh
#plt.plot(tsh,y,'k-',color=[0.75,0.75,0.75],lw=2)
#plt.plot(dGY['Age'],dGY['Stemwood biomass (tC/ha)'],'k-',color=[0.5,0.5,0.5],lw=2)
gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Regression_BiomassVsAge_BySSC','png',900)

#%% Modelling by HRComp1

ikp=np.where( (zB['Data']>0) & (zLCC1['Data']==1) & (zH['Data']==0) )
mu_ref=np.mean(zB['Data'][ikp])

ikp=np.where( (zB['Data']>0) & (zH_Buf['Data']>0) & (zAc['Data']-zH['Data']>=10) & (zAc['Data']-zH['Data']<50) )

ivl=5
B=zB['Data'][ikp][0::ivl]
tsh=zAc['Data'][ikp][0::ivl]-zH['Data'][ikp][0::ivl]
ss_id=zR['Data'][ikp][0::ivl]
ss=np.array(['No Class' for _ in range(ss_id.size)],dtype=object)
idx=gu.IndicesFromUniqueArrayValues(ss_id)
for k in idx.keys():
	if k==0:
		continue
	ss[idx[k]]=u1ha.lut_n2s(meta['LUT']['Derived']['harvret1'],k)[0]
u,N=gu.CountByCategories(ss,'Percent')

#x=np.flip(np.array([B,tsh]).T,axis=0)
df=pd.DataFrame({'B':B,'TSH':tsh,'SS':ss})
print(df.describe())

# Use the formula method to perform a regression
md=smf.ols("B~TSH+C(SS)+TSH:C(SS)",data=df)
rs=md.fit(maxiter=100)
print(rs.summary())

tsh=np.arange(1,61)
lw=1.25
plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(14,10));
plt.plot([0,100],[mu_ref,mu_ref],'k-',color=[0.75,0.75,0.75],lw=2.5,label='Reference level\n(landscape mean biomass)')
y=rs.params['Intercept']+tsh*rs.params['TSH']
plt.plot(tsh,y,'g-',lw=lw,label='Dispersed')
y=rs.params['Intercept']+tsh*rs.params['TSH']+rs.params['C(SS)[T.Harvested with no reserves]']+rs.params['TSH:C(SS)[T.Harvested with no reserves]']*tsh
plt.plot(tsh,y,'k-',lw=lw,label='Harvest with no reserves')
y=rs.params['Intercept']+tsh*rs.params['TSH']+rs.params['C(SS)[T.Group]']+rs.params['TSH:C(SS)[T.Group]']*tsh
plt.plot(tsh,y,'c-',lw=lw,label='Group retention')
y=rs.params['Intercept']+tsh*rs.params['TSH']+rs.params['C(SS)[T.Wildlife trees]']+rs.params['TSH:C(SS)[T.Wildlife trees]']*tsh
plt.plot(tsh,y,'y-',lw=lw,color=[0.5,0,1],label='Wildlife trees')
#y=rs.params['Intercept']+tsh*rs.params['TSH']+rs.params['C(SS)[T.Riparian]']+rs.params['TSH:C(SS)[T.Riparian]']*tsh
#plt.plot(tsh,y,'r--',lw=lw,label='Riparian reserves')
#y=rs.params['Intercept']+tsh*rs.params['TSH']+rs.params['C(SS)[T.Un]']+rs.params['TSH:C(SS)[T.Riparian]']*tsh
#plt.plot(tsh,y,'r--',lw=lw,label='Riparian reserves')
y=rs.params['Intercept']+tsh*rs.params['TSH']+rs.params['C(SS)[T.Other]']+rs.params['TSH:C(SS)[T.Other]']*tsh
plt.plot(tsh,y,'y--',lw=lw,label='Other reserves')
#plt.plot(dGY['Age'],dGY['Stemwood biomass (tC/ha)'],'k-',color=[0.75,0.75,0.75],lw=2)
ax.set(xticks=np.arange(0,100,5),yticks=np.arange(0,400,5),ylabel='Stemwood biomass (tC/ha)',xlabel='Time since harvest, years',xlim=[0,60],ylim=[0,62])
ax.legend(loc='lower right',facecolor=[1,1,1],frameon=False);
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
plt.tight_layout()
gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Regression_BiomassVsAge_ByHRC','png',900)

#%% Relationship between VRI volume and PFI volume

# ikp=np.where( (zV['Data']>1) & (zVL['Data']>1) )
# x=zVL['Data'][ikp][0::30]
# x2=zVL['Data'][ikp][0::30]+zVD['Data'][ikp][0::30]
# y=zV['Data'][ikp][0::30]

# plt.close('all')
# plt.plot([0,1200],[0,1200],'k-')
# plt.plot(x,y,'.')
# rs,txt=gu.GetRegStats(x,y)
# plt.plot(rs['xhat Line'],rs['yhat Line'],'r-')
# rs,txt=gu.GetRegStats(x2,y)
# plt.plot(rs['xhat Line'],rs['yhat Line'],'r--')

#%% Age responses 
ikp=np.where( (zB['Data']>0) & (zH['Data']>0) )
bw=1; bin=np.arange(1,60,bw)
#x=zA_VRI['Data'][ikp][0::10]
x=zAc['Data'][ikp][0::10]-zH['Data'][ikp][0::10]

plt.close('all'); 
fig,ax=plt.subplots(1,figsize=gu.cm2inch(11,8)); 
#plt.plot([0,300],[0,300],'k-')
y=zB['Data'][ikp][0::10]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
plt.plot(bin,mu,'o',label='Predictive Forest Inventory')
#y=zB_VRI['Data'][ikp][0::10]
#N,mu,med,sig,se=gu.discres(x,y,bw,bin)
#plt.plot(bin,mu,'s',label='Vegetation Resource Inventory')
ax.set(ylabel='Stemwood biomass (tC/ha)',xlabel='Age, years',ylim=[0,80],xlim=[0,65])
ax.legend(loc='lower right',frameon=False)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=1.5)    
plt.tight_layout()

#%%

# ikp=np.where( (zD['Data']>0) & (zV['Data']>1) & (zVL['Data']>1) & (zA['Data']>0) )
# bw=10; bin=np.arange(0,260,10)
# x=zA['Data'][ikp][0::10]

# plt.close('all')
# y=0.45*zV['Data'][ikp][0::10]
# N,mu,med,sig,se=gu.discres(x,y,bw,bin)
# plt.plot(bin,mu,'o')
# bw=10; bin=np.arange(0,210,10)
# y=zB['Data'][ikp][0::10]
# N,mu,med,sig,se=gu.discres(x,y,bw,bin)
# plt.plot(bin,mu,'s')



#%% Plot chronosequence
a=zH['Data'][zH['Data']>0].flatten()[0::5]
plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(14,10))
ax.hist(a,np.arange(1950,2025,1),facecolor=[0.5,0.75,1])
ax.set(ylabel='Area harvested (hectares/year)',xlabel='Time, years')
ax.plot([2018,2018],[0,30000],'g--',color=[0.5,0,0],lw=3)
ax.legend(loc='lower right',frameon=False)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=1.5)    
plt.tight_layout()

#%% Plot chronosequence
fig,ax=plt.subplots(1,figsize=gu.cm2inch(14,10))
ax.hist(2018-a,np.arange(0,60,1),facecolor=[0.5,0.75,1])
ax.set(ylabel='Area harvested (hectares/year)',xlabel='Age at time of LiDAR (Years)')
ax.plot([0,0],[0,30000],'g--',color=[0.5,0,0],lw=3)

#%% Time since harvest
ikp=np.where( (zV['Data']>1) & (zVL['Data']>1) & (zH['Data']>0) )
bw=1; bin=np.arange(0,260,bw)
x=2018-zH['Data'][ikp][0::10]

plt.close('all');
fig,ax=plt.subplots(1,figsize=gu.cm2inch(11,8)); 
y=zB['Data'][ikp][0::10]
N,mu,med,sig,se=gu.discres(x,y,bw,bin)
plt.plot(bin,mu,'o',label='Predictive Forest Inventory')
y=zB_VRI['Data'][ikp][0::10]
#N,mu,med,sig,se=gu.discres(x,y,bw,bin)
#plt.plot(bin,mu,'s',label='Vegetation Resource Inventory')
#plt.plot(bin,1.*bin,'-',label='Model')
plt.plot(dGY['Age'],dGY['Stemwood biomass (tC/ha)'],'-',label='Model')
ax.set(ylabel='Stemwood biomass (tC/ha)',xlabel='Time since harvest, years',ylim=[0,90],xlim=[0,65])
ax.legend(loc='lower right',frameon=False)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=1.5)    
plt.tight_layout()

