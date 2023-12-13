'''
TERRAIN
'''

#%% Import modules

import os
import numpy as np
import time
import matplotlib.pyplot as plt
import matplotlib.colors
import geopandas as gpd
import pandas as pd
import copy
import scipy.io as spio
import fiona
from rasterio import features
from shapely import geometry
from scipy.interpolate import griddata
import cv2
import rasterio
import os
from shapely.geometry import box
from matplotlib import cm
from mpl_toolkits import mplot3d

import fcgadgets.macgyver.util_general as gu
import fcgadgets.macgyver.util_gis as gis
import fcgadgets.macgyver.util_query_gdb as qgdb
import fcgadgets.bc1ha.bc1ha_util as u1ha
import fcgadgets.cbrunner.cbrun_util as cbu
gp=gu.SetGraphics('Manuscript')

#%% Import parameters

# Save ROI grid to file
pth0=r'C:\Users\rhember\Documents\Data\Terrain'
gis.SaveGeoTiff(roi['grd'],pth0 + '\\roi.tif')

# Revise to be 20 m res
gis.ResampleRaster(pth0 + '\\roi.tif',5)
zMask=gis.OpenGeoTiff(pth0 + '\\roi.tif')

#%% Mosaic
def Mosaic(pth0,zMask):
    elev=copy.deepcopy(zMask)
    slope=copy.deepcopy(zMask)
    aspect=copy.deepcopy(zMask)
    floacc=copy.deepcopy(zMask)
    floacc['Data']=floacc['Data'].astype('float32')
    twi=copy.deepcopy(zMask)
    twi['Data']=twi['Data'].astype('float32')
    pth1=r'C:\Users\rhember\Documents\Data\Terrain\CDEM'
    fn=os.listdir(pth1 + '\\Projected')
    for fni in fn:        
        pin=pth1 + '\\Projected\\' + fni
        ds=rasterio.open(pin)
        b=box(*ds.bounds)
        if roi['gdf']['bound'].intersects(b)[0]==True:
            
            pout=pth0 + '\\tmp.tif'
            gis.ClipToRaster_ByFile(pin,pout,pth0 + '\\roi.tif')
    
            tmp=gis.OpenGeoTiff(pth0 + '\\tmp.tif')
            ind=np.where(tmp['Data']>0)
            elev['Data'][ind]=tmp['Data'][ind]
            
            pin2=r'C:\Users\rhember\Documents\Data\Terrain\CDEM\Derived\Aspect' + '\\Aspect' + fni[3:]
            gis.ClipToRaster_ByFile(pin2,pout,pth0 + '\\roi.tif')
            tmp=gis.OpenGeoTiff(pth0 + '\\tmp.tif')
            aspect['Data'][ind]=tmp['Data'][ind]
            
            pin2=r'C:\Users\rhember\Documents\Data\Terrain\CDEM\Derived\Slope' + '\\Slope' + fni[3:]
            gis.ClipToRaster_ByFile(pin2,pout,pth0 + '\\roi.tif')
            tmp=gis.OpenGeoTiff(pth0 + '\\tmp.tif')
            slope['Data'][ind]=tmp['Data'][ind]
            
            pin2=r'C:\Users\rhember\Documents\Data\Terrain\CDEM\Derived\FlowAcc' + '\\FloAcc' + fni[3:]
            gis.ClipToRaster_ByFile(pin2,pout,pth0 + '\\roi.tif')
            tmp=gis.OpenGeoTiff(pth0 + '\\tmp.tif')
            floacc['Data'][ind]=tmp['Data'][ind]
            floacc['Data']=(floacc['Data']+1)*(20**2)
        
    twi['Data']=np.log(np.minimum(1e6,np.maximum(0.01,floacc['Data']/np.tan(np.radians(slope['Data'].astype('float32'))))))
    return elev,slope,aspect,twi

elev,slope,aspect,twi=Mosaic(pth0,zMask)

#%% Fill elevation

p=np.percentile(elev['Data'],0.25)
d=[-1,0,1]
for dx in d:
    for dy in d:
        e=gis.imshift(elev['Data'],dx,dy)
        ind=np.where( (elev['Data']<p) & (e!=0) )
        elev['Data'][ind]=e[ind]

#%%

#plt.matshow(elev['Data'],clim=[200,2000])
#plt.matshow(twi['Data'])
#plt.matshow(aspect['Data'])

#%% Plot TWI

cm0=np.vstack( ((0.32,0.19,0.19,1),(0.41,0.28,0.27,1),(0.49,0.38,0.36,1),(0.58,0.47,0.44,1),(0.66,0.57,0.53,1),(0.75,0.66,0.61,1),(0.83,0.76,0.7,1),(0.92,0.85,0.78,1),(1,0.95,0.87,1),(0.83,0.87,0.85,1),(0.67,0.78,0.82,1),(0.5,0.7,0.8,1),(0.33,0.61,0.78,1),(0.17,0.53,0.76,1),(0,0.45,0.74,1)) )
cm0=matplotlib.colors.LinearSegmentedColormap.from_list('twi',cm0,N=30)
#cmap=cm0(twi['Data']/np.amax(twi['Data']))
cmap=cm0( (twi['Data']-np.amin(twi['Data']))/(np.amax(twi['Data'])-np.amin(twi['Data'])) )

#e=20; a=305
e=20; a=360-2*45
#e=15; a=45
plt.close('all'); fig=plt.figure(figsize=gu.cm2inch(12,12))  
ax=fig.add_subplot(111, projection='3d')  
#ax=plt.axes(projection='3d')
ax.plot_surface(zMask['X'],zMask['Y'],elev['Data'],rstride=2,cstride=2,facecolors=cmap,linewidth=0,antialiased=False,shade=False) # ,
ax.set(position=[0,0,1,1],xlim=[zMask['xmin'],zMask['xmax']],ylim=[zMask['ymin'],zMask['ymax']])
ax.set_zlim(1000,9000)
ax.set_axis_off()
ax.view_init(elev=e,azim=a)
plt.tight_layout()
gu.PrintFig(r'C:\Users\rhember\Documents\Data\Terrain\test_' + str(e) + '_' + str(a),'png',900)
plt.close('all')

#%%


