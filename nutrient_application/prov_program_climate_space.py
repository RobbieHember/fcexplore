'''

CLIMATE SPACE OF THE FERTILIZATION PROGRAM

'''

#%% Import modules

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import fcgadgets.macgyver.utilities_general as gu
import fcgadgets.macgyver.utilities_gis as gis
import fcgadgets.bc1ha.bc1ha_utilities as u1ha

#%% Import parameters

meta=u1ha.Init()
meta=u1ha.ImportLUTs(meta)
gp=gu.SetGraphics('Manuscript')

#%% Import data

vList=['refg','feca_yr','tmean_ann_n','prcp_ann_n']
meta['Graphics']['Map']['RGSF']=5
z=u1ha.Import_Raster(meta,[],vList)

z['tmean_ann_n']['Data']=z['tmean_ann_n']['Data'].astype('float')/10

#%%

#ikp=np.where( (z['refg']['Data']==1) & (z['tmean_ann_n']['Data']>-50) & (z['prcp_ann_n']['Data']>=0) )
#plt.hist(z['tmean_ann_n']['Data'][ikp].flatten()[0::100])
#plt.hist(z['prcp_ann_n']['Data'][ikp].flatten()[0::100])

ikp=np.where( (z['refg']['Data']==1) & (z['tmean_ann_n']['Data']>-50) & (z['prcp_ann_n']['Data']>=0) & (z['feca_yr']['Data']>0) )
xmin=np.min(z['tmean_ann_n']['Data'][ikp])
xmax=np.max(z['tmean_ann_n']['Data'][ikp])
ymin=np.min(z['prcp_ann_n']['Data'][ikp])
ymax=3000#np.max(z['prcp_ann_n']['Data'][ikp])

m1=z['tmean_ann_n']['Data'][ikp]
m2=z['prcp_ann_n']['Data'][ikp]
X,Y=np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
positions=np.vstack([X.ravel(),Y.ravel()])
values=np.vstack([m1,m2])
kernel=stats.gaussian_kde(values)
Z=np.reshape(kernel(positions).T, X.shape)

#plt.hist(z['tmean_ann_n']['Data'][ikp])
#plt.hist(z['prcp_ann_n']['Data'][ikp])

#%%

plt.close('all')
#plt.matshow(z['tmean_ann_n']['Data'],clim=[-1,10]); plt.colorbar()
#plt.matshow(z['prcp_ann_n']['Data'],clim=[0,2000]); plt.colorbar()
fig=plt.figure()
ax=fig.add_subplot(111)
ax.matshow(np.rot90(Z).T,cmap=plt.cm.gist_earth,extent=[xmin,xmax,ymin,ymax],aspect='auto')
#ax.plot(m1, m2, 'k.', markersize=2)
ax.set(xlabel='Mean annual temperature (degC)',ylabel='Mean annual precipitation (mm/yr)',xticks=np.arange(-100,100,1),yticks=np.arange(0,10000,200),xlim=[xmin,xmax],ylim=[ymin,ymax])
plt.gca().xaxis.tick_bottom()
plt.tight_layout()
plt.show()
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Fertilization\ClimateSpace\FECA_ClimateSpace','png',300)

