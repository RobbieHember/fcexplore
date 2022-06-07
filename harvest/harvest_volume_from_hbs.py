#%% Import modules

from os import listdir
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from fcgadgets.macgyver import utilities_general as gu
import scipy.io
import copy
import geopandas as gpd
import fiona
import scipy.io as spio
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats, linalg
from fcgadgets.macgyver import utilities_inventory as invu
from fcgadgets.cbrunner import cbrun_utilities as cbu

#%% Graphics parameters

fs=7
params={'font.sans-serif':'Arial','font.size':fs,'axes.edgecolor':'black','axes.labelsize':fs,'axes.labelcolor':'black','axes.titlesize':fs,'axes.linewidth':0.5,'lines.linewidth':0.5,
        'text.color':'black','xtick.color':'black','xtick.labelsize':fs,'xtick.major.width':0.5,'xtick.major.size':3,'xtick.direction':'in','ytick.color':'black','ytick.labelsize':fs,
        'ytick.major.width':0.5,'ytick.major.size':3,'ytick.direction':'in','legend.fontsize':fs,'savefig.dpi':300,'savefig.transparent':True}
plt.rcParams.update(params)
   
#%% Import HBS 
# The system appears to be down often
# -> Date of Invoice, "Billing processed" only 
# Can only do 12months at a time

tv=np.arange(2007,2022,1)

d={}
d['V All m3']=np.zeros(tv.size)
d['V NP m3']=np.zeros(tv.size)
d['V NP With Negatives m3']=np.zeros(tv.size)
d['V NP Logs Only m3']=np.zeros(tv.size)

for iT in range(tv.size):
    print(tv[iT])
    df=pd.read_csv(r'C:\Users\rhember\Documents\Data\Harvest\HBS\HBS ' + str(tv[iT]) + '.csv',header=None)
    
    dHB={}
    for i in range(len(df.columns)):
        dHB[i]=df.iloc[:,i].to_numpy()
    
    Type=dHB[9]
    Material=dHB[6]
    Vol=dHB[11]
    
    indAll=np.where( (Vol>0) )[0]
    ind1=np.where( (Type=='NP') & (Vol>0) )[0]
    ind1b=np.where( (Type=='NP') )[0]
    ind2=np.where( (Type=='NP') & (Vol>0) & (Material=='Logs') )[0]
    
    d['V All m3'][iT]=np.sum(Vol[indAll])
    d['V NP m3'][iT]=np.sum(Vol[ind1])
    d['V NP With Negatives m3'][iT]=np.sum(Vol[ind1b])
    d['V NP Logs Only m3'][iT]=np.sum(Vol[ind2])

d['Year']=tv

#%% Save

gu.opickle(r'C:\Users\rhember\Documents\Data\Harvest\HBS\AnnualSum.pkl',d)

#%% 

plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(9,6)); 
plt.plot(tv,d['V NP m3']/1e6,'-bo')
ax.set(position=[0.085,0.125,0.88,0.84],xlim=[2006.5,2021.5],xticks=np.arange(1800,2120,1), \
       yticks=np.arange(0,100,10),ylabel='Volume removed (Million m$^3$ yr$^-$$^1$)',xlabel='Time, years')
#ax.legend(loc='upper left',facecolor=[1,1,1],frameon=False)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=1.5)
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Harvest\Harvest Volume HBS\VolumeRemoved_HBS','png',900)


#%%

