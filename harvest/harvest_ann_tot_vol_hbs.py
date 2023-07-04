
'''

SUMMARIZE ANNUAL VOLUME FROM HBS
*** SEE OTHER SCRIPT FOR VOLUME AND WASTE BY TIMBERMARK ***

'''

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

#%% Parameters

tv_hbs=np.arange(2007,2023,1)

gradeL=[' ','1','2','3','4','5','6','7','8','B','C','D','E','F','G','H','I','J','K','L','M','U','W','X','Y','Z']

gp=gu.SetGraphics('Manuscript')

#%% Import HBS and quantify annual totals

# The system appears to be down often
# -> Date of Invoice, "Billing processed" only
# Can only do 12months at a time

d={}
d['Year']=tv_hbs
d['V All Abs (Mm3/yr)']=np.zeros(tv_hbs.size)
d['V All Logs Abs (Mm3/yr)']=np.zeros(tv_hbs.size)
d['V All NonLogs Abs (Mm3/yr)']=np.zeros(tv_hbs.size)
d['V All Hogged Abs (Mm3/yr)']=np.zeros(tv_hbs.size)
d['V All Chips Abs (Mm3/yr)']=np.zeros(tv_hbs.size)
for Grade in gradeL:
    d['V All Logs Grade ' + Grade + ' (Mm3/yr)']=np.zeros(tv_hbs.size)
d['V Waste All Abs (Mm3/yr)']=np.zeros(tv_hbs.size)

#d['V All (Mm3/yr)']=np.zeros(tv_hbs.size)
#d['V All Logs (Mm3/yr)']=np.zeros(tv_hbs.size)
#d['V All NonLogs (Mm3/yr)']=np.zeros(tv_hbs.size)
#d['V All Hogged (Mm3/yr)']=np.zeros(tv_hbs.size)
#d['V All Chips (Mm3/yr)']=np.zeros(tv_hbs.size)
#d['V Waste All (Mm3/yr)']=np.zeros(tv_hbs.size)
#d['V Waste Logs (Mm3/yr)']=np.zeros(tv_hbs.size)
#d['V Waste NonLogs (Mm3/yr)']=np.zeros(tv_hbs.size)
#d['V Waste Hogged (Mm3/yr)']=np.zeros(tv_hbs.size)
#d['V Waste Chips (Mm3/yr)']=np.zeros(tv_hbs.size)

#d['V NP (Mm3/yr)']=np.zeros(tv_hbs.size)
#d['V NP With Negatives (Mm3/yr)']=np.zeros(tv_hbs.size)
#d['V NP Logs Only (Mm3/yr)']=np.zeros(tv_hbs.size)

for iT in range(tv_hbs.size):

    print(tv_hbs[iT])

    df=pd.read_csv(r'C:\Users\rhember\Documents\Data\Harvest\HBS\HBS ' + str(tv_hbs[iT]) + '.csv',header=None)

    dHB={}
    for i in range(len(df.columns)):
        dHB[i]=df.iloc[:,i].to_numpy()

    Type=dHB[9]
    Material=dHB[6]
    Vol=dHB[11]
    Grade=dHB[7]
    TenureType=dHB[19]
    #print(np.unique(Grade))

    iAllLogs=np.where( (Material=='Logs') )[0]
    iAllNonLogs=np.where( (Material!='Logs') )[0]
    iAllHogged=np.where( (Material=='Hogged Tree Material') )[0]
    iAllChips=np.where( (Material=='Woodchips') )[0]

    iWasteAll=np.where( (Type=='WA') | (Type=='WU') )[0]
    #iWasteLogs=np.where( (Material=='Logs') & (Type=='WA') | (Material=='Logs') & (Type=='WU') )[0]
    #iWasteNonLogs=np.where( (Material!='Logs') & (Type=='WA') | (Material!='Logs') & (Type=='WU') )[0]
    #iWasteHogged=np.where( (Material=='Hogged Tree Material') & (Type=='WA') | (Material=='Hogged Tree Material') & (Type=='WU') )[0]
    #iWasteChips=np.where( (Material=='Woodchips') & (Type=='WA') | (Material=='Woodchips') & (Type=='WU') )[0]

    #ind1=np.where( (Type=='NP') & (Vol>0) )[0]
    #ind1b=np.where( (Type=='NP') )[0]
    #ind2=np.where( (Type=='NP') & (Vol>0) & (Material=='Logs') )[0]

    #d['V All (Mm3/yr)'][iT]=np.sum(Vol)/1e6;
    d['V All Abs (Mm3/yr)'][iT]=np.sum(np.abs(Vol))/1e6
    #d['V All Logs (Mm3/yr)'][iT]=np.sum(Vol[iAllLogs])/1e6;
    d['V All Logs Abs (Mm3/yr)'][iT]=np.sum(np.abs(Vol[iAllLogs]))/1e6
    #d['V All NonLogs (Mm3/yr)'][iT]=np.sum(Vol[iAllNonLogs])/1e6;
    d['V All NonLogs Abs (Mm3/yr)'][iT]=np.sum(np.abs(Vol[iAllNonLogs]))/1e6
    #d['V All Hogged (Mm3/yr)'][iT]=np.sum(Vol[iAllHogged])/1e6;
    d['V All Hogged Abs (Mm3/yr)'][iT]=np.sum(np.abs(Vol[iAllHogged]))/1e6
    #d['V All Chips (Mm3/yr)'][iT]=np.sum(Vol[iAllChips])/1e6;
    d['V All Chips Abs (Mm3/yr)'][iT]=np.sum(np.abs(Vol[iAllChips]))/1e6

    d['V Waste All Abs (Mm3/yr)'][iT]=np.sum(np.abs(Vol[iWasteAll]))/1e6
    for iGrade in range(len(gradeL)):
        indG=np.where( (Grade==gradeL[iGrade])  )[0]
        if indG.size>0:
            d['V All Logs Grade ' + gradeL[iGrade] + ' (Mm3/yr)'][iT]=np.sum(np.abs(Vol[indG]))/1e6

    #d['V NP (Mm3/yr)'][iT]=np.sum(Vol[ind1])/1e6
    #d['V NP With Negatives (Mm3/yr)'][iT]=np.sum(Vol[ind1b])/1e6
    #d['V NP Logs Only (Mm3/yr)'][iT]=np.sum(Vol[ind2])/1e6

#%% Save annual

gu.opickle(r'C:\Users\rhember\Documents\Data\Harvest\HBS\HBS_AnnualSummary.pkl',d)

df=pd.DataFrame(d)
df.to_excel(r'C:\Users\rhember\Documents\Data\Harvest\HBS\HBS_AnnualSummary.xlsx',index=False)

#%% Plot

d=gu.ipickle(r'C:\Users\rhember\Documents\Data\Harvest\HBS\HBS_AnnualSummary.pkl')

plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(14,6));
plt.plot(d['Year'],d['V All Abs (Mm3/yr)'],'-bo')
#plt.plot(d['Year'],d['V All Abs (Mm3/yr)']/1e6,'-rs')
#plt.plot(d['Year'],d['V NP (Mm3/yr)']/1e6,'-g^')
ax.set(position=[0.085,0.125,0.88,0.84],xticks=np.arange(1800,2120,1), \
       yticks=np.arange(0,100,10),ylabel='Volume removed (Million m$^3$ yr$^-$$^1$)',xlabel='Time, years',xlim=[2006.5,2021.5])
#ax.legend(loc='upper left',facecolor=[1,1,1],frameon=False)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=1.5)
#gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Harvest\Harvest Volume HBS\VolumeRemoved_HBS','png',900)


#%%

