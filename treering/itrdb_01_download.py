#%% Import modules
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import urllib
from urllib.request import urlretrieve
from bs4 import BeautifulSoup
import requests
import os
import time
import fcgadgets.macgyver.util_general as gu

#%% Download raw tree ring data from ITRDB
def get_files(url):
   page = requests.get(url).text  
   soup = BeautifulSoup(page, 'html.parser')
   return [url + '/' + node.get('href') for node in soup.find_all('a') if 
           node.get('href').endswith('.rwl')]

# North America
natL=['canada','usa','mexico']
for nat in natL:
    url=r'https://www.ncei.noaa.gov/pub/data/paleo/treering/measurements/northamerica' + '//' + nat + '//'
    fl=get_files(url)
    for f in fl:
        if f[-4:]=='.rwl':
            nam=f[f.rfind('//')+2:]
            if nam.rfind('noaa')>=0:
                continue
            fout=r'C:\Data\TreeRings\ITRDB\2023\Given' + '\\' + nam
            urlretrieve(f,fout)

# Europe, Asia
natL=['europe','asia']
for nat in natL:
    url=r'https://www.ncei.noaa.gov/pub/data/paleo/treering/measurements' + '//' + nat + '//'
    fl=get_files(url)
    for f in fl:
        if f[-4:]=='.rwl':
            nam=f[f.rfind('//')+2:]
            if nam.rfind('noaa')>=0:
                continue
            fout=r'C:\Data\TreeRings\ITRDB\2023\Given' + '\\' + nam
            urlretrieve(f,fout)

#%% Import raw data
# Only ring width, density and early/late width excluded

f=os.listdir(r'C:\Data\TreeRings\ITRDB\2023\Given')
bw=[6,2,4,1,6,6,6,6,6,6,6,6,6,6]
n=50000000
d={'ID_Site':np.array(['' for _ in range(n)],dtype='object'),
   'ID_Tree':np.array(['' for _ in range(n)],dtype='object'),
   'ID_Tree_Unique':np.zeros(n,dtype='int32'),
   'ID_Species':np.zeros(n,dtype='int16'),
   'Lat':np.zeros(n,dtype='float'),
   'Lon':np.zeros(n,dtype='float'),
   'Year':np.zeros(n,dtype='int16'),
   'Age':np.zeros(n,dtype='float'),
   'Width':np.zeros(n,dtype='float')}
cnt0=0
for fn in f:    
    if (fn.rfind('e.rwl')>0) | (fn.rfind('i.rwl')>0) | (fn.rfind('l.rwl')>0) | (fn.rfind('n.rwl')>0) | (fn.rfind('p.rwl')>0) | (fn.rfind('t.rwl')>0) | (fn.rfind('x.rwl')>0):
        print('Skipping')
        continue
    
    df=pd.read_fwf(r'C:\Data\TreeRings\ITRDB\2023\Given' + '\\' + fn,widths=bw,header=2)
    try:
        W=df.iloc[:,4:].values
        W=np.char.replace(W.astype(str),' -','')
        W=np.char.replace(W.astype(str),'-','')
        W=np.char.replace(W.astype(str),'.','')    
        W=np.array(W,dtype='float64')
        W=np.reshape(W,(W.size,1))[:,0]
    except:
        print('Not working')
        continue
    m=W.size
    ID_Site=np.array([fn[0:fn.rfind('.')] for _ in range(m)],dtype=object)
    ID_Tree=np.array([' ' for _ in range(m)],dtype=object)
    Year=np.zeros(m)    
    cnt=0
    for i in range(len(df)):
        ID_Tree[cnt:cnt+10]=str(df.iloc[i,0])
        yr=df.iloc[i,2]
        if (yr<=0) | (np.isnan(yr)==True):
            continue
        Year[cnt:cnt+10]=np.arange(yr,yr+10,1)
        cnt=cnt+10
    d['ID_Site'][cnt0:cnt0+m]=ID_Site
    d['ID_Tree'][cnt0:cnt0+m]=ID_Tree
    d['Year'][cnt0:cnt0+m]=Year
    d['Width'][cnt0:cnt0+m]=W
    cnt0=cnt0+m
    print(cnt0)
 
# Truncate
for k in d.keys():
    d[k]=d[k][0:cnt0+m]

# Save
gu.opickle(r'C:\Data\TreeRings\ITRDB\2023\Processed\ITRDB_Data_L0.pkl',d)

#%% Calculate age and diameter

d=gu.ipickle(r'C:\Data\TreeRings\ITRDB\2023\Processed\ITRDB_Data_L0.pkl')

dM=gu.ReadExcel(r'C:\Data\TreeRings\ITRDB\2023\ITRDBmetadata12January2022.xlsx')
dSpc=gu.ReadExcel(r'C:\Data\TreeRings\ITRDB\2023\SpeciesCrosswalk.xlsx')

uS=np.unique(d['ID_Site'])

d['Dib']=np.nan*np.ones(d['Width'].size)
d['Dob']=np.nan*np.ones(d['Width'].size)
d['Bsw']=np.nan*np.ones(d['Width'].size)
d['Gsw']=np.nan*np.ones(d['Width'].size)
d['YearFirst']=np.nan*np.ones(d['Width'].size)
d['YearLast']=np.nan*np.ones(d['Width'].size)
cnt_TreeU=0
for iS in range(uS.size):
    ind0=np.where(d['ID_Site']==uS[iS])[0]
    ind1=np.where(dM['ITRDB_Code']==np.char.upper(uS[iS]))[0]
    if ind1.size==0:
        continue    
    d['Lat'][ind0]=dM['Latitude'][ind1]
    d['Lon'][ind0]=dM['Longitude'][ind1]
    indSpc=np.where(dSpc['Code']==dM['Species'][ind1])[0]
    d['ID_Species'][ind0]=dSpc['ID'][indSpc]
    idT=d['ID_Tree'][ind0]
    yrT=d['Year'][ind0]
    uT=np.unique(idT)
    
    d2b1=dSpc['d2b1'][indSpc]
    d2b2=dSpc['d2b2'][indSpc]
    
    #break
    for iT in range(uT.size):
        ind4=np.where(idT==uT[iT])[0]
        W=d['Width'][ind0[ind4]]
        
        # Unique tree ID
        d['ID_Tree_Unique'][ind0[ind4]]=cnt_TreeU        
        
        # Remove outliers and interpolate missing years
        mu=np.nanmean(W)
        sd=np.nanstd(W)
        z=(W-mu)/sd
        #plt.plot(z,'.b-')
        iBad=np.where(np.abs(z)>3)[0]
        W[iBad]=np.nan
        #plt.plot(W,'.b-')
        iGood=np.where( (np.isnan(W)==False) & (np.abs(z)<3) )[0]
        if iGood.size==0:
            continue
        A=np.arange(1,W.size+1,1)
        Wi=np.interp(A,A[iGood],W[iGood])        
        #plt.plot(Wi,'o')
        
        # Calculate inside-bark diameter (preliminary)
        Dib=np.cumsum(2*Wi/10) # cm
        
        # Apply scale factor to get (mm/yr)
        mx=np.nanmax(Dib)
        if (mx>=10000000):
            Wi=0.000001*Wi
        elif (mx>=1000000):
            Wi=0.00001*Wi
        elif (mx>=100000):
            Wi=0.0001*Wi
        elif (mx>=10000) & (mx<=100000):
            Wi=0.001*Wi
        elif (mx>=1000) & (mx<=10000):
            Wi=0.01*Wi
        elif (mx>=100) & (mx<=1000):
            Wi=0.1*Wi
        else:
            Wi=1.0*Wi
        
        # Calculate inside-bark diameter
        Dib=np.cumsum(2*Wi/10) # cm
        #print(np.nanmean(Dib))
        #plt.close('all'); plt.plot(A,Dib,'-b.')
        d['Age'][ind0[ind4]]=A
        d['Dib'][ind0[ind4]]=Dib
        
        Dob=(1/0.944)*Dib
        Bsw=d2b1*Dob**d2b2
        d['Dob'][ind0[ind4]]=Dob
        d['Bsw'][ind0[ind4]]=Bsw
        d['Gsw'][ind0[ind4]]=np.append(np.nan,np.diff(Bsw))
        
        d['YearFirst'][ind0[ind4]]=np.nanmin(yrT[ind4])
        d['YearLast'][ind0[ind4]]=np.nanmax(yrT[ind4])
        
        # Update counter
        cnt_TreeU=cnt_TreeU+1
        
    print(iS)

# Save
gu.opickle(r'C:\Data\TreeRings\ITRDB\2023\Processed\ITRDB_Data_L1.pkl',d)
   
#%% Calculate time series for each sample

d=gu.ipickle(r'C:\Data\TreeRings\ITRDB\2023\Processed\ITRDB_Data.pkl')

uSP=np.unique(d['ID_Species'][d['ID_Species']>0])
dT={}
dT['Year']=np.arange(1800,2021,1)
dT['Gsw']=np.zeros((dT['Year'].size,uSP.size))
for iS in range(uSP.size):
    print(iS)
    ind0=np.where( (d['ID_Species']==uSP[iS]) &
                 (d['Gsw']>0) & (d['Gsw']<100) &
                 (d['Age']>0) & (d['Age']<300) &
                 (d['YearFirst']<=1800) & (d['YearLast']>=2000) )[0]
    if ind0.size==0:
        continue
    Year0=d['Year'][ind0]
    Gsw0=d['Gsw'][ind0]
    indT=gu.IndicesFromUniqueArrayValues(Year0)
    for k in indT.keys():
        iT=np.where(dT['Year']==k)[0]
        dT['Gsw'][iT,iS]=np.nanmean(Gsw0[indT[k]])

#
iSP=2
ind=np.where(dSpc['ID']==uSP[iSP])[0]
print(dSpc['Code'][ind][0])
plt.close('all')
plt.plot(dT['Year'],dT['Gsw'][:,iSP],'-ko')

plt.close('all')
plt.plot(dT['Year'],np.nanmean(dT['Gsw'],axis=1),'-ko')


