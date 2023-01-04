'''

CONVERT CHARCOAL RECORD TO ANNUAL PROBABILITY OF WILDFIRE

'''

#%% Import Python modules

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import scipy.stats as stats 
from fcgadgets.macgyver import utilities_general as gu

#%% Import charcoal record digitized from Marlon et al. (2012)
# This must already be a moving average of the z-score as it does not exceed abs(1).

#pth=r'C:\Users\rhember\Documents\Data\Wildfire\WildfireReconstructionMarlon2012\WildfireReconstructionMarlon2012.xlsx'
#dfC=pd.read_excel(pth)
#
#tvC=np.arange(500,2000,1)
#yC=np.interp(tvC,dfC['Time'].values,dfC['Percent area burned'].values)
#plt.plot(tvC,yC,'.')
#
## Adjust so that it is relative mean during observed wildfire record
#it_mu=np.where(tvC>=1920)[0]
#yC_adj=(yC-np.mean(yC[it_mu]))*1
#
#plt.plot(tvC,yC_adj,'.')
#plt.grid('on')


pth=r'C:\Users\rhember\Documents\Data\DigitizedFigures\Marlonetal2012_Fig2.xlsx'
dfC=pd.read_excel(pth)

tvC=np.arange(500,2000,1)
yC=np.interp(tvC,dfC['Time'].values,dfC['Charcoal influx zscore'].values)
plt.close('all')
plt.plot(tvC,yC,'b-',lw=1.25)

# Adjust so that it is relative mean during observed wildfire record
it_mu=np.where(tvC>=1920)[0]
yC_adj=(yC-np.mean(yC[it_mu]))*1

plt.plot(tvC,yC_adj,'g-',lw=1.25)
plt.grid('on')

#%% Import wildfire disturbances

#pth=r'G:\My Drive\Data\FCI_RollupBySparseGrid\Outputs\Scenario0001\Metadata.pkl'
#pth=r'C:\Users\rhember\Documents\Data\FCI_Projects\FCI_RollupBySparseGrid\Outputs\Scenario0001\Metadata.pkl'
#fin=open(pth,'rb')
#meta=pickle.load(fin); fin.close()
meta=gu.ipickle(r'C:\Users\rhember\Documents\Data\FCI_Projects\FCI_SparseGrid\Outputs\Scenario0001\Metadata.pkl')

iScn=0; iEns=0; iBat=0
        
# Import disturbance history
dh=gu.ipickle(r'C:\Users\rhember\Documents\Data\FCI_Projects\FCI_SparseGrid\Inputs\Scenario0001\Disturbance_Ens0001_Bat0001.pkl')
#dh.append(gu.ipickle(r'C:\Users\rhember\Documents\Data\FCI_Projects\FCI_SparseGrid\Inputs\Scenario0001\Disturbance_Ens0001_Bat0002.pkl'))
#dh.append(gu.ipickle(r'C:\Users\rhember\Documents\Data\FCI_Projects\FCI_SparseGrid\Inputs\Scenario0001\Disturbance_Ens0001_Bat0003.pkl'))
    
DistTypes=list(meta['LUT Dist'].keys())
DistID=np.array(list(meta['LUT Dist'].values()))        
        
# Calculate area disturbed by disturbance type (hectares)
A=np.zeros((meta['Time'].size,len(DistTypes)))
for i in range(len(dh)):        
    for j in range(dh[i]['Year'].size):            
        iT=np.where(meta['Time']==dh[i]['Year'][j])[0]            
        iD=np.where(DistID==dh[i]['ID_Type'][j])[0]
        A[iT,iD]=A[iT,iD]+1        
AR=A/meta['N Stand']*100 # Convert to percent
AR_ma=gu.movingave(AR,10,'Historical')

it=np.where( (meta['Time']>=1920) & (meta['Time']<=2017) )[0]
tvF=meta['Time'][it]
yF=AR[it,0]
yF_ma=AR_ma[it,0]

#%% Fit distribution model to annual area burned

plt.close('all')
bins=np.arange(0,1.6,0.1)
h=plt.hist(yF,bins=bins,density=True)

# Exponential fit
loc,scale=stats.expon.fit(yF)
x=np.linspace(0,1.5,1000)
fx=stats.expon.pdf(x,loc,scale)
plt.close('all')
h=plt.hist(yF,bins=bins,density=True,facecolor=[0.8,0.8,0.8])
plt.plot(x,fx,'b-',linewidth=1.5,color=[0,0,0.8])
pval=stats.kstest(yF,"expon",args=(loc,scale))

# GEV
loc,scale,size=stats.genextreme.fit(yF)
x=np.linspace(0,1.5,1000)
fx=stats.genextreme.pdf(x,loc,scale,size)
plt.close('all')
h=plt.hist(yF,bins=bins,density=True,facecolor=[0.8,0.8,0.8])
plt.plot(x,fx,'b-',linewidth=1.5,color=[0,0,0.8])
pval=stats.kstest(yF,"genextreme",args=(loc,scale,size))

# Weibull
shape,loc,scale=stats.weibull_min.fit(yF)
x=np.linspace(0,1.5,1000)
fx=stats.weibull_min.pdf(x,shape,loc,scale)
plt.close('all')
h=plt.hist(yF,bins=bins,density=True,facecolor=[0.8,0.8,0.8])
plt.plot(x,fx,'b-',linewidth=1.5,color=[0,0,0.8])
pval=stats.kstest(yF,"weibull_min",args=(shape,loc,scale))

#%% Create final best-availabe estimate

tvFC=np.arange(-2000,2018,1)
yFC=np.zeros(tvFC.size)
yFC_adj=np.zeros(tvFC.size)

it=np.where((tvFC>=tvF[0]) & (tvFC<=tvF[-1]) )[0]
yFC[it]=yF

it=np.where(tvFC<1920)[0]
yFC[it]=stats.expon.rvs(loc=loc,scale=scale,size=it.size)

yFC_ma=gu.movingave(yFC,30,'Centre')

# Create an adjusted charcoal record (forced to match)
d_ma=np.zeros(tvFC.size)

# Add charcoal period
it=np.where( (tvFC>=tvC[0]) & (tvFC<=tvC[-1]) )[0]
d_ma=yFC_ma[it]-yC_adj
yFC_adj[it]=np.maximum(0,yFC[it]-d_ma)

# Add modern period
it=np.where(tvFC>=1920)[0]
yFC_adj[it]=yF

# Add pre-charcoal period
it1=np.where(tvFC<tvC[0])[0]
it2=np.where(tvFC>=tvC[0])[0]
yFC_adj[it1]=yFC[it1]+np.mean(yFC[it1])-np.mean(yFC[it2])

# Calculate moving average
yFC_adj_ma=gu.movingave(yFC_adj,30,'Centre')

# Plot
itPlot=np.where(tvFC>=500)
plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,6))
ax.plot(tvFC[itPlot],yFC[itPlot],'k-',color=[0.7,0.7,0.7])
ax.plot(tvFC[itPlot],yFC_ma[itPlot],'k-',linewidth=1.5)
ax.plot(tvFC[itPlot],yFC_adj[itPlot],'c-',color=[0.7,1,1],alpha=0.5)
ax.plot(tvFC[itPlot],yFC_adj_ma[itPlot],'c-',color=[0,0.8,0.8],linewidth=1.5,alpha=0.5)
ax.set(position=[0.07,0.13,0.9,0.84],xlim=[500,2000],ylim=[0,1.6],ylabel='Area burned (%)',xlabel='Time, years')
ax.set_xticks(np.arange(500,2000,100))
ax.grid()

#%% Save to file

df=pd.DataFrame(data=np.column_stack((tvFC,yFC_adj)),columns=['Time','Percent area burned'])
df.to_excel(r'I:\My Drive\Data\WildfireReconstructionMarlon2012\WildfireReconstructionMarlon2012.xlsx')

# Compare final reconstruction (zscore) againzed original charcoal zscore
yFC_adj_z=(yFC_adj-np.mean(yFC_adj))/np.std(yFC_adj)
yFC_adj_z_ma=bf.movingave(yFC_adj_z,25,'Centre')
plt.figure(2)
plt.plot(tvC,yC,'k-')
plt.plot(tvFC,yFC_adj_z_ma,'c-',color=[0,0.8,0.8],linewidth=1.5,alpha=0.5)




#%% Look at exponential area burned

tv=np.arange(0,2021,1)
loc=0.00001
scale=0.007
y=np.minimum(1,stats.expon.rvs(loc=loc,scale=scale,size=tv.size))
print(100*np.median(y))
print(100*np.mean(y))

plt.close('all')
plt.bar(tv,100*y,1,facecolor=[0.7,0,0])

#%% Look at Weibull area burned

tv=np.arange(1000,2021,1)

shape=0.0001
loc=0.00001
scale=0.000001
#y0=stats.weibull_min.pdf(tv,shape,loc,scale)
y0=stats.weibull_min.rvs(shape,loc=loc,scale=scale,size=tv.size)
y=np.minimum(1,y0)
print(100*np.median(y))
print(100*np.mean(y))

plt.close('all')
plt.bar(tv,100*y,1,facecolor=[0.7,0,0])


#%% Look at Generalize Pareto area burned

import scipy.stats as stats

tv=np.arange(0,2021,1)

shape=1.1
loc=0.06
scale=0.02
#y0=stats.weibull_min.pdf(tv,shape,loc,scale)
y=stats.genpareto.rvs(shape,loc=loc,scale=scale,size=tv.size)
print(np.median(y))
print(np.mean(y))

plt.close('all')
plt.bar(tv,y,1,facecolor=[0.7,0,0])




