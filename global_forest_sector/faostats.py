'''
Global Forest sector - Mine FAOSTATS for info
'''
#%% Prepare session
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd
import openpyxl
import copy
import gc as garc
import time
from scipy import stats
from fcgadgets.macgyver import util_general as gu
gp=gu.SetGraphics('Presentation Light')

#%% Import data

# Crosswalk between countries and regions
dfCW=gu.ReadExcel(r'C:\Data\FAOSTAT\CountryRegionCrosswalk.xlsx')

# Forest area
dFA=gu.ReadExcel(r'C:\Data\FAOSTAT\FAOSTAT Summary.xlsx','Area Forest Raw')

# Forest area primary forest
dFAP=gu.ReadExcel(r'C:\Data\FAOSTAT\FAOSTAT Summary.xlsx','Area Primary Forest')

# Forest conversion area
dFAC=gu.ReadExcel(r'C:\Data\FAOSTAT\FAOSTAT Summary.xlsx','Forest Conversion Raw')

# Forest products
dFP=gu.ReadExcel(r'C:\Data\FAOSTAT\FAOSTAT Summary.xlsx','Products Raw')

#%%

u=np.unique(dFA['Area'])
sFA={}
for iu in u:
	#break
	ind=np.where( (dFA['Area']==iu) & (dFA['Year']==2017) )[0]
	sFA[iu]=dFA['Value'][ind]/1000 # Mha
sFAP={}
for iu in u:
	ind=np.where( (dFAP['Area']==iu) & (dFAP['Year']==2017) & (dFAP['Element']=='Area') )[0]
	sFAP[iu]=dFAP['Value'][ind]/1000 # Mha
sFP={} # Converted to tonnes
for iu in u:
	ind=np.where( (dFP['Area']==iu) & (dFP['Item']=='Industrial roundwood') & (dFP['Element']=='Production') )[0]
	sFP[iu]=0.5*dFP['Y2017'][ind]/1e6 # Mt
sWF={} # Converted to tonnes
for iu in u:
	ind=np.where( (dFP['Area']==iu) & (dFP['Item']=='Wood fuel') & (dFP['Element']=='Production') )[0]
	sWF[iu]=0.5*dFP['Y2017'][ind]/1e6 # Mt

#%%
d={}
for iu in u:
	if sFP[iu].size==0:
		continue
	if sWF[iu].size==0:
		continue
	d[iu]=np.column_stack( (sFA[iu][0],sFP[iu][0],sWF[iu][0]) )
print(d['Germany'])
print(d['Canada'])

#%% Actual forest area and efficiency

j='World'
Target=(4.5*1e9)/(sFA[j]*1e6)
Current=sFP[j][0]/1000

cL=['Canada','Argentina','United States of America','Finland','Poland','Turkey','Brazil','Chile','China, mainland','Germany','Pakistan','Indonesia','Ethiopia','Russian Federation'];
y=np.zeros(len(cL))
x=np.zeros(len(cL))
for i in range(len(cL)):
	j=cL[i]
	y[i]=sFP[j][0]/sFA[j][0]
	if j=='Canada':
		x[i]=sFA[j][0]-62
	else:
		x[i]=sFA[j][0]
y=np.append(y,(0.5*74)/62) # British Columbia
x=np.append(x,62) # British Columbia
cL.append('British\nColumbia')

cL2=['Canada','Argentina','USA','Finland','Poland','Turkey','Brazil','Chile','China','Germany','Pakistan','Indonesia','Ethiopia','Russia','British Columbia','Target'];

plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(19.5,10))
plt.plot(x,y,'ko')
for i in range(y.size):
	plt.text(x[i]+10,y[i],cL2[i],va='center',ha='left',fontsize=10,color='k')

xhat=np.arange(0,500,1)
yhat=np.zeros(xhat.size)
for i in range(xhat.size):
	yhat[i]=(4.5*1e9)/(xhat[i]*1e6)
#plt.plot(xhat,yhat,'-',lw=2,color=[1,0.7,0.3])

ax.set(xlabel='Forest area',ylabel='Wood building materials produced\n(tonnes per hectare per year)',xticks=np.arange(0,1200,100),xlim=[0,900])
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=3,direction='out')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
#ax.get_xaxis().set_ticks([])
#ax.get_yaxis().set_ticks([])
plt.tight_layout()
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\GFS_Scatter1','png',900)

#%% Area and efficiency (with conservation)

j='World'
E_Target=(4.5*1e9)/(sFA[j]*1e6)
Yield_Target=4.5
Yield_Current=sFP[j][0]/1000

fTHLB=0.30

cL=['Finland','Poland','Germany','Brazil','Ethiopia','Canada','Argentina','United States of America','Turkey','Chile','China, mainland','Pakistan','Indonesia','Russian Federation'];
Yield=np.zeros(len(cL))
Area=np.zeros(len(cL))
Area_new=np.zeros(len(cL))
E=np.zeros(len(cL))
E_new=np.zeros(len(cL))
for i in range(len(cL)):
	j=cL[i]
	Yield[i]=sFP[j][0]/1000
	E[i]=sFP[j][0]/sFA[j][0]
	if j=='Canada':
		Area[i]=sFA[j][0]-62
		Area_new[i]=fTHLB*(sFA[j][0]-62)
	else:
		Area[i]=sFA[j][0]
		Area_new[i]=fTHLB*sFA[j][0]
E=np.append(E,(0.5*74)/62) # British Columbia
Area=np.append(Area,62) # British Columbia
Area_new=np.append(Area_new,fTHLB*62) # British Columbia
cL.append('British\nColumbia')
cL2=['Finland','Poland','Germany','Brazil','Ethiopia','Canada','Argentina','USA','Turkey','Chile','China','Pakistan','Indonesia','Russia','British Columbia']

adj=E[1]/E
E_new=adj*E
# E_new=E.copy()
# E_new[0]=1.25*E_new[0]
# E_new[1]=1*E_new[1]
# E_new[2]=1.25*E_new[2] # Germany
# E_new[3]=7*E_new[3]
# E_new[4]=7*E_new[4]
# E_new[5]=7*E[5] # Canada
# E_new[6]=4*E[6]
# E_new[7]=4*E[7]
# E_new[8]=4*E[8]
# E_new[9]=2*E[9] # Chile
# E_new[10]=6*E[10] # China
# E_new[11]=5*E[11]
# E_new[12]=5*E[12] # India
# E_new[13]=14*E[13] # Russia
# E_new[14]=4*E[14] # British Columbia

#Yield_Target/(np.sum(E*Area)/1000)
Yield_new=E_new*Area
Yield_new2=E_new*Area_new

plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(19.5,10))
plt.plot(Area,E,'ko',color=[0.8,0.8,0.8])
plt.plot(Area_new,E,'ro',color=[0.7,0,0])
for i in range(E.size):
	plt.plot([Area_new[i],Area[i]],[E[i],E[i]],'r-',color=[0.8,0.8,0.8])
	plt.text(Area[i]+10,E[i],cL2[i],va='center',ha='left',fontsize=10,color='k')

ax.set(xlabel='Forest area',ylabel='Wood building materials produced\n(tonnes per hectare per year)',xticks=np.arange(0,1200,100),xlim=[0,900])
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=3,direction='out')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\GFS_Scatter2','png',900)

# Summary
print(np.sum(Yield_new)/1000)
print(np.sum(Yield_new2)/1000)

#%% Area and efficiency (with conservation and increased efficiency)

plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(19.5,10))
plt.plot(Area_new,E,'ko',color=[0.8,0.8,0.8])
plt.plot(Area_new[E_new-E!=0],E_new[E_new-E!=0],'ro',color=[0.7,0,0])
for i in range(E.size):
	plt.plot([Area_new[i],Area_new[i]],[E[i],E_new[i]],'r-',color=[0.8,0.8,0.8])
	#plt.text(Area_new[i]+10,E_new[i],cL2[i],va='center',ha='left',fontsize=10,color='k')

ax.set(xlabel='Forest area',ylabel='Wood building materials produced\n(tonnes per hectare per year)',xticks=np.arange(0,1200,100),xlim=[0,900])
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=3,direction='out')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\GFS_Scatter3','png',900)

# Summary
print(np.sum(Yield_new)/1000)
print(np.sum(Yield_new2)/1000)

#%% Area and efficiency (with conservation and increased efficiency)

plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(19.5,10))
plt.plot(Area,E,'ko',color=[0.8,0.8,0.8])
plt.plot(Area_new,E_new,'ro',color=[0.7,0,0])
for i in range(E.size):
	plt.plot([Area[i],Area_new[i]],[E[i],E_new[i]],'r-',color=[0.8,0.8,0.8])
	#plt.text(Area_new[i]+10,E_new[i],cL2[i],va='center',ha='left',fontsize=10,color='k')

ax.set(xlabel='Forest area',ylabel='Wood building materials produced\n(tonnes per hectare per year)',xticks=np.arange(0,1200,100),xlim=[0,900])
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=3,direction='out')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\GFS_Scatter4','png',900)

# Summary
print(np.sum(Yield_new)/1000)
print(np.sum(Yield_new2)/1000)

#%%

j='World'
y0=sFP[j]
y1=sFP[j]/(f_thlb*sFA[j])
Target=(4.5*1e9)/(sFA[j]*1e6)

cL=['Africa','Australia and New Zealand','Northern America','South America','Europe','Asia'];
y=np.zeros(len(cL))
x=np.zeros(len(cL))
for i in range(len(cL)):
	j=cL[i]
	y[i]=sFP[j][0]/sFA[j][0]
	if j=='Canada':
		x[i]=sFA[j][0]-62
	else:
		x[i]=sFA[j][0]
#y=np.append(y,(0.5*74)/62) # British Columbia
#x=np.append(x,62) # British Columbia
#cL.append('British\nColumbia')

cL2=['Africa','Australia and New Zealand','Northern America','South America','Europe','Asia','Target'];

plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(19.5,10))
plt.plot(x,y,'ko')
for i in range(y.size):
	plt.text(x[i]+10,y[i],cL2[i],va='center',ha='left',fontsize=10,color='k')

xhat=np.arange(0,500,1)
yhat=np.zeros(xhat.size)
for i in range(xhat.size):
	yhat[i]=(4.5*1e9)/(xhat[i]*1e6)
#plt.plot(xhat,yhat,'-',lw=2,color=[1,0.7,0.3])

ax.set(xlabel='Forest area',ylabel='Wood building materials produced\n(tonnes per hectare per year)',xticks=np.arange(0,1600,100),xlim=[0,1600])
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=3,direction='out')
#ax.spines['top'].set_visible(False)
#ax.spines['right'].set_visible(False)
#ax.get_xaxis().set_ticks([])
#ax.get_yaxis().set_ticks([])
plt.tight_layout()

#%% Barchart
f_thlb=20/62

#j='United States of America'
j='World'
y0=sFP[j]
y1=sFP[j]/(f_thlb*sFA[j])
Target=(4.5*1e9)/(f_thlb*sFA[j]*1e6)

#plt.close('all')
#plt.bar(0,y1)

cL=['World','Canada','United States of America','Finland','Turkey','Brazil','China, mainland','Germany','Pakistan','Ethiopia','Russian Federation'];
y=np.zeros(len(cL))
x=np.zeros(len(cL))
for i in range(len(cL)):
	j=cL[i]
	y[i]=sFP[j][0]/(f_thlb*sFA[j][0])
	x[i]=f_thlb*sFA[j][0]
y=np.append(y,(0.5*74)/22) # British Columbia
x=np.append(x,22) # British Columbia
cL.append('British\nColumbia')

plt.plot(x,y,'ko')

cL2=['Global','Canada','USA','Finland','Turkey','Brazil','China','Germany','Pakistan','Ethiopia','Russia','British\nColumbia','Target'];

ord=np.argsort(y)
y=y[ord]
cL2=np.array(cL2)[ord]

plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(19.5,10))
plt.bar(np.arange(y.size),y,0.9)
plt.bar(y.size,Target,0.9,fc=[1,0.7,0.3])
ax.set(xlim=[-0.5,y.size+0.5],xticks=np.arange(y.size),ylabel='Tonnes per hectare of forest',xticklabels=cL2)
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=3,direction='out')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
#ax.get_xaxis().set_ticks([])
#ax.get_yaxis().set_ticks([])
plt.tight_layout()
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\GlobalIndustrialRoundwoodPerHaForest','png',900)


#%% Compile regions
rL=['Europe','Northern America','South America','Asia','Oceania','Africa']
d={}
for r in rL:
	d[r]={'Area Forest':0.0,'Area Primary Forest':0.0,'Area Deforestation':0.0,'Roundwood Prod':0.0,'Industrial Roundwood Prod':0.0,'Sawnwood Prod':0.0}

# Add forest area (ha)
for i in range(len(dFA)):
	if dFA.loc[i,'Year']!=2020:
		continue
	ind=np.where(dfCW['name']==dFA.loc[i,'Area'])[0]
	if ind.size==0:
		continue
	if np.isnan(dFA.loc[i,'Value'])==True:
		continue
	reg=dfCW.loc[ind,'region'].values[0]
	sreg=dfCW.loc[ind,'sub-region'].values[0]
	ireg=dfCW.loc[ind,'intermediate-region'].values[0]
	if reg not in rL:
		if sreg not in rL:
			if ireg not in rL:
				continue
	if (sreg=='Northern America'):
		reg=sreg
	elif (ireg=='South America'):
		reg=ireg
	d[reg]['Area Forest']=d[reg]['Area Forest']+dFA.loc[i,'Value']*1000
	
# Add primary forest area (2017 is the last year with estimates)
for i in range(len(dFAP)):
	if dFAP.loc[i,'Year']!=2017:
		continue
	ind=np.where(dfCW['name']==dFAP.loc[i,'Area'])[0]
	if ind.size==0:
		continue
	if np.isnan(dFAP.loc[i,'Value'])==True:
		continue
	reg=dfCW.loc[ind,'region'].values[0]
	sreg=dfCW.loc[ind,'sub-region'].values[0]
	ireg=dfCW.loc[ind,'intermediate-region'].values[0]
	if reg not in rL:
		if sreg not in rL:
			if ireg not in rL:
				continue
	if (sreg=='Northern America'):
		reg=sreg
	elif (ireg=='South America'):
		reg=ireg
	d[reg]['Area Primary Forest']=d[reg]['Area Primary Forest']+dFAP.loc[i,'Value']*1000
 
# Add area deforestation
for i in range(len(dFCA)):
	if (dFCA.loc[i,'Year']!=2020) | (dFCA.loc[i,'Element']!='Area'):
		continue
	ind=np.where(dfCW['name']==dFCA.loc[i,'Area'])[0]
	if ind.size==0:
		continue
	if np.isnan(dFCA.loc[i,'Value'])==True:
		continue
	reg=dfCW.loc[ind,'region'].values[0]
	sreg=dfCW.loc[ind,'sub-region'].values[0]
	ireg=dfCW.loc[ind,'intermediate-region'].values[0]
	if reg not in rL:
		if sreg not in rL:
			if ireg not in rL:
				continue
	if (sreg=='Northern America'):
		reg=sreg
	elif (ireg=='South America'):
		reg=ireg
	d[reg]['Area Deforestation']=d[reg]['Area Deforestation']+dFCA.loc[i,'Value']*1000
	  
# Add roundwood production (m3/yr)
for i in range(len(dFP)):
	if (dFP.loc[i,'Item']!='Roundwood') | (dFP.loc[i,'Element']!='Production'):
		continue
	ind=np.where(dfCW['name']==dFP.loc[i,'Area'])[0]
	if ind.size==0:
		continue
	reg=dfCW.loc[ind,'region'].values[0]
	sreg=dfCW.loc[ind,'sub-region'].values[0]
	ireg=dfCW.loc[ind,'intermediate-region'].values[0]
	if reg not in rL:
		if sreg not in rL:
			if ireg not in rL:
				continue
	if (sreg=='Northern America'):
		reg=sreg
	elif (ireg=='South America'):
		reg=ireg
	d[reg]['Roundwood Prod']=d[reg]['Roundwood Prod']+dFP.loc[i,'Y2020']

# Add industrial roundwood production (m3/yr)
for i in range(len(dFP)):
	if (dFP.loc[i,'Item']!='Industrial roundwood') | (dFP.loc[i,'Element']!='Production'):
		continue
	ind=np.where(dfCW['name']==dFP.loc[i,'Area'])[0]
	if ind.size==0:
		continue
	reg=dfCW.loc[ind,'region'].values[0]
	sreg=dfCW.loc[ind,'sub-region'].values[0]
	ireg=dfCW.loc[ind,'intermediate-region'].values[0]
	if reg not in rL:
		if sreg not in rL:
			if ireg not in rL:
				continue
	if (sreg=='Northern America'):
		reg=sreg
	elif (ireg=='South America'):
		reg=ireg
	d[reg]['Industrial Roundwood Prod']=d[reg]['Industrial Roundwood Prod']+dFP.loc[i,'Y2020']
	   
# Add sawnwood production (m3/yr)
for i in range(len(dFP)):
	if (dFP.loc[i,'Item']!='Sawnwood') | (dFP.loc[i,'Element']!='Production'):
		continue
	ind=np.where(dfCW['name']==dFP.loc[i,'Area'])[0]
	if ind.size==0:
		continue
	reg=dfCW.loc[ind,'region'].values[0]
	sreg=dfCW.loc[ind,'sub-region'].values[0]
	ireg=dfCW.loc[ind,'intermediate-region'].values[0]
	if reg not in rL:
		if sreg not in rL:
			if ireg not in rL:
				continue
	if (sreg=='Northern America'):
		reg=sreg
	elif (ireg=='South America'):
		reg=ireg
	d[reg]['Sawnwood Prod']=d[reg]['Sawnwood Prod']+dFP.loc[i,'Y2020']

# Add fuelwood
for reg in d.keys():
	d[reg]['Fuelwood']=d[reg]['Roundwood Prod']-d[reg]['Industrial Roundwood Prod']

# Save to spreadsheet
df=pd.DataFrame(d)
df.to_excel(r'C:\Users\rhember\Documents\Data\FAOSTAT\SummaryStats.xlsx')



#%% Plot

fw=np.array([])
irw=np.array([])
for reg in d.keys():
	fw=np.append(fw,d[reg]['Fuelwood'])
	irw=np.append(irw,d[reg]['Industrial Roundwood Prod'])

plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,6))
ax.bar(np.arange(fw.size),fw/1e6,fc=[0.7,0.8,1])
ax.bar(np.arange(fw.size),irw/1e6,bottom=fw/1e6,fc=[0.8,1,0.7])
ax.set(position=[0.09,0.12,0.88,0.82],xticks=np.arange(fw.size),xticklabels=np.array(list(d.keys())),ylabel='Million m3/year')
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=1.5)
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Global Forest Sector\RegionalSummaryHarvest','png',500)
	
