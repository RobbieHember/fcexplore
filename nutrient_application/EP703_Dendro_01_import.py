'''
EP703 RETROSPECTIVE MONITORING STUDY - IMPORT DATA
Notes:
	The version of this used in Hember et al. 2023 was frozen
'''
#%% Import modules
from os import listdir
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.io
import copy
from scipy.optimize import curve_fit
import scipy.io as spio
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats, linalg
from matplotlib.patches import Rectangle
from fcgadgets.macgyver import util_general as gu

#%% Graphics parameters
gp=gu.SetGraphics('Manuscript')

#%% Preapre project
meta={}
meta['Paths']={}
meta['Paths']['Figures']=r'C:\Users\rhember\OneDrive - Government of BC\Figures\EP703\Retrospective Monitoring Study'

#%% Import raw data

# Import tree-level observations
dTL=gu.ReadExcel(r'C:\Users\rhember\OneDrive - Government of BC\Shared\EP703\Retrospective Monitoring Study\Data\Given\Core Data\Received 20211105 from JBerg\EP703_RM_Master_stacked_pith_est.xlsx')

# Import tree summary file
dTS=gu.ReadExcel(r'C:\Users\rhember\OneDrive - Government of BC\Shared\EP703\Retrospective Monitoring Study\Data\Given\Field Data\Received 20210125 from JBapty\RM Tree Summary.xlsx')

# Import plot summary file
dPS=gu.ReadExcel(r'C:\Users\rhember\OneDrive - Government of BC\Shared\EP703\Retrospective Monitoring Study\Data\Given\Field Data\Received 20210125 from JBapty\RM Plot Summary.xlsx')

# Climate data
clm=gu.ImportMat(r'C:\Users\rhember\OneDrive - Government of BC\Shared\EP703\Retrospective Monitoring Study\Data\Given\Climate\Received 20211024 from RHember\EP703_TR_Climate.mat','clm')

# Import EP703 tree-level data
dEPT=gu.ipickle(r'C:\Data\EP703\Processed\EP703_TL.pkl')

#%% Custom functions

# Function (tree height as a function of basal area)
def funcH(x,a,b):
	return a*x**b

#%% Add derived variables

# Add unique tree ID
dTL['ID_Tree_Unique']=np.zeros(dTL['TRW'].size,dtype=int)
cnt=0
u=np.unique(np.column_stack((dTL['ID_Inst'],dTL['ID_Plot'],dTL['ID_Tree'])),axis=0)
for iU in range(u.shape[0]):
	ind=np.where( (dTL['ID_Inst']==u[iU,0]) & (dTL['ID_Plot']==u[iU,1]) & (dTL['ID_Tree']==u[iU,2]) )[0]
	dTL['ID_Tree_Unique'][ind]=cnt
	cnt=cnt+1

# Add unique plot ID
dTL['ID_Plot_Unique']=np.zeros(dTL['TRW'].size,dtype=int)
cnt=0
u=np.unique(np.column_stack( (dTL['ID_Inst'],dTL['ID_Plot']) ),axis=0)
for iU in range(u.shape[0]):
	ind=np.where( (dTL['ID_Inst']==u[iU,0]) & (dTL['ID_Plot']==u[iU,1]) )[0]
	dTL['ID_Plot_Unique'][ind]=cnt
	cnt=cnt+1

# Adjust units of BAI (cm2/yr)
dTL['BAI']=dTL['BAI']/100

# Add info from site and tree summary files and add some derived variables
dTL['BGC SS']=np.zeros(dTL['TRW'].size)
dTL['Dob 2020']=np.zeros(dTL['TRW'].size)
dTL['H 2020']=np.zeros(dTL['TRW'].size)
dTL['Year First']=np.zeros(dTL['TRW'].size)
dTL['Year Last']=np.zeros(dTL['TRW'].size)
dTL['Age']=np.zeros(dTL['TRW'].size)
dTL['Dib']=np.zeros(dTL['TRW'].size)
dTL['Dib Last']=np.zeros(dTL['TRW'].size)
dTL['H obs EP']=np.zeros(dTL['TRW'].size)
dTL['Dob obs EP']=np.zeros(dTL['TRW'].size)
dTL['Bsw']=np.zeros(dTL['TRW'].size)
dTL['Gsw']=np.zeros(dTL['TRW'].size)
dTL['RGR']=np.zeros(dTL['TRW'].size)
dTL['TRW RR']=np.zeros(dTL['TRW'].size)
dTL['BAI RR']=np.zeros(dTL['TRW'].size)
dTL['Bsw RR']=np.zeros(dTL['TRW'].size)
dTL['Gsw AD']=np.zeros(dTL['TRW'].size)
dTL['Gsw RR']=np.zeros(dTL['TRW'].size)
dTL['RGR RR']=np.zeros(dTL['TRW'].size)
dTL['Lat']=np.zeros(dTL['TRW'].size)
dTL['Lon']=np.zeros(dTL['TRW'].size)

uTID=np.unique(dTL['ID_Tree_Unique'])
for iU in range(uTID.size):

	# Index to unique tree ID
	indUT=np.where(dTL['ID_Tree_Unique']==uTID[iU])[0]

	# Define first and last year of measurement for each tree
	iWithData=np.where(dTL['Dib'][indUT]>-1.0)[0]
	dTL['Year First'][indUT]=dTL['Year'][indUT[iWithData[0]]]
	dTL['Year Last'][indUT]=dTL['Year'][indUT[iWithData[-1]]]
	dTL['Age'][indUT[iWithData]]=np.arange(1,iWithData.size+1,1)

	# Index to site summary info
	indPS=np.where( (dPS['ID_Inst']==dTL['ID_Inst'][indUT[0]]) & (dPS['ID_Plot']==dTL['ID_Plot'][indUT[0]]) )[0]

	# Site series
	dTL['BGC SS'][indUT]=dPS['BGC site series'][indPS]
	dTL['Lat'][indUT]=dPS['Lat'][indPS]
	dTL['Lon'][indUT]=dPS['Lon'][indPS]

	# Index to tree summary file
	indTS=np.where( (dTS['ID_Inst']==dTL['ID_Inst'][indUT[0]]) & (dTS['ID_Plot']==dTL['ID_Plot'][indUT[0]]) & (dTS['ID_Tree']==dTL['ID_Tree'][indUT[0]]) )[0]

	if indTS.size>0:
		try:
			dTL['Dob 2020'][indUT]=dTS['DBH_2020 (cm)'][indTS].astype('float')
		except:
			print('Something odd about DBH 2020')
		dTL['H 2020'][indUT]=dTS['Ht_2020 (m)'][indTS][0]
	else:
		print('No crosswalk found with tree summary')
		continue

	# Inside-bark diameter (cm)
	dTL['Dib'][indUT]=0.1*np.cumsum(2*np.nan_to_num(dTL['TRW'][indUT]))

	# Inside-bark diameter for last year (cm)
	# *** Only populated if last year >= 2019 (for comparison with 2020 tape measurement)
	if dTL['Year'][indUT[iWithData[-1]]]>=2019:
		dTL['Dib Last'][indUT[iWithData[-1]]]=dTL['Dib'][indUT[iWithData[-1]]]

	# Stemwood biomass from allometric function (kgDM) - DBH only (Ung et al. 2008)
	dTL['Bsw'][indUT]=0.0204*dTL['Dib'][indUT]**2.6974

	# Stemwood biomass growth (kgDM/yr)
	dTL['Gsw'][indUT[1:]]=np.diff(dTL['Bsw'][indUT])

	# Relative growth rate
	dTL['RGR'][indUT[1:]]=np.log(dTL['Bsw'][indUT[1:]])-np.log(dTL['Bsw'][indUT[0:-1]])

	# Standardized growth relative to mean growth during a period leading up to N application
	ind0=np.where( (dTL['ID_Tree_Unique']==uTID[iU]) & (dTL['Year']>=1971-5) & (dTL['Year']<=1970) )[0]
	dTL['TRW RR'][indUT]=dTL['TRW'][indUT]/np.mean(dTL['TRW'][ind0])
	dTL['BAI RR'][indUT]=dTL['BAI'][indUT]/np.mean(dTL['BAI'][ind0])
	dTL['Bsw RR'][indUT]=dTL['Bsw'][indUT]/np.mean(dTL['Bsw'][ind0])
	dTL['Gsw AD'][indUT]=dTL['Gsw'][indUT]-np.mean(dTL['Gsw'][ind0])
	dTL['Gsw RR'][indUT]=dTL['Gsw'][indUT]/np.maximum(0.1,np.mean(dTL['Gsw'][ind0]))
	dTL['RGR RR'][indUT]=dTL['RGR'][indUT]/np.mean(dTL['RGR'][ind0])

# Change name of treatment
dTL['Treatment']=dTL['Treatment.EP703']

# Delete fields
del dTL['Unnamed: 0'],dTL['Treatment.EP703']

# Summarize QA flags
dTL['QA Summary']=np.zeros(dTL['TRW'].size)

# Cross dating
ind=np.where( (dTL['QC_xdate_flags']=='Good_correlation') | (dTL['QC_xdate_flags']=='Weak_correlation') )[0]
dTL['QA Summary'][ind]=1

# Unrealistic size removed
ind=np.where( (dTL['Dib Last']/dTL['Dob 2020']>1.0) )[0]
dTL['QA Summary'][ind]=0

print('% retained = ' + str(100*np.where(dTL['QA Summary']==1)[0].size/dTL['QA Summary'].size))

# Define unique installations
uInst=np.unique(dTL['ID_Inst'])

# Define unique site series
uSS=np.unique(dTL['BGC SS'])

# Combine SS 5 and 6 (may not be a good idea for all analysis!)
dTL['BGC SS Comb']=dTL['BGC SS'].copy()
ind=np.where( (dTL['BGC SS']==5) | (dTL['BGC SS']==6) )[0]
dTL['BGC SS Comb'][ind]=99

#%% Add detrended ring width (from J Axelson in dplR)

# Path to data
pthin=r'C:\Users\rhember\OneDrive - Government of BC\Shared\EP703\Retrospective Monitoring Study\Data\Given\Core Data\Received 20211222 from JAxelson'

# Initialize
dTL['TRW S NE DPLR']=np.zeros(dTL['Year'].size)

# Loop through installations
for iI in range(uInst.size):

	if uInst[iI]<10:
		nI='0' + str(uInst[iI])
	else:
		nI=str(uInst[iI])

	# Determine the letter associated with this installation
	letL=['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p']
	for let in letL:
		try:
			d0=gu.ReadExcel(pthin + '\\703_' + nI + '\\Series\\703_' + nI + let + '_rwi_NE.xlsx')
			break
		except:
			pass

	# Pull year out
	Year=d0['Unnamed: 0']; del d0['Unnamed: 0']

	# Loop through trees
	for k in d0.keys():

		# Get tree and plot IDs
		ind=0
		for i in range(len(k)):
			if k[i]=='_':
				ind=i
				break
		ID_Plot=int(k[0:ind])
		ID_Tree=int(k[ind+1:])

		# Index to tree-level structure
		ind0=np.where( (dTL['ID_Inst']==uInst[iI]) & (dTL['ID_Plot']==ID_Plot) & (dTL['ID_Tree']==ID_Tree) & (dTL['Year']>=Year[0]) & (dTL['Year']<=Year[-1]) )[0]
		if ind0.size==0:
			print('Tree not found in DB')
			continue
		#ind0=np.where( (dTL['ID_Inst']==uInst[iI]) & (dTL['ID_Plot']==ID_Plot) )[0]
		#np.unique(dTL['ID_Tree'][ind0])
		#ind0.size

		# Populate
		dTL['TRW S NE DPLR'][ind0]=d0[k]

# Standardize the detrended time series
dTL['TRW S NE DPLR RR']=np.zeros(dTL['TRW S NE DPLR'].size)
for iU in range(uTID.size):

	# Index to unique tree ID
	indUT=np.where(dTL['ID_Tree_Unique']==uTID[iU])[0]

	# Standardized growth relative to mean growth during a period leading up to N application
	ind0=np.where( (dTL['ID_Tree_Unique']==uTID[iU]) & (dTL['Year']>=1971-5) & (dTL['Year']<=1970) )[0]
	dTL['TRW S NE DPLR RR'][indUT]=dTL['TRW S NE DPLR'][indUT]/np.mean(dTL['TRW S NE DPLR'][ind0])

#%% Add climate data

# Define unique coordinates
uLL=np.unique(np.column_stack((dTL['ID_Inst'],dTL['ID_Plot'],dTL['Lat'],dTL['Lon'])),axis=0)

# Export to process climate data
flg=0
if flg==1:
	df=pd.DataFrame(uLL)
	df.to_excel(r'G:\My Drive\EP703 TR coordinates.xlsx')

# Define base period for normals
iTn_ann=np.where( (clm['tv'][:,0]>=1971) & (clm['tv'][:,0]<=2000) )[0]
iTn_ws=np.where( (clm['tv'][:,0]>=1971) & (clm['tv'][:,0]<=2000) & (clm['tv'][:,1]>=5) & (clm['tv'][:,0]<=9) )[0]

yrs=np.arange(1950,2022,1)

dTL['tmean_ann_n']=np.zeros(dTL['TRW'].size)
dTL['prcp_ann_n']=np.zeros(dTL['TRW'].size)
dTL['tmean_gs_r']=np.zeros(dTL['TRW'].size)
dTL['prcp_gs_r']=np.zeros(dTL['TRW'].size)
dTL['cwd_gs_r']=np.zeros(dTL['TRW'].size)
dTL['ws_gs_r']=np.zeros(dTL['TRW'].size)
for iLL in range(uLL.shape[0]):
	ind=np.where( (dTL['Lat']==uLL[iLL,2]) & (dTL['Lon']==uLL[iLL,3]) )[0]
	dTL['tmean_ann_n'][ind]=np.mean(clm['tmean'][iTn_ann,iLL])
	dTL['prcp_ann_n'][ind]=np.sum(clm['prcp'][iTn_ann,iLL])/(iTn_ann.size/12)
	for iY in range(yrs.size):
		iT=np.where( (clm['tv'][:,0]==yrs[iY]) & (clm['tv'][:,1]>=4) & (clm['tv'][:,1]<=9) )[0]
		iT2=np.where(dTL['Year'][ind]==yrs[iY])[0]
		dTL['tmean_gs_r'][ind[iT2]]=np.mean(clm['tmean'][iT,iLL])
		dTL['prcp_gs_r'][ind[iT2]]=np.mean(clm['prcp'][iT,iLL])
		dTL['cwd_gs_r'][ind[iT2]]=np.mean(clm['etp_tmw'][iT,iLL]-clm['eta_tmw'][iT,iLL])
		dTL['ws_gs_r'][ind[iT2]]=np.mean(clm['ws_tmw'][iT,iLL])

uC=np.unique(np.column_stack((dTL['tmean_ann_n'],dTL['prcp_ann_n'])),axis=0)

#plt.plot(uC[:,0],uC[:,1],'.')

#%% Get stand-level information from EP703 dataset
sobs=gu.ipickle(r'C:\Data\EP703\Processed\EP703_SL.pkl')
sobs['RGR']=(np.log(sobs['Bsw_t1'])-np.log(sobs['Bsw_t0']))/sobs['DT']

dByInst={}
dByInst['DA_BA_G']=np.zeros(uInst.size)
dByInst['DR_BA_G']=np.zeros(uInst.size)
dByInst['DA_Hobs_G']=np.zeros(uInst.size)
dByInst['DR_Hobs_G']=np.zeros(uInst.size)
dByInst['DA_Bsw_G']=np.zeros(uInst.size)
dByInst['DR_Bsw_G']=np.zeros(uInst.size)
dByInst['DA_RGR']=np.zeros(uInst.size)
dByInst['DR_RGR']=np.zeros(uInst.size)
for iInst in range(uInst.size):
	iC=np.where( (sobs['ID_Site']==uInst[iInst]) & (sobs['N_Dose']==0) & (sobs['TSF_t0']>=0) & (sobs['TSF_t1']<=9) )[0]
	iF=np.where( (sobs['ID_Site']==uInst[iInst]) & (sobs['N_Dose']==225) & (sobs['TSF_t0']>=0) & (sobs['TSF_t1']<=9) )[0]
	dByInst['DA_BA_G'][iInst]=np.mean(sobs['BA_G'][iF]-sobs['BA_G'][iC])
	dByInst['DR_BA_G'][iInst]=np.mean( (sobs['BA_G'][iF]-sobs['BA_G'][iC])/sobs['BA_G'][iC]*100 )
	dByInst['DA_Hobs_G'][iInst]=np.nanmean(sobs['H_obs_G'][iF]-sobs['H_obs_G'][iC])
	dByInst['DR_Hobs_G'][iInst]=np.nanmean( (sobs['H_obs_G'][iF]-sobs['H_obs_G'][iC])/sobs['H_obs_G'][iC]*100 )
	dByInst['DA_Bsw_G'][iInst]=np.mean(sobs['Bsw_G'][iF]-sobs['Bsw_G'][iC])
	dByInst['DR_Bsw_G'][iInst]=np.mean( (sobs['Bsw_G'][iF]-sobs['Bsw_G'][iC])/sobs['Bsw_G'][iC]*100 )
	dByInst['DA_RGR'][iInst]=np.median(sobs['RGR'][iF]-sobs['RGR'][iC])
	dByInst['DR_RGR'][iInst]=np.median( (sobs['RGR'][iF]-sobs['RGR'][iC])/sobs['RGR'][iC]*100 )

	#plt.close('all')
	#plt.plot(sobs['TSF_t0'][iC],sobs['N_t0'][iC],'bo')
	#plt.plot(sobs['TSF_t0'][iF],sobs['N_t0'][iF],'rs')

	#plt.plot(sobs['TSF_t0'][iC],sobs['Bsw_G'][iC],'bo')
	#plt.plot(sobs['TSF_t0'][iF],sobs['Bsw_G'][iF],'rs')

#%% Get tree-level EP703 data
dTL['H obs EP']=np.zeros(dTL['TRW'].size)
dTL['Dob obs EP']=np.zeros(dTL['TRW'].size)
for iInst in range(uInst.size):
	ind0=np.where( (dTL['ID_Inst']==uInst[iInst]) & (dTL['Treatment']=='C') | (dTL['ID_Inst']==uInst[iInst]) & (dTL['Treatment']=='F1') )[0]
	uPlot=np.unique(dTL['ID_Plot'][ind0])
	for iPlot in range(uPlot.size):
		ind1=np.where( (dTL['ID_Inst']==uInst[iInst]) & (dTL['ID_Plot']==uPlot[iPlot]) & (dTL['QA Summary']==1) )[0]
		uTree=np.unique(dTL['ID_Tree'][ind1])
		for iTree in range(uTree.size):
			ind2=np.where( (dEPT['ID_Site']==uInst[iInst]) & (dEPT['ID_Plot']==uPlot[iPlot]) & (dEPT['ID_Tree']==uTree[iTree]) & (dEPT['H_obs_t0']>0) )[0]
			for iY in range(ind2.size):
				ind3=np.where( (dTL['ID_Inst']==uInst[iInst]) & (dTL['ID_Plot']==uPlot[iPlot]) & (dTL['ID_Tree']==uTree[iTree]) & (dTL['Year']==dEPT['Year_t0'][ind2[iY]]) )[0]
				dTL['H obs EP'][ind3]=dEPT['H_obs_t0'][ind2[iY]]
				dTL['Dob obs EP'][ind3]=dEPT['D_t0'][ind2[iY]]

#%% Comparison between Dob and Dib from cores
def CompareDiameters(dTL,meta):

	# Isolate data
	ikp=np.where( (dTL['Dob 2020']>0) & (dTL['Dib Last']>0) & (dTL['QA Summary']==1) )[0]
	#ikp=np.where( (dTL['Dob 2020']>0) & (dTL['Dib Last']>0) & (dTL['Pith_status']=='Yes') )[0]

	# Fit linear best fit relationship
	y=dTL['Dib Last'][ikp]
	x=dTL['Dob 2020'][ikp]
	x1=sm.tools.tools.add_constant(x)
	md=sm.OLS(y,x1).fit()
	md.summary()
	xhat=np.linspace(np.min(x1[:,1]),np.max(x1[:,1]),10)
	yhat=md.predict(np.c_[np.ones(xhat.size),xhat])

	# Plot
	plt.close('all')
	fig,ax=plt.subplots(1,figsize=gu.cm2inch(8.5,8.5))
	ax.plot([0,90],[0,90],'k-',lw=2,color=[0.8,0.8,0.8])
	ax.plot(x,y,'o',ms=4,mec='w',mfc='k')
	ax.plot(xhat,yhat,'r-',lw=1.25)
	ax.set(xlim=[0,80],ylim=[0,80],xlabel='Tape-based Dob in 2020 (cm)',ylabel='Core-based Dib in 2019 (cm)')
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')
	plt.tight_layout()
	gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\EP703\TreeRings\DiameterComparison','png',150)
	return

#%% Save
df=pd.DataFrame(dTL)
df.to_excel(r'C:\Users\rhember\OneDrive - Government of BC\Shared\EP703\Retrospective Monitoring Study\Data\EP703_RM_TreeDataAnalysisReady.xlsx',index=False)

#%% Look at individual trees
def LookAtIndividualTrees():

#	ind=np.where( (dTL['ID_Tree_Unique']==25) )[0]
#
#	print(dTL['Est/mm'][ind[0]])
#
#	plt.close('all')
#	fig,ax=plt.subplots(1,figsize=gu.cm2inch(8.5,8.5))
#	#ax.plot([0,90],[0,90],'k-',lw=2,color=[0.8,0.8,0.8])
#	ax.plot(dTL['Age'][ind],dTL['TRW'][ind],'-ko',ms=3,mec='k',mfc='k')
#	#ax.plot(dTL['Year'][ind],dTL['TRW'][ind],'-ko',ms=3,mec='k',mfc='k')
#
#
#	fig,ax=plt.subplots(1,figsize=gu.cm2inch(8.5,8.5))
#	#ax.plot([0,90],[0,90],'k-',lw=2,color=[0.8,0.8,0.8])
#	ax.plot(dTL['Age'][ind],dTL['Dib'][ind],'-ko',ms=3,mec='k',mfc='k')
#	#ax.plot(dTL['Year'][ind],dTL['TRW'][ind],'-ko',ms=3,mec='k',mfc='k')

	u=np.unique(dGF['ID_Tree_Unique'])
	for iTree in range(u.size):
		ind=np.where( (dGF['ID_Tree_Unique']==u[iTree]) )[0]
		if ind.size==0:
			continue
		plt.close('all')
		fig,ax=plt.subplots(1,figsize=gu.cm2inch(8.5,8.5))
		ax.plot(dGF['Year'][ind],dGF['Dib'][ind],'-b')
		gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\EP703\TreeRings\IndividualTree\Tree_' + str(int(u[iTree])),'png',300)

	return

#%% Gap-fill early missing rings using Weibul function
# *** Exploritory, not used in RM study. ***

# Weibul distribution function
def func(x,a,b,c,d):
	y=a*(1+((b*(x/c)**d-1)/(np.exp(x/c))))
	return y

def GapFill():

	uTID=np.unique(dTL['ID_Tree_Unique'])

	# Initialize new structure
	dGF={}
	ExclusionList=['TRW S NE DPLR', 'TRW S NE DPLR RR', 'tmean_ann_n', 'prcp_ann_n', 'tmean_gs_r', 'prcp_gs_r', 'cwd_gs_r', 'ws_gs_r']
	for k in dTL.keys():
		if np.isin(k,ExclusionList)==True:
			continue
		if (dTL[k].dtype=='float') | (dTL[k].dtype=='int32'):
			dGF[k]=np.zeros(200000)
		else:
			dGF[k]=np.array(['' for _ in range(200000)],dtype=dTL[k].dtype)

	cnt=0
	for iU in range(uTID.size):

		# Index to unique tree ID
		iTree=np.where( (dTL['ID_Tree_Unique']==uTID[iU]) & (dTL['QA Summary']==1) )[0]
		if iTree.size==0:
			continue

		TRW=dTL['TRW'][iTree]
		Age=dTL['Age'][iTree]
		Year=dTL['Year'][iTree]
		D_ob20=dTL['Dob 2020'][iTree]
		iBad=np.where( (np.isnan(TRW)==True) )[0]

		flg=0
		if flg==1:
			plt.close('all')
			fig,ax=plt.subplots(1)
			ax.plot(Age,TRW,'.b-')
			p0=np.array([1.2,4,13,2.2])
			TRW_hat=func(x,p0[0],p0[1],p0[2],p0[3])
			ax.plot(x,TRW_hat,'g--',lw=1.5)

		iFit=np.where( (np.isnan(TRW)==False) )[0]
		x=np.append(np.arange(0,5,1),Age[iFit])
		y=np.append(np.zeros(5),TRW[iFit])
		p0=np.array([1.2,4,13,2.2])
		#p,pcov=curve_fit(func,x,y,p0)
		try:
			p,pcov=curve_fit(func,x,y,p0)
		except:
			print('Fit failed, skipping')
			continue
		TRW_hat=func(Age,p[0],p[1],p[2],p[3])
		#ax.plot(Age,TRW_hat,'g--',lw=1.5)

		TRW_gf=TRW.copy()
		TRW_gf[iBad]=TRW_hat[iBad]
		#fig,ax=plt.subplots(1)
		#ax.plot(Age,TRW_gf,'b-')

		# Diameter from core
		Dg_core_gf=0.1*2*TRW_gf
		D_core_gf=np.cumsum(Dg_core_gf)
		#fig,ax=plt.subplots(1)
		#ax.plot(Age,D_core_gf,'b-')

		# Force to match measurement of Dob in 2020
		rat=D_ob20[-1]/D_core_gf[-1]
		D_core_gf=rat*D_core_gf
		Dg_core_gf=np.append(0,np.diff(D_core_gf))

		#fig,ax=plt.subplots(1)
		#ax.plot(Age,D_core_gf,'b-')
		#ax.plot(Age[-1],D_ob20[-1],'ks')

		# Index to final data dictionary
		ind=np.arange(cnt,cnt+Year.size,1)

		# Populate scalar variables
		for k in dGF.keys():
			dGF[k][ind]=dTL[k][ iTree[0] ]

		#Filler=np.ones(Year_full.size)

		dGF['Year'][ind]=Year
		dGF['Age'][ind]=Age
		dGF['TRW'][ind]=TRW_gf
		dGF['Dib'][ind]=D_core_gf
		dGF['Dib Last'][ind]=np.max(D_core_gf)

		BA=np.pi*(dGF['Dib'][ind]/2)**2
		dGF['BAI'][ind]=np.append(0,np.diff(BA))

		# Stemwood biomass from allometric function (kgDM) - DBH only (Ung et al. 2008)
		dGF['Bsw'][ind]=0.0204*dGF['Dib'][ind]**2.6974

		# Stemwood biomass growth (kgDM/yr)
		dGF['Gsw'][ind]=np.append(0,np.diff(dGF['Bsw'][ind]))

		# Relative growth rate
		dGF['RGR'][ind]=np.append(0,np.log(dGF['Bsw'][ind[1:]])-np.log(dGF['Bsw'][ind[0:-1]]))

		# Standardized growth relative to mean growth during a period leading up to N application
		ind0=np.where( (dGF['ID_Tree_Unique']==dGF['ID_Tree_Unique'][ind[0]]) & (dGF['Year']>=1971-5) & (dGF['Year']<=1970) )[0]
		dGF['TRW RR'][ind]=dGF['TRW'][ind]/np.mean(dGF['TRW'][ind0])
		dGF['BAI RR'][ind]=dGF['BAI'][ind]/np.mean(dGF['BAI'][ind0])
		dGF['Bsw RR'][ind]=dGF['Bsw'][ind]/np.mean(dGF['Bsw'][ind0])
		dGF['Gsw AD'][ind]=dGF['Gsw'][ind]-np.mean(dGF['Gsw'][ind0])
		dGF['Gsw RR'][ind]=dGF['Gsw'][ind]/np.maximum(0.1,np.mean(dGF['Gsw'][ind0]))
		dGF['RGR RR'][ind]=dGF['RGR'][ind]/np.mean(dGF['RGR'][ind0])

		# Update counter
		cnt=cnt+Year.size

	# Remove excess data
	for k in dGF.keys():
		dGF[k]=dGF[k][0:cnt]

	return dGF

dGF=GapFill()

#%% Comparison between Dob and Dib from cores (after gap-filling)
def CompareDiametersAfterGF():
	# Isolate data
	ikp=np.where( (dGF['Dob 2020']>0) & (dGF['Dib Last']>0) & (dGF['QA Summary']==1) )[0]
	#ikp=np.where( (dGF['Dob 2020']>0) & (dGF['Dib Last']>0) & (dGF['Pith_status']=='Yes') )[0]

	# Fit linear best fit relationship
	y=dGF['Dib Last'][ikp]
	x=dGF['Dob 2020'][ikp]
	x1=sm.tools.tools.add_constant(x)
	md=sm.OLS(y,x1).fit()
	md.summary()
	xhat=np.linspace(np.min(x1[:,1]),np.max(x1[:,1]),10)
	yhat=md.predict(np.c_[np.ones(xhat.size),xhat])

	# Plot
	#plt.close('all')
	fig,ax=plt.subplots(1,figsize=gu.cm2inch(8.5,8.5))
	ax.plot([0,90],[0,90],'k-',lw=2,color=[0.8,0.8,0.8])
	ax.plot(x,y,'o',ms=4,mec='w',mfc='k')
	ax.plot(xhat,yhat,'r-',lw=1.25)
	ax.set(xlim=[0,80],ylim=[0,80],xlabel='Tape-based Dob in 2020 (cm)',ylabel='Core-based Dib in 2019 (cm)')
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')
	plt.tight_layout()
	#gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\EP703\TreeRings\DiameterComparison_GF','png',150)
	return

#%% QA - Plot of each tree
plt.close('all')
fig,ax=plt.subplots(1,figsize=gu.cm2inch(8.5,8.5))
uT=np.unique(dTL['ID_Tree_Unique'])
for iT in range(uT.size):
	ind=np.where(dGF['ID_Tree']==uT[iT])[0]
	cl=np.random.random(3)
	ax.plot(dGF['Year'][ind],dGF['Age'][ind],'k-',color=cl)


