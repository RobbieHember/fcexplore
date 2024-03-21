''' 
EP703 UTILITIES
'''
#%% Import modules
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyproj
from scipy.optimize import curve_fit
import scipy.io as spio
import warnings
from fcgadgets.macgyver import util_general as gu
warnings.filterwarnings('ignore')
gp=gu.SetGraphics('Manuscript')

#%% Configuration
def Config():
	meta={}
	meta['Paths']={}
	meta['Paths']['Project']=r'C:\Data\EP703\Data'
	meta['Paths']['Raw Data']=meta['Paths']['Project'] + '\\Raw Data\\Received20240216'
	meta['Paths']['Outputs']=meta['Paths']['Project'] + '\\Processed'
	meta['Paths']['NACID']=r'E:\Data\Climate\NACID'

	meta['Param']={}
	meta['Param']['Allometry']={} # Allometry - DBH plus height (Ung et al. 2008)
	meta['Param']['Allometry']['FD']={'Stemwood':[0.0191,1.5365,1.3634],'Foliage':[0.0718,2.2935,-0.4744],'Branch':[0.0351,2.24811,0],'Bark':[0.0083,2.4811,0]}
	meta['Param']['Allometry']['HW']={'Stemwood':[0.0113,1.9332,1.1125],'Foliage':[0.2656,2.0107,-0.7963],'Branch':[0.0609,2.0021,0],'Bark':[0.0019,2.3356,0.6371]}
	meta['Param']['Allometry']['CW']={'Stemwood':[0.0188,1.3376,1.5293],'Foliage':[0.1097,1.553,0.0000],'Branch':[0.0611,1.9208,0],'Bark':[0.0002,2.4369,1.1315]}
	meta['Param']['Allometry']['BA']={'Stemwood':[0.0315,1.8297,0.8056],'Foliage':[0.0453,2.4867,-0.4982],'Branch':[0.0420,2.0313,0],'Bark':[0.0067,2.6970,-0.3105]}
	# Allometry - DBH only (Ung et al. 2008)
	# d2b={'Species':'FD','Stemwood':[0.0204,2.6974],'Foliage':[0.1233,1.6636],'Branch':[0.0404,2.1388],'Bark':[0.0069,2.5462]}
	# pA.append(d2b.copy())
	# d2b={'Species':'HW','Stemwood':[0.0141,2.8668],'Foliage':[0.1676,1.4339],'Branch':[0.0703,1.9547],'Bark':[0.0025,2.8062]}
	# pA.append(d2b.copy())
	# d2b={'Species':'CW','Stemwood':[0.0111,2.8027],'Foliage':[0.1233,1.5152],'Branch':[0.1158,1.7196],'Bark':[0.0003,3.2721]}
	# pA.append(d2b.copy())
	# d2b={'Species':'BA','Stemwood':[0.0424,2.4289],'Foliage':[0.0645,1.9400],'Branch':[0.0322,2.1313],'Bark':[0.0057,2.4786]}
	# pA.append(d2b.copy())
	
	meta['Param']['kg2Mg']=0.001 # Convert kg to Mg
	meta['Param']['dm2c']=0.5 # Convert dry mass to carbon
	meta['Param']['da_list']=['A', 'AB', 'AD', 'AE', 'AX', 'BT', 'CK', 'D', 'DDF', 'DDP', 'DDS','DM', 'DN', 'DR', 'DRA', 'DS', 'DT', 'FCK', 'FK', 'MT', 'N', 'NGC','NW', 'NWS', 'NX', 'NY', 'PC', 'SC', 'SCB', 'SCM', 'SCT'] # List of DBHage agents
	return meta

#%% Function (tree height as a function of basal area)
def funcH(x,a,b):
	return a*x**b

#%% Import data
def ImportData(meta):

	# Import site data
	dSite=gu.ReadExcel(meta['Paths']['Project'] + '\\EP703 Plot Level Summary.xlsx')

	# Unique sites
	uSite=np.unique(dSite['ID Installation'])

	# Gap-filled heights from FAIB
	#dGFH=gu.ReadExcel(meta['Paths']['Project'] + '\\Raw Data\\GapFilledHeightsFromFAIB\\EP703_tdf_ht_vol.xlsx')

	bpL=['Stemwood','Foliage','Branch','Bark']
	sL=['FD','HW','CW','BA']

	# Initialize stand-level dictionary
	varS=['ID_Site','ID_Plot','Lat','Lon','BGC','SI','EstabType','AreaPlot',
		'N_Init','pH_Min_Init','pH_Humus_Init','P_Min_Init','P_Humus_Init','Ca_Min_Init',
		'NA_Dose','NA_Type','NA_Method',
		'Thin_Level','Thin_Method','Thin_Format','Thin_TargetSPH',
		'NumTreesMeas','Num_H_G_NegOrMis',
		'Spc1_CD_t0','Spc1_Pct_t0','Spc2_CD_t0','Spc2_Pct_t0',
		'Year_t0','Year_t1',
		'DT',
		'Age_t0','Age_t1',
		'TSF_t0','TSF_t1',
		'N_t0','N_t1','N_Net',
		'DBH_t0','DBH_t1',
		'BA_t0','BA_t1','BA_G',
		'H_obs_t0','H_obs_t1','H_obs_G',
		'H_gf_t0','H_gf_t1','H_gf_G',
		'Cbk_t0','Cbk_t1','Cbk_G','Cbk_M','Cbk_R','Cbk_L','Cbk_Net',
		'Cbr_t0','Cbr_t1','Cbr_G','Cbr_M','Cbr_R','Cbr_L','Cbr_Net',
		'Cf_t0','Cf_t1','Cf_G','Cf_M','Cf_R','Cf_L','Cf_Net',
		'Csw_t0','Csw_t1',
		'Csw_G','Csw_M','Csw_L','Csw_R','Csw_Net',
		'Csw_G_cc1','Csw_G_cc2','Csw_G_cc3','Csw_G_cc4','Csw_G_cc5','Csw_G_cc6',
		'Csw_M_cc1','Csw_M_cc2','Csw_M_cc3','Csw_M_cc4','Csw_M_cc5','Csw_M_cc6',
		'Csw_Net_cc1','Csw_Net_cc2','Csw_Net_cc3','Csw_Net_cc4','Csw_Net_cc5','Csw_Net_cc6',
		'DA_A_t0','DA_AB_t0','DA_AD_t0','DA_AE_t0','DA_AX_t0','DA_BT_t0',
		'DA_CK_t0','DA_D_t0','DA_DDF_t0','DA_DDP_t0','DA_DDS_t0',
		'DA_DM_t0','DA_DN_t0','DA_DR_t0','DA_DRA_t0','DA_DS_t0','DA_DT_t0',
		'DA_FCK_t0','DA_FK_t0','DA_MT_t0','DA_N_t0','DA_NGC_t0','DA_NW_t0','DA_NWS_t0',
		'DA_NX_t0','DA_NY_t0','DA_PC_t0','DA_SC_t0','DA_SCB_t0','DA_SCM_t0','DA_SCT_t0',
		'DA_A_t1','DA_AB_t1','DA_AD_t1','DA_AE_t1','DA_AX_t1','DA_BT_t1',
		'DA_CK_t1','DA_D_t1','DA_DDF_t1','DA_DDP_t1','DA_DDS_t1',
		'DA_DM_t1','DA_DN_t1','DA_DR_t1','DA_DRA_t1','DA_DS_t1','DA_DT_t1',
		'DA_FCK_t1','DA_FK_t1','DA_MT_t1','DA_N_t1','DA_NGC_t1','DA_NW_t1','DA_NWS_t1',
		'DA_NX_t1','DA_NY_t1','DA_PC_t1','DA_SC_t1','DA_SCB_t1','DA_SCM_t1','DA_SCT_t1']

	varL_str=['BGC','EstabType','NA_Type','NA_Method','Spc1_CD_t0','Spc2_CD_t0','Thin_Method','Thin_Format','Thin_TargetSPH']
	sobs={}
	for v in varS:
		if np.isin(v,varL_str)==True:
			sobs[v]=np.array(['' for _ in range(25000)],dtype=object)
		else:
			sobs[v]=np.zeros(25000,dtype='float')

	# Initialize a counter for stand-level dictionary
	cntS=0
	#sobs['BGC']=sobs['BGC'].astype('U')
	#sobs['EstabType']=sobs['EstabType'].astype('U')
	#sobs['NA_Type']=sobs['NA_Type'].astype('U')
	#sobs['NA_Method']=sobs['NA_Method'].astype('U')
	#sobs['Spc1_CD_t0']=sobs['Spc1_CD_t0'].astype('U')
	#sobs['Spc2_CD_t0']=sobs['Spc2_CD_t0'].astype('U')
	
	# Initialize tree-level dictionary
	varT=['ID_Site','ID_Plot','ID_Tree','Lat','Lon','BGC','SI','AEF','N_Init','pH_Min_Init',
		'NA_Dose','Year_t0','Year_t1','DT',
		'TSF_t0','TSF_t1',
		'SA_t0','SA_t1',
		'SN_t0','SN_t1',
		'SCsw_t0','SCsw_t1','SCswLT_t0','SCswLT_t1',
		'ID_Spc_t0','Mort','Rec','TreeClass_t0','CrownClass_t0','TopDead_t0','TopDead_t1',
		'TopBroke_t0','TopBroke_t1','DA_t0','DA_t1','D_t0','D_t1','BA_t0','BA_t1','BA_G',
		'H_obs_t0','H_obs_t1','H_obs_G','H_mod_t0','H_mod_t1','H_mod_G','H_gf_t0','H_gf_t1','H_gf_G',
		'Csw_t0','Csw_t1','Csw_G','Csw_mod_t0','Csw_mod_t1','H_Fit_ITM_t0','H_Fit_ITM_t1']
	tobs={}
	for k in varT:
		tobs[k]=np.zeros(1500000)
	
	# Initialize a counter for tree-level dictionary
	cntT=0
	
	tobs['BGC']=tobs['BGC'].astype('U')
	tobs['DA_t0']=tobs['DA_t0'].astype('U')
	tobs['DA_t1']=tobs['DA_t1'].astype('U')

	# Loop through installations (40 sec)
	for iS in range(uSite.size):

		# Import tree-level data for site
		if uSite[iS]<10:
			fnam='0' + str(uSite[iS])
		else: 
			fnam=str(uSite[iS])
		pthin=meta['Paths']['Raw Data'] + '\\EP703_inst' + fnam + '.csv'
		try:
			dInst=gu.ReadCSV(pthin)
		except:
			dInst=gu.ReadExcel(meta['Paths']['Raw Data'] + '\\EP703_inst' + fnam + '.xlsx')

		#N=N+dInst['tree']'].unique().size

		# Calculate basal area (cm2)
		dInst['ba']=np.pi*(dInst['dbh']/2)**2

		# Vital status
		dInst['VS']=np.ones(dInst['tree'].size)
		dInst['VS'][dInst['mort']!='nan']=0

		# Initialize modelled and gap-filled tree height
		dInst['ht_mod']=np.nan*np.ones(dInst['tree'].size)
		dInst['ht_gf']=dInst['ht'].copy()

		# Indicator
		dInst['H_Fit_ITM']=np.zeros(dInst['tree'].size)

		# Index to site table
		indSite=np.where(dSite['ID Installation']==uSite[iS])[0]

		# Unique plots within this installation
		uPlot=dSite['ID Plot'][indSite]

		# Unique years
		#uYear=np.sort(np.unique(dInst['year']))

		# Area expansion factor
		Area=dSite['Plot Area Ha'][indSite[0]]
		AEF=1/Area

		# Age 
		SA=dSite['Initial Age'][indSite[0]]

		# Year of establishment
		YrEst=dSite['Year Established'][indSite[0]]

		# Geographic coordinates
		Lat=dSite['Lat Deg'][indSite[0]]+dSite['Lat Min'][indSite[0]]/60+dSite['Lat Sec'][indSite[0]]/(60*60)
		Lon=-dSite['Lon Deg'][indSite[0]]-dSite['Lon Min'][indSite[0]]/60-dSite['Lon Sec'][indSite[0]]/(60*60)

		# Gap fill heights
		dInst=GapFillHeights(meta,dInst)

		# Import FAIB cleaned tree data
		#dInst=ImportCleanedTreeDataFromFAIB(dInst):

		# Loop through plots and years
		for iP in range(len(uPlot)):

			# Index to plot
			indP=np.where( (dSite['ID Installation']==uSite[iS]) & (dSite['ID Plot']==uPlot[iP]) )[0][0]

			# Unique measurement years for this plot
			indP2=np.where( (dInst['inst']==uSite[iS]) & (dInst['plot']==uPlot[iP]) )[0]
			uYear=np.sort(np.unique(dInst['year'][indP2]))

			for iY in range(len(uYear)-1):
				print(str(uSite[iS]) + ' ' + str(uPlot[iP]) + ' ' + str(uYear[iY]))

				yr0=uYear[iY]
				yr1=uYear[iY+1]
				dt=yr1-yr0

				# Index to first measurement
				ind0=np.where((dInst['plot']==uPlot[iP]) & (dInst['year']==uYear[iY]))[0]

				# Index to second measurement
				ind1=np.where((dInst['plot']==uPlot[iP]) & (dInst['year']==uYear[iY+1]))[0]

				# First measurement
				id0_all=dInst['tree'][ind0]
				sp0_all=dInst['spp'][ind0]
				vs0_all=dInst['VS'][ind0]
				dbh0_all=dInst['dbh'][ind0]
				h_obs0_all=dInst['ht'][ind0]
				h_mod0_all=dInst['ht_mod'][ind0]
				h_gf0_all=dInst['ht_gf'][ind0]
				ba0_all=np.pi*(dbh0_all/2)**2
				TreeClass0_all=dInst['TrClass'][ind0]
				CrownClass0_all=dInst['CrClass'][ind0]
				btop0_all=dInst['bTop'][ind0]
				dtop0_all=dInst['dTop'][ind0]
				da0_all=dInst['DmAg1'][ind0]
				H_Fit_ITM0_all=dInst['H_Fit_ITM'][ind0]
				b0_all={}
				for v in bpL:
					b=meta['Param']['Allometry']['FD'][v]
					b0_all[v]=b[0]*dbh0_all**b[1]*h_gf0_all**b[2]
					for s in sL:
						b=meta['Param']['Allometry'][s][v]
						b0_all[v]=b[0]*dbh0_all**b[1]*h_gf0_all**b[2]

				# Second measurement
				id1_all=dInst['tree'][ind1]
				sp1_all=dInst['spp'][ind1]
				vs1_all=dInst['VS'][ind1]
				dbh1_all=dInst['dbh'][ind1]
				h_obs1_all=dInst['ht'][ind1]
				h_mod1_all=dInst['ht_mod'][ind1]
				h_gf1_all=dInst['ht_gf'][ind1]
				H_Fit_ITM1_all=dInst['H_Fit_ITM'][ind1]
				ba1_all=np.pi*(dbh1_all/2)**2
				TreeClass1_all=dInst['TrClass'][ind1]
				CrownClass1_all=dInst['CrClass'][ind1]
				btop1_all=dInst['bTop'][ind1]
				dtop1_all=dInst['dTop'][ind1]
				da1_all=dInst['DmAg1'][ind1]
				b1_all={}
				for v in bpL:
					b=meta['Param']['Allometry']['FD'][v]
					b1_all[v]=b[0]*dbh1_all**b[1]*h_gf1_all**b[2]
					for s in sL:
						b=meta['Param']['Allometry'][s][v]
						b1_all[v]=b[0]*dbh1_all**b[1]*h_gf1_all**b[2]

				# Indices to trees found at each measurement
				c,ia,ib=gu.intersect(id0_all,id1_all)

				# Find index to lost trees - trees that were counted in t0 and
				# disappear (ia2).
				ia2=np.where(np.isin(id0_all,id1_all)==False)[0]

				# Find trees that were Recruited - trees that were
				# not counted at t0, but were at t1.
				ib2=np.where(np.isin(id1_all,id0_all)==False)[0]

				# Index to all live trees at t0
				ind0_live=np.where(vs0_all==1)[0]

				# Index to all live trees at t1
				ind1_live=np.where(vs1_all==1)[0]

				# Index to trees that survived
				ind_surv=np.where((vs0_all[ia]==1) & (vs1_all[ib]==1))[0]

				ind_surv_cc1=np.where((vs0_all[ia]==1) & (vs1_all[ib]==1) & (CrownClass0_all[ia]==1))[0]
				ind_surv_cc2=np.where((vs0_all[ia]==1) & (vs1_all[ib]==1) & (CrownClass0_all[ia]==2))[0]
				ind_surv_cc3=np.where((vs0_all[ia]==1) & (vs1_all[ib]==1) & (CrownClass0_all[ia]==3))[0]
				ind_surv_cc4=np.where((vs0_all[ia]==1) & (vs1_all[ib]==1) & (CrownClass0_all[ia]==4))[0]
				ind_surv_cc5=np.where((vs0_all[ia]==1) & (vs1_all[ib]==1) & (CrownClass0_all[ia]==5))[0]
				ind_surv_cc6=np.where((vs0_all[ia]==1) & (vs1_all[ib]==1) & (CrownClass0_all[ia]==6))[0]

				# Index to trees that died
				ind_mort=np.where((vs0_all[ia]==1) & (vs1_all[ib]==0))[0]

				ind_mort_cc1=np.where((vs0_all[ia]==1) & (vs1_all[ib]==0) & (CrownClass0_all[ia]==1))[0]
				ind_mort_cc2=np.where((vs0_all[ia]==1) & (vs1_all[ib]==0) & (CrownClass0_all[ia]==2))[0]
				ind_mort_cc3=np.where((vs0_all[ia]==1) & (vs1_all[ib]==0) & (CrownClass0_all[ia]==3))[0]
				ind_mort_cc4=np.where((vs0_all[ia]==1) & (vs1_all[ib]==0) & (CrownClass0_all[ia]==4))[0]
				ind_mort_cc5=np.where((vs0_all[ia]==1) & (vs1_all[ib]==0) & (CrownClass0_all[ia]==5))[0]
				ind_mort_cc6=np.where((vs0_all[ia]==1) & (vs1_all[ib]==0) & (CrownClass0_all[ia]==6))[0]

				# Index to lost live trees
				ind_lost=np.where( (vs0_all[ia2]==1) )[0]
				#print(id0_all[ia2[ind_lost]])
				#print(dbh0_all[ia2[ind_lost]])
				#print(dbh1_all[ib2[ind_lost]])

				# Index to recruited trees
				ind_rec=np.where((vs1_all[ib2]==1))[0]
				#print(id0_all[ia2[ind_rec]])

				# Calculate stand-level variables

				# Assign each species with a numerical ID
				ID_Spc_all_t0=np.zeros((id0_all.shape))
				ind=np.where(sp0_all=='FD')[0]; ID_Spc_all_t0[ind]=1
				ind=np.where(sp0_all=='HW')[0]; ID_Spc_all_t0[ind]=2
				ind=np.where(sp0_all=='CW')[0]; ID_Spc_all_t0[ind]=3
				ind=np.where(sp0_all=='BA')[0]; ID_Spc_all_t0[ind]=4
				ID_Spc_all_t1=np.zeros((id1_all.shape))
				ind=np.where(sp1_all=='FD')[0]; ID_Spc_all_t1[ind]=1
				ind=np.where(sp1_all=='HW')[0]; ID_Spc_all_t1[ind]=2
				ind=np.where(sp1_all=='CW')[0]; ID_Spc_all_t1[ind]=3
				ind=np.where(sp1_all=='BA')[0]; ID_Spc_all_t1[ind]=4

				# Species composition
				uSpc=np.unique(sp0_all)
				nSpc=np.zeros((uSpc.shape[0],2))
				for s in range(len(uSpc)):
					ind=np.where(sp0_all==uSpc[s])[0]
					nSpc[s,0]=s
					nSpc[s,1]=ind.size
				nSpc=np.flip(nSpc[nSpc[:,1].argsort(),:],axis=0)
				Spc1_CD_t0=uSpc[int(nSpc[0,0])]
				Spc1_Pct_t0=nSpc[0,1]/ind0.size*100
				if uSpc.size>1:
					Spc2_CD_t0=uSpc[int(nSpc[1,0])]
					Spc2_Pct_t0=nSpc[1,1]/ind0.size*100
				else:
					Spc2_CD_t0=' '
					Spc2_Pct_t0=0

				# Number of height measurements
				Num_H_Obs_t0=np.where(h_gf0_all[ind0_live]>0)[0].size
				
				# Number of missing or negative height growth
				H_G=(h_gf1_all[ib[ind_surv]]-h_gf0_all[ia[ind_surv]])/dt
				ind=np.where((H_G<0) | (np.isnan(H_G==True)))[0]
				Num_H_G_NegOrMis=ind.size
				
				# Stand age
				SA_t0=SA+uYear[iY]-YrEst
				SA_t1=SA+uYear[iY+1]-YrEst
				
				# Stand density (stems ha-1)
				SN_t0=AEF*ind0_live.size
				SN_t1=AEF*ind1_live.size
				SN_Net=(SN_t1-SN_t0)/dt
				
				# Mean diameter (cm)
				SDBH_t0=np.nanmean(dbh0_all[ind0_live])
				SDBH_t1=np.nanmean(dbh1_all[ind1_live])
	
				# Basal area (m2/ha)
				SBA_t0=AEF*np.nansum(ba0_all[ind0_live])/10000
				SBA_t1=AEF*np.nansum(ba1_all[ind1_live])/10000
				SBA_G=(SBA_t1-SBA_t0)/dt
				
				# Mean tree height (m)
				SH_obs_t0=np.nanmean(h_obs0_all[ind0_live])
				SH_obs_t1=np.nanmean(h_obs1_all[ind1_live])
				SH_obs_G=(SH_obs_t1-SH_obs_t0)/dt
				
				# Mean tree height (m)
				SH_gf_t0=np.nanmean(h_gf0_all[ind0_live])
				SH_gf_t1=np.nanmean(h_gf1_all[ind1_live])
				SH_gf_G=(SH_gf_t1-SH_gf_t0)/dt
	
				# Biomass, live (MgC/ha)
				SCsw_t0=meta['Param']['dm2c']*meta['Param']['kg2Mg']*AEF*np.nansum(b0_all['Stemwood'][ind0_live])
				SCf_t0=meta['Param']['dm2c']*meta['Param']['kg2Mg']*AEF*np.nansum(b0_all['Foliage'][ind0_live])
				SCbr_t0=meta['Param']['dm2c']*meta['Param']['kg2Mg']*AEF*np.nansum(b0_all['Branch'][ind0_live])
				SCbk_t0=meta['Param']['dm2c']*meta['Param']['kg2Mg']*AEF*np.nansum(b0_all['Bark'][ind0_live])
				SCsw_t1=meta['Param']['dm2c']*meta['Param']['kg2Mg']*AEF*np.nansum(b1_all['Stemwood'][ind1_live])
				SCf_t1=meta['Param']['dm2c']*meta['Param']['kg2Mg']*AEF*np.nansum(b1_all['Foliage'][ind1_live])
				SCbr_t1=meta['Param']['dm2c']*meta['Param']['kg2Mg']*AEF*np.nansum(b1_all['Branch'][ind1_live])
				SCbk_t1=meta['Param']['dm2c']*meta['Param']['kg2Mg']*AEF*np.nansum(b1_all['Bark'][ind1_live])
				
				# Biomass growth (MgC/ha/yr) 
				tmp0=np.nansum(b0_all['Stemwood'][ia[ind_surv]].flatten())
				tmp1=np.nansum(b1_all['Stemwood'][ib[ind_surv]].flatten())
				SCsw_G=meta['Param']['dm2c']*meta['Param']['kg2Mg']*AEF*(tmp1-tmp0)/dt
				
				tmp0=np.nansum(b0_all['Stemwood'][ia[ind_surv_cc1]].flatten())
				tmp1=np.nansum(b1_all['Stemwood'][ib[ind_surv_cc1]].flatten())
				SCsw_G_cc1=meta['Param']['dm2c']*meta['Param']['kg2Mg']*AEF*(tmp1-tmp0)/dt
				
				tmp0=np.nansum(b0_all['Stemwood'][ia[ind_surv_cc2]].flatten())
				tmp1=np.nansum(b1_all['Stemwood'][ib[ind_surv_cc2]].flatten())
				SCsw_G_cc2=meta['Param']['dm2c']*meta['Param']['kg2Mg']*AEF*(tmp1-tmp0)/dt
				
				tmp0=np.nansum(b0_all['Stemwood'][ia[ind_surv_cc3]].flatten())
				tmp1=np.nansum(b1_all['Stemwood'][ib[ind_surv_cc3]].flatten())
				SCsw_G_cc3=meta['Param']['dm2c']*meta['Param']['kg2Mg']*AEF*(tmp1-tmp0)/dt
				
				tmp0=np.nansum(b0_all['Stemwood'][ia[ind_surv_cc4]].flatten())
				tmp1=np.nansum(b1_all['Stemwood'][ib[ind_surv_cc4]].flatten())
				SCsw_G_cc4=meta['Param']['dm2c']*meta['Param']['kg2Mg']*AEF*(tmp1-tmp0)/dt
				
				tmp0=np.nansum(b0_all['Stemwood'][ia[ind_surv_cc5]].flatten())
				tmp1=np.nansum(b1_all['Stemwood'][ib[ind_surv_cc5]].flatten())
				SCsw_G_cc5=meta['Param']['dm2c']*meta['Param']['kg2Mg']*AEF*(tmp1-tmp0)/dt
				
				tmp0=np.nansum(b0_all['Stemwood'][ia[ind_surv_cc6]].flatten())
				tmp1=np.nansum(b1_all['Stemwood'][ib[ind_surv_cc6]].flatten())
				SCsw_G_cc6=meta['Param']['dm2c']*meta['Param']['kg2Mg']*AEF*(tmp1-tmp0)/dt
				
				tmp0=np.nansum(b0_all['Foliage'][ia[ind_surv]].flatten())
				tmp1=np.nansum(b1_all['Foliage'][ib[ind_surv]].flatten())
				SCf_G=meta['Param']['dm2c']*meta['Param']['kg2Mg']*AEF*(tmp1-tmp0)/dt
				
				tmp0=np.nansum(b0_all['Bark'][ia[ind_surv]].flatten())
				tmp1=np.nansum(b1_all['Bark'][ib[ind_surv]].flatten())
				SCbk_G=meta['Param']['dm2c']*meta['Param']['kg2Mg']*AEF*(tmp1-tmp0)/dt
				
				tmp0=np.nansum(b0_all['Branch'][ia[ind_surv]].flatten())
				tmp1=np.nansum(b1_all['Branch'][ib[ind_surv]].flatten())
				SCbr_G=meta['Param']['dm2c']*meta['Param']['kg2Mg']*AEF*(tmp1-tmp0)/dt
				
				# Biomass mortality (MgC/ha/yr)
				SCsw_M=meta['Param']['dm2c']*meta['Param']['kg2Mg']*AEF*np.nansum(b1_all['Stemwood'][ib[ind_mort]].flatten())/dt
				SCf_M=meta['Param']['dm2c']*meta['Param']['kg2Mg']*AEF*np.nansum(b1_all['Foliage'][ib[ind_mort]].flatten())/dt
				SCbr_M=meta['Param']['dm2c']*meta['Param']['kg2Mg']*AEF*np.nansum(b1_all['Branch'][ib[ind_mort]].flatten())/dt
				SCbk_M=meta['Param']['dm2c']*meta['Param']['kg2Mg']*AEF*np.nansum(b1_all['Bark'][ib[ind_mort]].flatten())/dt
				
				SCsw_M_cc1=meta['Param']['dm2c']*meta['Param']['kg2Mg']*AEF*np.nansum(b1_all['Stemwood'][ib[ind_mort_cc1]].flatten())/dt
				SCsw_M_cc2=meta['Param']['dm2c']*meta['Param']['kg2Mg']*AEF*np.nansum(b1_all['Stemwood'][ib[ind_mort_cc2]].flatten())/dt
				SCsw_M_cc3=meta['Param']['dm2c']*meta['Param']['kg2Mg']*AEF*np.nansum(b1_all['Stemwood'][ib[ind_mort_cc3]].flatten())/dt
				SCsw_M_cc4=meta['Param']['dm2c']*meta['Param']['kg2Mg']*AEF*np.nansum(b1_all['Stemwood'][ib[ind_mort_cc4]].flatten())/dt
				SCsw_M_cc5=meta['Param']['dm2c']*meta['Param']['kg2Mg']*AEF*np.nansum(b1_all['Stemwood'][ib[ind_mort_cc5]].flatten())/dt
				SCsw_M_cc6=meta['Param']['dm2c']*meta['Param']['kg2Mg']*AEF*np.nansum(b1_all['Stemwood'][ib[ind_mort_cc6]].flatten())/dt
	
				# Biomass recruitment (MgC/ha/yr)
				if ind_rec.size!=0:
					SCsw_R=meta['Param']['dm2c']*meta['Param']['kg2Mg']*AEF*np.nansum(b1_all['Stemwood'][ib2[ind_rec]].flatten())/dt
					SCf_R=meta['Param']['dm2c']*meta['Param']['kg2Mg']*AEF*np.nansum(b1_all['Foliage'][ib2[ind_rec]].flatten())/dt
					SCbr_R=meta['Param']['dm2c']*meta['Param']['kg2Mg']*AEF*np.nansum(b1_all['Branch'][ib2[ind_rec]].flatten())/dt
					SCbk_R=meta['Param']['dm2c']*meta['Param']['kg2Mg']*AEF*np.nansum(b1_all['Bark'][ib2[ind_rec]].flatten())/dt
				else:
					SCsw_R=0.0
					SCf_R=0.0
					SCbr_R=0.0
					SCbk_R=0.0

				# Biomass lost (MgC/ha/yr)
				if ind_lost.size!=0:
					SCsw_L=meta['Param']['dm2c']*meta['Param']['kg2Mg']*AEF*np.nansum(b0_all['Stemwood'][ia2[ind_lost]].flatten())/dt
					SCf_L=meta['Param']['dm2c']*meta['Param']['kg2Mg']*AEF*np.nansum(b0_all['Foliage'][ia2[ind_lost]].flatten())/dt
					SCbr_L=meta['Param']['dm2c']*meta['Param']['kg2Mg']*AEF*np.nansum(b0_all['Branch'][ia2[ind_lost]].flatten())/dt
					SCbk_L=meta['Param']['dm2c']*meta['Param']['kg2Mg']*AEF*np.nansum(b0_all['Bark'][ia2[ind_lost]].flatten())/dt
				else:
					SCsw_L=0.0
					SCf_L=0.0
					SCbr_L=0.0
					SCbk_L=0.0

				# Net growth
				SCsw_Net=SCsw_R+SCsw_G-SCsw_M-SCsw_L
				SCsw_Net_cc1=SCsw_G_cc1-SCsw_M_cc1
				SCsw_Net_cc2=SCsw_G_cc2-SCsw_M_cc2
				SCsw_Net_cc3=SCsw_G_cc3-SCsw_M_cc3
				SCsw_Net_cc4=SCsw_G_cc4-SCsw_M_cc4
				SCsw_Net_cc5=SCsw_G_cc5-SCsw_M_cc5
				SCsw_Net_cc6=SCsw_G_cc6-SCsw_M_cc6

				# Damage agents
				da0_t={}
				da1_t={}
				for da in meta['Param']['da_list']:
					ind=np.where(da0_all==da)[0]
					da0_t['DA_' + str(da) + '_t0']=0
					if ind.size>0:
						da0_t['DA_' + str(da) + '_t0']=ind.size
					ind=np.where(da1_all==da)[0]
					da1_t['DA_' + str(da) + '_t1']=0
					if ind.size>0:
						da1_t['DA_' + str(da) + '_t1']=ind.size

				# Add to stand-level dataframe
				DataToAdd={'ID_Site':uSite[iS],
				   'ID_Plot':iP,
				   'Lat':Lat,
				   'Lon':Lon,
				   'BGC':dSite['BGC class'][indSite[0]],
				   'SI':dSite['Site Index'][indSite[0]],
				   'EstabType':dSite['Establishment Type'][indSite[0]],
				   'AreaPlot':Area,
				   'N_Init':dSite['Initial SPH'][indSite[0]].astype('float'),
				   'pH_Min_Init':dSite['pH Mineral'][indSite[0]],
				   'pH_Humus_Init':dSite['pH Humus'][indSite[0]],
				   'P_Min_Init':dSite['P Mineral'][indSite[0]],
				   'P_Humus_Init':dSite['P Humus'][indSite[0]],
				   'Ca_Min_Init':dSite['Ca Mineral'][indSite[0]],
				   'NA_Dose':dSite['N Dose'][indP],
				   'NA_Type':dSite['N Type'][indP],
				   'NA_Method':dSite['N Method'][indP],
				   'Thin_Level':dSite['Thinning Level'][indP],
				   'Thin_Method':dSite['Thinning Method'][indP],
				   'Thin_Format':dSite['Thinning Format'][indP],
				   'Thin_TargetSPH':dSite['Thinning Target SPH'][indP],
				   'NumTreesMeas':np.array(ia.size+ib.size).astype('float'),
				   'Num_H_G_NegOrMis':Num_H_G_NegOrMis,
				   'Spc1_CD_t0':Spc1_CD_t0,
				   'Spc1_Pct_t0':Spc1_Pct_t0,
				   'Spc2_CD_t0':Spc2_CD_t0,
				   'Spc2_Pct_t0':Spc2_Pct_t0,
				   'Year_t0':yr0,
				   'Year_t1':yr1,
				   'DT':dt,
				   'Age_t0':SA_t0.astype('float'),
				   'Age_t1':SA_t1.astype('float'),
				   'TSF_t0':yr0-YrEst,
				   'TSF_t1':yr1-YrEst,
				   'N_t0':SN_t0.astype('float'),
				   'N_t1':SN_t1.astype('float'),
				   'N_Net':SN_Net.astype('float'),
				   'DBH_t0':SDBH_t0,
				   'DBH_t1':SDBH_t1,
				   'BA_t0':SBA_t0,
				   'BA_t1':SBA_t1,
				   'BA_G':SBA_G,
				   'H_obs_t0':SH_obs_t0,
				   'H_obs_t1':SH_obs_t1,
				   'H_obs_G':SH_obs_G,
				   'H_gf_t0':SH_gf_t0,
				   'H_gf_t1':SH_gf_t1,
				   'H_gf_G':SH_gf_G,
				   'Csw_t0':SCsw_t0,
				   'Csw_t1':SCsw_t1,
				   'Csw_G':SCsw_G,
				   'Csw_M':SCsw_M,
				   'Csw_R':SCsw_R,
				   'Csw_L':SCsw_L,
				   'Csw_Net':SCsw_Net,
				   'Cf_t0':SCf_t0,
				   'Cf_t1':SCf_t1,
				   'Cf_G':SCf_G,
				   'Cf_M':SCf_M,
				   'Cf_R':SCf_R,
				   'Cf_L':SCf_L,
				   'Cf_Net':SCf_R+SCf_G-SCf_M,
				   'Cbr_t0':SCbr_t0,
				   'Cbr_t1':SCbr_t1,
				   'Cbr_G':SCbr_G,
				   'Cbr_M':SCbr_M,
				   'Cbr_R':SCbr_R,
				   'Cbr_L':SCbr_L,
				   'Cbr_Net':SCbr_R+SCbr_G-SCbr_M,
				   'Cbk_t0':SCbk_t0,
				   'Cbk_t1':SCbk_t1,
				   'Cbk_G':SCbk_G,
				   'Cbk_M':SCbk_M,
				   'Cbk_R':SCbk_R,
				   'Cbk_L':SCbk_L,
				   'Cbk_Net':SCbk_R+SCbk_G-SCbk_M,
				   'Csw_G_cc1':SCsw_G_cc1,
				   'Csw_G_cc2':SCsw_G_cc2,
				   'Csw_G_cc3':SCsw_G_cc3,
				   'Csw_G_cc4':SCsw_G_cc4,
				   'Csw_G_cc5':SCsw_G_cc5,
				   'Csw_G_cc6':SCsw_G_cc6,
				   'Csw_M_cc1':SCsw_M_cc1,
				   'Csw_M_cc2':SCsw_M_cc2,
				   'Csw_M_cc3':SCsw_M_cc3,
				   'Csw_M_cc4':SCsw_M_cc4,
				   'Csw_M_cc5':SCsw_M_cc5,
				   'Csw_M_cc6':SCsw_M_cc6,
				   'Csw_Net_cc1':SCsw_Net_cc1,
				   'Csw_Net_cc2':SCsw_Net_cc2,
				   'Csw_Net_cc3':SCsw_Net_cc3,
				   'Csw_Net_cc4':SCsw_Net_cc4,
				   'Csw_Net_cc5':SCsw_Net_cc5,
				   'Csw_Net_cc6':SCsw_Net_cc6,
				   'DA_A_t0':da0_t['DA_A_t0'],
				   'DA_AB_t0':da0_t['DA_AB_t0'],
				   'DA_AD_t0':da0_t['DA_AD_t0'],
				   'DA_AE_t0':da0_t['DA_AE_t0'],
				   'DA_AX_t0':da0_t['DA_AX_t0'],
				   'DA_BT_t0':da0_t['DA_BT_t0'],
				   'DA_CK_t0':da0_t['DA_CK_t0'],
				   'DA_D_t0':da0_t['DA_D_t0'],
				   'DA_DDF_t0':da0_t['DA_DDF_t0'],
				   'DA_DDP_t0':da0_t['DA_DDP_t0'],
				   'DA_DDS_t0':da0_t['DA_DDS_t0'],
				   'DA_DM_t0':da0_t['DA_DM_t0'],
				   'DA_DN_t0':da0_t['DA_DN_t0'],
				   'DA_DR_t0':da0_t['DA_DR_t0'],
				   'DA_DRA_t0':da0_t['DA_DRA_t0'],
				   'DA_DS_t0':da0_t['DA_DS_t0'],
				   'DA_DT_t0':da0_t['DA_DT_t0'],
				   'DA_FCK_t0':da0_t['DA_FCK_t0'],
				   'DA_FK_t0':da0_t['DA_FK_t0'],
				   'DA_MT_t0':da0_t['DA_MT_t0'],
				   'DA_N_t0':da0_t['DA_N_t0'],
				   'DA_NGC_t0':da0_t['DA_NGC_t0'],
				   'DA_NW_t0':da0_t['DA_NW_t0'],
				   'DA_NWS_t0':da0_t['DA_NWS_t0'],
				   'DA_NX_t0':da0_t['DA_NX_t0'],
				   'DA_NY_t0':da0_t['DA_NY_t0'],
				   'DA_PC_t0':da0_t['DA_PC_t0'],
				   'DA_SC_t0':da0_t['DA_SC_t0'],
				   'DA_SCB_t0':da0_t['DA_SCB_t0'],
				   'DA_SCM_t0':da0_t['DA_SCM_t0'],
				   'DA_SCT_t0':da0_t['DA_SCT_t0'],
				   'DA_A_t1':da1_t['DA_A_t1'],
				   'DA_AB_t1':da1_t['DA_AB_t1'],
				   'DA_AD_t1':da1_t['DA_AD_t1'],
				   'DA_AE_t1':da1_t['DA_AE_t1'],
				   'DA_AX_t1':da1_t['DA_AX_t1'],
				   'DA_BT_t1':da1_t['DA_BT_t1'],
				   'DA_CK_t1':da1_t['DA_CK_t1'],
				   'DA_D_t1':da1_t['DA_D_t1'],
				   'DA_DDF_t1':da1_t['DA_DDF_t1'],
				   'DA_DDP_t1':da1_t['DA_DDP_t1'],
				   'DA_DDS_t1':da1_t['DA_DDS_t1'],
				   'DA_DM_t1':da1_t['DA_DM_t1'],
				   'DA_DN_t1':da1_t['DA_DN_t1'],
				   'DA_DR_t1':da1_t['DA_DR_t1'],
				   'DA_DRA_t1':da1_t['DA_DRA_t1'],
				   'DA_DS_t1':da1_t['DA_DS_t1'],
				   'DA_DT_t1':da1_t['DA_DT_t1'],
				   'DA_FCK_t1':da1_t['DA_FCK_t1'],
				   'DA_FK_t1':da1_t['DA_FK_t1'],
				   'DA_MT_t1':da1_t['DA_MT_t1'],
				   'DA_N_t1':da1_t['DA_N_t1'],
				   'DA_NGC_t1':da1_t['DA_NGC_t1'],
				   'DA_NW_t1':da1_t['DA_NW_t1'],
				   'DA_NWS_t1':da1_t['DA_NWS_t1'],
				   'DA_NX_t1':da1_t['DA_NX_t1'],
				   'DA_NY_t1':da1_t['DA_NY_t1'],
				   'DA_PC_t1':da1_t['DA_PC_t1'],
				   'DA_SC_t1':da1_t['DA_SC_t1'],
				   'DA_SCB_t1':da1_t['DA_SCB_t1'],
				   'DA_SCM_t1':da1_t['DA_SCM_t1'],
				   'DA_SCT_t1':da1_t['DA_SCT_t1']}

				for k in DataToAdd:
					sobs[k][cntS]=DataToAdd[k]
				cntS=cntS+1
				#sobs=sobs.append(DataToAdd.copy(),ignore_index=True)
				
				# Add survivor data to tree-level dataframe
				ID_Spc_t0=ID_Spc_all_t0[ia[ind_surv]]
				ID_Tree=id0_all[ia[ind_surv]]
				VS_t0=vs0_all[ia[ind_surv]]
				D_t0=dbh0_all[ia[ind_surv]]
				D_t1=dbh1_all[ib[ind_surv]]
				BA_t0=ba0_all[ia[ind_surv]]
				BA_t1=ba1_all[ib[ind_surv]]
				BA_G=(BA_t1-BA_t0)/dt
				H_obs_t0=h_obs0_all[ia[ind_surv]]
				H_obs_t1=h_obs1_all[ib[ind_surv]]
				H_obs_G=(H_obs_t1-H_obs_t0)/dt
				H_mod_t0=h_mod0_all[ia[ind_surv]]
				H_mod_t1=h_mod1_all[ib[ind_surv]]
				H_mod_G=(H_mod_t1-H_mod_t0)/dt
				H_gf_t0=h_gf0_all[ia[ind_surv]]
				H_gf_t1=h_gf1_all[ib[ind_surv]]
				H_gf_G=(H_gf_t1-H_gf_t0)/dt
				H_Fit_ITM_t0=H_Fit_ITM0_all[ia[ind_surv]]
				H_Fit_ITM_t1=H_Fit_ITM1_all[ib[ind_surv]]
				Csw_t0=meta['Param']['dm2c']*b0_all['Stemwood'][ia[ind_surv]]
				Csw_t1=meta['Param']['dm2c']*b1_all['Stemwood'][ib[ind_surv]]
				Csw_G=meta['Param']['dm2c']*(b1_all['Stemwood'][ib[ind_surv]]-b0_all['Stemwood'][ia[ind_surv]])/dt
				Csw_mod_t0=meta['Param']['dm2c']*b0_all['Stemwood'][ia[ind_surv]]
				Csw_mod_t1=meta['Param']['dm2c']*b1_all['Stemwood'][ib[ind_surv]]
				TreeClass_t0=TreeClass0_all[ia[ind_surv]]
				TreeClass_t1=TreeClass1_all[ib[ind_surv]]
				CrownClass_t0=CrownClass0_all[ia[ind_surv]]
				CrownClass_t1=CrownClass1_all[ib[ind_surv]]
				TopDead_t0=dtop0_all[ia[ind_surv]]
				TopDead_t1=dtop1_all[ib[ind_surv]]
				TopBroke_t0=btop0_all[ia[ind_surv]]
				TopBroke_t1=btop1_all[ib[ind_surv]]
				DA_t0=da0_all[ia[ind_surv]]
				DA_t1=da1_all[ib[ind_surv]]
				
				# Stand biomass of larger trees
				tlist=np.zeros((ind_surv.size,3))
				tlist[:,0]=np.arange(0,ind_surv.size,1)
				tlist[:,1]=Csw_t0
				tlist=np.nan_to_num(tlist)
				tlist=np.flip(tlist[tlist[:,1].argsort(),:],axis=0)
				tlist[:,2]=np.cumsum(tlist[:,1]*meta['Param']['kg2Mg']*AEF)
				SCswLT_t0=np.zeros(ind_surv.shape)
				SCswLT_t0[tlist[:,0].astype(int)]=tlist[:,2]
				#plt.plot(Csw_t0,SCswLT,'.')
				
				for v in range(ind_surv.size):
					DataToAdd={'ID_Site':uSite[iS],
					   'ID_Plot':iP,
					   'ID_Tree':ID_Tree[v],
					   'Lat':Lat,
					   'Lon':Lon,'BGC':dSite['BGC class'][indSite[0]],
					   'SI':dSite['Site Index'][indSite[0]],
					   'AEF':AEF,'N_Init':dSite['Initial SPH'][indSite[0]],
					   'pH_Min_Init':dSite['pH Mineral'][indSite[0]],
					   'NA_Dose':dSite['N Dose'][indP].astype('float'),
					   'Year_t0':yr0,'Year_t1':yr1,
					   'DT':dt,
					   'TSF_t0':yr0-YrEst,
					   'TSF_t1':yr1-YrEst,
					   'SA_t0':SA_t0.astype('float'),
					   'SA_t1':SA_t1.astype('float'),
					   'SN_t0':SN_t0,
					   'SN_t1':SN_t1,
					   'SCsw_t0':SCsw_t0,
					   'SCsw_t1':SCsw_t1,
					   'SCswLT_t0':SCswLT_t0[v],
					   'SCswLT_t1':0,
					   'ID_Spc_t0':ID_Spc_t0[v],
					   'Mort':0,'Rec':0,
					   'TreeClass_t0':TreeClass_t0[v],
					   'CrownClass_t0':CrownClass_t0[v],
					   'TopDead_t0':TopDead_t0[v],'TopDead_t1':TopDead_t1[v],
					   'TopBroke_t0':TopBroke_t0[v],'TopBroke_t1':TopBroke_t1[v],
					   'DA_t0':DA_t0[v],'DA_t1':DA_t1[v],
					   'D_t0':D_t0[v],'D_t1':D_t1[v],
					   'BA_t0':BA_t0[v],'BA_t1':BA_t1[v],'BA_G':BA_G[v],
					   'H_obs_t0':H_obs_t0[v],'H_obs_t1':H_obs_t1[v],'H_obs_G':H_obs_G[v],
					   'H_mod_t0':H_mod_t0[v],'H_mod_t1':H_mod_t1[v],'H_mod_G':H_mod_G[v],
					   'H_gf_t0':H_gf_t0[v],'H_gf_t1':H_gf_t1[v],'H_gf_G':H_gf_G[v],
					   'Csw_t0':Csw_t0[v],'Csw_t1':Csw_t1[v],'Csw_G':Csw_G[v],
					   'Csw_mod_t0':Csw_mod_t0[v],'Csw_mod_t1':Csw_mod_t1[v],
					   'H_Fit_ITM_t0':H_Fit_ITM_t0[v],'H_Fit_ITM_t1':H_Fit_ITM_t1[v]}
					for k in DataToAdd:
						tobs[k][cntT]=DataToAdd[k]
					cntT=cntT+1
	
				# Add mortality data to tree-level dataframe
				if ind_mort.size>0:
					ID_Spc_t0=ID_Spc_all_t0[ia[ind_mort]]
				
					ID_Tree=id0_all[ia[ind_mort]]
				
					VS_t0=vs0_all[ia[ind_mort]]
					  
					D_t0=dbh0_all[ia[ind_mort]]
					D_t1=dbh1_all[ib[ind_mort]]
					
					BA_t0=ba0_all[ia[ind_mort]]
					BA_t1=ba1_all[ib[ind_mort]]
					BA_G=(BA_t1-BA_t0)/dt
				
					H_obs_t0=h_obs0_all[ia[ind_mort]]
					H_obs_t1=h_obs1_all[ib[ind_mort]]
					H_obs_G=np.nan*np.ones(ind_mort.size,dtype=float)
					
					H_mod_t0=h_mod0_all[ia[ind_mort]]
					H_mod_t1=h_mod1_all[ib[ind_mort]]
					H_mod_G=np.nan*np.ones(ind_mort.size,dtype=float)
					
					H_gf_t0=h_gf0_all[ia[ind_mort]]
					H_gf_t1=h_gf1_all[ib[ind_mort]]
					H_gf_G=np.nan*np.ones(ind_mort.size,dtype=float)
	
					H_Fit_ITM_t0=H_Fit_ITM0_all[ia[ind_mort]]
					H_Fit_ITM_t1=H_Fit_ITM1_all[ib[ind_mort]]
				
					Csw_t0=meta['Param']['dm2c']*b0_all['Stemwood'][ia[ind_mort]]
					Csw_t1=meta['Param']['dm2c']*b1_all['Stemwood'][ib[ind_mort]]
					Csw_G=meta['Param']['dm2c']*(b1_all['Stemwood'][ib[ind_mort]]-b0_all['Stemwood'][ia[ind_mort]])/dt
					
					TreeClass_t0=TreeClass0_all[ia[ind_mort]]
					TreeClass_t1=TreeClass1_all[ib[ind_mort]]
				
					CrownClass_t0=CrownClass0_all[ia[ind_mort]]
					CrownClass_t1=CrownClass1_all[ib[ind_mort]]
				
					TopDead_t0=dtop0_all[ia[ind_mort]]
					TopDead_t1=dtop1_all[ib[ind_mort]]
				
					TopBroke_t0=btop0_all[ia[ind_mort]]
					TopBroke_t1=btop1_all[ib[ind_mort]]
				
					DA_t0=da0_all[ia[ind_mort]]
					DA_t1=da1_all[ib[ind_mort]]
					
					# Stand biomass of larger trees
					tlist=np.zeros((ind_mort.size,3))
					tlist[:,0]=np.arange(0,ind_mort.size,1)
					tlist[:,1]=Csw_t0
					tlist=np.nan_to_num(tlist)
					tlist=np.flip(tlist[tlist[:,1].argsort(),:],axis=0)
					tlist[:,2]=np.cumsum(tlist[:,1]*meta['Param']['kg2Mg']*AEF)
					SCswLT_t0=np.zeros(ind_mort.shape)
					SCswLT_t0[tlist[:,0].astype(int)]=tlist[:,2]
					#plt.plot(Csw_t0,SCswLT,'.')

					for v in range(ind_mort.size):
						DataToAdd={'ID_Site':uSite[iS],
								   'ID_Plot':iP,
								   'ID_Tree':ID_Tree[v],
								   'Lat':Lat,
								   'Lon':Lon,'BGC':dSite['BGC class'][indSite[0]],'SI':dSite['Site Index'][indSite[0]],
								   'AEF':AEF,'N_Init':dSite['Initial SPH'][indSite[0]],
								   'pH_Min_Init':dSite['pH Mineral'][indSite[0]],
								   'NA_Dose':dSite['N Dose'][indP].astype('float'),
								   'Year_t0':yr0,'Year_t1':yr1,'DT':dt,
								   'TSF_t0':yr0-YrEst,'TSF_t1':yr1-YrEst,
								   'SA_t0':SA_t0.astype('float'),'SA_t1':SA_t1.astype('float'),
								   'SN_t0':SN_t0,'SN_t1':SN_t1,
								   'SCsw_t0':SCsw_t0,'SCsw_t1':SCsw_t1,'SCswLT_t0':SCswLT_t0[v],'SCswLT_t1':0,
								   'ID_Spc_t0':ID_Spc_t0[v],'Mort':1,'Rec':0,'TreeClass_t0':TreeClass_t0[v],
								   'CrownClass_t0':CrownClass_t0[v],'TopDead_t0':TopDead_t0[v],'TopDead_t1':TopDead_t1[v],
								   'TopBroke_t0':TopBroke_t0[v],'TopBroke_t1':TopBroke_t1[v],'DA_t0':DA_t0[v],'DA_t1':DA_t1[v],
								   'D_t0':D_t0[v],'D_t1':D_t1[v],
								   'BA_t0':BA_t0[v],'BA_t1':BA_t1[v],'BA_G':BA_G[v],
								   'H_obs_t0':H_obs_t0[v],'H_obs_t1':H_obs_t1[v],'H_obs_G':H_obs_G[v],
								   'H_mod_t0':H_mod_t0[v],'H_mod_t1':H_mod_t1[v],'H_mod_G':H_mod_G[v],
								   'H_gf_t0':H_gf_t0[v],'H_gf_t1':H_gf_t1[v],'H_gf_G':H_gf_G[v],
								   'Csw_t0':Csw_t0[v],'Csw_t1':Csw_t1[v],'Csw_G':Csw_G[v],
								   'Csw_mod_t0':np.nan,'Csw_mod_t1':np.nan,
								   'H_Fit_ITM_t0':H_Fit_ITM_t0[v],'H_Fit_ITM_t1':H_Fit_ITM_t1[v]}
						for k in DataToAdd:
							tobs[k][cntT]=DataToAdd[k]
						cntT=cntT+1 

				# Add recruitment data to tree-level dataframe
				if ib2.size>0:
					n=ib2[ind_rec].size
					ID_Spc_t1=ID_Spc_all_t1[ib2[ind_rec]]
					ID_Tree=id1_all[ib2[ind_rec]]
					VS_t1=vs1_all[ib2[ind_rec]]
						 
					D_t0=np.nan*np.ones(n,dtype=float)
					D_t1=dbh1_all[ib2[ind_rec]]
					
					BA_t0=np.nan*np.ones(n,dtype=float)
					BA_t1=ba1_all[ib2[ind_rec]]
					BA_G=np.nan*np.ones(n,dtype=float)
				
					H_obs_t0=np.nan*np.ones(n,dtype=float)
					H_obs_t1=h_obs1_all[ib2[ind_rec]]
					H_obs_G=np.nan*np.ones(n,dtype=float)
				
					H_mod_t0=np.nan*np.ones(n,dtype=float)
					H_mod_t1=h_mod1_all[ib2[ind_rec]]
					H_mod_G=np.nan*np.ones(n,dtype=float)
					
					H_gf_t0=np.nan*np.ones(n,dtype=float)
					H_gf_t1=h_gf1_all[ib2[ind_rec]]
					H_gf_G=np.nan*np.ones(n,dtype=float)
				
					H_Fit_ITM_t0=np.nan*np.ones(n,dtype=float)
					H_Fit_ITM_t1=H_Fit_ITM1_all[ib2[ind_rec]]
				
					Csw_t0=np.nan*np.ones(n,dtype=float)
					Csw_t1=meta['Param']['dm2c']*b1_all['Stemwood'][ib2[ind_rec]]
					Csw_G=np.nan*np.ones(n,dtype=float)
					
					TreeClass_t0=np.nan*np.ones(n,dtype=float)
					TreeClass_t1=TreeClass1_all[ib2[ind_rec]]
				
					CrownClass_t0=np.nan*np.ones(n,dtype=float)
					CrownClass_t1=CrownClass1_all[ib2[ind_rec]]
				
					TopDead_t0=np.nan*np.ones(n,dtype=float)
					TopDead_t1=dtop1_all[ib2[ind_rec]]
				
					TopBroke_t0=np.nan*np.ones(n,dtype=float)
					TopBroke_t1=btop1_all[ib2[ind_rec]]
				
					DA_t0=np.nan*np.ones(n,dtype=float)
					DA_t1=da1_all[ib2[ind_rec]]
					
					# Stand biomass of larger trees
					tlist=np.zeros((ind_rec.size,3))
					tlist[:,0]=np.arange(0,ind_rec.size,1)
					tlist[:,1]=Csw_t1
					tlist=np.nan_to_num(tlist)
					tlist=np.flip(tlist[tlist[:,1].argsort(),:],axis=0)
					tlist[:,2]=np.cumsum(tlist[:,1]*meta['Param']['kg2Mg']*AEF)
					SCswLT_t1=np.zeros(ind_rec.shape)
					SCswLT_t1[tlist[:,0].astype(int)]=tlist[:,2]
					#plt.plot(Csw_t0,SCswLT,'.')

					for v in range(ind_rec.size):
						DataToAdd={'ID_Site':uSite[iS],
								   'ID_Plot':iP,
								   'ID_Tree':ID_Tree[v],
								   'Lat':Lat,
								   'Lon':Lon,'BGC':dSite['BGC class'][indSite[0]],'SI':dSite['Site Index'][indSite[0]],
								   'AEF':AEF,'N_Init':dSite['Initial SPH'][indSite[0]],
								   'pH_Min_Init':dSite['pH Mineral'][indSite[0]],
								   'NA_Dose':dSite['N Dose'][indP].astype('float'),
								   'Year_t0':yr0,'Year_t1':yr1,'DT':dt,
								   'TSF_t0':yr0-YrEst,'TSF_t1':yr1-YrEst,
								   'SA_t0':SA_t0.astype('float'),'SA_t1':SA_t1.astype('float'),
								   'SN_t0':SN_t0,'SN_t1':SN_t1,
								   'SCsw_t0':SCsw_t0,'SCsw_t1':SCsw_t1,'SCswLT_t0':0,'SCswLT_t1':SCswLT_t1[v],
								   'ID_Spc_t0':ID_Spc_t1[v],'Mort':0,'Rec':1,'TreeClass_t0':TreeClass_t0[v],
								   'CrownClass_t0':CrownClass_t0[v],'TopDead_t0':TopDead_t0[v],'TopDead_t1':TopDead_t1[v],
								   'TopBroke_t0':TopBroke_t0[v],'TopBroke_t1':TopBroke_t1[v],'DA_t0':DA_t0[v],'DA_t1':DA_t1[v],
								   'D_t0':D_t0[v],'D_t1':D_t1[v],
								   'BA_t0':BA_t0[v],'BA_t1':BA_t1[v],'BA_G':BA_G[v],
								   'H_obs_t0':H_obs_t0[v],'H_obs_t1':H_obs_t1[v],'H_obs_G':H_obs_G[v],
								   'H_mod_t0':H_mod_t0[v],'H_mod_t1':H_mod_t1[v],'H_mod_G':H_mod_G[v],
								   'H_gf_t0':H_gf_t0[v],'H_gf_t1':H_gf_t1[v],'H_gf_G':H_gf_G[v],
								   'Csw_t0':Csw_t0[v],'Csw_t1':Csw_t1[v],'Csw_G':Csw_G[v],
								   'Csw_mod_t0':np.nan,'Csw_mod_t1':np.nan,
								   'H_Fit_ITM_t0':H_Fit_ITM_t0[v],'H_Fit_ITM_t1':H_Fit_ITM_t1[v]}
						for k in DataToAdd:
							tobs[k][cntT]=DataToAdd[k]
						cntT=cntT+1

	# Remove excess zeros
	for k in sobs:
		sobs[k]=sobs[k][0:cntS]
	for k in tobs:
		tobs[k]=tobs[k][0:cntT]
	
	# Convert to float
	tobs['SA_t0']=tobs['SA_t0'].astype('float')
	tobs['SCswLT_t0']=tobs['SCswLT_t0'].astype('float')
	tobs['TSF_t0']=tobs['TSF_t0'].astype('float')
	tobs['TSF_t1']=tobs['TSF_t1'].astype('float')
	tobs['NA_Dose']=tobs['NA_Dose'].astype('float')
	
	# Save to file
	gu.opickle(meta['Paths']['Outputs'] + '\\EP703_SL.pkl',sobs)
	gu.opickle(meta['Paths']['Outputs'] + '\\EP703_TL.pkl',tobs)
	
	return sobs

#%% Add environmental data
def ImportEnvironmentalData(meta):
	
	# Import data
	sobs=gu.ipickle(meta['Paths']['Project']['Outputs'] + '\\EP703_SL.pkl')
	tobs=gu.ipickle(meta['Paths']['Project']['Outputs'] + '\\EP703_TL.pkl')
	
	# Unique location of installations
	uLL=np.unique(np.column_stack([sobs['Lat'],sobs['Lon']]),axis=0)
	
	# Spatial reference system of NACID
	SRS_NACID=pyproj.Proj('+proj=lcc +lat_1=49 +lat_2=77 +lat_0=40 +lon_0=-100 +x_0=0 +y_0=0 +no_defs +a=6378137 +rf=298.257222101 +to_meter=1')
	
	# Map of plots
	plt.close('all')
	fig,ax=plt.subplots(figsize=(9,9))
	#gdf1.plot(ax=ax,color=[0.5,0.5,0.5],edgecolor='k',label='Land')
	plt.plot(uLL[:,1],uLL[:,0],'go')
	ax.grid(color='k', linestyle='-', linewidth=0.25)
	
	# Calculate projected coordinates
	x,y=SRS_NACID(uLL[:,1],uLL[:,0])
	
	# Import NACID X and Y grids
	mat=spio.loadmat(meta['Paths']['NACID'] + '\\grid.mat',squeeze_me=True)
	ix=np.where(np.asarray(mat['grd'].dtype.names)=='xL')[0][0]
	iy=np.where(np.asarray(mat['grd'].dtype.names)=='yL')[0][0]
	X=mat['grd'][()][ix].astype('float')[0,:]
	Y=mat['grd'][()][iy].astype('float')[:,0]
	
	# Index to NACID grid
	ixy=np.zeros((x.shape[0],2),dtype=int)
	for i in range(len(ixy)):
		adx=np.abs(x[i]-X)
		ady=np.abs(y[i]-Y)
		ixy[i,0]=int(np.where(adx==np.min(adx))[0][0])
		ixy[i,1]=int(np.where(ady==np.min(ady))[0][0])
	
	# Check that we are indexing the right location
	#plt.close('all')
	#plt.imshow(tmean0,vmin=0,vmax=20)
	#for i in range(len(x)):
	#	plt.plot(ixy[i,0],ixy[i,1],'yo')
	
	# Define a time period for the annual data
	tv=np.arange(1970,2016,1)
	
	# Import annual total nitrogen deposition
	ndep=np.zeros((len(tv),len(x)))
	for i in range(len(tv)):
		if tv[i]>2013:
			continue # No data past 2013
		pthin=meta['Paths']['NACID'] + '\\Grids\\NACID_ndep_tot_ann_abso_comp_hist_v1\\NACID_ndep_tot_ann_abso_comp_hist_v1_' + str(tv[i]) + '.mat'
		z=spio.loadmat(pthin,squeeze_me=True)
		ndep0=z['z'][()][1].astype('float')*z['z'][()][0]
		for j in range(len(x)):
			print(str(i) + ' ' + str(j))
			ndep[i,j]=ndep0[ixy[j,1],ixy[j,0]]
	#plt.plot(tv,np.mean(ndep,axis=1))
	
	# Import mean annual temperature normal
	tmean_ann_n=np.zeros((1,len(x)))
	for i in range(12):
		pthin=meta['Paths']['NACID'] + '\\Grids\\NACID_tmean_mon_norm_1971to2000_si_hist_v1\\NACID_tmean_mon_norm_1971to2000_si_hist_v1_' + str(i+1) + '.mat'
		z=spio.loadmat(pthin,squeeze_me=True)
		idat=np.where(np.asarray(z['z'].dtype.names)=='Data')[0][0]
		iSF=np.where(np.asarray(z['z'].dtype.names)=='ScaleFactor')[0][0]
		tmean0=z['z'][()][idat].astype('float')*z['z'][()][iSF]	
		for j in range(len(x)):
			tmean_ann_n[0,j]=tmean_ann_n[0,j]+tmean0[ixy[j,1],ixy[j,0]]
	tmean_ann_n=np.tile(tmean_ann_n/12,(tv.size,1))
	 
	# Import mean annual precipitation normal
	prcp_ann_n=np.zeros((1,len(x)))
	for i in range(12):
		pthin=meta['Paths']['NACID'] + '\\Grids\\NACID_prcp_mon_norm_1971to2000_si_hist_v1\\NACID_prcp_mon_norm_1971to2000_si_hist_v1_' + str(i+1) + '.mat'
		z=spio.loadmat(pthin,squeeze_me=True)
		idat=np.where(np.asarray(z['z'].dtype.names)=='Data')[0][0]
		iSF=np.where(np.asarray(z['z'].dtype.names)=='ScaleFactor')[0][0]
		prcp0=z['z'][()][idat].astype('float')*z['z'][()][iSF]	
		for j in range(len(x)):
			prcp_ann_n[0,j]=prcp_ann_n[0,j]+prcp0[ixy[j,1],ixy[j,0]]
	prcp_ann_n=np.tile(prcp_ann_n,(tv.size,1))	   
	
	# Import mean annual temperature actual
	tmean_ann_r=np.zeros((len(tv),len(x)))
	for i in range(len(tv)):	
		pthin=meta['Paths']['NACID'] + '\\Grids\\NACID_tmean_ann_abso_si_hist_v1\\NACID_tmean_ann_abso_si_hist_v1_' + str(tv[i]) + '.mat'
		z=spio.loadmat(pthin,squeeze_me=True)
		idat=np.where(np.asarray(z['z'].dtype.names)=='Data')[0][0]
		iSF=np.where(np.asarray(z['z'].dtype.names)=='ScaleFactor')[0][0]
		tmean0=z['z'][()][idat].astype('float')*z['z'][()][iSF]
		for j in range(len(x)):
			print(str(i) + ' ' + str(j))
			tmean_ann_r[i,j]=tmean0[ixy[j,1],ixy[j,0]]
			
	# Import annual precipitation actual
	prcp_ann_r=np.zeros((len(tv),len(x)))
	for i in range(len(tv)):	
		pthin=meta['Paths']['NACID'] + '\\Grids\\NACID_prcp_ann_abso_si_hist_v1\\NACID_prcp_ann_abso_si_hist_v1_' + str(tv[i]) + '.mat'
		z=spio.loadmat(pthin,squeeze_me=True)
		idat=np.where(np.asarray(z['z'].dtype.names)=='Data')[0][0]
		iSF=np.where(np.asarray(z['z'].dtype.names)=='ScaleFactor')[0][0]
		prcp0=z['z'][()][idat].astype('float')*z['z'][()][iSF]
		for j in range(len(x)):
			print(str(i) + ' ' + str(j))
			prcp_ann_r[i,j]=prcp0[ixy[j,1],ixy[j,0]]
	
	# Import warm-season mean temperature
	tmean_ws_r=np.zeros((len(tv),len(x)))
	for i in range(len(tv)):	
		pthin=meta['Paths']['NACID'] + '\\Grids\\NACID_tmean_gs_abso_si_hist_v1\\NACID_tmean_gs_abso_si_hist_v1_' + str(tv[i]) + '.mat'
		z=spio.loadmat(pthin,squeeze_me=True)
		idat=np.where(np.asarray(z['z'].dtype.names)=='Data')[0][0]
		iSF=np.where(np.asarray(z['z'].dtype.names)=='ScaleFactor')[0][0]
		tmean0=z['z'][()][idat].astype('float')*z['z'][()][iSF]
		for j in range(len(x)):
			print(str(i) + ' ' + str(j))
			tmean_ws_r[i,j]=tmean0[ixy[j,1],ixy[j,0]]
	
	# Import warm-season mean soil water content
	ws_ws_r=np.zeros((len(tv),len(x)))
	for i in range(len(tv)):	
		pthin=meta['Paths']['NACID'] + '\\Grids\\NACID_ws_tmw_gs_abso_comp_hist_v1\\NACID_ws_tmw_gs_abso_comp_hist_v1_' + str(tv[i]) + '.mat'
		z=spio.loadmat(pthin,squeeze_me=True)
		idat=np.where(np.asarray(z['z'].dtype.names)=='Data')[0][0]
		iSF=np.where(np.asarray(z['z'].dtype.names)=='ScaleFactor')[0][0]
		ws0=z['z'][()][idat].astype('float')*z['z'][()][iSF]
		for j in range(len(x)):
			print(str(i) + ' ' + str(j))
			ws_ws_r[i,j]=ws0[ixy[j,1],ixy[j,0]]
			
	# Import warm-season mean climatic water deficit
	cwd_ws_r=np.zeros((len(tv),len(x)))
	for i in range(len(tv)):	
		pthin=meta['Paths']['NACID'] + '\\Grids\\NACID_cwd_tmw_gs_abso_comp_hist_v1\\NACID_cwd_tmw_gs_abso_comp_hist_v1_' + str(tv[i]) + '.mat'
		z=spio.loadmat(pthin,squeeze_me=True)
		idat=np.where(np.asarray(z['z'].dtype.names)=='Data')[0][0]
		iSF=np.where(np.asarray(z['z'].dtype.names)=='ScaleFactor')[0][0]
		cwd0=z['z'][()][idat].astype('float')*z['z'][()][iSF]
		for j in range(len(x)):
			print(str(i) + ' ' + str(j))
			cwd_ws_r[i,j]=cwd0[ixy[j,1],ixy[j,0]] 
	
	# Import warm-season mean climatic water deficit
	etp_ws_r=np.zeros((len(tv),len(x)))
	for i in range(len(tv)):	
		pthin=meta['Paths']['NACID'] + '\\Grids\\NACID_etp_tmw_gs_abso_comp_hist_v1\\NACID_etp_tmw_gs_abso_comp_hist_v1_' + str(tv[i]) + '.mat'
		z=spio.loadmat(pthin,squeeze_me=True)
		idat=np.where(np.asarray(z['z'].dtype.names)=='Data')[0][0]
		iSF=np.where(np.asarray(z['z'].dtype.names)=='ScaleFactor')[0][0]
		etp0=z['z'][()][idat].astype('float')*z['z'][()][iSF]
		for j in range(len(x)):
			print(str(i) + ' ' + str(j))
			etp_ws_r[i,j]=etp0[ixy[j,1],ixy[j,0]] 
	
	# Add to SOBS structure
	sobs['cwd_mjjas_r']=np.zeros(sobs['ID_Site'].size)
	sobs['etp_mjjas_r']=np.zeros(sobs['ID_Site'].size)
	sobs['ndep_ann_r']=np.zeros(sobs['ID_Site'].size)
	sobs['prcp_ann_r']=np.zeros(sobs['ID_Site'].size)
	sobs['prcp_ann_n']=np.zeros(sobs['ID_Site'].size)
	sobs['tmean_ann_r']=np.zeros(sobs['ID_Site'].size)
	sobs['tmean_ann_n']=np.zeros(sobs['ID_Site'].size)
	sobs['tmean_mjjas_r']=np.zeros(sobs['ID_Site'].size)
	sobs['ws_mjjas_r']=np.zeros(sobs['ID_Site'].size)
	for i in range(len(x)):
		indSite=np.where( (sobs['Lat']==uLL[i,0]) & (sobs['Lon']==uLL[i,1]) )[0]
		for j in range(len(indSite)):
			it=np.where( (tv>=sobs['Year_t0'][indSite[j]]) & (tv<=sobs['Year_t1'][indSite[j]]) )[0]
			sobs['ndep_ann_r'][indSite[j]]=np.mean(ndep[it,i])
			sobs['cwd_mjjas_r'][indSite[j]]=np.mean(cwd_ws_r[it,i])
			sobs['etp_mjjas_r'][indSite[j]]=np.mean(etp_ws_r[it,i])
			sobs['prcp_ann_n'][indSite[j]]=np.mean(prcp_ann_n[it,i])
			sobs['prcp_ann_r'][indSite[j]]=np.mean(prcp_ann_r[it,i])
			sobs['tmean_ann_n'][indSite[j]]=np.mean(tmean_ann_n[it,i])
			sobs['tmean_ann_r'][indSite[j]]=np.mean(tmean_ann_r[it,i])
			sobs['tmean_mjjas_r'][indSite[j]]=np.mean(tmean_ws_r[it,i])
			sobs['ws_mjjas_r'][indSite[j]]=np.mean(ws_ws_r[it,i])

	# Save to file
	gu.opickle(meta['Paths']['Project']['Outputs'] + '\\EP703_SL.pkl',sobs)
	gu.opickle(meta['Paths']['Project']['Outputs'] + '\\EP703_TL.pkl',tobs)

	return

#%%
def GapFillHeights(meta,dInst):
	uPlotS=np.unique(dInst['plot'])
	for iPlotS in range(len(uPlotS)):
		
		# Index to plot
		indP=np.where((dInst['plot']==uPlotS[iPlotS]))[0]
		idP_all=dInst['tree'][indP]
		hP_all=dInst['ht'][indP]
		baP_all=dInst['ba'][indP]
		
		# Global fit
		indP_Fit=np.where((np.isnan(hP_all+baP_all)==False))[0]
		x=baP_all[indP_Fit].astype('float')
		y=hP_all[indP_Fit].astype('float')
		try:
			poptG,pcovG=curve_fit(funcH,x,y,[5,0.25])
		except:
			continue

		flg=0
		if flg==1:
			yhat=funcH(x,poptG[0],poptG[1])
			N_Fit=indP_Fit.size
			b_Perf=np.polyfit(y,yhat,1)
			r_Perf=np.corrcoef(y,yhat)[0,1]
			plt.close('all')
			fig,ax=plt.subplots(1,2)
			ax[0].plot(x,y,'.')
			ax[0].plot(x,yhat,'.')
			ax[1].plot(y,yhat,'.')
		
			#x2=np.c_[np.ones(x.size),np.log(x)]
			#md=sm.OLS(np.log(y),x2).fit()
			#md.summary()
		
		indGF=np.where( (dInst['plot']==uPlotS[iPlotS]) & (np.isnan(dInst['ba'])==False) )[0]
		dInst['ht_mod'][indGF]=funcH(dInst['ba'][indGF],poptG[0],poptG[1])
		
		# Individual tree fit
		uTree=np.unique(idP_all[indP_Fit])
		for k in range(uTree.size):
			indCal=np.where( (idP_all==uTree[k]) & (np.isnan(hP_all+baP_all)==False) )[0]
			if indCal.size>2:
				x=baP_all[indCal].astype('float')
				y=hP_all[indCal].astype('float')
				try:
					j=1
					poptI,pcovI=curve_fit(funcH,x,y,[5,0.25])
					yhat=funcH(x,poptI[0],poptI[1])
					ax[0].plot(x,yhat,'-')
					indGF=np.where( (dInst['plot']==uPlotS[j]) & (dInst['tree']==uTree[k]) & (np.isnan(dInst['ba'])==False) )[0]
					dInst['ht_mod'][indGF]=funcH(dInst['ba'][indGF],poptI[0],poptI[1])
					dInst['H_Fit_ITM'][indGF]=1
				except:
					continue

	# Create gap-filled variable
	ind=np.where( (np.isnan(dInst['ht'])==True) & (np.isnan(dInst['ht_mod'])==False) )[0]
	dInst['ht_gf'][ind]=dInst['ht_mod'][ind]

	return dInst

#%%
def ReviseFAIBSummaryByPlot():
	d=gu.ReadExcel(r'C:\Data\EP703\Data\Raw Data\Plot Summary From FAIB\703 plot vs trtmt by inst.xlsx')
	u=np.unique(np.column_stack((d['install'],d['plot'])),axis=0)
	
	z=np.zeros((u.shape[0],4))
	for i in range(u.shape[0]):
		ind=np.where( (d['install']==u[i,0]) & (d['plot']==u[i,1]) )[0]
		z[i,0]=u[i,0]
		z[i,1]=u[i,1]
		if d['PCode'][ind[0]]=='T0F0':
			z[i,2]=0; z[i,3]=0
		elif d['PCode'][ind[0]]=='T0F1':
			z[i,2]=0; z[i,3]=225
		elif d['PCode'][ind[0]]=='T0F1':
			z[i,2]=0; z[i,3]=225
		elif d['PCode'][ind[0]]=='T0F2':
			z[i,2]=0; z[i,3]=450
		elif d['PCode'][ind[0]]=='T0F3':
			z[i,2]=0; z[i,3]=675
		elif d['PCode'][ind[0]]=='T0F4':
			z[i,2]=0; z[i,3]=900
		elif d['PCode'][ind[0]]=='T1F0':
			z[i,2]=20; z[i,3]=0
		elif d['PCode'][ind[0]]=='T1F1':
			z[i,2]=20; z[i,3]=225
		elif d['PCode'][ind[0]]=='T1F2':
			z[i,2]=20; z[i,3]=450
		elif d['PCode'][ind[0]]=='T1F3':
			z[i,2]=20; z[i,3]=675
		elif d['PCode'][ind[0]]=='T1F4':
			z[i,2]=20; z[i,3]=900
		elif d['PCode'][ind[0]]=='T2F0':
			z[i,2]=35; z[i,3]=0
		elif d['PCode'][ind[0]]=='T2F1':
			z[i,2]=35; z[i,3]=225
		elif d['PCode'][ind[0]]=='T2F2':
			z[i,2]=35; z[i,3]=450
		elif d['PCode'][ind[0]]=='T2F3':
			z[i,2]=35; z[i,3]=675
		elif d['PCode'][ind[0]]=='T2F4':
			z[i,2]=35; z[i,3]=900
		elif d['PCode'][ind[0]]=='T3F0':
			z[i,2]=50; z[i,3]=0
		elif d['PCode'][ind[0]]=='T3F1':
			z[i,2]=50; z[i,3]=225
		elif d['PCode'][ind[0]]=='T3F2':
			z[i,2]=50; z[i,3]=450
		elif d['PCode'][ind[0]]=='T3F3':
			z[i,2]=50; z[i,3]=675
		elif d['PCode'][ind[0]]=='T3F4':
			z[i,2]=50; z[i,3]=900
		elif d['PCode'][ind[0]]=='T0F5':
			z[i,2]=0; z[i,3]=d['Namout'][ind[0]]
		elif d['PCode'][ind[0]]=='T0F6':
			z[i,2]=0; z[i,3]=d['Namout'][ind[0]]
	df=pd.DataFrame(z)
	df.to_excel(r'C:\Data\EP703\Data\Summ.xlsx')
	return

#%%
def ImportCleanedTreeDataFromFAIB(dInst):
	# 		Nm=0
	# 		for i in range(dInst['dbh'].size):
	# 			ind=np.where( (dGFH['install']==uSite[iS]) & (dGFH['plot']==dInst['plot'][i]) & (dGFH['treeno']==dInst['tree'][i]) & (dGFH['meas']-1==dInst['meas'][i]) )[0]
	# 			if ind.size==0:
	# 				Nm=Nm+1
	# 				print('missing')
	# 				break
	# 				continue
	# 			dInst['dbh'][i]=dGFH['dbh'][ind]
	# 			dInst['ht'][i]=dGFH['ht'][ind]
	# 			dInst['ht_gf'][i]=dGFH['ht'][ind]
	return