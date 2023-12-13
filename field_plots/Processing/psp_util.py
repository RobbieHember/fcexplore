#%% Import modules
import numpy as np
from scipy.io import loadmat
from fcgadgets.macgyver import util_general as gu
import fcgadgets.bc1ha.bc1ha_util as u1ha

#%% Import data
def ImportGroundPlotData(meta,**kwargs):

    # Define paths
    if 'Paths' not in meta:
        meta['Paths']={}
    meta['Paths']['GP']={}
    meta['Paths']['GP']['DB']=r'C:\Data\GroundPlots\PSP-NADB2'
    meta['Paths']['GP']['Figures']=r'C:\Users\rhember\OneDrive - Government of BC\Figures\Ground Plots'

    # Import LUTs and parameters
    meta=ImportParameters(meta)

    if kwargs['type']=='Stand':
        data=gu.ipickle(meta['Paths']['GP']['DB'] + '\\Processed\\L2\\L2_BC.pkl')['sobs']
    elif kwargs['type']=='Tree':
        data=gu.ipickle(meta['Paths']['GP']['DB'] + '\\Processed\\L2\\L2_BC.pkl')['tobs']
    elif kwargs['type']=='Just Parameters':
        data=[]
    else:
        data=[]
        print('Type not recognized')

    #if kwargs['include_soil']=='True':
    soc=gu.ipickle(r'C:\Data\Soils\Shaw et al 2018 Database\SITES.pkl')

    return meta,data,soc

#%% Import LUTs
def ImportParameters(meta):

    if 'Param' not in meta.keys():
        meta['Param']={}
    if 'LUT' not in meta.keys():
        meta['LUT']={}

    # Parameters
    meta['Param']['GP']={}
    meta['Param']['GP']['Raw Tables']={}
    tabL=['Source','Territory','General','Species','Species Crosswalk','Biomass Allometry','Plot Type','Damage Agent','DA Crosswalk','Vital Status','Crown Class','Crown Class Crosswalk','Stature','Management','PFT','Ecozone BC L1','Ecozone BC L2','Ecozone CA L1']
    for tab in tabL:
        meta['Param']['GP']['Raw Tables'][tab]=gu.ReadExcel(meta['Paths']['GP']['DB'] + '\\Parameters\\PSP-NADB2 Parameters.xlsx',tab)
    
    for i in range(meta['Param']['GP']['Raw Tables']['General']['Name'].size):
        meta['Param']['GP'][ meta['Param']['GP']['Raw Tables']['General']['Name'][i] ]=meta['Param']['GP']['Raw Tables']['General']['Value'][i]

    # Look up tables
    meta['LUT']['GP']={}
    vL=['Source','Territory','Species','Plot Type','Vital Status','Damage Agent','Crown Class','Stature','PFT','Ecozone BC L1','Ecozone BC L2','Ecozone CA L1']
    for v in vL:
        meta['LUT']['GP'][v]={}
        for i in range(meta['Param']['GP']['Raw Tables'][v]['ID'].size):
                meta['LUT']['GP'][v][ meta['Param']['GP']['Raw Tables'][v]['Code'][i] ]=meta['Param']['GP']['Raw Tables'][v]['ID'][i]
    
    # # List of jurisdictions
    # meta['LUT']['GP']['List of Jurisidictions']={}
    # meta['LUT']['GP']['List of Jurisidictions']['CA']=['BC','AB','SK','MB','ON','QC','NL','NS','NB']
    # meta['LUT']['GP']['List of Jurisidictions']['US']=['AK','AL','AR','AZ','CA','CO','CT','DE','FL','GA','IA','ID','IL',
    #     'IN','KS','KY','LA','MA','MD','ME','MI','MN','MO','MS','MT','NC','ND',
    #     'NE','NH','NJ','NM','NV','NY','OK','OH','OR','PA','RI','SC','SD',
    #     'TN','TX','UT','VA','VT','WA','WI','WV','WY']

    return meta

#%% Get code from LUT ID
def lut_id2cd(lut,id):
    for k in lut.keys():
        if lut[k]==id:
            cd=k
    return cd

#%% Read tree level data by species
def Read_L3_TL_BySpc(spp):

    #spp=['SW']

    pthin=r'E:\Data\ForestInventory\PSP-NADB\Data\03_IntegratedDatabases\TreeLevel\BySpecies'

    for iS in range(len(spp)):

        z=loadmat(pthin + '\\FI_PSP_L3_TL_BySpc_' + spp[iS] + '.mat')

        d={}
        for iV in range(z['z']['Name'][0].size):
            d[z['z']['Name'][0][iV][0]]=z['z']['Data'][0][iV].flatten().astype(float)*z['z']['ScaleFactor'][0][iV][0][0]

    #--------------------------------------------------------------------------
    # Fix stand age variable in BC
    #--------------------------------------------------------------------------

#    ind0=np.where(d['ID_DB']==1)[0]
#
#    ID_Plot=d['ID_Plot'][ind0]
#    ID_Tree=d['ID_Tree'][ind0]
#    Age_t0_Stand=d['Age_t0_Stand'][ind0]
#    Age_t1_Stand=d['Age_t1_Stand'][ind0]
#    DT=d['DT'][ind0]
#
#    u=np.unique( np.column_stack((ID_Plot,ID_Tree)),axis=1 )
#
#    for i in range(u.shape[0]):
#        ind=np.where( (ID_Plot==u[i,0]) & (ID_Tree==u[i,1]) )[0]
#        if ind.size==0:
#            continue
#        if ind.shape[0]==1:
#            continue
#        Age_t0_Stand[ind]=Age_t0_Stand[ind[0]]+np.append(0,np.cumsum(DT[ind[0:-1]]))
#        Age_t1_Stand[ind]=Age_t0_Stand[ind[0]]+np.cumsum(DT[ind])
#
#    d['Age_t0_Stand'][ind0]=Age_t0_Stand
#    d['Age_t1_Stand'][ind0]=Age_t1_Stand

    return d

#%%
def CalcStatsByBGCZone(meta,gpt,soc):
    vL=['Ctot L t0','Ctot D t0','Cbk L t0','Cbr L t0','Cf L t0','Cr L t0','Csw L t0',
        'Ctot G Surv','Ctot G Recr','Ctot Mort Nat','Ctot Mort Harv','Ctot Net',
        'Vws G Surv','Vws G Recr','Vws Mort Nat','Vws Mort Harv','Vws Net','Area']
    
    uZ=np.unique(gpt['Ecozone BC L1'])
    uZ=uZ[uZ>0]
    lab=np.array(['' for _ in range(uZ.size)],dtype=object)

    d={'CN':{},'CNV':{},'YSM':{},'VRI':{}}
    for k in d.keys():
        for v in vL:
            d[k][v]={'N':np.zeros(uZ.size),'mu':np.zeros(uZ.size),'sd':np.zeros(uZ.size),'se':np.zeros(uZ.size),'sum':np.zeros(uZ.size)}
        for iU in range(uZ.size):
            lab[iU]=lut_id2cd(meta['LUT']['GP']['Ecozone BC L1'],uZ[iU])
            for v in vL:
                if np.isin(v,['Ctot L t0','Cbk L t0','Cbr L t0','Cf L t0','Cr L t0','Csw L t0'])==True:
                    ind=np.where( (gpt['Ecozone BC L1']==uZ[iU]) & (gpt['PTF ' + k]==1) &
                             (gpt['Cbk L t0']>=0) & (gpt['Cbk L t0']<2000) &
                             (gpt['Cbr L t0']>=0) & (gpt['Cbr L t0']<2000) &
                             (gpt['Cf L t0']>=0) & (gpt['Cf L t0']<2000) &
                             (gpt['Cr L t0']>=0) & (gpt['Cr L t0']<2000) &
                             (gpt['Csw L t0']>=0) & (gpt['Csw L t0']<2000) &
                             (gpt['Ctot L t0']>=0) & (gpt['Ctot L t0']<2000))[0]
                elif v=='Ctot D t0':
                    ind=np.where( (gpt['Ecozone BC L1']==uZ[iU]) & (gpt['PTF ' + k]==1) & (gpt['Ctot D t0']>=0) )[0]
                else:
                    ind=np.where( (gpt['Ecozone BC L1']==uZ[iU]) & (gpt['PTF ' + k]==1) & (gpt['Ctot Net']>=-1000) & (gpt['Ctot Net']<1000) )[0]
                d[k][v]['N'][iU]=ind.size
                d[k][v]['mu'][iU]=np.nanmean(gpt[v][ind])
                d[k][v]['sd'][iU]=np.nanstd(gpt[v][ind])
                d[k][v]['se'][iU]=np.nanstd(gpt[v][ind])/np.sqrt(ind.size)
            ind2=np.where(meta['Param']['BE']['BGC Zone Averages']['Name']==lab[iU])[0]
            d[k]['Area']['sum'][iU]=1e6*meta['Param']['BE']['BGC Zone Averages']['Area Treed (Mha)'][ind2]
    
    d['CN']['Ctot L t0']['mu']
    
    # Add soils
    vL_Soil=['TOT_C_THA','MIN_C_THA','ORG_C_THA']
    d['Soil']={}
    for v in vL_Soil:
        d['Soil'][v]={'N':np.zeros(uZ.size),'mu':np.zeros(uZ.size),'sd':np.zeros(uZ.size),'se':np.zeros(uZ.size),'sum':np.zeros(uZ.size)}
    #d['Soil']['Area']={'sum':np.zeros(uZ.size)}
    for iU in range(uZ.size):
        cd=lut_id2cd(meta['LUT']['GP']['Ecozone BC L1'],uZ[iU])
        id=meta['LUT']['BEC_BIOGEOCLIMATIC_POLY']['ZONE'][cd]
        for v in vL_Soil:
            ind=np.where( (soc['bgcz']==id) )[0]
            d['Soil'][v]['N'][iU]=ind.size
            d['Soil'][v]['mu'][iU]=np.nanmean(soc[v][ind])
            d['Soil'][v]['sd'][iU]=np.nanstd(soc[v][ind])
            d['Soil'][v]['se'][iU]=np.nanstd(soc[v][ind])/np.sqrt(ind.size)
        #ind2=np.where(meta['Param']['BE']['BGC Zone Averages']['Name']==lab[iU])[0]
        #d['Soil']['Area']['sum'][iU]=1e6*meta['Param']['BE']['BGC Zone Averages']['Area Treed (Mha)'][ind2]
    
    for k in d.keys():
        for v in d[k].keys():
            if v=='Area':
                continue
            d[k][v]['sum']=d[k][v]['mu']*d['CN']['Area']['sum']/1e9
    
    # Add Total ecosystem carbon
    d['Soil']['TEC']={'N':np.zeros(uZ.size),'mu':np.zeros(uZ.size),'sd':np.zeros(uZ.size),'se':np.zeros(uZ.size),'sum':np.zeros(uZ.size)}
    d['Soil']['TEC']['mu']=d['CN']['Ctot D t0']['mu']+d['CN']['Cbk L t0']['mu']+d['CN']['Cbr L t0']['mu']+d['CN']['Cf L t0']['mu']+d['CN']['Cr L t0']['mu']+d['CN']['Csw L t0']['mu']+d['Soil']['TOT_C_THA']['mu']
    d['Soil']['TEC']['sum']=(d['CN']['Ctot D t0']['mu']+d['CN']['Cbk L t0']['mu']+d['CN']['Cbr L t0']['mu']+d['CN']['Cf L t0']['mu']+d['CN']['Cr L t0']['mu']+d['CN']['Csw L t0']['mu']+d['Soil']['TOT_C_THA']['mu'])*d['CN']['Area']['sum']/1e9
    d1={'id':uZ,'code':lab,'data':d}
    return d1

#%%
def CalcStatsByAgeClass(meta,gpt):
	bw=25; bin=np.arange(bw,250+bw,bw)
	ptfL=['PTF CN','PTF CNY','PTF YSM','PTF VRI']
	vL=['Ctot L t0','Ctot G Surv','Ctot G Recr','Ctot Mort','Ctot Net','Ctot Mort Nat','Ctot Mort Harv']
	d={}
	for ptf in ptfL:
		d[ptf]={'Coast':{},'Interior':{}}
		for v in vL:
			reg='Coast'
			ind=np.where( (gpt[ptf]==1) & (gpt['Ecozone BC L1']==meta['LUT']['GP']['Ecozone BC L1']['CWH']) | (gpt[ptf]==1) & (gpt['Ecozone BC L1']==meta['LUT']['GP']['Ecozone BC L1']['MH']) )[0]
			x=gpt['Age Med t0'][ind]
			y=gpt[v][ind]
			N,mu,med,sig,se=gu.discres(x,y,bw,bin)
			d[ptf][reg][v]={'N':N,'mu':mu,'se':se}
	
			reg='Interior'
			ind=np.where( (gpt[ptf]==1) & (gpt['Ecozone BC L1']!=meta['LUT']['GP']['Ecozone BC L1']['CWH']) & (gpt['Ecozone BC L1']!=meta['LUT']['GP']['Ecozone BC L1']['MH']) & (gpt['Ecozone BC L1']!=meta['LUT']['GP']['Ecozone BC L1']['ICH']) )[0]
			x=gpt['Age Med t0'][ind]
			y=gpt[v][ind]
			N,mu,med,sig,se=gu.discres(x,y,bw,bin)
			d[ptf][reg][v]={'N':N,'mu':mu,'se':se}
	d1={'bin':bin,'bw':bw,'data':d}
	return d1

