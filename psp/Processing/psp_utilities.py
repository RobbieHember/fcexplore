
#%% Import modules

import numpy as np
from scipy.io import loadmat
from fcgadgets.macgyver import utilities_general as gu

#%% Import data

def ImportPSPs(**kwargs):
    meta={}
    meta['Paths']={}
    meta['Paths']['DB']=r'C:\Users\rhember\Documents\Data\GroundPlots\PSP-NADB2'
    meta['Paths']['Figs']=r'C:\Users\rhember\OneDrive - Government of BC\Figures\Ground Plots'
    meta=ImportParameters(meta)
    d=gu.ipickle(meta['Paths']['DB'] + '\\Processed\\L2\\L2_BC.pkl')
    #d=gu.ipickle(meta['Paths']['DB'] + '\\Processed\\L2\\L2_BC_WithLID.pkl')
    if kwargs['type']=='Stand':
        data=d['sobs']
    elif kwargs['type']=='Tree':
        data=d['tobs']
    elif kwargs['type']=='Just Parameters':
        data=[]
    else:
        data=[]
        print('Type not recognized')

    return meta,data

#%% Import LUTs

def ImportParameters(meta):

    # Allometry
    meta['Allo B']=gu.ReadExcel(meta['Paths']['DB'] + '\\Parameters\\Parameters Allometry Biomass.xlsx')
    meta['Allo V Tot Nigh16']=gu.ReadExcel(meta['Paths']['DB'] + '\\Parameters\\Parameters Allometry Volume Total BC Nigh 2016.xlsx')

    # Species kyes
    meta['Species']={}
    meta['Species']['BC']=gu.ReadExcel(meta['Paths']['DB'] + '\\Parameters\\Parameters Species Codes.xlsx','Key_BC')

    meta['LUT Tables']={}
    tabL=['Source','Jurisdiction','Plot Type','Damage Agents','DA BC','Vital Status','Management','Insect Codes','Pathogen Codes','PFT',
          'Crown Class','Ecozone BC L1','Ecozone BC L2','Ecozone CA L1']
    for tab in tabL:
        meta['LUT Tables'][tab]=gu.ReadExcel(meta['Paths']['DB'] + '\\Parameters\\Parameters LUTs.xlsx',tab)

    meta['LUT']={}

    meta['LUT']['Source']={}
    for i in range(meta['LUT Tables']['Source']['ID'].size):
        meta['LUT']['Source'][ meta['LUT Tables']['Source']['Value'][i] ]=meta['LUT Tables']['Source']['ID'][i]

    meta['LUT']['Species']={}
    for i in range(meta['Allo B']['ID'].size):
            meta['LUT']['Species'][ meta['Allo B']['Code'][i] ]=meta['Allo B']['ID'][i]

    meta['LUT']['Plot Type BC']={}
    for i in range(meta['LUT Tables']['Plot Type']['Given'].size):
        if meta['LUT Tables']['Plot Type']['Jurisdiction'][i]=='BC':
            meta['LUT']['Plot Type BC'][ meta['LUT Tables']['Plot Type']['Given'][i] ]=meta['LUT Tables']['Plot Type']['Value'][i]

    meta['LUT']['Vital Status']={}
    for i in range(meta['LUT Tables']['Vital Status']['ID'].size):
        meta['LUT']['Vital Status'][ meta['LUT Tables']['Vital Status']['Value'][i] ]=meta['LUT Tables']['Vital Status']['ID'][i]

    meta['LUT']['Damage Agents']={}
    for i in range(meta['LUT Tables']['Damage Agents']['ID'].size):
        meta['LUT']['Damage Agents'][ meta['LUT Tables']['Damage Agents']['Value'][i] ]=meta['LUT Tables']['Damage Agents']['ID'][i]

    meta['LUT']['Ecozone BC L1']={}
    for i in range(meta['LUT Tables']['Ecozone BC L1']['ID'].size):
        meta['LUT']['Ecozone BC L1'][ meta['LUT Tables']['Ecozone BC L1']['Value'][i] ]=meta['LUT Tables']['Ecozone BC L1']['ID'][i]

    meta['LUT']['Ecozone BC L2']={}
    for i in range(meta['LUT Tables']['Ecozone BC L2']['ID'].size):
        meta['LUT']['Ecozone BC L2'][ meta['LUT Tables']['Ecozone BC L2']['Value'][i] ]=meta['LUT Tables']['Ecozone BC L2']['ID'][i]

    meta['LUT']['Ecozone CA L1']={}
    for i in range(meta['LUT Tables']['Ecozone CA L1']['ID'].size):
        meta['LUT']['Ecozone CA L1'][ meta['LUT Tables']['Ecozone CA L1']['Value'][i] ]=meta['LUT Tables']['Ecozone CA L1']['ID'][i]

    meta['LUT']['Insect Codes']={}
    for i in range(meta['LUT Tables']['Insect Codes']['ID'].size):
        meta['LUT']['Insect Codes'][ meta['LUT Tables']['Insect Codes']['Value'][i] ]=meta['LUT Tables']['Insect Codes']['ID'][i]

    meta['LUT']['Pathogen Codes']={}
    for i in range(meta['LUT Tables']['Pathogen Codes']['ID'].size):
        meta['LUT']['Pathogen Codes'][ meta['LUT Tables']['Pathogen Codes']['Value'][i] ]=meta['LUT Tables']['Pathogen Codes']['ID'][i]

    meta['LUT']['Management']={}
    for i in range(meta['LUT Tables']['Management']['ID'].size):
        meta['LUT']['Management'][ meta['LUT Tables']['Management']['Value'][i] ]=meta['LUT Tables']['Management']['ID'][i]

    meta['LUT']['Crown Class']={}
    for i in range(meta['LUT Tables']['Crown Class']['ID'].size):
        meta['LUT']['Crown Class'][ meta['LUT Tables']['Crown Class']['Value'][i] ]=meta['LUT Tables']['Crown Class']['ID'][i]

    meta['LUT']['PFT']={}
    for i in range(meta['LUT Tables']['PFT']['ID'].size):
        meta['LUT']['PFT'][ meta['LUT Tables']['PFT']['Value'][i] ]=meta['LUT Tables']['PFT']['ID'][i]

    meta['LUT']['Stature']={}
    meta['LUT']['Stature']['Standing']=1
    meta['LUT']['Stature']['Fallen']=2

    meta['LUT']['ClimateClassCondensed']={}
    labL=['Hyper humid','Humid','Subhumid','Semi arid','Arid','Outside Boundary']
    for iLab in range(len(labL)):
        meta['LUT']['ClimateClassCondensed'][labL[iLab]]=iLab+1

    # List of jurisdictions
    meta['List of Jurisidictions']={}
    meta['List of Jurisidictions']['CA']=['BC','AB','SK','MB','ON','QC','NL','NS','NB']
    meta['List of Jurisidictions']['US']=['AK','AL','AR','AZ','CA','CO','CT','DE','FL','GA','IA','ID','IL',
        'IN','KS','KY','LA','MA','MD','ME','MI','MN','MO','MS','MT','NC','ND',
        'NE','NH','NJ','NM','NV','NY','OK','OH','OR','PA','RI','SC','SD',
        'TN','TX','UT','VA','VT','WA','WI','WV','WY']

    return meta

#%% Get code from LUT ID

def lut_id2cd(meta,nam_lut,id):
    for k in meta['LUT'][nam_lut].keys():
        if meta['LUT'][nam_lut][k]==id:
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