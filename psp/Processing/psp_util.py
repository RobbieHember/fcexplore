
#%% Import modules

import numpy as np
from scipy.io import loadmat
from fcgadgets.macgyver import util_general as gu

#%% Import data

def ImportPlotData(meta,**kwargs):

    # Define paths
    if 'Paths' not in meta:
        meta['Paths']={}
    meta['Paths']['GP']={}
    meta['Paths']['GP']['DB']=r'C:\Users\rhember\Documents\Data\GroundPlots\PSP-NADB2'
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

    return meta,data

#%% Import LUTs

def ImportParameters(meta):

    if 'Param' not in meta.keys():
        meta['Param']={}
    if 'LUT' not in meta.keys():
        meta['LUT']={}

    # Parameters
    meta['Param']['GP']={}

    meta['Param']['GP']['Carbon Content']=0.5

    # Allometry
    meta['Param']['GP']['Allo B']=gu.ReadExcel(meta['Paths']['GP']['DB'] + '\\Parameters\\Parameters Allometry Biomass.xlsx')
    meta['Param']['GP']['Allo V Tot Nigh16']=gu.ReadExcel(meta['Paths']['GP']['DB'] + '\\Parameters\\Parameters Allometry Volume Total BC Nigh 2016.xlsx')

    # Look up tables
    meta['LUT']['GP']={}

    # Raw tables
    meta['LUT']['GP']['Raw Tables']={}
    tabL=['Source','Jurisdiction','Plot Type','Damage Agents','DA BC','Vital Status','Management','Insect Codes','Pathogen Codes','PFT',
          'Crown Class','Ecozone BC L1','Ecozone BC L2','Ecozone CA L1']
    for tab in tabL:
        meta['LUT']['GP']['Raw Tables'][tab]=gu.ReadExcel(meta['Paths']['GP']['DB'] + '\\Parameters\\Parameters LUTs.xlsx',tab)

    # Species kyes
    meta['LUT']['GP']['Species Given']={}
    meta['LUT']['GP']['Species Given']['BC']=gu.ReadExcel(meta['Paths']['GP']['DB'] + '\\Parameters\\Parameters Species Codes.xlsx','Key_BC')

    meta['LUT']['GP']['Species']={}
    for i in range(meta['Param']['GP']['Allo B']['ID'].size):
            meta['LUT']['GP']['Species'][ meta['Param']['GP']['Allo B']['Code'][i] ]=meta['Param']['GP']['Allo B']['ID'][i]

    meta['LUT']['GP']['Source']={}
    for i in range(meta['LUT']['GP']['Raw Tables']['Source']['ID'].size):
        meta['LUT']['GP']['Source'][ meta['LUT']['GP']['Raw Tables']['Source']['Value'][i] ]=meta['LUT']['GP']['Raw Tables']['Source']['ID'][i]

    meta['LUT']['GP']['Plot Type BC']={}
    for i in range(meta['LUT']['GP']['Raw Tables']['Plot Type']['Given'].size):
        if meta['LUT']['GP']['Raw Tables']['Plot Type']['Jurisdiction'][i]=='BC':
            meta['LUT']['GP']['Plot Type BC'][ meta['LUT']['GP']['Raw Tables']['Plot Type']['Given'][i] ]=meta['LUT']['GP']['Raw Tables']['Plot Type']['Value'][i]

    meta['LUT']['GP']['Vital Status']={}
    for i in range(meta['LUT']['GP']['Raw Tables']['Vital Status']['ID'].size):
        meta['LUT']['GP']['Vital Status'][ meta['LUT']['GP']['Raw Tables']['Vital Status']['Value'][i] ]=meta['LUT']['GP']['Raw Tables']['Vital Status']['ID'][i]

    meta['LUT']['GP']['Damage Agents']={}
    for i in range(meta['LUT']['GP']['Raw Tables']['Damage Agents']['ID'].size):
        if meta['LUT']['GP']['Raw Tables']['Damage Agents']['Value'][i]=='nan':
            continue
        meta['LUT']['GP']['Damage Agents'][ meta['LUT']['GP']['Raw Tables']['Damage Agents']['Value'][i] ]=meta['LUT']['GP']['Raw Tables']['Damage Agents']['ID'][i]

    meta['LUT']['GP']['Ecozone BC L1']={}
    for i in range(meta['LUT']['GP']['Raw Tables']['Ecozone BC L1']['ID'].size):
        meta['LUT']['GP']['Ecozone BC L1'][ meta['LUT']['GP']['Raw Tables']['Ecozone BC L1']['Value'][i] ]=meta['LUT']['GP']['Raw Tables']['Ecozone BC L1']['ID'][i]

    meta['LUT']['GP']['Ecozone BC L2']={}
    for i in range(meta['LUT']['GP']['Raw Tables']['Ecozone BC L2']['ID'].size):
        meta['LUT']['GP']['Ecozone BC L2'][ meta['LUT']['GP']['Raw Tables']['Ecozone BC L2']['Value'][i] ]=meta['LUT']['GP']['Raw Tables']['Ecozone BC L2']['ID'][i]

    meta['LUT']['GP']['Ecozone CA L1']={}
    for i in range(meta['LUT']['GP']['Raw Tables']['Ecozone CA L1']['ID'].size):
        meta['LUT']['GP']['Ecozone CA L1'][ meta['LUT']['GP']['Raw Tables']['Ecozone CA L1']['Value'][i] ]=meta['LUT']['GP']['Raw Tables']['Ecozone CA L1']['ID'][i]

    meta['LUT']['GP']['Insect Codes']={}
    for i in range(meta['LUT']['GP']['Raw Tables']['Insect Codes']['ID'].size):
        meta['LUT']['GP']['Insect Codes'][ meta['LUT']['GP']['Raw Tables']['Insect Codes']['Value'][i] ]=meta['LUT']['GP']['Raw Tables']['Insect Codes']['ID'][i]

    meta['LUT']['GP']['Pathogen Codes']={}
    for i in range(meta['LUT']['GP']['Raw Tables']['Pathogen Codes']['ID'].size):
        meta['LUT']['GP']['Pathogen Codes'][ meta['LUT']['GP']['Raw Tables']['Pathogen Codes']['Value'][i] ]=meta['LUT']['GP']['Raw Tables']['Pathogen Codes']['ID'][i]

    meta['LUT']['GP']['Management']={}
    for i in range(meta['LUT']['GP']['Raw Tables']['Management']['ID'].size):
        meta['LUT']['GP']['Management'][ meta['LUT']['GP']['Raw Tables']['Management']['Value'][i] ]=meta['LUT']['GP']['Raw Tables']['Management']['ID'][i]

    meta['LUT']['GP']['Crown Class']={}
    for i in range(meta['LUT']['GP']['Raw Tables']['Crown Class']['ID'].size):
        meta['LUT']['GP']['Crown Class'][ meta['LUT']['GP']['Raw Tables']['Crown Class']['Value'][i] ]=meta['LUT']['GP']['Raw Tables']['Crown Class']['ID'][i]

    meta['LUT']['GP']['PFT']={}
    for i in range(meta['LUT']['GP']['Raw Tables']['PFT']['ID'].size):
        meta['LUT']['GP']['PFT'][ meta['LUT']['GP']['Raw Tables']['PFT']['Value'][i] ]=meta['LUT']['GP']['Raw Tables']['PFT']['ID'][i]

    meta['LUT']['GP']['Stature']={}
    meta['LUT']['GP']['Stature']['Standing']=1
    meta['LUT']['GP']['Stature']['Fallen']=2

    meta['LUT']['GP']['ClimateClassCondensed']={}
    labL=['Hyper humid','Humid','Subhumid','Semi arid','Arid','Outside Boundary']
    for iLab in range(len(labL)):
        meta['LUT']['GP']['ClimateClassCondensed'][labL[iLab]]=iLab+1

    # List of jurisdictions
    meta['LUT']['GP']['List of Jurisidictions']={}
    meta['LUT']['GP']['List of Jurisidictions']['CA']=['BC','AB','SK','MB','ON','QC','NL','NS','NB']
    meta['LUT']['GP']['List of Jurisidictions']['US']=['AK','AL','AR','AZ','CA','CO','CT','DE','FL','GA','IA','ID','IL',
        'IN','KS','KY','LA','MA','MD','ME','MI','MN','MO','MS','MT','NC','ND',
        'NE','NH','NJ','NM','NV','NY','OK','OH','OR','PA','RI','SC','SD',
        'TN','TX','UT','VA','VT','WA','WI','WV','WY']

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