
#%% Import modules

import numpy as np
from scipy.io import loadmat

#%% Read tree level data by species

def Read_L3_TL_BySpc(spp):
    
    #spp=['AT']
    
    pthin=r'D:\Data\ForestInventory\PSP-NADB\Data\03_IntegratedDatabases\TreeLevel\BySpecies'
    
    for iS in range(len(spp)):
        
        z=loadmat(pthin + '\\FI_PSP_L3_TL_BySpc_' + spp[iS] + '.mat')
    
        d={}
        for iV in range(z['z']['Name'][0].size):
            d[z['z']['Name'][0][iV][0]]=z['z']['Data'][0][iV].flatten().astype(float)*z['z']['ScaleFactor'][0][iV][0][0]
    
    
    return


function qtl=PSP_UTL_read_L3_TL_BySpc(spp);


for s=1:length(spp)
  
  if ischar(spp(s))==1
  
    eval(['load E:\Data\ForestInventory\PSP-NADB\Data\03_IntegratedDatabases\TreeLevel\BySpecies\FI_PSP_L3_TL_BySpc_' spp '.mat']);
    
    for i=1:length(z)
      tmp=double(z(i).Data).*z(i).ScaleFactor;
      eval(['qtl.' z(i).Name '=tmp;']);
    end
    
  else
    
    eval(['load E:\Data\ForestInventory\PSP-NADB\Data\03_IntegratedDatabases\TreeLevel\BySpecies\FI_PSP_L3_TL_BySpc_' spp{s} '.mat']);
    
    if s==1
      
      fn=[]; for i=1:length(z); fn{i}=z(i).Name; end
      
      for i=1:length(z)
        
        tmp=double(z(i).Data).*z(i).ScaleFactor;
        eval(['qtl.' z(i).Name '=tmp;']);
        
      end
      
    else
      
      for i=1:length(fn)
        
        tmp=double(z(i).Data).*z(i).ScaleFactor;
        
        eval(['qtl.' z(i).Name '=[qtl.' z(i).Name '; tmp];']);
        
      end
      
    end
    
  end
  
end

return