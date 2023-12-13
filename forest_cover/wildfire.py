
    
#%% QA

#meta['LUT']['RSLT_FOREST_COVER_INV_SVW']['STOCKING_TYPE_CODE'].keys()
#Out[35]: dict_keys(['ART', 'BR', 'FOR', 'NAT', 'NPL', 'PL', 'RD', 'RHR', 'UNN'])

#meta['LUT']['RSLT_FOREST_COVER_INV_SVW']['STOCKING_STATUS_CODE'].keys()
#Out[36]: dict_keys(['A', 'AF', 'C', 'G', 'IMM', 'L', 'M', 'MAT', 'NC', 'NF', 'NP', 'NSR', 'OR', 'R', 'RES', 'S', 'U'])

z=u1ha.Import_Raster(meta,[],['refg','fc_yr','fc_stc','fc_ssc','fire_yr'])

z=u1ha.Import_Raster(meta,[],['fc_yr','fire_yr'])
for k in z.keys(): z[k]=z[k]['Data']

yd=z['fc_yr']-z['fire_yr']

ind=np.where( (z['fire_yr']>0) & (yd>-100) )

plt.hist(yd[ind][0::20],np.arange(-30,30))


    