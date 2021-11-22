'''
EP703 RETROSPECTIVE MONITORING STUDY - ANALYSIS
'''

#%% Import data

# Run import script

flag_savefigs='On'



#%% Calculate responses following Ballard and Majid (1985) and Brockley (2015)

# Custom functions
def GetCIsForRR(lo0,hi0,lo1,hi1):
    tmp=np.array([lo1/lo0,lo1/hi0,hi1/lo0,hi1/hi0]).T
    lo=np.min(tmp,axis=0)
    hi=np.max(tmp,axis=0)    
    return lo,hi

def GetCIsForDR(lo0,hi0,lo1,hi1):
    a=(lo1-lo0)/lo0*100
    b=(lo1-hi0)/hi0*100 
    c=(hi1-lo0)/lo0*100 
    d=(hi1-hi0)/hi0*100
    tmp=np.array([ a,b,c,d ]).T
    lo=np.min(tmp,axis=1)
    hi=np.max(tmp,axis=1)    
    return lo,hi

# Define reference periods before and after N application
t_ref0=[1970-4,1970]
t_ref1=[1971,1971+9]

# Initialize
vList=['TRW','BAI','Gsw','RGR']
tList=['C','F1','SS3','SS5','SS6']
d1={}
for vn in vList:
    d1[vn]={}
    for tn in tList:
        d1[vn][tn]={}
        d1[vn][tn]['Mean']=np.zeros(uInst.size)
        d1[vn][tn]['Median']=np.zeros(uInst.size)
        #d1[vn][tn]['SD']=np.zeros(uInst.size)
        d1[vn][tn]['CI Lo']=np.zeros(uInst.size)
        d1[vn][tn]['CI Hi']=np.zeros(uInst.size)

# Populate
for iI in range(uInst.size):
    for vn in vList:
        
        # Indices
        indC0=np.where( (dTL['ID_Inst']==uInst[iI]) & (dTL['Year']>=t_ref0[0]) & (dTL['Year']<=t_ref0[1]) & (dTL['Treatment']=='C') & (dTL[vn]>-1) & (dTL[vn]<5000) & (dTL['QA Summary']==1) )[0]
        indC1=np.where( (dTL['ID_Inst']==uInst[iI]) & (dTL['Year']>=t_ref1[0]) & (dTL['Year']<=t_ref1[1]) & (dTL['Treatment']=='C') & (dTL[vn]>-1) & (dTL[vn]<5000) & (dTL['QA Summary']==1) )[0]
        indF0=np.where( (dTL['ID_Inst']==uInst[iI]) & (dTL['Year']>=t_ref0[0]) & (dTL['Year']<=t_ref0[1]) & (dTL['Treatment']=='F1') & (dTL[vn]>-1) & (dTL[vn]<5000) & (dTL['QA Summary']==1) )[0]
        indF1=np.where( (dTL['ID_Inst']==uInst[iI]) & (dTL['Year']>=t_ref1[0]) & (dTL['Year']<=t_ref1[1]) & (dTL['Treatment']=='F1') & (dTL[vn]>-1) & (dTL[vn]<5000) & (dTL['QA Summary']==1) )[0]
        indSS30=np.where( (dTL['ID_Inst']==uInst[iI]) & (dTL['Year']>=t_ref0[0]) & (dTL['Year']<=t_ref0[1]) & (dTL['BGC SS']==3) & (dTL[vn]>-1) & (dTL[vn]<5000) & (dTL['Treatment']!='C') & (dTL['Treatment']!='F1') & (dTL['QA Summary']==1) )[0]
        indSS31=np.where( (dTL['ID_Inst']==uInst[iI]) & (dTL['Year']>=t_ref1[0]) & (dTL['Year']<=t_ref1[1]) & (dTL['BGC SS']==3) & (dTL[vn]>-1) & (dTL[vn]<5000) & (dTL['Treatment']!='C') & (dTL['Treatment']!='F1') & (dTL['QA Summary']==1) )[0]
        indSS50=np.where( (dTL['ID_Inst']==uInst[iI]) & (dTL['Year']>=t_ref0[0]) & (dTL['Year']<=t_ref0[1]) & (dTL['BGC SS']==5) & (dTL[vn]>-1) & (dTL[vn]<5000) & (dTL['Treatment']!='C') & (dTL['Treatment']!='F1') & (dTL['QA Summary']==1) )[0]
        indSS51=np.where( (dTL['ID_Inst']==uInst[iI]) & (dTL['Year']>=t_ref1[0]) & (dTL['Year']<=t_ref1[1]) & (dTL['BGC SS']==5) & (dTL[vn]>-1) & (dTL[vn]<5000) & (dTL['Treatment']!='C') & (dTL['Treatment']!='F1') & (dTL['QA Summary']==1) )[0]
        indSS60=np.where( (dTL['ID_Inst']==uInst[iI]) & (dTL['Year']>=t_ref0[0]) & (dTL['Year']<=t_ref0[1]) & (dTL['BGC SS']==6) & (dTL[vn]>-1) & (dTL[vn]<5000) & (dTL['Treatment']!='C') & (dTL['Treatment']!='F1') & (dTL['QA Summary']==1) )[0]
        indSS61=np.where( (dTL['ID_Inst']==uInst[iI]) & (dTL['Year']>=t_ref1[0]) & (dTL['Year']<=t_ref1[1]) & (dTL['BGC SS']==6) & (dTL[vn]>-1) & (dTL[vn]<5000) & (dTL['Treatment']!='C') & (dTL['Treatment']!='F1') & (dTL['QA Summary']==1) )[0]
    
        # Median 
        yC0=np.nanmedian(dTL[vn][indC0]); yC1=np.nanmedian(dTL[vn][indC1]); yF0=np.nanmedian(dTL[vn][indF0]); yF1=np.nanmedian(dTL[vn][indF1]); ySS30=np.nanmedian(dTL[vn][indSS30]); ySS31=np.nanmedian(dTL[vn][indSS31]); ySS50=np.nanmedian(dTL[vn][indSS50]); ySS51=np.nanmedian(dTL[vn][indSS51]); ySS60=np.nanmedian(dTL[vn][indSS60]); ySS61=np.nanmedian(dTL[vn][indSS61]);
        d1[vn]['C']['Median'][iI]=yC1/yC0
        d1[vn]['F1']['Median'][iI]=yF1/yF0
        d1[vn]['SS3']['Median'][iI]=ySS31/ySS30
        d1[vn]['SS5']['Median'][iI]=ySS51/ySS50
        d1[vn]['SS6']['Median'][iI]=ySS61/ySS60
        
        # Mean
        muC0=np.nanmean(dTL[vn][indC0]); muC1=np.nanmean(dTL[vn][indC1]); muF0=np.nanmean(dTL[vn][indF0]); muF1=np.nanmean(dTL[vn][indF1]); muSS30=np.nanmean(dTL[vn][indSS30]); muSS31=np.nanmean(dTL[vn][indSS31]); muSS50=np.nanmean(dTL[vn][indSS50]); muSS51=np.nanmean(dTL[vn][indSS51]); muSS60=np.nanmean(dTL[vn][indSS60]); muSS61=np.nanmean(dTL[vn][indSS61]);
        d1[vn]['C']['Mean'][iI]=muC1/muC0
        d1[vn]['F1']['Mean'][iI]=muF1/muF0
        d1[vn]['SS3']['Mean'][iI]=muSS31/muSS30
        d1[vn]['SS5']['Mean'][iI]=muSS51/muSS50
        d1[vn]['SS6']['Mean'][iI]=muSS61/muSS60
        
        # Error (1 x S.E.)
        multip=1
        eC0=multip*np.nanstd(dTL[vn][indC0])/np.sqrt(indC0.size); 
        eC1=multip*np.nanstd(dTL[vn][indC1])/np.sqrt(indC1.size);         
        loC0=muC0-eC0
        hiC0=muC0+eC0
        loC1=muC1-eC1
        hiC1=muC1+eC1
        d1[vn]['C']['CI Lo'][iI],d1[vn]['C']['CI Hi'][iI]=GetCIsForRR(loC0,hiC0,loC1,hiC1)
        
        eF0=multip*np.nanstd(dTL[vn][indF0])/np.sqrt(indF0.size); 
        eF1=multip*np.nanstd(dTL[vn][indF1])/np.sqrt(indF1.size); 
        loF0=muF0-eF0
        hiF0=muF0+eF0
        loF1=muF1-eF1
        hiF1=muF1+eF1
        d1[vn]['F1']['CI Lo'][iI],d1[vn]['F1']['CI Hi'][iI]=GetCIsForRR(loF0,hiF0,loF1,hiF1)
        
        eSS30=multip*np.nanstd(dTL[vn][indSS30])/np.sqrt(indSS30.size); 
        eSS31=multip*np.nanstd(dTL[vn][indSS31])/np.sqrt(indSS31.size); 
        loSS30=muSS30-eSS30
        hiSS30=muSS30+eSS30
        loSS31=muSS31-eSS31
        hiSS31=muSS31+eSS31
        d1[vn]['SS3']['CI Lo'][iI],d1[vn]['SS3']['CI Hi'][iI]=GetCIsForRR(loSS30,hiSS30,loSS31,hiSS31)
        
        eSS50=multip*np.nanstd(dTL[vn][indSS50])/np.sqrt(indSS50.size); 
        eSS51=multip*np.nanstd(dTL[vn][indSS51])/np.sqrt(indSS51.size); 
        loSS50=muSS50-eSS50
        hiSS50=muSS50+eSS50
        loSS51=muSS51-eSS51
        hiSS51=muSS51+eSS51
        d1[vn]['SS5']['CI Lo'][iI],d1[vn]['SS5']['CI Hi'][iI]=GetCIsForRR(loSS50,hiSS50,loSS51,hiSS51)
        
        eSS60=multip*np.nanstd(dTL[vn][indSS60])/np.sqrt(indSS60.size); 
        eSS61=multip*np.nanstd(dTL[vn][indSS61])/np.sqrt(indSS61.size); 
        loSS60=muSS60-eSS60
        hiSS60=muSS60+eSS60
        loSS61=muSS61-eSS61
        hiSS61=muSS61+eSS61
        d1[vn]['SS6']['CI Lo'][iI],d1[vn]['SS6']['CI Hi'][iI]=GetCIsForRR(loSS60,hiSS60,loSS61,hiSS61)
    
       
#%% Plot change in C and F1 plots (between pre- and post-application reference periods)

vn='BAI'
stat='Mean'

cl1=[0.8,0.8,0.8]; cl2=[0.55,0.55,0.55]

plt.close('all'); 
fig,ax=plt.subplots(3,1,figsize=gu.cm2inch(8,10)); 

#------------------------------------------------------------------------------
# Both C and F
#------------------------------------------------------------------------------

ax[0].bar(np.arange(uInst.size)-0.15,d1[vn]['C'][stat],0.3,fc=cl1,label='Control');
ax[0].bar(np.arange(uInst.size)+0.15,d1[vn]['F1'][stat],0.3,fc=cl2,label='Fertilized');
ax[0].errorbar(np.arange(uInst.size)-0.15,d1[vn]['C'][stat],yerr=[ d1[vn]['C'][stat]-d1[vn]['C']['CI Lo'], d1[vn]['C']['CI Hi']-d1[vn]['C'][stat] ],color='k',ls='',capsize=2)
ax[0].errorbar(np.arange(uInst.size)+0.15,d1[vn]['F1'][stat],yerr=[ d1[vn]['F1'][stat]-d1[vn]['F1']['CI Lo'], d1[vn]['F1']['CI Hi']-d1[vn]['C'][stat] ],color='k',ls='',capsize=2)
#ax[0].yaxis.set_ticks_position('both'); ax[0].xaxis.set_ticks_position('both')
ax[0].set(position=[0.16,0.7,0.8,0.29],xlim=[-0.5,uInst.size-0.5],xticks=np.arange(uInst.size),xticklabels='',ylabel='BAI (cm yr$^-$$^1$)')

#------------------------------------------------------------------------------
# Absolute differences
#------------------------------------------------------------------------------

da_be=d1[vn]['F1'][stat]-d1[vn]['C'][stat]
da_lo,da_hi=gu.GetCIsFromDifference(d1[vn]['C']['CI Lo'],d1[vn]['C']['CI Hi'],d1[vn]['F1']['CI Lo'],d1[vn]['F1']['CI Hi'])

da_be_mu=np.mean(da_be)
da_be_se=1*np.nanstd(da_be)/np.sqrt(da_be_mu.size)

da_be=np.append(da_be,da_be_mu)
da_lo=np.append(da_lo,da_be_mu-da_be_se)
da_hi=np.append(da_hi,da_be_mu+da_be_se)

#------------------------------------------------------------------------------
# Relative differences
#------------------------------------------------------------------------------

#plt.close('all'); 
#fig,ax=plt.subplots(1,figsize=gu.cm2inch(14.5,5)); 
ax[1].plot([-1,20],[0,0],'k-')
ax[1].bar(np.arange(uInst.size+1),da_be,0.5,fc=cl2)
ax[1].errorbar(np.arange(uInst.size+1),da_be,yerr=[ da_be-da_lo, da_hi-da_be ],color='k',ls='',capsize=2)
#ax[1].yaxis.set_ticks_position('both'); ax[1].xaxis.set_ticks_position('both')
ax[1].set(position=[0.16,0.39,0.8,0.29],xlim=[-0.5,uInst.size+0.5],xticks=np.arange(uInst.size+1),xticklabels='',ylabel='$\Delta$ BAI (cm yr$^-$$^1$)')

dr_be=(d1[vn]['F1'][stat]-d1[vn]['C'][stat])/d1[vn]['C'][stat]*100
dr_lo,dr_hi=GetCIsForDR(d1[vn]['C']['CI Lo'],d1[vn]['C']['CI Hi'],d1[vn]['F1']['CI Lo'],d1[vn]['F1']['CI Hi'])
dr_be_mu=np.mean(dr_be)
dr_be_se=1*np.nanstd(dr_be)/np.sqrt(dr_be_mu.size)
dr_be=np.append(dr_be,dr_be_mu)
dr_lo=np.append(dr_lo,dr_be_mu-dr_be_se)
dr_hi=np.append(dr_hi,dr_be_mu+dr_be_se)

#plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(14.5,5)); 
ax[2].plot([-1,20],[0,0],'k-')
ax[2].bar(np.arange(uInst.size+1),dr_be,0.5,fc=cl2)
ax[2].errorbar(np.arange(uInst.size+1),dr_be,yerr=[ dr_be-dr_lo, dr_hi-dr_be ],color='k',ls='',capsize=2)
#ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')
ax[2].set(position=[0.16,0.08,0.8,0.29],xlim=[-0.5,uInst.size+0.5],xticks=np.arange(uInst.size+1),xticklabels=np.append(uInst,'All'),ylabel='$\Delta$ BAI (%)',xlabel='Installation ID')
gu.axletters(ax,plt,0.015,0.89,LetterStyle='Caps',FontWeight='Bold')
if flag_savefigs=='On':
    gu.PrintFig(meta['Paths']['Figures'] + '\\Response_' + vn,'png',900)
    fig.savefig(meta['Paths']['Figures'] + '\\Response_' + vn + '.svg',format='svg',dpi=1200)


#%% Evaluate pseudocontrols

vn='Gsw'

dr_be=(d1[vn]['F1'][stat]-d1[vn]['C'][stat])/d1[vn]['C'][stat]*100
dr_lo,dr_hi=GetCIsForDR(d1[vn]['C']['CI Lo'],d1[vn]['C']['CI Hi'],d1[vn]['F1']['CI Lo'],d1[vn]['F1']['CI Hi'])
dr_be_mu=np.mean(dr_be)
dr_be_se=1*np.nanstd(dr_be)/np.sqrt(dr_be_mu.size)

dr_be=(d1[vn]['F1'][stat]-d1[vn]['SS3'][stat])/d1[vn]['SS3'][stat]*100
dr_lo,dr_hi=GetCIsForDR(d1[vn]['SS3']['CI Lo'],d1[vn]['SS3']['CI Hi'],d1[vn]['F1']['CI Lo'],d1[vn]['F1']['CI Hi'])
dr_be_mu=np.append(dr_be_mu,np.nanmean(dr_be))
dr_be_se=np.append(dr_be_se,1*np.nanstd(dr_be)/np.sqrt(dr_be_mu.size))

dr_be=(d1[vn]['F1'][stat]-d1[vn]['SS5'][stat])/d1[vn]['SS5'][stat]*100
dr_lo,dr_hi=GetCIsForDR(d1[vn]['SS5']['CI Lo'],d1[vn]['SS5']['CI Hi'],d1[vn]['F1']['CI Lo'],d1[vn]['F1']['CI Hi'])
dr_be_mu=np.append(dr_be_mu,np.nanmean(dr_be))
dr_be_se=np.append(dr_be_se,1*np.nanstd(dr_be)/np.sqrt(dr_be_mu.size))

plt.close('all'); 
fig,ax=plt.subplots(1,figsize=gu.cm2inch(12,8)); 
ax.bar([1,2,3],dr_be_mu,0.8)

#%%
ax.errorbar(bin,be,yerr=[ be-lo, hi-be ],color=[0.8,0.8,0.8],ls='',capsize=0.5,elinewidth=2)
ax.annotate('Application',xy=(1971,1),xytext=(1971,35),
    arrowprops={'color':'red','arrowstyle':'->'},ha='center',color='r');
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')
ax.set(position=[0.09,0.11,0.89,0.86],ylim=[-17,38],yticks=np.arange(-15,40,5),xlim=[bin[0]-0.5,bin[-1]+0.5],xticks=np.arange(1970,2025,5),ylabel='Standardized BAI',xlabel='Time, years')

if flag_savefigs=='On':
    gu.PrintFig(meta['Paths']['Figures'] + '\\Pseudocontrols test','png',900)
    fig.savefig(meta['Paths']['Figures'] + '\\Pseudocontrols test.svg',format='svg',dpi=1200)


  
#%% Time response

#vd='TRW Standardized RR'
vd='BAI Standardized RR'
#vd='Gsw Standardized RR'
#vd='RGR Standardized RR'

bin=np.arange(1964,2022,1)

muC=np.zeros(bin.size)
muF=np.zeros(bin.size)
muSS3=np.zeros(bin.size)
muSS5=np.zeros(bin.size)
seC=np.zeros(bin.size)
seF=np.zeros(bin.size)
seSS3=np.zeros(bin.size)
seSS5=np.zeros(bin.size)

for i in range(bin.size):
    ind=np.where( (dTL['Year']==bin[i]) & (dTL['Treatment']=='C') & (dTL[vd]>0) & (dTL[vd]<2000) & (dTL['QA Summary']==1) )[0]
    muC[i]=np.nanmean(dTL[vd][ind])
    seC[i]=np.nanstd(dTL[vd][ind])/np.sqrt(ind.size)
    ind=np.where( (dTL['Year']==bin[i]) & (dTL['Treatment']=='F1') & (dTL[vd]>0) & (dTL[vd]<2000) & (dTL['QA Summary']==1) )[0]
    muF[i]=np.nanmean(dTL[vd][ind])    
    seF[i]=np.nanstd(dTL[vd][ind])/np.sqrt(ind.size)
    ind=np.where( (dTL['Year']==bin[i]) & (dTL['BGC SS']==3) & (dTL[vd]>0) & (dTL[vd]<2000) & (dTL['Treatment']!='C') & (dTL['Treatment']!='F1') & (dTL['QA Summary']==1) )[0]
    muSS3[i]=np.nanmean(dTL[vd][ind])
    ind=np.where( (dTL['Year']==bin[i]) & (dTL['BGC SS']==5) & (dTL[vd]>0) & (dTL[vd]<2000) & (dTL['Treatment']!='C') & (dTL['Treatment']!='F1') & (dTL['QA Summary']==1) )[0]
    muSS5[i]=np.nanmean(dTL[vd][ind])

be=(muF-muC)/muC*100
lo,hi=GetCIsForDR(muC-seC,muC+seC,muF-seF,muF+seF)

plt.close('all'); 
fig,ax=plt.subplots(1,figsize=gu.cm2inch(12,8)); 
ax.plot(bin,np.zeros(bin.size),'k-',lw=1,color='k')
ax.plot(bin,be,'-ko',ms=3.5,mfc='w',mec='k',lw=1,color='k')
ax.errorbar(bin,be,yerr=[ be-lo, hi-be ],color=[0.8,0.8,0.8],ls='',capsize=0.5,elinewidth=2)
ax.annotate('Application',xy=(1971,1),xytext=(1971,35),
    arrowprops={'color':'red','arrowstyle':'->'},ha='center',color='r');
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')
ax.set(position=[0.09,0.11,0.89,0.86],ylim=[-17,38],yticks=np.arange(-15,40,5),xlim=[bin[0]-0.5,bin[-1]+0.5],xticks=np.arange(1970,2025,5),ylabel='Standardized BAI',xlabel='Time, years')

if flag_savefigs=='On':
    gu.PrintFig(meta['Paths']['Figures'] + '\\' + vd + '_vs_Time','png',900)
    fig.savefig(meta['Paths']['Figures'] + '\\' + vd + '_vs_Time.svg',format='svg',dpi=1200)

#%% Plot time series of standardized responses using pesudo controls

#vd='TRW Standardized RR'
vd='BAI Standardized RR'
#vd='Gsw Standardized RR'
#vd='RGR Standardized RR'

plt.close('all'); 
fig,ax=plt.subplots(1,figsize=gu.cm2inch(12,8)); 
ax.plot(bin,np.zeros(bin.size),'k-',lw=1,color='k')
be=(muF-muC)/muC*100
#ax.plot(bin,be,'-ko',ms=3.5,mfc='w',mec='k',lw=1,color='k')
be=(muF-muSS3)/muSS3*100
ax.plot(bin,be,'-ko',ms=3.5,mfc='w',mec='k',lw=1,color='k')
be=(muF-muSS5)/muSS5*100
#ax.plot(bin,be,'-ko',ms=3.5,mfc='w',mec='k',lw=1,color='k')


##%% Time response
#
##vd='TRW Standardized RR'
#vd='RGR Standardized RR'
#vd='Gsw'
##vd='Gsw Standardized RR'
#
#id=uInst[2]
#
#ind=np.where( (dTL['ID_Inst']==id) )[0]
#print(np.mean(dTL['tmean_ann_n'][ind]))
#
#bin=np.arange(1920,2022,1)
#
#yC_mu=np.zeros(bin.size)
#yF_mu=np.zeros(bin.size)
#ySS3_mu=np.zeros(bin.size)
#ySS5_mu=np.zeros(bin.size)
#ySS6_mu=np.zeros(bin.size)
#
#yC_med=np.zeros(bin.size)
#yF_med=np.zeros(bin.size)
#ySS3_med=np.zeros(bin.size)
#ySS5_med=np.zeros(bin.size)
#ySS6_med=np.zeros(bin.size)
#for i in range(bin.size):
#    ind=np.where( (dTL['ID_Inst']==id) & (dTL['Year']==bin[i]) & (dTL['Treatment']=='C') & (dTL[vd]>0) & (dTL[vd]<2000) & (dTL['QA Summary']==1) )[0]
#    yC_mu[i]=np.nanmean(dTL[vd][ind])
#    yC_med[i]=np.nanmedian(dTL[vd][ind])
#    ind=np.where( (dTL['ID_Inst']==id) & (dTL['Year']==bin[i]) & (dTL['BGC SS']==3) & (dTL[vd]>0) & (dTL[vd]<2000) & (dTL['QA Summary']==1) )[0]
#    ySS3_mu[i]=np.nanmean(dTL[vd][ind])
#    ySS3_med[i]=np.nanmedian(dTL[vd][ind])
#    ind=np.where( (dTL['ID_Inst']==id) & (dTL['Year']==bin[i]) & (dTL['BGC SS']==5) & (dTL[vd]>0) & (dTL[vd]<2000) & (dTL['QA Summary']==1) )[0]
#    ySS5_mu[i]=np.nanmean(dTL[vd][ind])
#    ySS5_med[i]=np.nanmedian(dTL[vd][ind])
#    ind=np.where( (dTL['ID_Inst']==id) & (dTL['Year']==bin[i]) & (dTL['BGC SS']==6) & (dTL[vd]>0) & (dTL[vd]<2000) & (dTL['QA Summary']==1) )[0]
#    ySS6_mu[i]=np.nanmean(dTL[vd][ind])
#    ySS6_med[i]=np.nanmedian(dTL[vd][ind])
#    ind=np.where( (dTL['ID_Inst']==id) & (dTL['Year']==bin[i]) & (dTL['Treatment']=='F1') & (dTL[vd]>0) & (dTL[vd]<2000) & (dTL['QA Summary']==1) )[0]
#    yF_mu[i]=np.nanmean(dTL[vd][ind])
#    yF_med[i]=np.nanmedian(dTL[vd][ind])
#
#plt.close('all'); ms=3
#plt.plot(bin,yC_mu,'-ko',ms=ms)
#plt.plot(bin,yF_mu,'-go',ms=ms)
#plt.plot(bin,ySS3_mu,'-rs',ms=ms)
#plt.plot(bin,ySS5_mu,'-bs',ms=ms)
#plt.plot(bin,ySS6_mu,'-cs',ms=ms)
#
#plt.figure()
#plt.plot(bin,yC_med,'-ko',ms=ms)
#plt.plot(bin,yF_med,'-go',ms=ms)
#plt.plot(bin,ySS3_med,'-rs',ms=ms)
#plt.plot(bin,ySS5_med,'-bs',ms=ms)
#plt.plot(bin,ySS6_med,'-cs',ms=ms)
#
#plt.figure()
#plt.plot(bin,yF_mu/yC_mu,'-bo')
##plt.plot(bin,yF_mu/ySS3_mu,'-ro')
##plt.plot(bin,yF_mu/ySS5_mu,'-go')
#
#
#
#
##%% Time response (by tree)
#
##vd='TRW Standardized RR'
##vd='RGR Standardized RR'
#vd='Gsw'
##vd='Gsw Standardized RR'
#
#ind0=np.where( (dTL['Treatment']=='F1') & (dTL[vd]>0) & (dTL[vd]<5000) & (dTL['QA Summary']==1) )[0]
#
#u=np.unique(dTL['ID_Tree Unique'][ind0])
#
#ind=np.where( (dTL['ID_Tree Unique']==u[7]) & (dTL[vd]>0) & (dTL[vd]<5000) & (dTL['QA Summary']==1) )[0]
#print(dTL['ID_Inst'][ind[0]])
#
#plt.close('all')
#plt.plot(dTL['Year'][ind],dTL[vd][ind],'-ko')
#
#plt.figure()
#plt.plot(dTL['Year'][ind],dTL['ws_gs_r'][ind],'-ko')
#
#
#
#
##%% Age response
#
#vd='TRW Standardized RR'
#vd='Gsw Standardized RR'
#
#bin=np.arange(1,80,1)
#yC=np.zeros(bin.size)
#yCP=np.zeros(bin.size)
#yF=np.zeros(bin.size)
#for i in range(bin.size):
#    ind=np.where( (dTL['Time since first ring']==bin[i]) & (dTL['Treatment']=='C') )[0]
#    yC[i]=np.nanmedian(dTL[vd][ind])
#    #ind=np.where( (dTL['Age']==bin[i]) & (dTL['Treatment']!='C') & (dTL['Treatment']!='F1') )[0]
#    ind=np.where( (dTL['Time since first ring']==bin[i]) & (dTL['BGC SS']==3) )[0]
#    yCP[i]=np.nanmedian(dTL[vd][ind])
#    ind=np.where( (dTL['Time since first ring']==bin[i]) & (dTL['Treatment']=='F1') )[0]
#    yF[i]=np.nanmedian(dTL[vd][ind])
#
#plt.close('all')
#plt.plot(bin,yC,'-bo')
##plt.plot(bin,yCP,'-cd')
#plt.plot(bin,yF,'-gs')
#
##%% Age response
#
#id=uInst[2]
#
#bin=np.arange(1,80,1)
#yC=np.zeros(bin.size)
#yCP=np.zeros(bin.size)
#yF=np.zeros(bin.size)
#for i in range(bin.size):
#    ind=np.where( (dTL['ID_Inst']==id) & (dTL['Time since first ring']==bin[i]) & (dTL['Treatment']=='C') )[0]
#    yC[i]=np.nanmean(dTL['Gsw Standardized RR'][ind])
#    #ind=np.where( (dTL['ID_Inst']==id) & (dTL['Age']==bin[i]) & (dTL['Treatment']!='C') & (dTL['Treatment']!='F1') )[0]
#    ind=np.where( (dTL['ID_Inst']==id) & (dTL['Time since first ring']==bin[i]) & (dTL['BGC SS']==3) )[0]
#    yCP[i]=np.nanmean(dTL['Gsw Standardized RR'][ind])
#    ind=np.where( (dTL['ID_Inst']==id) & (dTL['Time since first ring']==bin[i]) & (dTL['Treatment']=='F1') )[0]
#    yF[i]=np.nanmean(dTL['Gsw Standardized RR'][ind])
#
#plt.close('all')
#plt.plot(bin,yC,'-bo')
##plt.plot(bin,yCP,'-cd')
#plt.plot(bin,yF,'-gs')
#
##%% Biomass response
#
#id=uInst[1]
#
#bin=np.arange(5,150,5)
#yC=np.zeros(bin.size)
#yCP=np.zeros(bin.size)
#yF=np.zeros(bin.size)
#for i in range(bin.size):
#    ind=np.where( (dTL['ID_Inst']==id) & (np.abs(dTL['Bsw']-bin[i])<=2) & (dTL['Treatment']=='C') )[0]
#    yC[i]=np.nanmean(dTL['Gsw'][ind])
#    ind=np.where( (dTL['ID_Inst']==id) & (np.abs(dTL['Bsw']-bin[i])<=2) & (dTL['BGC SS']==3) & (dTL['Treatment']!='C') & (dTL['Treatment']!='F1') )[0]
#    yCP[i]=np.nanmean(dTL['Gsw'][ind])
#    ind=np.where( (dTL['ID_Inst']==id) & (np.abs(dTL['Bsw']-bin[i])<=2) & (dTL['Treatment']=='F1') )[0]
#    yF[i]=np.nanmean(dTL['Gsw'][ind])
#
#plt.close('all')
#plt.plot(bin,yC,'-bo')
#plt.plot(bin,yCP,'-cd')
#plt.plot(bin,yF,'-gs')
#
##%% Comparison of growth trends across site series
#
#vd='Gsw'
#tv=np.arange(1920,2022,1)
#y=np.zeros((tv.size,uSS.size))
#for iSS in range(uSS.size):
#    for iT in range(tv.size):
#        ind=np.where( (dTL['Year First']<1966) & (dTL['ID_Inst']==uInst[0]) & (dTL['Year']==tv[iT]) & (dTL['BGC SS']==uSS[iSS]) & (dTL['Treatment']!='F1') & (dTL[vd]>0) & (dTL[vd]<2000) )[0]
#        y[iT,iSS]=np.nanmedian(dTL[vd][ind])
#    
#plt.close('all')
#plt.plot(tv,y[:,0],'-ko')
#plt.plot(tv,y[:,1],'-ro')
#plt.plot(tv,y[:,2],'-bo')
#plt.plot(tv,y[:,3],'-gs')
#
##%%
#
#vd='BAI Standardized RR'
#vd='Gsw'
#tv=np.arange(1950,2022,1)
#y=np.zeros((tv.size,uSS.size))
#for iSS in range(uSS.size):
#    for iT in range(tv.size):
#        ind=np.where( (dTL['Year First']<1966) & (dTL['Year']==tv[iT]) & (dTL['BGC SS']==uSS[iSS]) & (dTL['Treatment']!='F1') & (dTL[vd]>0) & (dTL[vd]<2000) )[0]
#        y[iT,iSS]=np.nanmedian(dTL[vd][ind])
#    y[:,iSS]=y[:,iSS]-gu.movingave(y[:,iSS],5,'Centre')
#    
#plt.close('all')
#plt.plot(tv,y[:,0],'-ko')
#plt.plot(tv,y[:,1],'-ro')
#plt.plot(tv,y[:,2],'-bo')
#plt.plot(tv,y[:,3],'-gs')
#
#plt.close('all')
#plt.plot(tv,np.zeros(tv.size),lw=3,color=[0.8,0.8,0.8])
#plt.plot(tv,gu.movingave(y[:,0],10,'center'),'-ko')
#
#plt.close('all')
#plt.plot(tv,np.zeros(tv.size),lw=3,color=[0.8,0.8,0.8])
#plt.plot(tv,y[:,0]-gu.movingave(y[:,0],10,'center'),'-ko')
#plt.plot(tv,y[:,1]-gu.movingave(y[:,1],10,'center'),'-ro')
#plt.plot(tv,y[:,2]-gu.movingave(y[:,2],10,'center'),'-bo')
#plt.plot(tv,y[:,3]-gu.movingave(y[:,3],10,'center'),'-gs')
#
##%%
#
#tv=np.arange(1,80,1)
#y=np.zeros((tv.size,uSS.size))
#for iSS in range(uSS.size):
#    for iT in range(tv.size):
#        ind=np.where( (dTL['Year First']<1966) & (dTL['Time since first ring']==tv[iT]) & (dTL['BGC SS']==uSS[iSS]) & (dTL['Treatment']!='F1') )[0]
#        y[iT,iSS]=np.nanmedian(dTL['Gsw Standardized RR'][ind])
#    
#plt.close('all')
#plt.plot(tv,y[:,0],'-ko')
#plt.plot(tv,y[:,1],'-ro')
#plt.plot(tv,y[:,2],'-bo')
#plt.plot(tv,y[:,3],'-gs')
#
#
##%% Comparison between Dob and H in 2020
#
## Isolate data
#iC=np.where( (dTL['Dob 2020']>0) & (dTL['H 2020']>0) & (dTL['Treatment']=='C') )[0]
#iF1=np.where( (dTL['Dob 2020']>0) & (dTL['H 2020']>0) & (dTL['Treatment']=='F1') )[0]
#
## Fit linear best fit relationship
#yC=dTL['H 2020'][iC]
#xC=dTL['Dob 2020'][iC]
#x1C=sm.tools.tools.add_constant(xC)
#mdC=sm.OLS(yC,x1C).fit()
#mdC.summary()
#xhat=np.linspace(np.min(x1C[:,1]),np.max(x1C[:,1]),10)
#yhatC=mdC.predict(np.c_[np.ones(xhat.size),xhat])
#
#yF1=dTL['H 2020'][iF1]
#xF1=dTL['Dob 2020'][iF1]
#x1F1=sm.tools.tools.add_constant(xF1)
#mdF1=sm.OLS(yF1,x1F1).fit()
#mdF1.summary()
#yhatF1=mdF1.predict(np.c_[np.ones(xhat.size),xhat])
#
#
## Plot
#plt.close('all')
#fig,ax=plt.subplots(1,figsize=gu.cm2inch(8.5,8.5))
##ax.plot([0,90],[0,90],'k-',lw=2,color=[0.8,0.8,0.8])
#ax.plot(xC,yC,'o',ms=4,mec='w',mfc='b')
#ax.plot(xF1,yF1,'o',ms=4,mec='w',mfc='r')
#
#ax.plot(xhat,yhatC,'b-',lw=1.25)
#ax.plot(xhat,yhatF1,'r--',lw=1.25)
#ax.set(xlim=[0,80],ylim=[0,80],xlabel='Inside-bark diameter from cores 2019 (cm)',ylabel='Outside-bark diameter from 2020 (cm)')
#ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both')
#plt.tight_layout()