
#%% Import modules


tva=np.arange(1920,2022,1)

vr='Gsw'

y0=np.zeros(tva.size)
for iT in range(tva.size):
    ind=np.where( (dTL['Year']==tva[iT]) & (dTL['Treatment']!='F1') & (dTL[vr]>0) & (dTL[vr]<5000) & (dTL['QA Summary']==1) )[0]
    y0[iT]=np.mean(dTL[vr][ind])

plt.close('all')
plt.plot(tva,y0,'-bo')

y1=np.zeros(tva.size)
for iT in range(tva.size):
    ind=np.where( (dTL['Year']==tva[iT]) & (dTL['Treatment']!='F1') & (dTL[vr]>0) & (dTL[vr]<5000) & (dTL['QA Summary']==1) )[0]
    y1[iT]=np.mean(dTL[vr][ind])
    
plt.plot(tva,y1,'-gs')    


plt.close('all')
plt.plot(dTL['Year'],dTL['Age'],'ko')

#%% Time response (by installations)

tva=np.arange(1920,2022,1)

datT=[]
for iI in range(uInst.size):
    ListSS=[]
    for iSS in range(uSS.size): 
        d0={}
        for k in dGF.keys():            
            d0[k]=np.zeros(tva.size)
            for iT in range(tva.size):
                try:
                    ind=np.where( (dGF['ID_Inst']==uInst[iI]) & (dGF['Year']==tva[iT]) & (dGF['BGC SS']==uSS[iSS]) & (dGF['Treatment']!='F1') & (dGF[k]>0) & (dGF[k]<5000) & (dGF['QA Summary']==1) )[0]
                except:
                    continue
                if ind.size>0:
                    try:
                        d0[k][iT]=np.nanmean(dGF[k][ind])
                    except:
                        pass
        ListSS.append(copy.deepcopy(d0))
    datT.append(copy.deepcopy(ListSS))

#%%
 
vn='Age'
iI=0
plt.close('all'); ms=3
plt.plot(tva,datT[iI][0][vn],'-ko',ms=ms)
    
vn='Gsw'
iI=7
plt.close('all'); ms=3
plt.plot(tva,datT[iI][0][vn],'-ko',ms=ms)
plt.plot(tva,datT[iI][1][vn],'-rs',ms=ms)
plt.plot(tva,datT[iI][2][vn],'-bs',ms=ms)
plt.plot(tva,datT[iI][3][vn],'-gs',ms=ms)

plt.figure()
plt.plot(tva,datT[iI][0]['tmean_gs_r'],'-ko',ms=ms)

#%%

iI=6
plt.close('all'); ms=3
plt.plot(dat[iI][0]['ws_gs_r'],dat[iI][0][vn],'ko',ms=ms)



#%% Age response (by installations)

Age=np.arange(1,100,1)

datA=[]
for iI in range(uInst.size):
    ListSS=[]
    for iSS in range(uSS.size): 
        d0={}
        for k in dGF.keys():            
            d0[k]=np.zeros(Age.size)
            for iT in range(Age.size):
                try:
                    ind=np.where( (dGF['ID_Inst']==uInst[iI]) & (dGF['Age']==Age[iT]) & (dGF['BGC SS']==uSS[iSS]) & (dGF['Treatment']!='F1') & (dGF[k]>0) & (dGF[k]<5000) & (dGF['QA Summary']==1) )[0]
                except:
                    continue
                if ind.size>0:
                    try:
                        d0[k][iT]=np.nanmean(dGF[k][ind])
                    except:
                        pass
        ListSS.append(copy.deepcopy(d0))
    datA.append(copy.deepcopy(ListSS))

#%%
    
vn='Dib'
iI=7
plt.close('all'); ms=3
plt.plot(Age,datA[iI][0][vn],'-ko',ms=ms)
plt.plot(Age,datA[iI][1][vn],'-rs',ms=ms)
plt.plot(Age,datA[iI][2][vn],'-bs',ms=ms)
plt.plot(Age,datA[iI][3][vn],'-gs',ms=ms)   



#%% Time response (by tree)

vd='Dib'

ind0=np.where( (dTL['Treatment']=='C') & (dTL[vd]>0) & (dTL[vd]<5000) & (dTL['QA Summary']==1) )[0]

u=np.unique(dTL['ID_Tree_Unique'][ind0])

ind=np.where( (dTL['ID_Tree_Unique']==u[33]) & (dTL[vd]>0) & (dTL[vd]<5000) & (dTL['QA Summary']==1) )[0]
Inst=dTL['ID_Inst'][ind[0]]
print(Inst)
print(dTL['Pith_status'][ind[0]])
print(dTL['Est/mm'][ind[0]]/10)
iInst=np.where( (dTL['ID_Inst']==Inst) & (dTL[vd]>0) & (dTL[vd]<5000) & (dTL['QA Summary']==1) )[0]

plt.close('all')
plt.plot(dTL['Year'][iInst],dTL[vd][iInst],'k.',mec=[0.88,0.88,0.88],mfc=[0.88,0.88,0.88])
plt.plot(dTL['Year'][ind],dTL[vd][ind],'-ko')
plt.plot(2020,dTL['Dob 2020'][ind[0]],'gs')


#%% Mixed model

dTL['SS']=dTL['BGC SS']
dTL['Age']=dTL['Age']
dTL['Age2']=dTL['Age']**2
dTL['LnGsw']=np.log(dTL['Gsw'])
dTL['LnBsw']=np.log(dTL['Bsw'])

d0=copy.deepcopy(dTL)
ikp=np.where( (d0['Year']>1970) & (d0['QA Summary']==1) & \
             (d0['BAI']>0) & (d0['BAI']<5000) & \
             (d0['LnGsw']>-2) & (d0['LnGsw']<20) & (d0['Age']>0) & (d0['Age']<200) & \
             (d0['ws_gs_r']>=0) & (d0['ws_gs_r']<=200) )[0]
for k in d0.keys():
    d0[k]=d0[k][ikp]

data=pd.DataFrame(d0)

mstr='LnGsw ~ LnBsw + Age + cwd_gs_r + C(SS) + C(SS):cwd_gs_r'

md=smf.mixedlm(mstr,data,groups=data['ID_Tree'],re_formula='~LnBsw')
mdf=md.fit(method=["lbfgs"])

print(mdf.summary())


