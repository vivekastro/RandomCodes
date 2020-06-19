import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import os
import glob
import scipy as sp
from pydl.pydl.pydlutils import yanny
import pydl as pydl
import pydl.pyspherematch as pshm
from pydl.pydl.pydlutils.spheregroup import *


def init_plotting():
    plt.rcParams['font.size'] = 15
    plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['axes.labelsize'] = 1.5*plt.rcParams['font.size']
    plt.rcParams['axes.titlesize'] = 1.5*plt.rcParams['font.size']
    plt.rcParams['legend.fontsize'] = 1.25*plt.rcParams['font.size']
    plt.rcParams['xtick.labelsize'] = 1.25*plt.rcParams['font.size']
    plt.rcParams['ytick.labelsize'] = 1.25*plt.rcParams['font.size']
    plt.rcParams['xtick.major.size'] = 3
    plt.rcParams['xtick.minor.size'] = 3
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['xtick.minor.width'] = 1
    plt.rcParams['ytick.major.size'] = 3
    plt.rcParams['ytick.minor.size'] = 3
    plt.rcParams['ytick.major.width'] = 1
    plt.rcParams['ytick.minor.width'] = 1
    plt.rcParams['legend.frameon'] = False
    plt.rcParams['legend.loc'] = 'center left'
    plt.rcParams['axes.linewidth'] = 2

    #plt.gca().spines['right'].set_color('none')
    #plt.gca().spines['top'].set_color('none')
    plt.gca().xaxis.set_ticks_position('bottom')
    plt.gca().yaxis.set_ticks_position('left')

def negativeRAs(ralist):
    newralist=[]
    for ra in ralist:
        if ra >= 300 :
            t=ra - 360.0
            ra = t
        newralist.append(ra)
    return newralist

ra1,dec1=np.loadtxt('boss_survey_outline_points_sorted_ra-60_300_NGC.txt').T
ra2,dec2=np.loadtxt('boss_survey_outline_points_sorted_ra-60_300_SGC.txt').T

tolerance_arcsec=1.5
tolerance_deg = tolerance_arcsec/3600.
sdsscount = 0
bosscount = 0
sdssbosscount = 0

#Readin Hewett & Wild 2010 Redshifts 
hra,hdec,hz,hzerr,first,ahz,zmethod,hplate,hmjd,hfiber = np.loadtxt('/Users/vzm83/Downloads/1045471_Supplementary_Data/mnras0408-2302-SD1.txt',usecols=(1,2,3,4,5,6,7,8,9,10)).T
hplate=[int(x) for x in hplate] ; hmjd=[int(x) for x in hmjd] ; hfiber = [int(x) for x in hfiber]


#For Old file used for SDSS-III containing 2109 tatgets
#baltargets = yanny.read_table_yanny(filename='master-BAL-targets-yanny-format1.dat.txt',tablename='TARGET')

baltargets = yanny.read_table_yanny(filename='green01-TDSS_FES_VARBALmaster1.par.txt',tablename='TARGET')
newtargets=yanny.read_table_yanny('targeting13-explained_more_TDSS_FES_VARBAL_201605.dat',tablename='TARGET')
print len(baltargets)
print baltargets['ra'],baltargets['dec']
#baltargetsra = baltargets['ra']#np.concatenate((baltargets['ra'],newtargets['ra']))
#baltargetsdec = baltargets['dec']#np.concatenate((baltargets['dec'],newtargets['dec']))

baltargetsra = np.concatenate((baltargets['ra'],newtargets['ra']))
baltargetsdec = np.concatenate((baltargets['dec'],newtargets['dec']))

#dr14q = fits.open('../DR14QSO/DR12Q.fits')[1].data
dr14q = fits.open('../DR14QSO/DR14Q_v4_4.fits')[1].data

spAllfile = 'spAll-v5_10_7.fits'
spAll = fits.open(spAllfile)[1].data

sequels = fits.open('spAll-v5_10_0.fits')[1].data
smatch = np.where(sequels['PROGRAMNAME'] == 'sequels')[0]
nsequels = sequels[smatch]  

finaltargets = list()
for ii in range(len(baltargetsra)):
    if (baltargetsra[ii],baltargetsdec[ii]) not in finaltargets :
        finaltargets.append((baltargetsra[ii],baltargetsdec[ii]))

print len(finaltargets)

index1,index2,dist = spherematch(baltargetsra,baltargetsdec,dr14q['RA'],dr14q['DEC'],tolerance_deg,maxmatch=0)
index3,index4,dist = spherematch(baltargetsra,baltargetsdec,spAll['RA'],spAll['DEC'],tolerance_deg,maxmatch=0)
index5,index6,dist = spherematch(baltargetsra,baltargetsdec,sequels['RA'],sequels['DEC'],tolerance_deg,maxmatch=0)
indexs5,indexs6,dist = spherematch(baltargetsra,baltargetsdec,nsequels['RA'],nsequels['DEC'],tolerance_deg,maxmatch=0)
#outputfiles
out = open('Surveybasedcounts_ntotal_3028.txt','w')
out1=open('Threeepochdata_Input_ntotal_3028.txt','w')
completera = []; completedec = [];completeredshift = [];  completeredshift_tag = []
complete_sdsspmf = []; complete_bosspmf = [];complete_ebosspmf = []
complete_sdssmjd = []; complete_bossmjd = [];complete_ebossmjd = []
for i in range(len(index1)):
    match2 = np.where(index3 == index1[i])[0]
    match3 = np.where(index5 == index1[i])[0]
    print index1[i],index2[i],index3[match2],index4[match2]
    print 'TARget RA {0:10.5f}\t{1}\t{2:10.5f}\t{3}'.format( baltargetsra[index1[i]],baltargetsra[index3[match2]],dr14q['RA'][index2[i]],spAll['RA'][index4[match2]])
    print 'TARget DEC {0:10.5f}\t{1}\t{2:10.5f}\t{3}'.format( baltargetsdec[index1[i]],baltargetsdec[index3[match2]],dr14q['DEC'][index2[i]],spAll['DEC'][index4[match2]])
    print 'DR14 PLATE-MJD-FIBER',dr14q['PLATE'][index2[i]],dr14q['MJD'][index2[i]],dr14q['FIBERID'][index2[i]]
    #print 'DR14 SDSS PLATE-MJD-FIBER',dr14q['PLATE_DR7'][index2[i]],dr14q['MJD_DR7'][index2[i]],dr14q['FIBERID_DR7'][index2[i]]
    print 'DR14 DUP PLATE-MJD-FIBER',dr14q['PLATE_DUPLICATE'][index2[i]],dr14q['MJD_DUPLICATE'][index2[i]],dr14q['FIBERID_DUPLICATE'][index2[i]]
    print 'spALL PLATE-MJD_FIBER',spAll['PLATE'][index4[match2]],spAll['MJD'][[index4[match2]]],spAll['FIBERID'][index4[match2]]
    #print '{0:10.5f}\t{1:10.5f}\t{2}\t\t{3}\t{4}'.format(baltargetsra[index1[i]],baltargetsdec[index1[i]],match2,dr14q['N_SPEC_SDSS'][index2[i]],dr14q['N_SPEC_BOSS'][index2[i]])
    sdupplate=[]
    dupplate = dr14q['PLATE_DUPLICATE'][index2[i]]
    dupmjd = dr14q['MJD_DUPLICATE'][index2[i]]
    dupfiber = dr14q['FIBERID_DUPLICATE'][index2[i]]
    dr14z= dr14q['Z'][index2[i]]
    spAllplate = spAll['PLATE'][index4[match2]]
    spAllmjd = spAll['MJD'][[index4[match2]]]
    spAllfiber = spAll['FIBERID'][index4[match2]]
    spAllz = spAll['Z'][index4[match2]]
    sequelsplate = sequels['PLATE'][index6[match3]]
    sequelsmjd = sequels['MJD'][[index6[match3]]]
    sequelsfiber = sequels['FIBERID'][index6[match3]]
    sequelsz = sequels['Z'][index6[match3]]
    for j in range(len(dupplate)):
        xx=np.where(dupplate > 1)[0]
        if len(xx) > 0:
            sdupplate = dupplate[xx]
            sdupmjd = dupmjd[xx]
            sdupfiber = dupfiber[xx]
    finalplate =list();finalmjd =list();finalfiber=list()        
    allpmfstrings =list() 
    allmjdlist = list()
    #create pmf strings
    dr14pmf ='{}-{}-{}'.format(dr14q['PLATE'][index2[i]],dr14q['MJD'][index2[i]],dr14q['FIBERID'][index2[i]])
    allpmfstrings.append(dr14pmf)
    allmjdlist.append(dr14q['MJD'][index2[i]])
    for k in range(len(sdupplate)):
        dupstring = '{0:04d}-{1:05d}-{2:04d}'.format(sdupplate[k],sdupmjd[k],sdupfiber[k])
        allpmfstrings.append(dupstring)
        allmjdlist.append(sdupmjd[k])
    for l in range(len(spAllplate)):
        spAllstring = '{0:04d}-{1:05d}-{2:04d}'.format(spAllplate[l],spAllmjd[l],spAllfiber[l])
        allpmfstrings.append(spAllstring)
        allmjdlist.append(spAllmjd[l])
    for ll in range(len(sequelsplate)):
        sequelsstring = '{0:04d}-{1:05d}-{2:04d}'.format(sequelsplate[ll],sequelsmjd[ll],sequelsfiber[ll])
        allpmfstrings.append(sequelsstring)
        allmjdlist.append(sequelsmjd[ll])
    if len(index2) == 1:
        piperedshift = dr14z
    elif len(spAllz) == 1:
        piperedshift = spAllz[0]
    elif len(sequelsz)==1 :
        piperedshift = sequelsz[0]
    else:
        piperedshift = dr14z
        
    print 'redshift',dr14z,spAllz,sequelsz,piperedshift
    mra=baltargetsra[index1[i]];mdec=baltargetsdec[index1[i]]
    a,b,ds= pshm.spherematch(mra,mdec,hra,hdec,tol=tolerance_deg)
    if len(a) > 0:
        print 'Hewett & Wild',hra[b],hdec[b], hz[b]
        redshift = hz[b]
        redshift_tag = "HW10"
    else:
        redshift = piperedshift 
        redshift_tag = "PIPE"
       # print 'No match in Hewett & Wild'


    print 'redshift',redshift,redshift_tag
    finalpmfstrings =list()
    finalmjd = list()
    for m in range(len(allpmfstrings)):
        if allpmfstrings[m] not in finalpmfstrings:
            finalpmfstrings.append(allpmfstrings[m])
            finalmjd.append(allmjdlist[m])
    print finalpmfstrings,finalmjd
    eboss = 0; boss=0;sdss=0
    sdsspmf=[];bosspmf=[];ebosspmf=[]
    sdssmjd=[];bossmjd=[];ebossmjd=[]
    for jj in range(len(finalmjd)):
        if finalmjd[jj] < 55050:
            sdss+=1
            sdsspmf.append(finalpmfstrings[jj])
            sdssmjd.append(finalmjd[jj])
        elif np.logical_and(finalmjd[jj] >= 55050., finalmjd[jj] <= 56660.):
            boss +=1
            bosspmf.append(finalpmfstrings[jj])
            bossmjd.append(finalmjd[jj])
        else :
            eboss +=1
            ebosspmf.append(finalpmfstrings[jj])
            ebossmjd.append(finalmjd[jj])
    sdsspmf=np.array(sdsspmf);sdssmjd=np.array(sdssmjd)
    bosspmf=np.array(bosspmf);bossmjd=np.array(bossmjd)
    ebosspmf=np.array(ebosspmf);ebossmjd=np.array(ebossmjd)
    sdssinds = np.array(sdssmjd).argsort()
    sortedsdsspmf = sdsspmf[sdssinds];sortedsdssmjd = sdssmjd[sdssinds]
    bossinds = np.array(bossmjd).argsort()
    sortedbosspmf = bosspmf[bossinds];sortedbossmjd = bossmjd[bossinds]
    ebossinds = np.array(ebossmjd).argsort()
    sortedebosspmf = ebosspmf[ebossinds];sortedebossmjd = ebossmjd[ebossinds]
    completera.append(mra);completedec.append(mdec);completeredshift.append(redshift);completeredshift_tag.append(redshift_tag)
    complete_sdsspmf.append(sortedsdsspmf);complete_sdssmjd.append(sortedsdssmjd)
    complete_bosspmf.append(sortedbosspmf);complete_bossmjd.append(sortedbossmjd)
    complete_ebosspmf.append(sortedebosspmf);complete_ebossmjd.append(sortedebossmjd)
    print '--'*101
    print>>out,'{0}\t{1}\t{2}\t{3}\t{4}'.format(baltargetsra[index1[i]],baltargetsdec[index1[i]],sdss,boss,eboss)
        #if  dr14q['N_SPEC_SDSS'][index2[i]] > 0 :
    #    sdsscount += 1
    #if  dr14q['N_SPEC_BOSS'][index2[i]] > 0 :
    #    bosscount += 1
    #if  (dr14q['N_SPEC_SDSS'][index2[i]] > 0) & (dr14q['N_SPEC_BOSS'][index2[i]] > 0) :
    #    sdssbosscount +=1
    #khfkah=raw_input()
#print sdsscount,bosscount,sdssbosscount
    if ((len(sortedsdsspmf) > 0) & (len(sortedbosspmf) > 0) & (len(sortedebosspmf) > 0) ):
        if redshift_tag=='HW10':
            print>>out1,'{0:10.5f}\t{1:10.5f}\t{2:10.5f}\t{3}\t{4}\t{5}'.format(mra,mdec,redshift[0],sortedsdsspmf[-1],sortedbosspmf[-1],sortedebosspmf[-1])
        else:
            print>>out1,'{0:10.5f}\t{1:10.5f}\t{2:10.5f}\t{3}\t{4}\t{5}'.format(mra,mdec,redshift,sortedsdsspmf[-1],sortedbosspmf[-1],sortedebosspmf[-1])
    #fasdfkl=raw_input()

out.close()
out1.close()
np.savez('Compile_PMF_info_ntotal_30128.npz',
        ra = completera ,
        dec = completedec ,
        redshift = completeredshift ,
        redshift_tag = completeredshift_tag ,
        sdsspmf = complete_sdsspmf ,
        bosspmf = complete_bosspmf ,
        ebosspmf = complete_ebosspmf ,
        sdssmjd = complete_sdssmjd ,
        bossmjd = complete_bossmjd ,
        ebossmjd = complete_ebossmjd)

data = np.genfromtxt('Surveybasedcounts_ntotal_3028.txt',names=['ra','dec','sdss','boss','eboss'],dtype=(float,float,int,int,int))
triple_epochs = 0

sdss_boss=0

sdss_eboss=0

boss_eboss = 0

for i in range(len(data)):
    if ((data['sdss'][i] >0) & (data['boss'][i] > 0) & (data['eboss'][i]> 0)): 
        triple_epochs+=1
    if ((data['sdss'][i] >0) & (data['boss'][i] > 0)):
        sdss_boss+=1
    if ((data['sdss'][i] >0) & (data['eboss'][i] > 0)):
        sdss_eboss+=1
    if ((data['boss'][i] > 0) & (data['eboss'][i]> 0)):
        boss_eboss+=1
print '--'*51
print 'Summary of BAL QSO observations'
print '--'*51
print 'No of BAL QSOs with SDSS epochs  ',len(data[np.where(data['sdss']>0)])
print 'No of BAL QSOs with BOSS epochs  ',len(data[np.where(data['boss']>0)])
print 'No of BAL QSOs with eBOSS epochs  ',len(data[np.where(data['eboss']>0)])
print '--'*51
print 'No of BAL QSOs with SDSS and BOSS epochs  ',sdss_boss
print 'No of BAL QSOs with SDSS and eBOSS epochs  ',sdss_eboss
print 'No of BAL QSOs with BOSS and eBOSS epochs  ',boss_eboss
print '--'*51
print 'No of BAL QSOs with SDSS, BOSS and eBOSS epochs  ',triple_epochs
print '--'*51

x1 = np.where(data['sdss'] > 0)[0]
x2 = np.where(data['boss'] > 0)[0]
x3 = np.where(data['eboss'] > 0)[0]

x12 = np.where((data['sdss'] > 0) & (data['boss'] > 0))[0]
x13 = np.where((data['sdss'] > 0) & (data['eboss'] > 0))[0]
x23 = np.where((data['boss'] > 0) & (data['eboss'] > 0))[0]

x123 = np.where((data['sdss'] > 0) & (data['boss'] > 0)&  (data['eboss'] > 0))[0]
fig,(ax,ax1)=plt.subplots(1,2,figsize=(25,10))
init_plotting()
ax.plot(ra1,dec1,'-',color='black',alpha=0.5)
ax.plot(ra2,dec2,'-',color='black',alpha=0.5)
ax1.plot(ra1,dec1,'-',color='black',alpha=0.5)
ax1.plot(ra2,dec2,'-',color='black',alpha=0.5)


ax.plot(negativeRAs(baltargetsra),baltargetsdec,'.',markersize=2,color='black',label='Parent Sample'+'(#'+str(len(baltargetsra))+')')
ax.plot(negativeRAs(data['ra'][x1]),data['dec'][x1],'+',markersize=2,color='red',label='SDSS Sample'+'(#'+str(len(x1))+')')
ax.plot(negativeRAs(data['ra'][x2]),data['dec'][x2],'o',markerfacecolor='none',markersize=3,color='blue',label='BOSS Sample'+'(#'+str(len(x2))+')')
ax.plot(negativeRAs(data['ra'][x3]),data['dec'][x3],'s',markerfacecolor='none',markersize=4,color='magenta',label='eBOSS Sample'+'(#'+str(len(x3))+')')


ax1.plot(negativeRAs(baltargetsra),baltargetsdec,'.',markersize=2,color='black',label='Parent Sample'+'(#'+str(len(baltargetsra))+')')
ax1.plot(negativeRAs(data['ra'][x12]),data['dec'][x12],'x',markersize=2,color='red',label='SDSS + BOSS Sample'+'(#'+str(len(x12))+')')
ax1.plot(negativeRAs(data['ra'][x13]),data['dec'][x13],'o',markerfacecolor='none',markersize=3,lw=2,color='blue',label='SDSS + eBOSS Sample'+'(#'+str(len(x13))+')')
ax1.plot(negativeRAs(data['ra'][x23]),data['dec'][x23],'s',markerfacecolor='none',markersize=4,color='magenta',label='BOSS + eBOSS Sample'+'(#'+str(len(x23))+')')
#ax1.plot(negativeRAs(baltargetsra[indexs5]),baltargetsdec[indexs5],'d',markersize=2,color='green',label='SEQUELS :'+ str(len(indexs5)))
ax1.plot(negativeRAs(data['ra'][x123]),data['dec'][x123],'v',markerfacecolor='none',markersize=5,color='cyan',label='SDSS + BOSS + eBOSS Sample'+'(#'+str(len(x123))+')')

ax.set_xlabel(r'RA',fontsize=20)
ax1.set_xlabel(r'RA',fontsize=20)
ax.set_ylabel(r'DEC',fontsize=20)
ax1.set_ylabel(r'DEC',fontsize=20)
ax1.legend(loc=2)
ax.legend(loc=2)
fig.suptitle('BAL Targets on Sky')

fig.tight_layout()
fig.savefig('SDSSIV_BALQSO_Observation_status_multiple_epochs.jpeg')
fig.savefig('SDSSIV_BALQSO_Observation_status_multiple_epochs.pdf')
plt.show()

