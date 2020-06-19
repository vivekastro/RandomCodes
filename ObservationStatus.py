import numpy as np
from astropy.io import fits
import scipy as sp
import matplotlib.pyplot as plt
#from pyspherematch import *
from pydl.pydl.pydlutils import yanny
from pydl.pydl.pydlutils.spheregroup import *
import os
from astropy.time import Time

'''
This program plots the cumulative histogram and footprint of BAL quasars in the SDSS-IV
The program mainly does a spherematch between the spAll file and the targets in the initial
BAL catalogue
SEQUELS data is included in the plot. If not required, just uncomment the MJD based filtering
'''

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

#For Old file used for SDSS-III containing 2109 tatgets
#baltargets = yanny.read_table_yanny(filename='master-BAL-targets-yanny-format1.dat.txt',tablename='TARGET')

#For SDSS-IV USe the following file containing 2958 sources
baltargets = yanny.read_table_yanny(filename='green01-TDSS_FES_VARBALmaster1.par.txt',tablename='TARGET')
newtargets=yanny.read_table_yanny('targeting13-explained_more_TDSS_FES_VARBAL_201605.dat',tablename='TARGET')
print len(baltargets)
print baltargets['ra'],baltargets['dec']
baltargetsra = np.concatenate((baltargets['ra'],newtargets['ra']))
baltargetsdec = np.concatenate((baltargets['dec'],newtargets['dec']))


spAllfile = 'spAll-v5_10_10.fits'

spAll = fits.open(spAllfile)[1].data

tolerance_arcsec=1.5
tolerance_deg = tolerance_arcsec/3600.

#eb = np.where(spAll['MJD'] >= 56890)[0]
#print len(spAll),len(eb)
eboss = spAll
index1,index2,dist = spherematch(baltargetsra,baltargetsdec,eboss['RA'],eboss['DEC'],tolerance_deg,maxmatch=0)
out=open('spherematch_results_all_tol_1.5_v5_10_10.txt','w')
for i in range(len(index1)):
    print>>out, '{0:10.5f}\t{1:10.5f}\t{2:10.5f}\t{3:10.5f}\t{4}\t{5}'.format(baltargetsra[index1[i]],eboss['RA'][index2[i]],baltargetsdec[index1[i]],eboss['DEC'][index2[i]],eboss['MJD'][index2[i]],eboss['PLATE'][index2[i]])
#
print len(index1)
out.close()

pra,sra,pdec,sdec,mjd,plate=np.loadtxt('spherematch_results_all_tol_1.5_v5_10_10.txt').T
print mjd

radec = zip(pra,pdec)
unq, unq_inv, unq_cnt = np.unique(pra, return_inverse=True, return_counts=True)
oldlist = [];newlist=[];un_mjd=[]

for k in range(len(radec)):
    if radec[k] not in oldlist:
        oldlist.append(radec[k])
        newlist.append(radec[k])
        un_mjd.append(mjd[k])
    else:
        pass

ura = [];udec=[]
n_epochs=[]
for j in range(len(newlist)):
    mra=newlist[j][0];mdec=newlist[j][1]
    ep =  np.where((pra == mra) & (pdec == mdec))[0]
    n_epochs.append(len(ep))
    ura.append(mra);udec.append(mdec)

n_epochs =np.array(n_epochs)
ura =np.array(ura);udec=np.array(udec)

print len(oldlist),len(newlist),len(un_mjd)
n_bins = 100
fig, (ax1,ax) = plt.subplots(1,2,figsize=(20, 10))
init_plotting()
n, bins, patches = ax.hist(mjd, n_bins,color='blue', normed=0, histtype='step',
                           cumulative=True, label='Spectra')
n1, bins1, patches1 = ax.hist(un_mjd, n_bins,color='red', normed=0, histtype='step',
                           cumulative=True, label='QSO')

ax.set_title('Observations vs MJD')
ax.set_xlabel('MJD',fontsize=20)
ax.set_ylabel(r'Cumulative number of BAL quasars & Spectra',fontsize=20)
ax.legend(loc=2)
ylim = ax.get_ylim()
xlim = ax.get_xlim()
ax.text(xlim[1]-0.1*(xlim[1] - xlim[0]),ylim[1]-0.051*(ylim[1]-ylim[0]),str(len(mjd)),color='blue',fontsize=20)
ax.text(xlim[1]-0.15*(xlim[1] - xlim[0]),ylim[1]-0.25*(ylim[1]-ylim[0]),str(len(un_mjd)),color='red',fontsize=20)
for label in (ax.get_xticklabels() + ax.get_yticklabels()):
    label.set_fontsize(13)



#fig1,ax1=plt.subplots(figsize=(15, 8))
init_plotting()
x1 = np.where(n_epochs == 1)[0]
x2 = np.where(n_epochs >= 2)[0]
ax1.plot(ra1,dec1,'-',color='black',alpha=0.5)
ax1.plot(ra2,dec2,'-',color='black',alpha=0.5)
ax1.set_xlim(-55,300)
ax1.set_ylim(-15,100)
ax1.plot(negativeRAs(baltargetsra),baltargetsdec,'.',markersize=3,color='black',label='Parent Sample'+'(#'+str(len(baltargetsra))+')')
ax1.plot(negativeRAs(ura[x1]),udec[x1],'o',color='red',markersize=3,label='1 epoch'+'(#'+str(len(x1))+')')
ax1.plot(negativeRAs(ura[x2]),udec[x2],'s',color='blue',markersize=3,label='2 or more epochs'+'(#'+str(len(x2))+')')
ax1.set_xlabel(r'RA',fontsize=20)
ax1.set_ylabel(r'DEC',fontsize=20)
ax1.set_title('Targets on Sky')
ylim = ax1.get_ylim()
xlim = ax1.get_xlim()
throughdate=Time(np.max(mjd),format='mjd')
throughdate.format = 'fits'
print throughdate.value[0:10]
ax1.text(xlim[0]+0.05*(xlim[1] - xlim[0]),ylim[1]-0.05*(ylim[1] - ylim[0]),'Through '+str(throughdate.value[0:10]),fontsize=20)
ax1.legend(loc=1)
for label1 in (ax1.get_xticklabels() + ax1.get_yticklabels()):
    label1.set_fontsize(13)

fig.tight_layout()
fig.savefig('SDSSIV_BALQSO_Observation_status_v5_10_10_all_new_tol_1.5_sky.jpeg')
#plt.show()

