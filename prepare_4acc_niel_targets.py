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


spAllfile = 'spAll-v5_10_7.fits'

spAll = fits.open(spAllfile)[1].data

tolerance_arcsec=1.5
tolerance_deg = tolerance_arcsec/3600.

#eb = np.where(spAll['MJD'] >= 56890)[0]
#print len(spAll),len(eb)
eboss = spAll
out=open('spherematch_results_allspAll_tol_1.5.txt','w')

index1,index2,dist = spherematch(baltargetsra,baltargetsdec,eboss['RA'],eboss['DEC'],tolerance_deg,maxmatch=0)
for i in range(len(index1)):
    print>>out, '{0:10.5f}\t{1:10.5f}\t{2:10.5f}\t{3:10.5f}\t{4}\t{5}\t{6}'.format(baltargetsra[index1[i]],eboss['RA'][index2[i]],baltargetsdec[index1[i]],eboss['DEC'][index2[i]],eboss['PLATE'][index2[i]],eboss['MJD'][index2[i]],eboss['FIBERID'][index2[i]])
#
out.close()
print len(index1)


spAllfile = 'spAll-v5_10_0.fits'

spAll = fits.open(spAllfile)[1].data
eboss = spAll
out1=open('spherematch_results_allsequels_tol_1.5.txt','w')

index1,index2,dist = spherematch(baltargetsra,baltargetsdec,eboss['RA'],eboss['DEC'],tolerance_deg,maxmatch=0)
for i in range(len(index1)):
    print>>out1, '{0:10.5f}\t{1:10.5f}\t{2:10.5f}\t{3:10.5f}\t{4}\t{5}\t{6}'.format(baltargetsra[index1[i]],eboss['RA'][index2[i]],baltargetsdec[index1[i]],eboss['DEC'][index2[i]],eboss['PLATE'][index2[i]],eboss['MJD'][index2[i]],eboss['FIBERID'][index2[i]])
#
out1.close()


spAll = fits.open('../DR14QSO/DR14Q_v4_4.fits')[1].data
eboss = spAll
out2=open('spherematch_results_alldr14q_tol_1.5.txt','w')

index1,index2,dist = spherematch(baltargetsra,baltargetsdec,eboss['RA'],eboss['DEC'],tolerance_deg,maxmatch=0)
for i in range(len(index1)):
    print>>out2, '{0:10.5f}\t{1:10.5f}\t{2:10.5f}\t{3:10.5f}\t{4}\t{5}\t{6}'.format(baltargetsra[index1[i]],eboss['RA'][index2[i]],baltargetsdec[index1[i]],eboss['DEC'][index2[i]],eboss['PLATE'][index2[i]],eboss['MJD'][index2[i]],eboss['FIBERID'][index2[i]])
#
out2.close()


radec = zip(baltargetsra,baltargetsdec)

