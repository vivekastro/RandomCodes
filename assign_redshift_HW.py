import numpy as np
from astropy.io import fits
import scipy as sp
import matplotlib.pyplot as plt
from pyspherematch import *
from pydl.pydl.pydlutils import yanny
#from pydl.pydl.pydlutils.spheregroup import *
import os
from astropy.time import Time

tolerance_arcsec=1.5
tolerance_deg = tolerance_arcsec/3600.

data= np.genfromtxt('Targets_with3epochdata_niels_input.txt',names=['ra','dec','pmf1','pmf2','pmf3'],dtype=(float,float,'|S25','|S25','|S25'))

hra,hdec,hz,hzerr,first,ahz,zmethod,hplate,hmjd,hfiber = np.loadtxt('/Users/vzm83/Downloads/1045471_Supplementary_Data/mnras0408-2302-SD1.txt',usecols=(1,2,3,4,5,6,7,8,9,10)).T
hplate=[int(x) for x in hplate] ; hmjd=[int(x) for x in hmjd] ; hfiber = [int(x) for x in hfiber]
completera = []; completedec = [];completeredshift = []; complete_epochs = []; completeredshift_tag = []

dr14q = fits.open('../DR14QSO/DR14Q_v4_4.fits')[1].data
out = open('Targets_with3epochdata_niels_input_with_redshift.txt','w')

for j in range(len(data)):
    mra=data['ra'][j];mdec=data['dec'][j]
    a,b,ds= spherematch(mra,mdec,hra,hdec,tol=tolerance_deg)
    if len(a) > 0:
        #print 'Hewett & Wild',hra[b],hdec[b], hz[b]
        redshift = hz[b[0]]
        redshift_tag = "HW10"
    else:
        aa,bb,ds= spherematch(mra,mdec,dr14q['RA'],dr14q['DEC'],tol=tolerance_deg)
        if len(aa) > 0:
            #print 'ZPIPE',dr14q['RA'][bb],dr14q['DEC'][bb], dr14q['Z'][bb]
            redshift = dr14q['Z'][bb[0]] 
            redshift_tag = "PIPE"
        else :
            print '******************************************************notfound'
            redshift = -999
            redshift_tag = "None"
            print 'Redshift: ',redshift, 'Redshift_Tag:',redshift_tag
    print>>out, '{0:10.5f}\t{1:10.5f}\t{2:10.5f}\t{3}\t{4}\t{5}'.format(data['ra'][j],data['dec'][j],redshift,data['pmf1'][j],data['pmf2'][j],data['pmf3'][j])

out.close()
