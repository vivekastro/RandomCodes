import numpy as np
import matplotlib.pyplot as plt
import glob 
import astropy.io.fits as fits
import os
import urllib2

def download_spectra(plate, mjd, fiber, dirname='.'):
    '''  Downloads SDSS spectra from DR14 and puts it in dirname
         Change the SDSS URL to download from a different location
    '''
    FITS_FILENAME = 'spec-%(plate)04i-%(mjd)05i-%(fiber)04i.fits'
    try :
        SDSS_URL = ('https://data.sdss.org/sas/dr14/eboss/spectro/redux/v5_10_0/spectra/%(plate)04i/'
                'spec-%(plate)04i-%(mjd)05i-%(fiber)04i.fits')
        urllib2.urlopen(SDSS_URL % dict(plate=plate,mjd=mjd,fiber=fiber))
        print 'Downloadin from dr14'
    except urllib2.HTTPError as err:
        if err.code == 404 :
            SDSS_URL = ('https://data.sdss.org/sas/dr8/sdss/spectro/redux/26/spectra/%(plate)04i/'
            		'spec-%(plate)04i-%(mjd)05i-%(fiber)04i.fits')
            print 'Downloadin from dr8'
    print SDSS_URL % dict(plate=plate,mjd=mjd,fiber=fiber)
    download_url = 'wget   '+SDSS_URL % dict(plate=plate,mjd=mjd,fiber=fiber)
    print download_url
    os.system(download_url)
    mv_cmd='mv '+FITS_FILENAME % dict(plate=plate,mjd=mjd,fiber=fiber) + ' '+dirname+'/.'
    #print mv_cmd
    os.system(mv_cmd)



#Reading the Sky server matched list corresponding to DR7, DR12 and DR14
data7=np.genfromtxt('dr7_result.csv',names=['Name','objID','ra','dec','run','rerun','camcol','field','typei','modelMag_u','modelMag_g','modelMag_r','modelMag_i','modelMag_z','specObjID','plate','mjd','fiberID','z','zConf'],dtype=('|S15',int,float,float,int,int,int,int,'|S15',float,float,float,float,float,int,int,int,int,float,float),delimiter=',')
data12 = fits.open('Skyserver_CrossIDWhitedwarfs_DR12.fits')[1].data  
data14 = fits.open('Skyserver_CrossIDWhitedwarfs_DR14.fits')[1].data

out = open('Whitedwarfs_in_SDSS_with_3epochs.txt','w')
tdata=np.genfromtxt('Whitedwarfs_in_spAll10_10.txt',dtype= ('|S20',float,float),names=['name','ra','dec'],skip_header=1)
count = 0
for i in range(len(tdata)):
    sdsspmf = 'None' ; bosspmf = 'None' ; ebosspmf = 'None'
    xx=np.where((data14['Name'] == tdata['name'][i]) )[0]
    #print xx,data14['ra'][xx],tdata['ra'][i],data14['dec'][xx],tdata['dec'][i],data14['plate'][xx],data14['mjd'][xx],data14['fiberID'][xx]
    for j in range(len(xx)):
        if data14['mjd'][xx[j]] <= 55050 :
            sdsspmf = '{0:04d}-{1:05d}-{2:04d}'.format(data14['plate'][xx[j]],data14['mjd'][xx[j]],data14['fiberID'][xx[j]])
        if ((data14['mjd'][xx[j]] > 55050 ) & (data14['mjd'][xx[j]]  <= 56660 )) :
            bosspmf = '{0:04d}-{1:05d}-{2:04d}'.format(data14['plate'][xx[j]],data14['mjd'][xx[j]],data14['fiberID'][xx[j]])
        if data14['mjd'][xx[j]] > 56660 :
            ebosspmf = '{0:04d}-{1:05d}-{2:04d}'.format(data14['plate'][xx[j]],data14['mjd'][xx[j]],data14['fiberID'][xx[j]])
    if ( (sdsspmf != 'None') & (bosspmf != 'None') & (ebosspmf != 'None')):
        print>>out, '{0}\t{1:10.5f}\t{2:10.5f}\t{3}\t{4}\t{5}'.format(data14['name'][xx[0]],data14['ra'][xx[0]],data14['dec'][xx[0]],sdsspmf,bosspmf,ebosspmf)
        #download_spectra(int(sdsspmf.split('-')[0]),int(sdsspmf.split('-')[1]),int(sdsspmf.split('-')[2]),dirname='Whitedwarf_Spectra')
        #download_spectra(int(bosspmf.split('-')[0]),int(bosspmf.split('-')[1]),int(bosspmf.split('-')[2]),dirname='Whitedwarf_Spectra')
        #download_spectra(int(ebosspmf.split('-')[0]),int(ebosspmf.split('-')[1]),int(ebosspmf.split('-')[2]),dirname='Whitedwarf_Spectra')
        count +=1

print count
out.close()
