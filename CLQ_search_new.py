import numpy as np
import astropy.io.fits as fits
import scipy as sp
import matplotlib.pyplot as plt
import astroML 
from matplotlib.ticker import NullFormatter
import pandas as pd
import os
import urllib2
from scipy import   stats
from pydl.pydl.pydlutils.spheregroup import *
from astroML.plotting import hist
from matplotlib.backends.backend_pdf import PdfPages
from scipy.interpolate import interp1d
from astropy.modeling.models import Voigt1D
from astropy import constants as const
from astropy import units as U
from astropy.coordinates import SkyCoord
from dustmaps.sfd import SFDQuery
from specutils import extinction 
from matplotlib.ticker import MaxNLocator
import seaborn as sns
from statsmodels.stats.outliers_influence import summary_table
import statsmodels.formula.api as smf
import statsmodels.api as sm
from pydl.pydl.pydlutils import yanny
from astropy.table import Table
from astropy.time import Time
from scipy.stats import kde
from scipy.ndimage.filters import gaussian_filter1d,uniform_filter1d

"""
Program to explore the search of CLQs in eBOSS quasar catalog
@author : Vivek M.
@date   : 12/April/2019
@version: 1.0
"""
spallversion='v5_13_0'
params = {
   'axes.labelsize': 18,
   'axes.linewidth': 1.5,
   #'text.fontsize': 8,
   'legend.fontsize': 15,
   'xtick.labelsize': 18,
   'ytick.labelsize': 18,
   'text.usetex': True,
   #'figure.figsize': [16, 5]
   'legend.frameon': False,
   'font.family': 'Times New Roman'
   }
plt.rcParams.update(params)

hfont = {'fontname':'Times New Roman'}

def computeSN1700(wave,flux,err):
    ww=np.where((wave >=1650) & (wave <= 1750))[0]
    return np.median(flux[ww])/np.median(err[ww])

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


def checkUniqueN():
    data=fits.open('CrossMatch_DR12Q_spAll_v5_13_0.fits')[1].data
    print 'Number of entries matching DR12 and spAll_v5-13-0:',len(data['SDSS_NAME'])
    print 'Number of quasars matching DR12 and spAll_v5-13-0:',len(np.unique(data['SDSS_NAME']))

def MultipleEpochQSOs():
    redux= '/uufs/chpc.utah.edu/common/home/sdss/ebosswork/eboss/spectro/redux/v5_13_0/spectra/lite'
    data=fits.open('CrossMatch_DR12Q_spAll_v5_13_0.fits')[1].data
    uname = np.unique(data['SDSS_NAME'])
    count = 0
    gcount = 0
    out=open('Candidate_CLQ_search_DR16.txt','w')
    name = [] ; ra=[] ; dec=[];zvi=[]
    umag1 =[] ;gmag1=[];rmag1=[];imag1=[];zmag1=[]
    umag2 =[] ;gmag2=[];rmag2=[];imag2=[];zmag2=[]
    #umag1err =[] ;gmag1err=[];rmag1err=[];imag1err=[];zmag1err=[]
    #umag2err =[] ;gmag2err=[];rmag2err=[];imag2err=[];zmag2err=[]

    plate1 = [] ; mjd1=[];fiber1=[]
    plate2 = [] ; mjd2=[];fiber2=[]
    gdiff = [] ; rdiff = [] ; idiff = []
    for i in range(len(uname)):
    #for i in range(150):
        xx=np.where(data['SDSS_NAME'] == uname[i])[0]
        if len(xx)>1:
            ndata = data[xx]
            #print uname[i],len(xx),xx,data['PLATE_1'][xx[0]],data['MJD_1'][xx[0]],data['FIBERID_1'][xx[0]],data['PLATE_2'][xx[0]],data['MJD_2'][xx[0]],data['FIBERID_2'][xx[0]],data['PLATE_2'][xx[-1]],data['MJD_2'][xx[-1]],data['FIBERID_2'][xx[-1]],data['FIBERMAG'][xx[0],1],data['FIBER2MAG'][xx[-1],1],data['FIBERFLUX'][xx[0],2],data['FIBERFLUX'][xx[-1],2]
            mjdl = ndata['MJD_2']
            maxmjd = max(mjdl)
            minmjd = min(mjdl)
            xmax = np.where(mjdl == maxmjd)[0][0]
            xmin = np.where(mjdl == minmjd)[0][0]
            print  mjdl,maxmjd,minmjd
            print xmax,xmin,ndata['MJD_2'][xmax],ndata['PLATE_2'][xmax],ndata['FIBERID_2'][xmax],ndata['MJD_2'][xmin],ndata['PLATE_2'][xmin],ndata['FIBERID_2'][xmin]
            #ksjhdf=raw_input()
            #print 'Check', ndata['MJD_2'],ndata['SDSS_NAME'],ndata['PLATE_2'],ndata['FIBERID_2']
            #plate1 = data['PLATE_2'][xx[0]] ; plate2 = data['PLATE_2'][xx[-1]]
            #pmf1 = '{0:04d}-{1:05d}-{2:04d}'.format(data['PLATE_2'][xx[0]],data['MJD_2'][xx[0]],data['FIBERID_2'][xx[0]])
            #pmf2 = '{0:04d}-{1:05d}-{2:04d}'.format(data['PLATE_2'][xx[-1]],data['MJD_2'][xx[-1]],data['FIBERID_2'][xx[-1]])
            #data1 = fits.open(os.path.join(redux,plate1,'spec-'+pmf1+'.fits'))[1].data
            gdiff.append((22.5 - 2.5*np.log10(ndata['SPECTROFLUX'][xmin,1])) - (22.5 - 2.5*np.log10(ndata['SPECTROFLUX'][xmax,1])))
            rdiff.append((22.5 - 2.5*np.log10(ndata['SPECTROFLUX'][xmin,2])) - (22.5 - 2.5*np.log10(ndata['SPECTROFLUX'][xmax,2])))
            idiff.append((22.5 - 2.5*np.log10(ndata['SPECTROFLUX'][xmin,3])) - (22.5 - 2.5*np.log10(ndata['SPECTROFLUX'][xmax,3])))
            if np.abs((22.5 - 2.5*np.log10(ndata['SPECTROFLUX'][xmin,1])) - (22.5 - 2.5*np.log10(ndata['SPECTROFLUX'][xmax,1]))) > 1 :
                print>>out,'{0}\t{1}\t{2:10.5f}\t{3:10.5f}\t{4:10.5f}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}'.format( uname[i],len(xx),ndata['RA_1'][0],ndata['DEC_1'][0],ndata['Z_VI'][0],ndata['PLATE_2'][xmin],ndata['MJD_2'][xmin],ndata['FIBERID_2'][xmin],ndata['PLATE_2'][xmax],ndata['MJD_2'][xmax],ndata['FIBERID_2'][xmax])
                gcount +=1
            name.append(ndata['SDSS_NAME'][0])
            ra.append(ndata['RA_1'][0])
            dec.append(ndata['DEC_1'][0])
            zvi.append(ndata['Z_VI'][0])
            plate1.append(ndata['PLATE_2'][xmin])
            mjd1.append(ndata['MJD_2'][xmin])
            fiber1.append(ndata['FIBERID_2'][xmin])
            plate2.append(ndata['PLATE_2'][xmax])
            mjd2.append(ndata['MJD_2'][xmax])
            fiber2.append(ndata['FIBERID_2'][xmax])
 
            umag1.append((22.5 - 2.5*np.log10(ndata['SPECTROFLUX'][xmin,0])))
            gmag1.append((22.5 - 2.5*np.log10(ndata['SPECTROFLUX'][xmin,1])))
            rmag1.append((22.5 - 2.5*np.log10(ndata['SPECTROFLUX'][xmin,2])))
            imag1.append((22.5 - 2.5*np.log10(ndata['SPECTROFLUX'][xmin,3])))
            zmag1.append((22.5 - 2.5*np.log10(ndata['SPECTROFLUX'][xmin,4])))

            umag2.append((22.5 - 2.5*np.log10(ndata['SPECTROFLUX'][xmax,0])))
            gmag2.append((22.5 - 2.5*np.log10(ndata['SPECTROFLUX'][xmax,1])))
            rmag2.append((22.5 - 2.5*np.log10(ndata['SPECTROFLUX'][xmax,2])))
            imag2.append((22.5 - 2.5*np.log10(ndata['SPECTROFLUX'][xmax,3])))
            zmag2.append((22.5 - 2.5*np.log10(ndata['SPECTROFLUX'][xmax,4])))
                
            count +=1

        print str(i+1)+'/'+str(len(uname))+' Running and Found candidates: '+str(gcount) 
    np.savez('CLQsearch_MasterList_Plate-MJD-Fiber.npz',
            name = np.array(name) ,
            ra = np.array(ra) ,
            dec = np.array(dec) ,
            zvi = np.array(zvi) ,
            plate1 = np.array(plate1) ,
            mjd1 = np.array(mjd1) ,
            fiber1 = np.array(fiber1) ,
            plate2 = np.array(plate2) ,
            mjd2 = np.array(mjd2) ,
            fiber2 = np.array(fiber2) ,
            umag1=np.array(umag1) ,
            gmag1=np.array(gmag1) ,
            rmag1=np.array(rmag1) ,
            imag1=np.array(imag1) ,
            zmag1=np.array(zmag1) ,
            umag2=np.array(umag2) ,
            gmag2=np.array(gmag2) ,
            rmag2=np.array(rmag2) ,
            imag2=np.array(imag2) ,
            zmag2=np.array(zmag2) ,
            )
    
    #print count
    #print gdiff

    out.close()
    gdiff=np.array(gdiff)
    rdiff=np.array(rdiff)
    idiff=np.array(idiff)
    yy=np.where(np.abs(gdiff) > 1)[0]

    fig,(ax,ax1,ax2)=plt.subplots(1,3,figsize=(10,5))
    
    

    ax.plot(gdiff,idiff,'.',color='black',label='gmag vs imag')
    nbins=20
    #k = kde.gaussian_kde((gdiff,idiff))
    #xi, yi = np.mgrid[gdiff.min():gdiff.max():nbins*1j, idiff.min():y.max():nbins*1j]
    #zi = k(np.vstack([xi.flatten(), yi.flatten()]))
    #ax.pcolormesh(xi, yi, zi.reshape(xi.shape), shading='gouraud', cmap=plt.cm.BuGn_r)
    #ax.contour(xi, yi, zi.reshape(xi.shape) )
    ax.plot(gdiff,idiff,'.',color='black',label='gmag vs imag')
    
    ax1.plot(gdiff,rdiff,'.',color='red',label='gmag vs rmag')
    
    ax2.plot(rdiff,idiff,'.',color='blue',label='rmag vs imag')
    ax.set(xlabel='$\Delta$g-mag',ylabel='$\Delta$i-mag')
    ax1.set(xlabel='$\Delta$g-mag',ylabel='$\Delta$r-mag')
    ax2.set(xlabel='$\Delta$r-mag',ylabel='$\Delta$i-mag')
    fig.tight_layout()
    fig.savefig('Candidate_CLQsearch_DR16_color-magn.jpg')
    plt.show()


def download_spectraCLQ(print_cmd=False):
    clqC=np.genfromtxt('Candidate_CLQ_search_DR16c.txt',names=['name','nepoch','ra','dec','zvi','plate1','mjd1','fiber1','plate2','mjd2','fiber2'],dtype=('|S30',int,float,float,float,int,int,int,int,int,int))
    if print_cmd :
        cptxt = open('CLQcopyspectrumfromBOSS_SPECTRO_REDUX_v5_13.txt','w')
    checkdownload=[]
    specdir = 'CLQsearch_spectra'
    for i in range(len(clqC['name'])):
    #for i in range(10):
        plates=[clqC['plate1'][i],clqC['plate2'][i]]
        mjds=[clqC['mjd1'][i],clqC['mjd2'][i]]
        fibers=[clqC['fiber1'][i],clqC['fiber2'][i]]
        for j in range(2): 
            if plates[j] >= 10000:
                FITS_FILENAME = 'spec-{0:05d}-{1:05d}-{2:04d}.fits'.format(plates[j],mjds[j],fibers[j])
            else:
                FITS_FILENAME = 'spec-{0:04d}-{1:05d}-{2:04d}.fits'.format(plates[j],mjds[j],fibers[j])
            if not print_cmd:
                if not ((os.path.isfile(FITS_FILENAME)) | (os.path.isfile(os.path.join(specdir,FITS_FILENAME)))):
                    download_spectra(plates[j],mjds[j],fibers[j],specdir)
                else : 
                    print 'Spectrum already downloaded'

            if print_cmd:
                print>>cptxt, 'cp   /uufs/chpc.utah.edu/common/home/sdss/ebosswork/eboss/spectro/redux/{0}/spectra/lite/{1}/{2}  ~/BOSS_BALDATA_CLQ/.'.format(spallversion,plates[j],FITS_FILENAME)        
            else :
                if not ((os.path.isfile(FITS_FILENAME)) & (os.path.isfile(os.path.join(specdir,FITS_FILENAME)))):
                    checkdownload.append(FITS_FILENAME)
            print '--'*31
   #         kfgh=raw_input()
    print checkdownload
    if print_cmd:
        cptxt.close()
def smoothCRTS(mjd,mag,magerr):
    mjd= np.array(mjd) ; mag = np.array(mag) ; magerr= np.array(magerr)
    minmjd =min(mjd)
    msort = np.argsort(mjd)
    mjd = mjd[msort] ; mag = mag[msort] ; magerr = magerr[msort]
    mjddiff = mjd[1:] - mjd[0:-1]

    gp = np.where(mjddiff > 100)[0]
    
    print type(gp)
    ngp = np.insert(np.insert(gp,0,0),len(gp)+1,len(mjd)-1)
    #fig,ax = plt.subplots(figsize=(10,5))
    #ax.plot(mjd,mag,'ok',alpha=0.2,label=name)
    #ax.set_xlabel('MJD')
    #ax.set_ylabel('CRTS V-mag')
    #ax.legend(loc=1)
    medmjd = []
    medmag = []
    medmagerr = []
    for ig,g in enumerate(ngp[0:-1]):
        #print mjd[g]
        if ig == 0:
            xg = np.where((mjd >= mjd[ngp[ig]]-10) & (mjd <=mjd[ngp[ig+1]]+10))[0]
            medmag.append(np.median(mag[xg]))
            medmjd.append(np.mean(mjd[xg]))
            medmagerr.append(np.std(mag[xg])/np.sqrt(len(xg)))
     #       ax.axvline(mjd[g]-10,ls='--',color='red')
        else:
            xg = np.where((mjd >= mjd[ngp[ig]]+10) & (mjd <=mjd[ngp[ig+1]]+10))[0]
            medmag.append(np.median(mag[xg]))
            medmjd.append(np.mean(mjd[xg]))
            medmagerr.append(np.std(mag[xg])/np.sqrt(len(xg)))
    return medmjd,medmag,medmagerr

def plot_spectra():
    text_font = {'fontname':'Times New Roman', 'size':'14'}
    pp = PdfPages('CLQsearches_plot_spectra_sn1700gt6.pdf') 
    clqC=np.genfromtxt('Candidate_CLQ_search_DR16c.txt',names=['name','nepoch','ra','dec','zvi','plate1','mjd1','fiber1','plate2','mjd2','fiber2'],dtype=('|S30',int,float,float,float,int,int,int,int,int,int))
    master =np.load('CLQsearch_MasterList_Plate-MJD-Fiber.npz')
    data=fits.open('CrossMatch_DR12Q_spAll_v5_13_0.fits')[1].data

    crts = pd.read_csv('CRTS_lc_CLQsearchSample_sn1700gt6.csv')
    specdir = 'CLQsearch_spectra'
    linelist = np.genfromtxt('/Users/vzm83/Proposals/linelist_speccy.txt',usecols=(0,1,2),dtype=('|S10',float,'|S5'),names=True)
    out = open('CLQsearch_deltamgs_sn1700gt6.txt','w') 
    cptxt = open('copyCLQadditionalspectra.txt','w')
    for i in range(len(clqC['name'])):
    #for i in range(100):
        print 'Working on ',i,' source'
        xx  = np.where((master['plate1'] == clqC['plate1'][i]) &(master['mjd1'] == clqC['mjd1'][i]) & (master['fiber1'] == clqC['fiber1'][i]))[0][0]
        xy =np.where(data['SDSS_NAME'] == clqC['name'][i])[0]
        ndata = data[xy]
        print xx
        if clqC['plate1'][i] >= 10000:
            FITS_FILENAME1 = 'spec-{0:05d}-{1:05d}-{2:04d}.fits'.format(clqC['plate1'][i],clqC['mjd1'][i],clqC['fiber1'][i])
        else:
            FITS_FILENAME1 = 'spec-{0:04d}-{1:05d}-{2:04d}.fits'.format(clqC['plate1'][i],clqC['mjd1'][i],clqC['fiber1'][i])
        if clqC['plate2'][i] >= 10000:
            FITS_FILENAME2 = 'spec-{0:05d}-{1:05d}-{2:04d}.fits'.format(clqC['plate2'][i],clqC['mjd2'][i],clqC['fiber2'][i])
        else:
            FITS_FILENAME2 = 'spec-{0:04d}-{1:05d}-{2:04d}.fits'.format(clqC['plate2'][i],clqC['mjd2'][i],clqC['fiber2'][i]) 

        data1 = fits.open(os.path.join(specdir,FITS_FILENAME1))[1].data
        data2 = fits.open(os.path.join(specdir,FITS_FILENAME2))[1].data
        gflux1 = (data1['flux']*(data1['and_mask'] == 0)).copy()
        gflux2 = (data2['flux']*(data2['and_mask'] == 0)).copy()
        zvi = clqC['zvi'][i]
        sn1 =computeSN1700(10**data1['loglam']/(1.0+zvi),data1['flux'],1.0/np.sqrt(data1['ivar']))
        sn2 =computeSN1700(10**data2['loglam']/(1.0+zvi),data2['flux'],1.0/np.sqrt(data2['ivar']))
        if ((np.median(gflux1) != 0) & (np.median(gflux2) != 0) & (sn1 > 6) & (sn2 > 6)) :
            fig=plt.figure(figsize=(18,8))
            ax=plt.subplot2grid((2, 3), (0, 0), colspan=2,rowspan=2)
            ax1 = plt.subplot2grid((2, 3), (0, 2))
            ax.plot(10**data1['loglam']/(1.0+zvi),gaussian_filter1d(data1['flux'],2),color='black',alpha=0.5,label=FITS_FILENAME1.split('.')[0][5:])
            ax.plot(10**data2['loglam']/(1.0+zvi),gaussian_filter1d(data2['flux'],2),color='red',alpha=0.5,label=FITS_FILENAME2.split('.')[0][5:])

            ax.plot(10**data1['loglam']/(1.0+zvi),1.0/np.sqrt(data1['ivar']),color='black',alpha=0.1)
            ax.plot(10**data2['loglam']/(1.0+zvi),1.0/np.sqrt(data2['ivar']),color='red',alpha=0.1)

            string1 = 'SDSS J{0}\tZ\_VI: {1:4.4f}\tN$\_{{spec}}$: {2}'.format(clqC['name'][i],zvi,clqC['nepoch'][i])
            
            string2 = 'RA: {0:4.4f}\tDEC: {1:4.4f}'.format(clqC['ra'][i],clqC['dec'][i])
            string3 = '{0:20}- {1:3.2f}  {2:3.2f}  {3:3.2f}  {4:3.2f}  {5:3.2f}'.format(FITS_FILENAME1.split('.')[0][5:], master['umag1'][xx], master['gmag1'][xx], master['rmag1'][xx], master['imag1'][xx], master['zmag1'][xx])
            string4 = '{0:20}- {1:3.2f}  {2:3.2f}  {3:3.2f}  {4:3.2f}  {5:3.2f}'.format(FITS_FILENAME2.split('.')[0][5:], master['umag2'][xx], master['gmag2'][xx], master['rmag2'][xx], master['imag2'][xx], master['zmag2'][xx])
            string5 = '{0:20}- {1:3.2f}  {2:3.2f}  {3:3.2f}  {4:3.2f}  {5:3.2f}'.format('$\Delta m(2-1)$', master['umag2'][xx] -master['umag1'][xx] , master['gmag2'][xx]-master['gmag1'][xx], master['rmag2'][xx]-master['rmag1'][xx], master['imag2'][xx]-master['imag1'][xx], master['zmag2'][xx] - master['zmag1'][xx])

            ax.set(xlabel='Rest wavelength ($\AA$)',ylabel='Flux',ylim=(-2,max(np.median(gflux1)+3*np.std(gflux1),np.median(gflux2)+3*np.std(gflux2) )))
            xlim,ylim=ax.get_xlim(),ax.get_ylim()
            print string1,xlim,ylim
            ax.text(xlim[0]+0.05*(xlim[1] - xlim[0]),ylim[1]-0.05*(ylim[1] - ylim[0]), string1,fontsize=18)
            ax.text(xlim[0]+0.05*(xlim[1] - xlim[0]),ylim[1]-0.09*(ylim[1] - ylim[0]), string2,fontsize=18)
            ax.text(xlim[0]+0.05*(xlim[1] - xlim[0]),ylim[1]-0.13*(ylim[1] - ylim[0]), string3,fontsize=18)
            ax.text(xlim[0]+0.05*(xlim[1] - xlim[0]),ylim[1]-0.17*(ylim[1] - ylim[0]), string4,fontsize=18)
            ax.text(xlim[0]+0.05*(xlim[1] - xlim[0]),ylim[1]-0.21*(ylim[1] - ylim[0]), string5,fontsize=18)


            obslambda = linelist['lambda']#*(1.+zvi)
            x = np.where((obslambda > xlim[0]) & (obslambda < xlim[1]))[0]
	    plotlambda = obslambda[x]
	    plotname = linelist['Name'][x]
	    plota_e = linelist['a_e'][x]
	    #print plotlambda
	    for k in range(len(plotlambda)):
	        if plota_e[k].strip() == 'Abs.' : 
	    	    ax.axvline(x=plotlambda[k], color='lawngreen', linestyle=':')
	    	    ax.text(plotlambda[k],ylim[0]+0.75*(ylim[1]-ylim[0]),plotname[k],color='Orange',ha='center',rotation=90,**text_font)
	        else :
	    	    ax.axvline(x=plotlambda[k], color='lightblue', linestyle=':')
	    	    ax.text(plotlambda[k],ylim[0]+0.75*(ylim[1]-ylim[0]),plotname[k],color='Brown',ha='center',rotation=90,**text_font)
            #Download and plot the other epoch data
            dupmjd=[]
            if clqC['nepoch'][i] > 2:
                xyz = np.where((ndata['MJD_2'] !=clqC['mjd1'][i]) & (ndata['MJD_2'] !=clqC['mjd2'][i]) )[0]
                for k in range(len(xyz)):
                    if ndata['PLATE_2'][xyz[k]] >=10000 :
                         FITS_FILENAME = 'spec-{0:05d}-{1:05d}-{2:04d}.fits'.format(ndata['PLATE_2'][xyz[k]],ndata['MJD_2'][xyz[k]],ndata['FIBERID_2'][xyz[k]])
                    else:
                         FITS_FILENAME = 'spec-{0:04d}-{1:05d}-{2:04d}.fits'.format(ndata['PLATE_2'][xyz[k]],ndata['MJD_2'][xyz[k]],ndata['FIBERID_2'][xyz[k]])

                    if not ((os.path.isfile(FITS_FILENAME)) | (os.path.isfile(os.path.join(specdir,FITS_FILENAME)))):
                        download_spectra(ndata['PLATE_2'][xyz[k]],ndata['MJD_2'][xyz[k]],ndata['FIBERID_2'][xyz[k]],specdir)
                    print>>cptxt, 'cp   /uufs/chpc.utah.edu/common/home/sdss/ebosswork/eboss/spectro/redux/{0}/spectra/lite/{1}/{2}  ~/BOSS_BALDATA_CLQ/.'.format(spallversion,ndata['PLATE_2'][xyz[k]],FITS_FILENAME)
                    data0 = fits.open(os.path.join(specdir,FITS_FILENAME))[1].data
                    ax.plot(10**data0['loglam']/(1.0+zvi),gaussian_filter1d(data0['flux'],2),color=plt.cm.RdYlBu(k*300),alpha=0.5,label=FITS_FILENAME.split('.')[0][5:])
                    ax.plot(10**data0['loglam']/(1.0+zvi),1.0/np.sqrt(data0['ivar']),color=plt.cm.RdYlBu(k*300),alpha=0.1)
                    dupmjd.append(ndata['MJD_2'][xyz[k]])
            ax.legend(loc=1)

            crm = np.where(crts['InputID'] ==clqC['name'][i] )[0]
            if len(crm) > 0 :
                CRTS_Vmag = crts['Mag'][crm]
                CRTS_Verr = crts['Magerr'][crm]
                CRTS_MJD = crts['MJD'][crm]
                ax1.errorbar(CRTS_MJD,CRTS_Vmag,yerr=CRTS_Verr,fmt='v',color='gold',label='CRTS',alpha=0.75)
                CRTS_medmjd,CRTS_medmag,CRTS_medmagerr = smoothCRTS(CRTS_MJD,CRTS_Vmag,CRTS_Verr) 
                ax1.errorbar(CRTS_medmjd,CRTS_medmag,yerr=CRTS_medmagerr,fmt='v',color='brown',alpha=0.75)
                ax1.set_ylim(ax1.get_ylim()[::-1])
                ax1.set(xlabel='MJD',ylabel='V-mag')
                ax1_ylim = ax1.get_ylim()
                ax1.legend(loc=1)
                ax1.axvline(clqC['mjd1'][i],ls='--',color='black',lw=3,zorder=-1,alpha=0.45)
                ax1.text(clqC['mjd1'][i], ax1_ylim[0]+0.2*(ax1_ylim[1] - ax1_ylim[0]),str(clqC['mjd1'][i]),fontsize=12,rotation='vertical')
                ax1.axvline(clqC['mjd2'][i],ls='--',color='red',lw=3,zorder=-1,alpha=0.45)
                ax1.text(clqC['mjd2'][i], ax1_ylim[0]+0.2*(ax1_ylim[1] - ax1_ylim[0]),str(clqC['mjd2'][i]),fontsize=12,rotation='vertical')
                for mm in dupmjd:
                    ax1.axvline(mm,ls='--',color='blue',lw=3,zorder=-1,alpha=0.45)
                    ax1.text(mm, ax1_ylim[0]+0.2*(ax1_ylim[1] - ax1_ylim[0]),str(mm),fontsize=12,rotation='vertical')

            fig.tight_layout()
            fig.savefig(pp,format='pdf')
            print>>out,  '{0}\t{1}\t{2:10.5f}\t{3:10.5f}\t{4:10.5f}\t{5:10.5f}\t{6:10.5f}\t{7:10.5f}\t{8:10.5f}\t{9:10.5f}'.format(clqC['name'][i],clqC['nepoch'][i],clqC['ra'][i],clqC['dec'][i],clqC['zvi'][i], master['umag2'][xx] -master['umag1'][xx] , master['gmag2'][xx]-master['gmag1'][xx], master['rmag2'][xx]-master['rmag1'][xx], master['imag2'][xx]-master['imag1'][xx], master['zmag2'][xx] - master['zmag1'][xx])
    out.close()
    pp.close()
def main():
    print '-------'
    #checkUniqueN()
    #MultipleEpochQSOs()
    #download_spectraCLQ()
    plot_spectra()
if __name__== "__main__":
    main()
