import numpy as np
import astropy.io.fits as fits
import scipy as sp
import matplotlib.pyplot as plt
import astroML 
from matplotlib.ticker import NullFormatter
import pandas as pd
import os,glob
import urllib2
from scipy import   stats
#from pydl.pydl.pydlutils.spheregroup import *
from astroML.plotting import hist
from matplotlib.backends.backend_pdf import PdfPages
from scipy.interpolate import interp1d
from astropy.modeling.models import Voigt1D
from astropy import constants as const
from astropy import units as U
from astropy.coordinates import SkyCoord
#from dustmaps.sfd import SFDQuery
#from specutils import extinction 
from matplotlib.ticker import MaxNLocator
import seaborn as sns
from statsmodels.stats.outliers_influence import summary_table
import statsmodels.formula.api as smf
import statsmodels.api as sm
#from pydl.pydl.pydlutils import yanny
from astropy.table import Table
from astropy.time import Time
from scipy.ndimage.filters import gaussian_filter1d,uniform_filter1d
from matplotlib.ticker import NullFormatter

"""
Program to explore the full eBOSS BALQSO sample. 
@author : Vivek M.
@date   : 22/Feb/2019
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
def plot_hline(y,**kwargs):
    data = kwargs.pop("data") #get the data frame from the kwargs
    plt.axhline(y=y, c='black',linestyle='--',zorder=-1) #zorder places the line underneath the other points

n_v = 1240.81
si_iv = 1393.755
c_iv = 1549.48
al_iii = 1857.4
mg_ii = 2799.117
p_v = 1122.9925#1128.008#1117.977 #
def negativeRAs(ralist):
    newralist=[]
    for ra in ralist:
        if ra >= 300 :
            t=ra - 360.0
            ra = t
        newralist.append(ra)
    return newralist

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


def onclick(event):
    '''
    Record the points by clicking on the plots.
    '''
    global ix, iy
    ix, iy = event.xdata, event.ydata

    print 'x = %d, y = %d'%(ix, iy)

    # assign global variable to access outside of function
    #global coords
    #coords.append((ix, iy))

    # Disconnect after 20 clicks
    #if len(coords) == 20:
   #     fig.canvas.mpl_disconnect(cid)
   #     plt.close(1)
    return (ix,iy)


def ray_tracing_method(x,y,poly):
    '''
    Check if the point is inside the polygon
    '''
    n = len(poly)
    inside = False

    p1x,p1y = poly[0]
    for i in range(n+1):
        p2x,p2y = poly[i % n]
        if y > min(p1y,p2y):
            if y <= max(p1y,p2y):
                if x <= max(p1x,p2x):
                    if p1y != p2y:
                        xints = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xints:
                        inside = not inside
        p1x,p1y = p2x,p2y

    return inside







def merge2initialSamples():
    
    #For Old file used for SDSS-III containing 2109 tatgets
    #baltargets = yanny.read_table_yanny(filename='master-BAL-targets-yanny-format1.dat.txt',tablename='TARGET')

    #For SDSS-IV USe the following file containing 2958 sources
    baltargets = yanny.read_table_yanny(filename='green01-TDSS_FES_VARBALmaster1.par.txt',tablename='TARGET')
    newtargets=yanny.read_table_yanny('targeting13-explained_more_TDSS_FES_VARBAL_201605.dat',tablename='TARGET')
    out = open('Master_initial_sample.txt','w') 
    baltargetsra = np.concatenate((baltargets['ra'],newtargets['ra']))
    baltargetsdec = np.concatenate((baltargets['dec'],newtargets['dec']))
    for i in range(len(baltargetsra)):
        print>>out, 'Target{0:04d}\t{1:10.5f}\t{2:10.5f}'.format(i+1,baltargetsra[i],baltargetsdec[i])
    print 'Initial sample has {} sources; {} from initial and {} from later'.format(len(baltargetsra),len(baltargets),len(newtargets))
    print 'Initial sample has unique {} sources; {} from initial and {} from later'.format(len(np.unique(baltargetsra)),len(np.unique(baltargets)),len(np.unique(newtargets)))
    
    out.close()

def getMulti_epochInfo():

    #odata = fits.open('Skyserver_CrossID_DR14_gibson.fits')[1].data
    #ndata = fits.open('topcatMatch_spAll_v5_13_0_gibson_BAL_sample.fits')[1].data
    #master = np.genfromtxt('Master_gibson_2005_targets_cor.txt',names=['name','ra','dec'],dtype=('|S15',float,float),skip_header=1)
    
    odata = fits.open('Skyserver_CrossID_DR14_master.fits')[1].data
    ndata = fits.open('topcatMatch_spAll_v5_13_0_master_BAL_sample_allmatches.fits')[1].data
    master = np.genfromtxt('Master_initial_sample.txt',names=['name','ra','dec'],dtype=('|S15',float,float),skip_header=1)

    masterpmf = [] ; masterplate = []; mastermjd = []; masterfiber=[];mastername = [] ; masterra=[]; masterdec=[]
    for i in range(len(master)):
    #for i in range(20):
        xx = np.where(ndata['col1'] == master['name'][i])[0]
        
        #yy = np.where(odata['col1'] == master['name'][i])[0]
        yy = np.where(odata['Name'] == master['name'][i])[0] #******CHANGE HERE for 3028***********


        plate = [] ; mjd = []; fiber = []
        plate_mjd_fiber = []

        if len(yy) > 0 :
            yodata = odata[yy]
            for j in range(len(yodata)):
                pmf = '{0:05d}-{1:05d}-{2:04d}'.format(yodata['plate'][j],yodata['mjd'][j],yodata['fiberID'][j])
                if pmf not in plate_mjd_fiber :
                    plate_mjd_fiber.append(pmf)
                    plate.append(yodata['plate'][j])
                    mjd.append(yodata['mjd'][j])
                    fiber.append(yodata['fiberID'][j])

        if len(xx) > 0 :
            xndata = ndata[xx]
            for k in range(len(xndata)):
                pmf = '{0:05d}-{1:05d}-{2:04d}'.format(xndata['PLATE'][k],xndata['MJD'][k],xndata['FIBERID'][k])
                if pmf not in plate_mjd_fiber :
                    plate_mjd_fiber.append(pmf)
                    plate.append(xndata['PLATE'][k])
                    mjd.append(xndata['MJD'][k])
                    fiber.append(xndata['FIBERID'][k])
        if ((master['name'][i] == 'nTarget0106') |(master['name'][i] ==  'Target0453')) :
             pmf = '{0:05d}-{1:05d}-{2:04d}'.format(707,52177,392)
             if pmf not in plate_mjd_fiber :
                plate_mjd_fiber.append(pmf)
                plate.append(707)
                mjd.append(52177)
                fiber.append(392)
        if ((master['name'][i] == 'nTarget0164') |(master['name'][i] ==  'Target0517')) :
             pmf = '{0:05d}-{1:05d}-{2:04d}'.format(431,51877,550)
             if pmf not in plate_mjd_fiber :
                plate_mjd_fiber.append(pmf)
                plate.append(431)
                mjd.append(51877)
                fiber.append(550)
        if ((master['name'][i] == 'nTarget1984') |(master['name'][i] ==  'Target2870')) :
            pmf = '{0:05d}-{1:05d}-{2:04d}'.format(645,52203,480)
            if pmf not in plate_mjd_fiber :
                plate_mjd_fiber.append(pmf)
                plate.append(645)
                mjd.append(52203)
                fiber.append(480)
        if ((master['name'][i] == 'nTarget0153') |(master['name'][i] ==  'Target0505')) :
            pmf = '{0:05d}-{1:05d}-{2:04d}'.format(1865,53312,261)
            if pmf not in plate_mjd_fiber :
                plate_mjd_fiber.append(pmf)
                plate.append(1865)
                mjd.append(53312)
                fiber.append(261)
        if ((master['name'][i] == 'nTarget0346') |(master['name'][i] ==  'Target0731')) :
            pmf = '{0:05d}-{1:05d}-{2:04d}'.format(467,51901,192)
            if pmf not in plate_mjd_fiber :
                plate_mjd_fiber.append(pmf)
                plate.append(467)
                mjd.append(51901)
                fiber.append(192)
        if ((master['name'][i] == 'nTarget0774') |(master['name'][i] ==  'Target1230')) :
            pmf = '{0:05d}-{1:05d}-{2:04d}'.format(1002,52646,297)
            if pmf not in plate_mjd_fiber :
                plate_mjd_fiber.append(pmf)
                plate.append(1002)
                mjd.append(52646)
                fiber.append(297)


        if len(mjd)>0:
            print '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(i+1,master['name'][i],len(xx),len(yy),min(mjd),max(mjd),max(plate),len(plate))
            masterpmf.append(plate_mjd_fiber)
            masterplate.append(plate)
            mastermjd.append(mjd)
            masterfiber.append(fiber)
            mastername.append(master['name'][i])
            masterra.append(master['ra'][i])
            masterdec.append(master['dec'][i])
        if len(plate)>0:
            print plate_mjd_fiber#,masterpmf
    data = {'name': np.array(master['name']),'pmf': np.array(masterpmf),'plate':np.array(masterplate),'mjd':np.array(mastermjd), 'fiber':np.array(masterfiber)}
    #pdata = Table(data)
    print len(data)
    #c1 = fits.Column(name='name', array=np.array(master['name']), format='15A')
    #c2 = fits.Column(name='pmf', array=np.array(masterpmf,dtype=np.object), format='PA(16)')
    #c3 = fits.Column(name='plate', array=np.array(masterplate,dtype=np.object), format='PI()')
    #c4 = fits.Column(name='mjd', array=np.array(mastermjd,dtype=np.object), format='PI()')
    #c5 = fits.Column(name='fiber', array=np.array(masterfiber,dtype=np.object), format='PI()')
    #tfits = fits.BinTableHDU.from_columns([c1,c2,  c3, c4, c5])
    #tfits.writeto('MasterList_Plate-MJD-Fiber.fits')
    #np.savez('MasterList_Plate-MJD-Fiber_2005.npz',
    #        name = np.array(mastername) ,
    #        ra = np.array(masterra) ,
    #        dec = np.array(masterdec) ,
    #        pmf = np.array(masterpmf) ,
    #        plate = np.array(masterplate) ,
    #        mjd = np.array(mastermjd) ,
    #        fiber = np.array(masterfiber) ,
    #        )
    #data = np.load('MasterList_Plate-MJD-Fiber_2005.npz')
    #out = open('Master_multi-epoch_information_numberofEpochs_v5_13_2005.txt','w')
    #*********** CHANGE HERE for 3028
    np.savez('MasterList_Plate-MJD-Fiber.npz',
            name = np.array(mastername) ,
            ra = np.array(masterra) ,
            dec = np.array(masterdec) ,
            pmf = np.array(masterpmf) ,
            plate = np.array(masterplate) ,
            mjd = np.array(mastermjd) ,
            fiber = np.array(masterfiber) ,
            )
#
    data = np.load('MasterList_Plate-MJD-Fiber.npz')
    out = open('Master_multi-epoch_information_numberofEpochs_v5_13.txt','w')
    for i in range(len(data['ra'])):
        print>>out, '{0}\t{1:10.5f}\t{2:10.5f}\t{3}'.format(data['name'][i],data['ra'][i],data['dec'][i],len(data['mjd'][i]))
    out.close()

def plot_nepochs():
    #master = np.genfromtxt('Master_gibson_2005_targets_cor.txt',names=['name','ra','dec'],dtype=('|S15',float,float),skip_header=1)
    #msum = np.genfromtxt('Master_multi-epoch_information_numberofEpochs_v5_13_2005.txt',names=['name','ra','dec','nepochs'])
    
    master = np.genfromtxt('Master_initial_sample.txt',names=['name','ra','dec'],dtype=('|S15',float,float),skip_header=1)
    msum = np.genfromtxt('Master_multi-epoch_information_numberofEpochs_v5_13.txt',names=['name','ra','dec','nepochs'])
    
    ra1,dec1=np.loadtxt('boss_survey_outline_points_sorted_ra-60_300_NGC.txt').T
    ra2,dec2=np.loadtxt('boss_survey_outline_points_sorted_ra-60_300_SGC.txt').T
    
    fig,ax = plt.subplots(figsize=(10,8))
    x1 = np.where(msum['nepochs'] == 1)[0]
    x2 = np.where(msum['nepochs'] == 2)[0]
    x3 = np.where(msum['nepochs'] >= 3)[0]
    print len(x1),len(x2)
    ax.plot(ra1,dec1,'-',color='black',alpha=0.5)
    ax.plot(ra2,dec2,'-',color='black',alpha=0.5)
    ax.plot(negativeRAs(master['ra']),master['dec'],'.',markersize=3,color='black',label='Parent Sample'+'(\#'+str(len(master['ra']))+')')#label='Parent Sample(\#2005)')#
    ax.plot(negativeRAs(msum['ra'][x1]),msum['dec'][x1],'o',color='red',markersize=3,label='1 epoch'+'(\#'+str(len(x1))+')')
    ax.plot(negativeRAs(msum['ra'][x2]),msum['dec'][x2],'s',color='blue',markersize=3,label='2  epochs'+'(\#'+str(len(x2))+')')
    ax.plot(negativeRAs(msum['ra'][x3]),msum['dec'][x3],'s',color='magenta',markersize=3,label='3 or more epochs'+'(\#'+str(len(x3))+')')
    ax.set_xlim(-55,300)
    ax.set_ylim(-15,100)
    ax.legend(loc=1)
    ax.grid()
    xlim,ylim=ax.get_xlim(),ax.get_ylim()
    throughdate=Time(58528,format='mjd')
    throughdate.format = 'fits'
    print throughdate.value[0:10]
    ax.text(xlim[0]+0.05*(xlim[1] - xlim[0]),ylim[1]-0.05*(ylim[1] - ylim[0]),'Through '+str(throughdate.value[0:10]),fontsize=20)

    ax.set(xlabel='RA',ylabel='DEC')
    fig.tight_layout()
    #fig.savefig('SDSS_IV_BALquasar_sample_nepochs_plot_v5_13_2005.jpg')
    fig.savefig('SDSS_IV_BALquasar_sample_nepochs_plot_v5_13.jpg')
   #plt.show()
    

def download_data(print_cmd=False):
    if print_cmd :
        #cptxt = open('copyspectrumfromBOSS_SPECTRO_REDUX_v5_13_2005.txt','w')
        cptxt = open('copyspectrumfromBOSS_SPECTRO_REDUX_v5_13.txt','w')
    else:
        checkdownload=[]
    #data = np.load('MasterList_Plate-MJD-Fiber_2005.npz')
    data = np.load('MasterList_Plate-MJD-Fiber.npz')
    specdir = 'SDSSIV_BALdata'
    for i in range(len(data['name'])):
    #for i in range(30):
        plates = data['plate'][i];mjds=data['mjd'][i];fibers=data['fiber'][i]
        for j in range(len(plates)):
            print data['name'][i], plates[j],mjds[j],fibers[j]
            plate = plates[j]; mjd = mjds[j]; fiber = fibers[j]
            if plate >= 10000:
                FITS_FILENAME = 'spec-{0:05d}-{1:05d}-{2:04d}.fits'.format(plate,mjd,fiber)
            else:
                FITS_FILENAME = 'spec-{0:04d}-{1:05d}-{2:04d}.fits'.format(plate,mjd,fiber)
            if not print_cmd:
                if not ((os.path.isfile(FITS_FILENAME)) | (os.path.isfile(os.path.join(specdir,FITS_FILENAME)))):
                    download_spectra(plates[j],mjds[j],fibers[j],specdir)
                else : 
                    print 'Spectrum already downloaded'

            if print_cmd:
                print>>cptxt, 'cp   /uufs/chpc.utah.edu/common/home/sdss/ebosswork/eboss/spectro/redux/{0}/spectra/lite/{1}/{2}  ~/BOSS_BALDATA/.'.format(spallversion,plate,FITS_FILENAME)        
            else :
                if not ((os.path.isfile(FITS_FILENAME)) & (os.path.isfile(os.path.join(specdir,FITS_FILENAME)))):
                    checkdownload.append(FITS_FILENAME)
            print '--'*31
    print checkdownload
    if print_cmd:
        cptxt.close()

def hist_multiepochs():
    #data = np.genfromtxt('Master_multi-epoch_information_numberofEpochs_v5_13_2005.txt',names=['name','ra','dec','nepochs'],dtype=('|S15',float,float,int))
    data = np.genfromtxt('Master_multi-epoch_information_numberofEpochs_v5_13.txt',names=['name','ra','dec','nepochs'],dtype=('|S15',float,float,int))
    fig,ax=plt.subplots(figsize=(5,5))
    bins=np.arange(1,11,1)
    hist(data['nepochs'],normed=False,bins=bins,histtype="step", linewidth= 2,alpha= 1, linestyle="-",color='red',align='mid', ax=ax)
    #ax.set(xlabel='Number of epochs',ylabel='Number of quasars',title='SDSS-IV: Initial sample 2005 sources')
    ax.set(xlabel='Number of epochs',ylabel='Number of quasars',title='SDSS-IV: Full sample 3028 sources')
    ax.grid()
    fig.tight_layout()
    #fig.savefig('eBOSS_fullsample_number_epochs_hist_v5_13_2005.jpg')
    fig.savefig('eBOSS_fullsample_number_epochs_hist_v5_13.jpg')
    #plt.show()

def cumulative_MJD():
    #data = np.load('MasterList_Plate-MJD-Fiber_2005.npz')
    data = np.load('MasterList_Plate-MJD-Fiber.npz')
    SDSSIV_MJDs = []; uSDSSIV_MJDs = []
    for i in range(len(data['name'])):
        tmjd = np.array(data['mjd'][i])
        xx=np.where(tmjd > 56777)[0]
        #print tmjd[xx]
        if len(xx) >0:
            uSDSSIV_MJDs.append(min(tmjd[xx]))
            for j in xx:
                print i,j,tmjd,xx
                SDSSIV_MJDs.append(tmjd[j])

    SDSSIV_MJDs = np.array(SDSSIV_MJDs)
    uSDSSIV_MJDs = np.array(uSDSSIV_MJDs)

    mjdbounds = [56870,57235,57601,57966,58331,58543]
    yr0 = np.where((SDSSIV_MJDs > 56777) & (SDSSIV_MJDs <= 56870))[0]
    yr1 = np.where((SDSSIV_MJDs > 56870) & (SDSSIV_MJDs <= 57235))[0]
    yr2 = np.where((SDSSIV_MJDs > 57235) & (SDSSIV_MJDs <= 57601))[0]
    yr3 = np.where((SDSSIV_MJDs > 57601) & (SDSSIV_MJDs <= 57966))[0]
    yr4 = np.where((SDSSIV_MJDs > 57966) & (SDSSIV_MJDs <= 58331))[0]
    yr5 = np.where((SDSSIV_MJDs > 58331) & (SDSSIV_MJDs <= 58543))[0]

    uyr0 = np.where((uSDSSIV_MJDs > 56777) & (uSDSSIV_MJDs <= 56870))[0]
    uyr1 = np.where((uSDSSIV_MJDs > 56870) & (uSDSSIV_MJDs <= 57235))[0]
    uyr2 = np.where((uSDSSIV_MJDs > 57235) & (uSDSSIV_MJDs <= 57601))[0]
    uyr3 = np.where((uSDSSIV_MJDs > 57601) & (uSDSSIV_MJDs <= 57966))[0]
    uyr4 = np.where((uSDSSIV_MJDs > 57966) & (uSDSSIV_MJDs <= 58331))[0]
    uyr5 = np.where((uSDSSIV_MJDs > 58331) & (uSDSSIV_MJDs <= 58543))[0]
    print "No. of quasars:",len(uyr1),len(uyr2),len(uyr3),len(uyr4),len(uyr5)
    print "No. of spectra:",len(yr1),len(yr2),len(yr3),len(yr4),len(yr5)
    print "No. of  quasars in sequels (56777-56870)",len(uyr0)
    print "No. of  spectra in sequels (56777-56870)",len(yr0)
    fig,ax=plt.subplots(figsize=(10,5))
    bins = np.arange(56820,58600)
    hist(SDSSIV_MJDs,normed=False,bins='knuth',histtype="step", linewidth= 2,alpha= 1, linestyle="-",color='red',align='mid', ax=ax)
    hist(SDSSIV_MJDs,normed=False,bins=bins,cumulative=True,histtype="step", linewidth= 3,alpha= 0.5, linestyle="--",color='blue',align='mid', ax=ax,label='spectra')
    hist(uSDSSIV_MJDs,normed=False,bins=bins,cumulative=True,histtype="step", linewidth= 3,alpha= 0.5, linestyle=":",color='magenta',align='mid', ax=ax,label='quasars')
    #ax.set(xlabel='MJDs',ylabel='Number of spectra',title='SDSS-IV rate of data collection: Initial 2005 sample')
    ax.set(xlabel='MJDs',ylabel='Number of spectra',title='SDSS-IV rate of data collection: Full 3028 sample')
    for mb in mjdbounds:
        ax.axvline(mb,ls=':',linewidth=3)
    xlim,ylim=ax.get_xlim(),ax.get_ylim()
    ax.text(56870+0.4*(57235-56870),ylim[1]-0.1*(ylim[1] - ylim[0]),'Year1')
    ax.text(57235+0.4*(57601-57235),ylim[1]-0.1*(ylim[1] - ylim[0]),'Year2')
    ax.text(57601+0.4*(57966-57601),ylim[1]-0.1*(ylim[1] - ylim[0]),'Year3')
    ax.text(57966+0.4*(58331-57966),ylim[1]-0.1*(ylim[1] - ylim[0]),'Year4')
    ax.text(58331+0.4*(58543-58331),ylim[1]-0.1*(ylim[1] - ylim[0]),'Year5')
    ax.text(56870+0.4*(57235-56870),ylim[1]-0.15*(ylim[1] - ylim[0]),str(len(yr1)))
    ax.text(57235+0.4*(57601-57235),ylim[1]-0.15*(ylim[1] - ylim[0]),str(len(yr2)))
    ax.text(57601+0.4*(57966-57601),ylim[1]-0.15*(ylim[1] - ylim[0]),str(len(yr3)))
    ax.text(57966+0.4*(58331-57966),ylim[1]-0.15*(ylim[1] - ylim[0]),str(len(yr4)))
    ax.text(58331+0.4*(58543-58331),ylim[1]-0.15*(ylim[1] - ylim[0]),str(len(yr5)))
    fig.legend(loc=3)
    fig.tight_layout()
    #fig.savefig('eBOSS_fullsample_rate_of_data_hist_2005.jpg')
    fig.savefig('eBOSS_fullsample_rate_of_data_hist.jpg')
    #plt.show()

def plot_spPlate():
    sp = fits.open('spAll-v5_13_0.fits')[1].data
    data = fits.open('spAll-v5_10_10.fits')[1].data
    yy=np.where(data['PROGRAMNAME'] == 'sequels')[0]
    sequels = data[yy]

    ra1,dec1=np.loadtxt('boss_survey_outline_points_sorted_ra-60_300_NGC.txt').T
    ra2,dec2=np.loadtxt('boss_survey_outline_points_sorted_ra-60_300_SGC.txt').T
    elgngc = np.where(sp['PROGRAMNAME'] == 'ELG_NGC')[0]
    elgsgc = np.where(sp['PROGRAMNAME'] == 'ELG_SGC')[0]
    eboss = np.where(sp['PROGRAMNAME'] == 'eboss')[0]
    speboss = sp[eboss]
    g = np.where(speboss['CLASS'] == 'GALAXY')[0]
    g1 = np.where(sequels['CLASS'] == 'GALAXY')[0]
    q = np.where(speboss['CLASS'] == 'QSO')[0]
    q1 = np.where(sequels['CLASS'] == 'QSO')[0]
    
    fig,(ax,ax1,ax2) = plt.subplots(1,3,figsize=(20,8))
    ax.plot(ra1,dec1,'-',color='black',alpha=0.5)
    ax.plot(ra2,dec2,'-',color='black',alpha=0.5)
    ax1.plot(ra1,dec1,'-',color='black',alpha=0.5)
    ax1.plot(ra2,dec2,'-',color='black',alpha=0.5)
    
    ax2.plot(ra1,dec1,'-',color='black',alpha=0.5)
    ax2.plot(ra2,dec2,'-',color='black',alpha=0.5)
    
 
    ax.plot(negativeRAs(speboss['RA'][g]),speboss['DEC'][g],'.',markersize=3,color='black',label='GALAXY(\#'+str(len(g)+len(g1))+')')
    ax.plot(negativeRAs(sequels['RA'][g1]),sequels['DEC'][g1],'.',markersize=3,color='black')
    ax1.plot(negativeRAs(speboss['RA'][q]),speboss['DEC'][q],'.',markersize=3,color='red',label='QSO(\#'+str(len(q)+len(q1))+')')
    ax1.plot(negativeRAs(sequels['RA'][q1]),sequels['DEC'][q1],'.',markersize=3,color='red')
    ax2.plot(negativeRAs(sp['RA'][elgngc]),sp['DEC'][elgngc],'.',markersize=3,color='blue',alpha=0.3,label='ELG(\#'+str(len(elgngc)+len(elgsgc))+')')
    ax2.plot(negativeRAs(sp['RA'][elgsgc]),sp['DEC'][elgsgc],'.',markersize=3,color='blue',alpha=0.3)
    ax.set(xlabel='RA',ylabel='DEC',title='eBOSS Galaxy footprint')
    ax.set_xlim(-55,300)
    ax.set_ylim(-15,100)
    ax.legend(loc=1)
    ax.grid()
    xlim,ylim=ax.get_xlim(),ax.get_ylim()
    throughdate=Time(58528,format='mjd')
    throughdate.format = 'fits'
    print throughdate.value[0:10]
    ax.text(xlim[0]+0.05*(xlim[1] - xlim[0]),ylim[1]-0.05*(ylim[1] - ylim[0]),'Through '+str(throughdate.value[0:10]),fontsize=20)
    
    ax1.set(xlabel='RA',ylabel='DEC',title='eBOSS QSO footprint')
    ax1.set_xlim(-55,300)
    ax1.set_ylim(-15,100)
    ax1.legend(loc=1)
    ax1.grid()
    xlim,ylim=ax1.get_xlim(),ax1.get_ylim()
    throughdate=Time(58528,format='mjd')
    throughdate.format = 'fits'
    print throughdate.value[0:10]
    ax1.text(xlim[0]+0.05*(xlim[1] - xlim[0]),ylim[1]-0.05*(ylim[1] - ylim[0]),'Through '+str(throughdate.value[0:10]),fontsize=20)

    ax2.set(xlabel='RA',ylabel='DEC',title='eBOSS ELG footprint')
    ax2.set_xlim(-55,300)
    ax2.set_ylim(-15,100)
    ax2.legend(loc=1)
    ax2.grid()
    xlim,ylim=ax2.get_xlim(),ax2.get_ylim()
    throughdate=Time(58528,format='mjd')
    throughdate.format = 'fits'
    print throughdate.value[0:10]
    ax2.text(xlim[0]+0.05*(xlim[1] - xlim[0]),ylim[1]-0.05*(ylim[1] - ylim[0]),'Through '+str(throughdate.value[0:10]),fontsize=20)

    fig.tight_layout()

    fig.savefig('SDSSIV_eBOSS_final_sample_footprint.jpg')
    #plt.show()

def three_epochs_count():
    data = np.load('MasterList_Plate-MJD-Fiber_2005.npz')
    out = open('Fullsample_3epoch_counts_2005.txt','w')
    #data = np.load('MasterList_Plate-MJD-Fiber.npz')
    #out = open('Fullsample_3epoch_counts.txt','w')
    for i in range(len(data['name'])):
        tmjd = np.array(data['mjd'][i])
        s1 = np.where((tmjd < 55176 ))[0]
        s2 = np.where((tmjd >= 55176 ) & (tmjd <= 56777))[0]
        s3 = np.where((tmjd > 56777 ))[0]
        
        if len(s1) > 0 :
            sdss1 = len(s1)
        else :
            sdss1 = 0
        if len(s2) > 0 :
            sdss2 = len(s2)
        else :
            sdss2 = 0
        if len(s3) > 0 :
            sdss3 = len(s3)
        else :
            sdss3 = 0

        print>>out,'{0}\t{1:10}\t{2:10}\t{3:10}\t{4:10}\t{5:10}'.format(data['name'][i],data['ra'][i],data['dec'][i],sdss1,sdss2,sdss3)
    out.close()

    mcount = np.genfromtxt('Fullsample_3epoch_counts_2005.txt',names=['name','ra','dec','sdss1','sdss2','sdss3'],dtype=('|S15',float,float,int,int,int))
    #mcount = np.genfromtxt('Fullsample_3epoch_counts.txt',names=['name','ra','dec','sdss1','sdss2','sdss3'],dtype=('|S15',float,float,int,int,int))
    xx111 = np.where((mcount['sdss1'] > 0) & (mcount['sdss2'] > 0) & (mcount['sdss3']>0))[0]
    xx110 = np.where((mcount['sdss1'] > 0) & (mcount['sdss2'] > 0) & (mcount['sdss3']==0))[0]
    xx101 = np.where((mcount['sdss1'] > 0) & (mcount['sdss2'] == 0) & (mcount['sdss3']>0))[0]
    xx011 = np.where((mcount['sdss1'] == 0) & (mcount['sdss2'] > 0) & (mcount['sdss3']>0))[0]
    print 'Sources with atleast 1 epoch each in SDSS-I/II, SDSS-III, SDSS-IV',len(xx111)
    print 'Sources with atleast 1 epoch in SDSS-I/II, SDSS-III, and no SDSS-IV',len(xx110)
    print 'Sources with atleast 1 epoch in SDSS-I/II, SDSS-IV, and no SDSS-III',len(xx101)
    print 'Sources with atleast 1 epoch in SDSS-III, SDSS-IV, and no SDSS-I/II',len(xx011)


def restFrame_timescales():
    fsm = np.genfromtxt('Master_full_sample_redshift.txt',names=['name','ra','dec','Z_VI','tag'],dtype=('|S25',float,float,float,'|S25'))#fits.open('Fullsample_3028_BAL_sources_crossmatch-DR12Q.fits')[1].data
    psm = np.genfromtxt('Master_gibson_sample_redshift.txt',names=['name','ra','dec','Z_VI','tag'],dtype=('|S25',float,float,float,'|S25'))#fits.open('Initialsample_2005_BAL_sources_crossmatch-DR12Q.fits')[1].data
    data = np.load('MasterList_Plate-MJD-Fiber.npz')    
    data1 = np.load('MasterList_Plate-MJD-Fiber_2005.npz')  

    Flbal = fits.open('LoBAL_sample_fullsample.fits')[1].data
    Ilbal = fits.open('LoBAL_sample_initialsample.fits')[1].data
    bins = np.arange(0,14,0.25)
    fmaxt = [] ;pmaxt=[] ;fmint=[] ; pmint=[] ;fzvi=[] ; pzvi=[] ;fhb=[];phb=[];flb=[];plb=[]
    for i in range(len(fsm)):
        xx=np.where(data['name'] == fsm['name'][i])[0]
        #fdata= data[xx]
        if len(xx) > 0:
            tmjd = sorted(np.array(data['mjd'][xx]))[0]
            if len(tmjd) > 1:
                maxtime = (max(tmjd) - min(tmjd))/(1.0+fsm['Z_VI'][i])/365.
                mintime = (tmjd[1] - tmjd[0])/(1.0+fsm['Z_VI'][i])/365.
                fmaxt.append(maxtime)
                fmint.append(mintime)
                fzvi.append(fsm['Z_VI'][i])
                if fsm['name'][i] in Flbal['col1']:
                    flb.append(maxtime)
                else:
                    fhb.append(maxtime)
                print 'Fmaxt', i

    for ii in range(len(psm)):
        pxx=np.where(data1['name'] == psm['name'][ii])[0]
        #pdata= data[pxx]
        if len(pxx) > 0:
            ttmjd = sorted(np.array(data1['mjd'][pxx]))[0]
            print ttmjd
            if len(ttmjd) > 1:
                pmaxtime = (max(ttmjd) - min(ttmjd))/(1.0+psm['Z_VI'][ii])/365
                pmintime = np.abs(ttmjd[1] - ttmjd[0])/(1.0+psm['Z_VI'][ii])/365
                pmaxt.append(pmaxtime)
                pmint.append(pmintime)
                pzvi.append(psm['Z_VI'][ii])
                if psm['name'][ii] in Ilbal['col1']:
                    plb.append(pmaxtime)
                else:
                    phb.append(pmaxtime)

                print 'Pmaxt', ii

    fmaxt = np.array(fmaxt) ;pmaxt=np.array(pmaxt) ;fmint=np.array(fmint) ; pmint=np.array(pmint) ;fzvi = np.array(fzvi);pzvi=np.array(pzvi)
    fhb=np.array(fhb)
    phb=np.array(phb)
    flb=np.array(flb)
    plb=np.array(plb)

    #fhb = np.where((fzvi> 1.6 ) & (fzvi < 5.6))[0]
    #phb = np.where((pzvi > 1.6 ) & (pzvi < 5.6))[0]
    #flb = np.where((fzvi> 0.25 ) & (fzvi < 3.95))[0]
    #plb = np.where((pzvi> 0.25 ) & (pzvi < 3.95))[0]
    print fmaxt
    fig,(ax,ax1)=plt.subplots(1,2,figsize=(10,5))
    hist(fmaxt,normed=False,bins=bins,histtype="step", linewidth= 2,alpha= 1, linestyle="-",color='black',align='mid',label='Full sample', ax=ax)
    hist(pmaxt,normed=False,bins=bins,histtype="step", linewidth= 2,alpha= 1, linestyle="-",color='red',align='mid',label='Initial sample', ax=ax)
    hist(fhb,normed=False,bins=bins,histtype="step", linewidth= 2,alpha= 1, linestyle="-",color='black',align='mid',label='Full sample-HiBAL', ax=ax1)
    hist(flb,normed=False,bins=bins,histtype="step", linewidth= 2,alpha= 1, linestyle="-",color='red',align='mid',label='Full sample-LoBAL', ax=ax1) 
    hist(phb,normed=False,bins=bins,histtype="step", linewidth= 2,alpha= 1, linestyle=":",color='black',align='mid',label='Initial sample-HiBAL', ax=ax1)
    hist(plb,normed=False,bins=bins,histtype="step", linewidth= 2,alpha= 1, linestyle=":",color='red',align='mid',label='Initial sample-LoBAL', ax=ax1) 
    ax.set(xlabel='Rest-frame time (years)', ylabel='Histogram',title='Maximum probed timescales')
    ax1.set(xlabel='Rest-frame time (years)', ylabel='Histogram',title='Maximum probed timescales')
    ax.legend(loc=1)
    ax1.legend(loc=1)
    ax.grid()
    ax1.grid()
    fig.tight_layout()
    fig.savefig('SDSSIV_BALQSO_restframetimescales.jpg')
    #plt.show()

def plot_redshift():
    fsm = np.genfromtxt('Master_full_sample_redshift.txt',names=['name','ra','dec','Z_VI','tag'],dtype=('|S25',float,float,float,'|S25'))#fits.open('Fullsample_3028_BAL_sources_crossmatch-DR12Q.fits')[1].data
    mg = np.where((fsm['Z_VI'] > 0.28) & (fsm['Z_VI'] < 2.28))[0] 
    siiv = np.where((fsm['Z_VI'] > 1.96) & (fsm['Z_VI'] < 5.55))[0] 
    civ = np.where((fsm['Z_VI'] > 1.68) & (fsm['Z_VI'] < 4.93))[0]
    al = np.where((fsm['Z_VI'] > 1.23) & (fsm['Z_VI'] < 3.93))[0]
    xx = np.where(fsm['Z_VI'] > 0)[0]
    bins = np.arange(0,6.5,0.25)

    fig,ax=plt.subplots(figsize=(5,5))
    hist(fsm['Z_VI'][xx],bins=bins,normed=0,histtype='step',color='black',lw=3,ax=ax)
    ax.set(xlabel='Z_VI',ylabel='Histogram',title='Redshift: full sample')
    xlim,ylim =ax.get_xlim(),ax.get_ylim()
    ax.text(xlim[0]+0.4*(xlim[1]-xlim[0]),ylim[1]-0.25*(ylim[1]-ylim[0]),'Z_VI between 1.96 and 5.55 :'+str(len(siiv)))
    ax.text(xlim[0]+0.4*(xlim[1]-xlim[0]),ylim[1]-0.2*(ylim[1]-ylim[0]),'Z_VI between 1.68 and 4.93 :'+str(len(civ)))
    ax.text(xlim[0]+0.4*(xlim[1]-xlim[0]),ylim[1]-0.1*(ylim[1]-ylim[0]),'Z_VI between 0.28 and 2.28 :'+str(len(mg)))
    ax.text(xlim[0]+0.4*(xlim[1]-xlim[0]),ylim[1]-0.15*(ylim[1]-ylim[0]),'Z_VI between 1.23 and 3.93 :'+str(len(al)))
    ax.grid()
    fig.tight_layout()
    fig.savefig('SDSSIV_BALs_fullsample_redshift.jpg')
    #plt.show()


def LoHIBAL_plot_redshift():
    ofsm = np.genfromtxt('Master_full_sample_redshift.txt',names=['name','ra','dec','Z_VI','tag'],dtype=('|S25',float,float,float,'|S25'))#fits.open('Fullsample_3028_BAL_sources_crossmatch-DR12Q.fits')[1].data
    y=np.where(ofsm['Z_VI'] > 0)[0]
    fsm= ofsm[y]
    psm = np.genfromtxt('Master_gibson_sample_redshift.txt',names=['name','ra','dec','Z_VI','tag'],dtype=('|S25',float,float,float,'|S25'))#fits.open('Initialsample_2005_BAL_sources_crossmatch-DR12Q.fits')[1].data
    Flbal = fits.open('LoBAL_sample_fullsample.fits')[1].data
    Ilbal = fits.open('LoBAL_sample_initialsample.fits')[1].data
    Flzvi = []; Fhzvi=[] ; Ilzvi=[] ; Ihzvi=[]
    bins = np.arange(0,6.5,0.25)
    for i in range(len(fsm)):
        if fsm['name'][i] in Flbal['col1']:
            Flzvi.append(fsm['Z_VI'][i])
        else:
            Fhzvi.append(fsm['Z_VI'][i])
    for ii in range(len(psm)):
        if psm['name'][ii] in Ilbal['col1']:
            Ilzvi.append(psm['Z_VI'][ii])
        else:
            Ihzvi.append(psm['Z_VI'][ii])
    Flzvi = np.array(Flzvi); Fhzvi=np.array(Fhzvi) ; Ilzvi=np.array(Ilzvi) ; Ihzvi=np.array(Ihzvi)
    mg = np.where((fsm['Z_VI'] > 0.28) & (fsm['Z_VI'] < 2.28))[0] 
    siiv = np.where((fsm['Z_VI'] > 1.96) & (fsm['Z_VI'] < 5.55))[0] 
    civ = np.where((fsm['Z_VI'] > 1.68) & (fsm['Z_VI'] < 4.93))[0]
    al = np.where((fsm['Z_VI'] > 1.23) & (fsm['Z_VI'] < 3.93))[0] 
    fig,ax=plt.subplots(figsize=(5,5))
    hist(Fhzvi,bins=bins,normed=0,histtype='step',color='black',lw=3,ax=ax,label='Full sample-HiBAL: '+str(len(Fhzvi)))
    hist(Ihzvi,bins=bins,normed=0,histtype='step',color='red',lw=3,ax=ax,label='Initial sample-HiBAL: '+str(len(Ihzvi)))
    hist(Flzvi,bins=bins,normed=0,histtype='step',color='black',ls=':',lw=3,ax=ax,label='Full sample-LoBAL: '+str(len(Flzvi)))
    hist(Ilzvi,bins=bins,normed=0,histtype='step',color='red',ls=':',lw=3,ax=ax,label='Initial sample-LoBAL: '+str(len(Ilzvi)))
    ax.set(xlabel='Z\_VI',ylabel='Histogram',title='Redshift distribution')
    xlim,ylim =ax.get_xlim(),ax.get_ylim()
    #ax.text(xlim[0]+0.5*(xlim[1]-xlim[0]),ylim[1]-0.25*(ylim[1]-ylim[0]),'Full sample-HiBAL :'+str(len(Fhzvi)))
    #ax.text(xlim[0]+0.5*(xlim[1]-xlim[0]),ylim[1]-0.2*(ylim[1]-ylim[0]),'Full sample-LoBAL :'+str(len(Flzvi)))
    #ax.text(xlim[0]+0.5*(xlim[1]-xlim[0]),ylim[1]-0.1*(ylim[1]-ylim[0]),'Initial sample-HiBAL :'+str(len(Ihzvi)))
    #ax.text(xlim[0]+0.5*(xlim[1]-xlim[0]),ylim[1]-0.15*(ylim[1]-ylim[0]),'Initial sample-LoBAL :'+str(len(Ilzvi)))
    #ax.grid()
    ax.legend(loc=1)
    ax.grid()
    fig.tight_layout()
    fig.savefig('SDSSIV_LoHiBALs_redshift.jpg')
    #plt.show()



def plot_spectra():
    text_font = {'fontname':'Times New Roman', 'size':'14'}
    pp = PdfPages('Fullanalysis_BAL_plot_spectra.pdf') 
    specdir = 'SDSSIV_BALdata'
    linelist = np.genfromtxt('/Users/vzm83/Proposals/linelist_speccy.txt',usecols=(0,1,2),dtype=('|S10',float,'|S5'),names=True)

    data = np.load('MasterList_Plate-MJD-Fiber.npz')
    #fsm = fits.open('Fullsample_3028_BAL_sources_crossmatch-DR12Q.fits')[1].data
    fsm = np.genfromtxt('Master_full_sample_redshift.txt',names=['name','ra','dec','Z_VI','tag'],dtype=('|S25',float,float,float,'|S25'))
    for i in range(len(data['name'])):
    #for i in range(10):
        plates = data['plate'][i]
        mjds = data['mjd'][i]
        fibers = data['fiber'][i]
        print data['name'][i]
        xx=np.where(fsm['name'] == data['name'][i])[0]
        if len(xx)>0:
            zvi = fsm['Z_VI'][xx[0]]
            print xx,data['name'][i],fsm['name'][xx]
            minflux =[] ;maxflux=[]
            fig,ax=plt.subplots(figsize=(15,8))
            for j in range(len(plates)):
                if plates[j] >=10000:
                    PMF = 'spec-{0:05d}-{1:05d}-{2:04d}.fits'.format(plates[j],mjds[j],fibers[j])
                else:
                    PMF = 'spec-{0:04d}-{1:05d}-{2:04d}.fits'.format(plates[j],mjds[j],fibers[j])
                print PMF
                if os.path.isfile(os.path.join(specdir,PMF)):
                    data1 = fits.open(os.path.join(specdir,PMF))[1].data
                    minflux.append(np.median((data1['flux']*(data1['and_mask'] == 0)).copy()))
                    maxflux.append(np.std((data1['flux']*(data1['and_mask'] == 0)).copy()))
                    ax.plot(10**data1['loglam']/(1.0+zvi),gaussian_filter1d(data1['flux'],2),color=plt.cm.RdYlBu(j*300),alpha=0.5,label=PMF.split('.')[0][5:])
                    ax.plot(10**data1['loglam']/(1.0+zvi),1.0/np.sqrt(data1['ivar']),color=plt.cm.RdYlBu(j*300),alpha=0.1)
                
                
                if j == len(plates)-1:
                    if (len(maxflux) > 0):
                        ax.set(xlabel='Rest wavelength ($\AA$)',ylabel='Flux',ylim=(-2,max(minflux)+3.0*max(maxflux)), title=data['name'][i])
                    xlim,ylim=ax.get_xlim(),ax.get_ylim()
                    #string1 = 'SDSS J{0}\tZ\_VI: {1:4.4f}\tN$\_{{spec}}$: {2}'.format(fsm['SDSS_NAME'][xx[0]],zvi,len(plates))
                    string1 = 'RA: {0:5.4f}\t DEC: {1:5.4f}  \tZ\_VI: {2:4.4f}\tN$\_{{spec}}$: {3}'.format(fsm['ra'][xx[0]],fsm['dec'][xx[0]],zvi,len(plates))
                    print string1,xlim,ylim
                    ax.text(xlim[0]+0.05*(xlim[1] - xlim[0]),ylim[1]-0.05*(ylim[1] - ylim[0]), string1,fontsize=18)
                
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
    
        ax.legend(loc=1)
        fig.tight_layout()
        fig.savefig(pp,format='pdf')
        #sdlfj=raw_input()
    pp.close()
                

def createBigFits():
    specdir = 'SDSSIV_BALdata'
    data = np.load('MasterList_Plate-MJD-Fiber.npz')
    fsm = fits.open('Fullsample_3028_BAL_sources_crossmatch-DR12Q.fits')[1].data
    balflux = []; balivar=[] ; balmask =[] ; balmap = []; balwave = []
    #for i in range(len(data['name'])):
    for i in range(10):
        plates = data['plate'][i]
        mjds = data['mjd'][i]
        fibers = data['fiber'][i]
        print data['name'][i]
        xx=np.where(fsm['col1'] == data['name'][i])[0]
        if len(xx)>0:
            zvi = fsm['Z_VI'][xx[0]]
        for j in range(len(plates)):
                if plates[j] >=10000:
                    PMF = 'spec-{0:05d}-{1:05d}-{2:04d}.fits'.format(plates[j],mjds[j],fibers[j])
                else:
                    PMF = 'spec-{0:04d}-{1:05d}-{2:04d}.fits'.format(plates[j],mjds[j],fibers[j])
                print PMF
                if os.path.isfile(os.path.join(specdir,PMF)):
                    data1 = fits.open(os.path.join(specdir,PMF))[1].data
                    balflux.append(data1['flux'])
                    balwave.append(10**data1['loglam'])
                    balivar.append(data1['ivar'])
                    balmask.append(data1['and_mask'])
                    tmpmap= fits.open(os.path.join(specdir,PMF))[2].data
                    tmpmap['z'] = zvi
                    balmap.append(tmpmap)
    #balflux = np.array(balflux)
    #balwave = np.array(balwave)
    #balivar = np.array(balivar)
    #balmask = np.array(balmask)
    #balmap = np.array(balmap)
    print balflux 
    bhdu0=fits.PrimaryHDU(balflux)
    bhdu1 = fits.ImageHDU(balivar)
    bhdu2 = fits.ImageHDU(balmask)
    bhdu3 = fits.BinTableHDU(data=balmap)
    bhdu4 = fits.ImageHDU(balwave)
    bhdulist = fits.HDUList([bhdu0,bhdu1,bhdu2,bhdu3,bhdu4])

    bhdulist.writeto('SDSSIV_BALQSO_BigFits.fits')

def plot_BALfootprint():
    print 'Begin'
    sp = fits.open('spAll-v5_13_0.fits')[1].data
    print 'Reading spAll-new'
    #data = fits.open('spAll-v5_10_10.fits')[1].data
    print 'Reading spAll-old'
    #yy=np.where(data['PROGRAMNAME'] == 'sequels')[0]
    #sequels = data[yy]
    elgngc = np.where(sp['PROGRAMNAME'] == 'ELG_NGC')[0]
    elgsgc = np.where(sp['PROGRAMNAME'] == 'ELG_SGC')[0]

    nn=pd.read_csv('NGC_bounds',names=['ra','dec'])
    dd=pd.read_csv('SGC_bounds',names=['ra','dec'])

    npolygon = [[ngcp[0],ngcp[1]] for ngcp in zip(nn['ra'],nn['dec'])]
    dpolygon = [[sgcp[0],sgcp[1]] for sgcp in zip(dd['ra'],dd['dec'])]
    
    
    


    radec = np.genfromtxt('Master_initial_sample.txt',usecols=(1,2),names=['ra','dec'],skip_header=1)
    #radec = np.genfromtxt('Master_gibson_2005_targets_cor.txt',usecols=(1,2),names=['ra','dec'],skip_header=1)
    ra=radec['ra'] ; dec = radec['dec']
    
    mcount = np.genfromtxt('Fullsample_3epoch_counts.txt',names=['name','ra','dec','sdss1','sdss2','sdss3'],dtype=('|S15',float,float,int,int,int))
    #mcount = np.genfromtxt('Fullsample_3epoch_counts_2005.txt',names=['name','ra','dec','sdss1','sdss2','sdss3'],dtype=('|S15',float,float,int,int,int))
    xx=np.where(mcount['sdss3']>0)[0]
    
    sqbe = fits.open('sequelsVARBALS_before_eBOSS56777_fullsample.fits')[1].data
    #sqbe = fits.open('sequelsVARBALS_before_eBOSS56777_initialsample.fits')[1].data
    mcountra = [ s for s in mcount['ra'][xx] ]
    mcountdec = [ s for s in mcount['dec'][xx] ]
    sqbera = [s for s in sqbe['RA']]
    sqbedec = [s for s in sqbe['dec']]
    
    ra1,dec1=np.loadtxt('boss_survey_outline_points_sorted_ra-60_300_NGC.txt').T
    ra2,dec2=np.loadtxt('boss_survey_outline_points_sorted_ra-60_300_SGC.txt').T
    eboss = np.where(sp['PROGRAMNAME'] == 'eboss')[0]
    speboss = sp[eboss]
    q = np.where(speboss['CLASS'] == 'QSO')[0]
    #q1 = np.where(sequels['CLASS'] == 'QSO')[0]
    print 'starting to plot'
    fig,ax1 = plt.subplots(figsize=(15,8))
    ax1.plot(ra1,dec1,'-',color='black',alpha=0.5)
    ax1.plot(ra2,dec2,'-',color='black',alpha=0.5)
    ax1.plot(nn['ra'],nn['dec'],'-',color='green',alpha=0.8,lw=3,label='eBOSS final Footprint')
    ax1.plot(dd['ra'],dd['dec'],'-',color='green',alpha=0.8,lw=3,label='')
    
 
    #ax1.plot(negativeRAs(speboss['RA'][q]),speboss['DEC'][q],'.',markersize=3,color='red',alpha=0.1)
    #ax1.plot(negativeRAs(sequels['RA'][q1]),sequels['DEC'][q1],'x',markersize=3,color='red',alpha=0.1)
    #ax1.plot(negativeRAs(sp['RA'][elgngc]),sp['DEC'][elgngc],'.',markersize=3,color='red',alpha=0.1)
    #ax1.plot(negativeRAs(sp['RA'][elgsgc]),sp['DEC'][elgsgc],'.',markersize=3,color='red',alpha=0.1)

    ax1.plot(negativeRAs(ra),dec,'x',markersize=2,color='blue',label='Targets')
    ax1.plot(negativeRAs(mcount['ra'][xx]),mcount['dec'][xx],'o',markersize=3,color='black',fillstyle='none',label='Observed')
    ax1.plot(negativeRAs(sqbera),sqbedec,'o',markersize=3,color='black',fillstyle='none',label='')
    ax1.set(xlabel='RA',ylabel='DEC',title='eBOSS QSO Full sample footprint')
    #ax1.set(xlabel='RA',ylabel='DEC',title='eBOSS QSO Initial sample footprint')
    ax1.set_xlim(-55,300)
    ax1.set_ylim(-15,100)
    ax1.legend(loc=1)
    ax1.grid()
    xlim,ylim=ax1.get_xlim(),ax1.get_ylim()
    ax1.text(xlim[0]+0.1*(xlim[1]-xlim[0]),ylim[1]-0.1*(ylim[1]-ylim[0]),'1908 targets inside final footprint',fontsize=18)
    ax1.text(xlim[0]+0.1*(xlim[1]-xlim[0]),ylim[1]-0.15*(ylim[1]-ylim[0]),'1757 observed inside final footprint',fontsize=18)
    #ax1.text(xlim[0]+0.1*(xlim[1]-xlim[0]),ylim[1]-0.1*(ylim[1]-ylim[0]),'1138 targets inside final footprint',fontsize=18)
    #ax1.text(xlim[0]+0.1*(xlim[1]-xlim[0]),ylim[1]-0.15*(ylim[1]-ylim[0]),'1039 observed inside final footprint',fontsize=18)

    fig.tight_layout()
    #print 'Interactive window should open now'
    coords = []
    NGC_targets = [ray_tracing_method(point[0], point[1], npolygon) for point in zip(negativeRAs(ra),dec)]
    SGC_targets = [ray_tracing_method(point[0], point[1], dpolygon) for point in zip(negativeRAs(ra),dec)]

    NGC_observed = [ray_tracing_method(point[0], point[1], npolygon) for point in zip(negativeRAs(mcountra),mcountdec)]
    NSQBE_observed = [ray_tracing_method(point[0], point[1], npolygon) for point in zip(negativeRAs(sqbera),sqbedec)]
    SGC_observed = [ray_tracing_method(point[0], point[1], dpolygon) for point in zip(negativeRAs(mcountra),mcountdec)]
    SSQBE_observed = [ray_tracing_method(point[0], point[1], dpolygon) for point in zip(negativeRAs(sqbera),sqbedec)]
    NGC_targets = np.array(NGC_targets)
    SGC_targets = np.array(SGC_targets)
    NGC_observed = np.array(NGC_observed)
    SGC_observed = np.array(SGC_observed)
    NSQBE_observed = np.array(NSQBE_observed)
    SSQBE_observed = np.array(SSQBE_observed)

    print NGC_targets
    print SGC_targets
    print NGC_observed
    print SGC_observed
    print len(np.where(NGC_targets == True)[0]) + len(np.where(SGC_targets == True)[0]), len(np.where(NGC_targets == True)[0]), len(np.where(SGC_targets == True)[0])
    print len(np.where(NGC_observed == True)[0]) + len(np.where(SGC_observed == True)[0])+len(np.where(NSQBE_observed == True)[0]) + len(np.where(SSQBE_observed == True)[0]), len(np.where(NGC_observed == True)[0])+len(np.where(NSQBE_observed == True)[0]) , len(np.where(SGC_observed == True)[0]) + +len(np.where(SSQBE_observed == True)[0]) 
    # Call click func
    #cid = fig.canvas.mpl_connect('button_press_event', onclick)

    #plt.show(1)
    
    fig.savefig('SDSSIV_eBOSS_VARBAL_footprint_fullsample.jpg')
    #fig.savefig('SDSSIV_eBOSS_VARBAL_footprint_initialsample.jpg')
    print len(ra)



    #plt.show()

def checkifSpectra():
    specdir = 'SDSSIV_BALdata'
    data = np.load('MasterList_Plate-MJD-Fiber.npz')
    count = 0
    for i in range(len(data['name'])):
        plates = data['plate'][i]
        mjds = data['mjd'][i]
        fibers = data['fiber'][i]
        for j in range(len(plates)):
            if plates[j] >=10000:
                PMF = 'spec-{0:05d}-{1:05d}-{2:04d}.fits'.format(plates[j],mjds[j],fibers[j])
            else:
                PMF = 'spec-{0:04d}-{1:05d}-{2:04d}.fits'.format(plates[j],mjds[j],fibers[j])
            
            filename = os.path.join(specdir,PMF)
            if not os.path.isfile(filename):
                print 'File not found {0:10.5f}\t{1:10.5f}\t{2}\t{3}'.format(data['ra'][i],data['dec'][i], PMF,i)
                count +=1
    print '--'*51
    print 'Total :', count


def VI_stats():
    files = glob.glob('/Users/vivekm/VI_out/*.npz')
    total = len(files)
    keys = ['_hibal_','_lobal_','_felobal_','_manynarrowsys_','_j0300analog_','_reddened_','_ironemitter_','_redshifted_','_accn_','_stable_','_emergence_','_disappear_','_lo2hibal_','_hi2lobal_','_lo2febal_','_fe2lobal_','_coordvar_','_xtreme_','_CIV_','_SiiV_','_AlIII_','_MgII_','_FeII_','_FeIII_','_AlII_']
    Llist =['Must get a follow up','May be interesting','Typical Source']
    
    print '--'*51
    for item in keys:
        count =0    
        for f in files:
            data = np.load(f)
            if data[item] == True:
                count +=1
        print '{0:30}\t\t{1}\t\t{2:3.4f}'.format(item,count,float(count)/total*100.)
    print '--'*51
    for ll in Llist:
        count =0    
        for f in files:
            data = np.load(f)
            if data['_LIST_'] == ll:
                count +=1
        print '{0:30}\t\t{1}\t\t{2:3.4f}'.format(ll,count,float(count)/total*100.)
    print '--'*51
def velocity(wave,ion_wave,z):
    #print wave
    #print z,ion_wave
    c = 299792.0# in km/s
    vel =np.zeros(len(wave))
    zabs = np.zeros(len(wave))
    for i in range(len(wave)):
        #print i,wave[i],zabs[i],vel[i]
        zabs[i] = (wave[i]/ion_wave) -1.0
        vel[i] = -((((1+z)/(1+zabs[i]))**2-1.0)/(((1+z)/(1+zabs[i]))**2+1))*c
    return vel

def plot_vel_felobals():
    n_v = 1240.81
    si_iv = 1399.8
    c_iv = 1549.48
    al_iii = 1857.4
    mg_ii = 2799.117
    p_v = 1122.9925#1128.008#1117.977 #
    c_ii = 1334.53
    c_iis = 1335.70
    si_iia = 1304.37
    si_iias = 1309.28
    si_iib = 1526.71
    si_iibs = 1533.43
    si_iic = 1808.01
    si_iics = 1816.93
    fe_ii = 2344.21
    fe_iis = 2396.36

    lines = {'n v' : 1240.81,
            'si iv': 1399.8,
            'c iv': 1549.48,
            'al iii': 1857.4,
            'mg ii' : 2799.117,
            'p v' : 1122.9925,
            'c ii' : 1334.53,
            'c ii*' : 1335.70,
            'si iia' : 1304.37,
            'si iia*' : 1309.28,
            'si iib' : 1526.71,
            'si iib*' :1533.43,
            'si iic' : 1808.01,
            'si iic*' : 1816.93,
            'fe ii' : 2344.21,
            'fe ii*' : 2396.3}

    felobals = pd.read_csv('felobals_list.txt',sep='\t',names=['name'])
    text_font = {'fontname':'Times New Roman', 'size':'14'}
    pp = PdfPages('Fullanalysis_BAL_plot_spectra.pdf') 
    specdir = 'BALQSOs_Spectra'#'SDSSIV_BALdata'
    linelist = np.genfromtxt('/Users/vivekm/Proposals/linelist_speccy.txt',usecols=(0,1,2),dtype=('|S10',float,'|S5'),names=True)

    data = np.load('MasterList_Plate-MJD-Fiber.npz')
    #fsm = fits.open('Fullsample_3028_BAL_sources_crossmatch-DR12Q.fits')[1].data
    fsm = np.genfromtxt('Master_full_sample_redshift.txt',names=['name','ra','dec','Z_VI','tag'],dtype=('|S25',float,float,float,'|S25'))
    #for i in range(len(felobals['name'])):
    for i in 11+np.arange(2):
        print '--'*21
        x0 = np.where(data['name'] == felobals['name'][i])[0][0]
        plates = data['plate'][x0]
        mjds = data['mjd'][x0]
        fibers = data['fiber'][x0]
        print felobals['name'][i],plates
        xx=np.where(fsm['name'] == felobals['name'][i])[0]
        print 'XX:',xx
        if len(xx)>0:
            zvi = fsm['Z_VI'][xx[0]]
            print 'Inside XX:',xx,felobals['name'][i],fsm['name'][xx]
            minflux =[] ;maxflux=[]
            fig=plt.figure( figsize=(15, 8))
            nullfmt = NullFormatter() 

            dim_ax = [0.05, 0.8, 0.9, 0.175]
            
            dim_vel0 = [0.1, 0.5, 0.25, 0.225]
            dim_vel0a = [0.1, 0.275, 0.25, 0.225]
            dim_vel0b = [0.1, 0.05, 0.25, 0.225]
    


            dim_vel1 = [0.4, 0.5, 0.25, 0.225]
            dim_vel1a = [0.4, 0.275, 0.25, 0.225]
            dim_vel1b = [0.4, 0.05, 0.25, 0.225]

            dim_vel2 = [0.7, 0.5, 0.25, 0.225]
            dim_vel3 = [0.7, 0.275, 0.25, 0.225]
            dim_vel4 = [0.7, 0.05, 0.25, 0.225]
    
            ax = plt.axes(dim_ax)

            
            ax_vel1 = plt.axes(dim_vel1)
            ax_vel1a = plt.axes(dim_vel1a)
            ax_vel1b = plt.axes(dim_vel1b)
            
            ax_vel0 = plt.axes(dim_vel0)
            ax_vel0a = plt.axes(dim_vel0a)
            ax_vel0b = plt.axes(dim_vel0b)

            ax_vel2 = plt.axes(dim_vel2)
            ax_vel3 = plt.axes(dim_vel3)
            ax_vel4 = plt.axes(dim_vel4)
            
            for j in range(len(plates)):
                print 'Inside jj:',plates[j]
                if plates[j] >=10000:
                    PMF = 'spec-{0:05d}-{1:05d}-{2:04d}.fits'.format(plates[j],mjds[j],fibers[j])
                else:
                    PMF = 'spec-{0:04d}-{1:05d}-{2:04d}.fits'.format(plates[j],mjds[j],fibers[j])
                print PMF
                if os.path.isfile(os.path.join(specdir,PMF)):

                    data1 = fits.open(os.path.join(specdir,PMF))[1].data
                    print 'Inside path:',min(data1['loglam']/(1.0+zvi))
                    minflux.append(np.median((data1['flux']*(data1['and_mask'] == 0)).copy()))
                    maxflux.append(np.std((data1['flux']*(data1['and_mask'] == 0)).copy()))
                    sc = np.median(gaussian_filter1d(data1['flux'],10)[np.where(gaussian_filter1d(data1['flux'],10)  > np.percentile(gaussian_filter1d(data1['flux'],10),80))])
                    ax.set(xlim=(1200,2000))
                    ax.plot(10**data1['loglam']/(1.0+zvi),gaussian_filter1d(data1['flux'],2),color=plt.cm.RdYlBu(j*300),alpha=0.5,label=PMF.split('.')[0][5:])
                    
                    ax_vel1.plot(velocity(10**data1['loglam']/(1.0+zvi),si_iv,0),gaussian_filter1d(data1['flux'],1)/sc,lw=2,color=plt.cm.RdYlBu(j*300),alpha=0.5)
                    ax_vel1a.plot(velocity(10**data1['loglam']/(1.0+zvi),si_iib,0),gaussian_filter1d(data1['flux'],1)/sc,lw=2,color=plt.cm.RdYlBu(j*300),alpha=0.5)
                    ax_vel1b.plot(velocity(10**data1['loglam']/(1.0+zvi),si_iibs,0),gaussian_filter1d(data1['flux'],1)/sc,lw=2,color=plt.cm.RdYlBu(j*300),alpha=0.5)

                    ax_vel0.plot(velocity(10**data1['loglam']/(1.0+zvi),c_ii,0),gaussian_filter1d(data1['flux'],1)/sc,lw=2,color=plt.cm.RdYlBu(j*300),alpha=0.5)
                    ax_vel0a.plot(velocity(10**data1['loglam']/(1.0+zvi),c_iis,0),gaussian_filter1d(data1['flux'],1)/sc,lw=2,color=plt.cm.RdYlBu(j*300),alpha=0.5)
                    ax_vel0b.plot(velocity(10**data1['loglam']/(1.0+zvi),c_iv,0),gaussian_filter1d(data1['flux'],1)/sc,lw=2,color=plt.cm.RdYlBu(j*300),alpha=0.5)



                    ax_vel2.plot(velocity(10**data1['loglam']/(1.0+zvi),si_iic,0),gaussian_filter1d(data1['flux'],1)/sc,lw=2,color=plt.cm.RdYlBu(j*300),alpha=0.5)
                    ax_vel3.plot(velocity(10**data1['loglam']/(1.0+zvi),si_iics,0),gaussian_filter1d(data1['flux'],1)/sc,lw=2,color=plt.cm.RdYlBu(j*300),alpha=0.5)
                    ax_vel4.plot(velocity(10**data1['loglam']/(1.0+zvi),fe_iis,0),gaussian_filter1d(data1['flux'],1)/sc,lw=2,color=plt.cm.RdYlBu(j*300),alpha=0.5)
                    ax_vel1.set_xlim(-5000,3000)    
                    ax_vel1a.set_xlim(-5000,3000)    
                    ax_vel1b.set_xlim(-5000,3000)    
                    ax_vel2.set_xlim(-5000,3000)    
                    
                    ax_vel0.set_xlim(-5000,3000)    
                    ax_vel0a.set_xlim(-5000,3000)    
                    ax_vel0b.set_xlim(-5000,3000)    
                    ax_vel3.set_xlim(-5000,3000)    
                    ax_vel4.set_xlim(-5000,3000)    
                    ax_vel0.set_ylabel('Norm Flux')
                    ax_vel0a.set_ylabel('Norm Flux')
                    ax_vel0b.set_ylabel('Norm Flux')
                    ax_vel1.set_ylabel('Norm Flux')
                    ax_vel1a.set_ylabel('Norm Flux')
                    ax_vel1b.set_ylabel('Norm Flux')
                    ax_vel2.set_ylabel('Norm Flux')
                    ax_vel3.set_ylabel('Norm Flux')
                    ax_vel4.set_ylabel('Norm Flux')
                    #for vl in [-4000,-3000,-2000,-1500,-1000,-500,0]:
                    for vl in [-2440,0]:
                        ax_vel0.axvline(vl,ls='--',lw=1,color='black')
                        ax_vel0a.axvline(vl,ls='--',lw=1,color='black')
                        ax_vel0b.axvline(vl,ls='--',lw=1,color='black')
                        ax_vel1.axvline(vl,ls='--',lw=1,color='black')
                        ax_vel1a.axvline(vl,ls='--',lw=1,color='black')
                        ax_vel1b.axvline(vl,ls='--',lw=1,color='black')
                        ax_vel2.axvline(vl,ls='--',lw=1,color='black')
                        ax_vel3.axvline(vl,ls='--',lw=1,color='black')
                        ax_vel4.axvline(vl,ls='--',lw=1,color='black')
                    ax_vel0.xaxis.set_major_formatter(nullfmt)
                    ax_vel0a.xaxis.set_major_formatter(nullfmt)
                    ax_vel1.xaxis.set_major_formatter(nullfmt)
                    ax_vel1a.xaxis.set_major_formatter(nullfmt)
                    ax_vel2.xaxis.set_major_formatter(nullfmt)
                    ax_vel3.xaxis.set_major_formatter(nullfmt)
                    ax_vel0.text(-4000,0.13,'C II  1334.53',fontsize=18)
                    ax_vel0a.text(-4000,0.13,'C II* 1335.70',fontsize=18)
                    ax_vel0b.text(-4000,0.13,'C IV 1550',fontsize=18)
                    ax_vel1.text(-4000,0.13,'Si IV 1393.',fontsize=18)
                    ax_vel1a.text(-4000,0.13,'Si IIb 1526.71',fontsize=18)
                    ax_vel1b.text(-4000,0.13,'Si IIb* 1533.43',fontsize=18)
                    ax_vel2.text(-4000,0.13,'Si IIc 1808.01',fontsize=18)
                    ax_vel3.text(-4000,0.13,'Si IIc* 1816.93',fontsize=18)
                    ax_vel4.text(-4000,0.13,'Fe II* 2396.36',fontsize=18)
                    ax_vel0.set_ylim(-0.1,1.25)
                    ax_vel0a.set_ylim(-0.1,1.25)
                    ax_vel0b.set_ylim(-0.1,1.25)
                    ax_vel1.set_ylim(-0.1,1.25)
                    ax_vel1a.set_ylim(-0.1,1.25)
                    ax_vel1b.set_ylim(-0.1,1.25)
                    ax_vel2.set_ylim(-0.1,1.25)
                    ax_vel3.set_ylim(-0.1,1.25)
                    ax_vel4.set_ylim(-0.1,1.25)


                   
                if j==0:
                    #if (len(maxflux) > 0):
                    ax.set(xlabel='Rest wavelength ($\AA$)',ylabel='Flux',ylim=(-2,max(minflux)+3.0*max(maxflux)), title=felobals['name'][i])
                    xlim,ylim=ax.get_xlim(),ax.get_ylim()
                    #string1 = 'SDSS J{0}\tZ\_VI: {1:4.4f}\tN$\_{{spec}}$: {2}'.format(fsm['SDSS_NAME'][xx[0]],zvi,len(plates))
                    string1 = 'RA: {0:5.4f}\t DEC: {1:5.4f}  \tZ\_VI: {2:4.4f}\tN$\_{{spec}}$: {3}'.format(fsm['ra'][xx[0]],fsm['dec'][xx[0]],zvi,len(plates))
                    print string1,xlim,ylim
                    ax.text(xlim[0]+0.05*(xlim[1] - xlim[0]),ylim[1]-0.05*(ylim[1] - ylim[0]), string1,fontsize=18)
                
                    obslambda = linelist['lambda']#*(1.+zvi)
                    x = np.where((obslambda > xlim[0]) & (obslambda < xlim[1]))[0]
	            plotlambda = obslambda[x]
	            plotname = linelist['Name'][x]
	            plota_e = linelist['a_e'][x]
	            #print plotlambda
	            #for k in range(len(plotlambda)):
	            #    if plota_e[k].strip() == 'Abs.' : 
	            #	    ax.axvline(x=plotlambda[k], color='lawngreen', linestyle=':')
	            #	    #ax.text(plotlambda[k],ylim[0]+0.75*(ylim[1]-ylim[0]),plotname[k],color='Orange',ha='center',rotation=90,**text_font)
	            #    else :
	            #	    ax.axvline(x=plotlambda[k], color='lightblue', linestyle=':')
	            	    #ax.text(plotlambda[k],ylim[0]+0.75*(ylim[1]-ylim[0]),plotname[k],color='Brown',ha='center',rotation=90,**text_font)
                    for k in lines:
                        ax.axvline(x=lines[k], color='lawngreen', linestyle=':')
	            	ax.text(lines[k],ylim[0]+0.75*(ylim[1]-ylim[0]),k,color='Brown',ha='center',rotation=90,**text_font)


                ax.legend(loc=1)


                fig.tight_layout()
                fig.savefig('FeloBALs_velocityplot_{}_new.jpg'.format(felobals['name'][i]),format='jpg')
            plt.show()

        

def checkMags_RADec_finestructure():
    finestr =  pd.read_csv('felobals_finestructure_list.txt',sep='\t',names=['name'])
    fsm = np.genfromtxt('Master_full_sample_redshift.txt',names=['name','ra','dec','Z_VI','tag'],dtype=('|S25',float,float,float,'|S25'))
    data = fits.open('Fullsample_3028_BAL_sources_crossmatch-DR12Q.fits')[1].data
    for i in range(len(finestr['name'])):
        x0 = np.where(data['col1'] == finestr['name'][i])[0]
        xx=np.where(fsm['name'] == finestr['name'][i])[0]
        print finestr['name'][i],fsm['ra'][xx[0]],fsm['dec'][xx[0]],fsm['Z_VI'][xx[0]],data['PSFMAG'][x0]

    

def main():     

    print '-------'
    #merge2initialSamples()
    #getMulti_epochInfo()
    #plot_nepochs()
    #download_data(print_cmd=False)
    #hist_multiepochs()
    #cumulative_MJD()
    #plot_spPlate()
    #three_epochs_count()
    #restFrame_timescales()
    #plot_spectra()
    #createBigFits()
    #plot_BALfootprint()
    #LoHIBAL_plot_redshift()
    #checkifSpectra()
    #VI_stats()
    plot_vel_felobals()
    #checkMags_RADec_finestructure()

if __name__== "__main__":
    main()
