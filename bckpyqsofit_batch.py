import scipy as sp
import glob, os,sys,timeit
import matplotlib
import numpy as np
from PyQSOFit import QSOFit
from astropy.io import fits
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

snfile = np.genfromtxt('SN1700_estimates_Master_full_sample_sources.txt',names=['name','plate','mjd','fiber','Z_VI','sn'],dtype=('|S25',int,int,int,float,float))

path1='/Users/vivekm/SDSSIV_BALQSO/pyqsofit/'
path2='/Users/vivekm/SDSSIV_BALQSO/pyqsofit/test/data/result/'
path3='/Users/vivekm/SDSSIV_BALQSO/pyqsofit/test/data/QA_other/'
path4 = '/Users/vivekm/Softwares/sfddata-master/'

specdir = '/Users/vivekm/SDSSIV_BALQSO/BALQSOs_Spectra/'
for i in range(len(snfile)):

    plate = snfile['plate'][i]
    mjd = snfile['mjd'][i]
    fiber = snfile['fiber'][i]
    if plate >=10000:
        PMF = 'spec-{0:05d}-{1:05d}-{2:04d}.fits'.format(plate,mjd,fiber)
    else:
        PMF = 'spec-{0:04d}-{1:05d}-{2:04d}.fits'.format(plate,mjd,fiber)
    
    data=fits.open(os.path.join(specdir,PMF)) 
    lam=10**data[1].data['loglam']        # OBS wavelength [A]
    flux=data[1].data['flux']             # OBS flux [erg/s/cm^2/A]
    err=1./np.sqrt(data[1].data['ivar'])  # 1 sigma error
    z=snfile['Z_VI'][i]                # Redshift

    #Optinal
    ra=data[0].header['plug_ra']          # RA 
    dec=data[0].header['plug_dec']        # DEC
    plateid = data[0].header['plateid']   # SDSS plate ID
    mjd = data[0].header['mjd']           # SDSS MJD
    fiberid = data[0].header['fiberid']   # SDSS fiber ID


    # get data prepared 
    q = QSOFit(lam, flux, err, z, ra = ra, dec = dec, plateid = plateid, mjd = mjd, fiberid = fiberid, path = path1)
    
    start = timeit.default_timer()
    # do the fitting
    q.Fit(name = None,nsmooth = 1, and_or_mask = False, deredden = True, reject_badpix = False, wave_range = [1270.,7000.],\
          wave_mask =None, decomposition_host = True, Mi = None, npca_gal = 5, npca_qso = 20, \
          Fe_uv_op = True, poly = True, BC = False, rej_abs = True, initial_guess = None, MC = True, \
          n_trails = 5, linefit = True, tie_lambda = True, tie_width = True, tie_flux_1 = True, tie_flux_2 = True,\
          save_result = True, plot_fig = True,save_fig = True, plot_line_name = True, plot_legend = True, \
          dustmap_path = path4, save_fig_path = path3, save_fits_path = path2,save_fits_name = None)
    
    end = timeit.default_timer()
    print 'Fitting {0} finished in : {1} s ' .format(PMF,str(np.round(end-start)))
    # grey shade on the top is the continuum windiows used to fit.
    
