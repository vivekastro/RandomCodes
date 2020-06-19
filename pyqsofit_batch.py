import scipy as sp
import glob, os,sys,timeit
import matplotlib
import numpy as np
from PyQSOFit import QSOFit
from astropy.io import fits
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
from multiprocessing import Pool




def runqsopar():
    path='/Users/vzm83/SDSSIV_BALQSO/pyqsofit/'
    newdata = np.rec.array([(6564.61,'Ha',6400.,6800.,'Ha_br',3,5e-3,0.004,0.05,0.015,0,0,0,0.05),\
                            (6564.61,'Ha',6400.,6800.,'Ha_na',1,1e-3,5e-4,0.0017,0.01,1,1,0,0.002),\
                            (6549.85,'Ha',6400.,6800.,'NII6549',1,1e-3,2.3e-4,0.0017,5e-3,1,1,1,0.001),\
                            (6585.28,'Ha',6400.,6800.,'NII6585',1,1e-3,2.3e-4,0.0017,5e-3,1,1,1,0.003),\
                            (6718.29,'Ha',6400.,6800.,'SII6718',1,1e-3,2.3e-4,0.0017,5e-3,1,1,2,0.001),\
                            (6732.67,'Ha',6400.,6800.,'SII6732',1,1e-3,2.3e-4,0.0017,5e-3,1,1,2,0.001),\
                            
                            (4862.68,'Hb',4640.,5100.,'Hb_br',1,5e-3,0.004,0.05,0.01,0,0,0,0.01),\
                            (4862.68,'Hb',4640.,5100.,'Hb_na',1,1e-3,2.3e-4,0.0017,0.01,1,1,0,0.002),\
                            (4960.30,'Hb',4640.,5100.,'OIII4959c',1,1e-3,2.3e-4,0.0017,0.01,1,1,0,0.002),\
                            (5008.24,'Hb',4640.,5100.,'OIII5007c',1,1e-3,2.3e-4,0.0017,0.01,1,1,0,0.004),\
                            #(4960.30,'Hb',4640.,5100.,'OIII4959w',1,3e-3,2.3e-4,0.004,0.01,1,1,0,0.001),\
                            #(5008.24,'Hb',4640.,5100.,'OIII5007w',1,3e-3,2.3e-4,0.004,0.01,1,1,0,0.002),\
                            #(4687.02,'Hb',4640.,5100.,'HeII4687_br',1,5e-3,0.004,0.05,0.005,0,0,0,0.001),\
                            #(4687.02,'Hb',4640.,5100.,'HeII4687_na',1,1e-3,2.3e-4,0.0017,0.005,1,1,0,0.001),\
                            
                            #(3934.78,'CaII',3900.,3960.,'CaII3934',2,1e-3,3.333e-4,0.0017,0.01,99,0,0,-0.001),\
                            
                            #(3728.48,'OII',3650.,3800.,'OII3728',1,1e-3,3.333e-4,0.0017,0.01,1,1,0,0.001),\
                            
                            #(3426.84,'NeV',3380.,3480.,'NeV3426',1,1e-3,3.333e-4,0.0017,0.01,0,0,0,0.001),\
                            #(3426.84,'NeV',3380.,3480.,'NeV3426_br',1,5e-3,0.0025,0.02,0.01,0,0,0,0.001),\
                            
                            (2798.75,'MgII',2700.,2900.,'MgII_br',1,5e-3,0.004,0.05,0.0017,0,0,0,0.05),\
                            (2798.75,'MgII',2700.,2900.,'MgII_na',2,1e-3,5e-4,0.0017,0.01,1,1,0,0.002),\
                            
                        
                            (1908.73,'CIII',1700.,1970.,'CIII_br',2,5e-3,0.004,0.05,0.015,99,0,0,0.01),\
                            #(1908.73,'CIII',1700.,1970.,'CIII_na',1,1e-3,5e-4,0.0017,0.01,1,1,0,0.002),\
                            (1892.03,'CIII',1700.,1970.,'SiIII1892',1,2e-3,0.001,0.015,0.003,1,1,0,0.005),\
                            (1857.40,'CIII',1700.,1970.,'AlIII1857',1,2e-3,0.001,0.015,0.003,1,1,0,0.005),\
                            #(1816.98,'CIII',1700.,1970.,'SiII1816',1,2e-3,0.001,0.015,0.01,1,1,0,0.0002),\
                            #(1786.7,'CIII',1700.,1970.,'FeII1787',1,2e-3,0.001,0.015,0.01,1,1,0,0.0002),\
                            #(1750.26,'CIII',1700.,1970.,'NIII1750',1,2e-3,0.001,0.015,0.01,1,1,0,0.001),\
                            #(1718.55,'CIII',1700.,1900.,'NIV1718',1,2e-3,0.001,0.015,0.01,1,1,0,0.001),\
                            
                            (1549.06,'CIV',1500.,1700.,'CIV_br',1,5e-3,0.004,0.05,0.015,0,0,0,0.05),\
                            (1549.06,'CIV',1500.,1700.,'CIV_na',1,1e-3,5e-4,0.0017,0.01,1,1,0,0.002),\
                            (1640.42,'CIV',1500.,1700.,'HeII1640',1,1e-3,5e-4,0.0017,0.008,1,1,0,0.002),\
                            #(1663.48,'CIV',1500.,1700.,'OIII1663',1,1e-3,5e-4,0.0017,0.008,1,1,0,0.002),\
                            #(1640.42,'CIV',1500.,1700.,'HeII1640_br',1,5e-3,0.0025,0.02,0.008,1,1,0,0.002),\
                            #(1663.48,'CIV',1500.,1700.,'OIII1663_br',1,5e-3,0.0025,0.02,0.008,1,1,0,0.002),\
                            
                            (1402.06,'SiIV',1290.,1450.,'SiIV_OIV1',1,5e-3,0.002,0.05,0.015,1,1,0,0.05),\
                            (1396.76,'SiIV',1290.,1450.,'SiIV_OIV2',1,5e-3,0.002,0.05,0.015,1,1,0,0.05),\
                            (1335.30,'SiIV',1290.,1450.,'CII1335',1,2e-3,0.001,0.015,0.01,1,1,0,0.001),\
                            (1304.35,'SiIV',1290.,1450.,'OI1304',1,2e-3,0.001,0.015,0.01,1,1,0,0.001),\
                            
                            (1215.67,'Lya',1150.,1290.,'Lya_br',1,5e-3,0.004,0.05,0.02,0,0,0,0.05),\
                            (1215.67,'Lya',1150.,1290.,'Lya_na',1,1e-3,5e-4,0.0017,0.01,0,0,0,0.002)\
                            ],\
                         formats='float32,a20,float32,float32,a20,float32,float32,float32,float32,\
                         float32,float32,float32,float32,float32',\
                         names='lambda,compname,minwav,maxwav,linename,ngauss,inisig,minsig,maxsig,voff,vindex,windex,findex,fvalue')
    #------header-----------------
    hdr = fits.Header()
    hdr['lambda'] = 'Vacuum Wavelength in Ang'
    hdr['minwav'] = 'Lower complex fitting wavelength range'
    hdr['maxwav'] = 'Upper complex fitting wavelength range'
    hdr['ngauss'] = 'Number of Gaussians for the line'
    hdr['inisig'] = 'Initial guess of linesigma [in lnlambda]'
    hdr['minsig'] = 'Lower range of line sigma [lnlambda]'  
    hdr['maxsig'] = 'Upper range of line sigma [lnlambda]'
    hdr['voff  '] = 'Limits on velocity offset from the central wavelength [lnlambda]'
    hdr['vindex'] = 'Entries w/ same NONZERO vindex constrained to have same velocity'
    hdr['windex'] = 'Entries w/ same NONZERO windex constrained to have same width'
    hdr['findex'] = 'Entries w/ same NONZERO findex have constrained flux ratios'
    hdr['fvalue'] = 'Relative scale factor for entries w/ same findex'
    #------save line info-----------
    hdu = fits.BinTableHDU(data=newdata,header=hdr,name='data')
    hdu.writeto(path+'qsopar.fits',overwrite=True)

def job(PMF):
        snfile = np.genfromtxt('/Users/vivekm/SDSSIV_BALQSO/SN1700_estimates_Master_full_sample_sources.txt',names=['name','plate','mjd','fiber','Z_VI','sn'],dtype=('|S25',int,int,int,float,float))
        path1='/Users/vivekm/SDSSIV_BALQSO/pyqsofit/'
        path2='/Users/vivekm/SDSSIV_BALQSO/pyqsofit/test/data/result/'
        path3='/Users/vivekm/SDSSIV_BALQSO/pyqsofit/test/data/QA_other/'
        path4 = '/Users/vivekm/Softwares/sfddata-master/'

        specdir = '/Users/vzm83/SDSSIV_BALQSO/BALQSOs_Spectra/'

        #Read in the data
        data=fits.open(PMF) 
        tlam=10**data[1].data['loglam']        # OBS wavelength [A]
        tflux=data[1].data['flux']             # OBS flux [erg/s/cm^2/A]
        terr=1./np.sqrt(data[1].data['ivar'])  # 1 sigma error
        and_mask = data[1].data['and_mask']   # and_mask

        
        #select only pixels with and_mask==0
        am=np.where(and_mask==0)[0]
        lam=tlam[am] ; flux=tflux[am] ; err = terr[am]
        
        #Optional
        ra=data[0].header['plug_ra']          # RA 
        dec=data[0].header['plug_dec']        # DEC
        plateid = data[0].header['plateid']   # SDSS plate ID
        mjd = data[0].header['mjd']           # SDSS MJD
        fiberid = data[0].header['fiberid']   # SDSS fiber ID
            
        vv=np.where((snfile['plate']==plateid) & (snfile['mjd']==mjd) & (snfile['fiber']==fiberid))[0][0]
        z=snfile['Z_VI'][vv]                # Redshift
        print plateid,mjd,fiberid
        # get data prepared 
        q = QSOFit(lam, flux, err, z, ra = ra, dec = dec, plateid = plateid, mjd = mjd, fiberid = fiberid, path = path1)
        print len(lam),len(flux),len(err),z 
        start = timeit.default_timer()
        # do the fitting
        print 'check'
        q.Fit(name = None,nsmooth = 1, and_or_mask = False, deredden = True, reject_badpix = False, wave_range = None,\
              wave_mask =None, decomposition_host = False, Mi = None, npca_gal = 5, npca_qso = 20, \
              Fe_uv_op = True, poly = True, BC = False, rej_abs = True, initial_guess = None, MC = False, \
              n_trails = 25, linefit = True, tie_lambda = True, tie_width = True, tie_flux_1 = True, tie_flux_2 = True,\
              save_result = True, plot_fig = False,save_fig = True, plot_line_name = True, plot_legend = True, \
              dustmap_path = path4, save_fig_path = path3, save_fits_path = path2,save_fits_name = None)
        
        outwave = q.wave
        outflux = q.flux_prereduced/(1.0+z)
        outerr = q.err_prereduced
        continuum = (q.Manygauss(np.log(q.wave),q.gauss_result)+q.f_conti_model)/(1.0+z)
        #q.stop()
        # MC run
        iteration = 50
        nlam = lam
        nerr = err
        contflux=np.zeros((iteration,len(q.wave)))
        for k in range(iteration):
            nflux = np.array([sp[0]+ np.random.randn()*sp[1] for sp in zip(flux,err)]) 
            q1 = QSOFit(nlam, nflux, nerr, z, ra = ra, dec = dec, plateid = plateid, mjd = mjd, fiberid = fiberid, path = path1)
            q1.Fit(name = None,nsmooth = 1, and_or_mask = False, deredden = True, reject_badpix = False, wave_range = None,\
              wave_mask =None, decomposition_host = False, Mi = None, npca_gal = 5, npca_qso = 20, \
              Fe_uv_op = True, poly = True, BC = False, rej_abs = True, initial_guess = None, MC = False, \
              n_trails = 25, linefit = True, tie_lambda = True, tie_width = True, tie_flux_1 = True, tie_flux_2 = True,\
              save_result = False, plot_fig = False,save_fig = False, plot_line_name = False, plot_legend = False, \
              dustmap_path = path4, save_fig_path = path3, save_fits_path = path2,save_fits_name = None)
            contflux[k,:] = q1.flux_prereduced
        good = [(~np.all(np.isnan(aa)) & (np.max(aa) < 1.25*np.max(q.flux_prereduced))) for aa in contflux]
        contfluxerr = np.std(contflux[good],axis=0)
        print 'Original Fit ',len(lam),len(q.wave),len(q.flux_prereduced),len(q.err_prereduced)
        print 'Continuum Flux err',len(contfluxerr)
        #Save the normalized flux
        normflux = outflux/continuum
        normerr=np.sqrt(outerr**2+contfluxerr**2)/continuum
        filename = os.path.join(path1,'test/data/norm/',os.path.basename(PMF).split('.')[0].replace('spec','norm')+'.txt')
        np.savetxt(filename,zip(outwave,normflux,normerr,outflux,outerr,continuum,contfluxerr), fmt='%10.5f')
        end = timeit.default_timer()
        print 'Fitting {0} finished in : {1} s ' .format(PMF,str(np.round(end-start)))
        #q1.stop()

snfile = np.genfromtxt('/Users/vivekm/SDSSIV_BALQSO/SN1700_estimates_Master_full_sample_sources.txt',names=['name','plate','mjd','fiber','Z_VI','sn'],dtype=('|S25',int,int,int,float,float))

path1='/Users/vivekm/SDSSIV_BALQSO/pyqsofit/'
path2='/Users/vivekm/SDSSIV_BALQSO/pyqsofit/test/data/oresult/'
path3='/Users/vivekm/SDSSIV_BALQSO/pyqsofit/test/data/oQA_other/'
path4 = '/Users/vivekm/Softwares/sfddata-master/'

specdir = '/Users/vivekm/SDSSIV_BALQSO/BALQSOs_Spectra/'

#if __name__ == '__main__':
start = timeit.default_timer()
count =0
tfiles = []
for path in zip(glob.glob(specdir+'*.fits')):
    tfiles.append(path[0])
print len(tfiles)
files=tfiles[0:4]
print files
pool = Pool()                         # Create a multiprocessing Pool
pool.map(job, files)
pool.close()
#except:
#    count +=1
#    unr=open('UnabletodoPyQSOFit_test.txt','a')
#    print>>unr,'Unable2fit-PyQSOFit: {0}'.format(count)
#    unr.close()

end = timeit.default_timer()
print('Fitting finished in : '+str(np.round(end-start))+'s')


##for i in range(len(snfile)):
#for i in range(350,1000):
#    try:
#        # Set the filename
#        plate = snfile['plate'][i]
#        mjd = snfile['mjd'][i]
#        fiber = snfile['fiber'][i]
#        if plate >=10000:
#            PMF = 'spec-{0:05d}-{1:05d}-{2:04d}.fits'.format(plate,mjd,fiber)
#        else:
#            PMF = 'spec-{0:04d}-{1:05d}-{2:04d}.fits'.format(plate,mjd,fiber)
#        
#        #Read in the data
#        data=fits.open(os.path.join(specdir,PMF)) 
#        tlam=10**data[1].data['loglam']        # OBS wavelength [A]
#        tflux=data[1].data['flux']             # OBS flux [erg/s/cm^2/A]
#        terr=1./np.sqrt(data[1].data['ivar'])  # 1 sigma error
#        and_mask = data[1].data['and_mask']   # and_mask
#        z=snfile['Z_VI'][i]                # Redshift
#        
#        #select only pixels with and_mask==0
#        am=np.where(and_mask==0)[0]
#        lam=tlam[am] ; flux=tflux[am] ; err = terr[am]
#        
#        #Optional
#        ra=data[0].header['plug_ra']          # RA 
#        dec=data[0].header['plug_dec']        # DEC
#        plateid = data[0].header['plateid']   # SDSS plate ID
#        mjd = data[0].header['mjd']           # SDSS MJD
#        fiberid = data[0].header['fiberid']   # SDSS fiber ID
#
#
#        # get data prepared 
#        q = QSOFit(lam, flux, err, z, ra = ra, dec = dec, plateid = plateid, mjd = mjd, fiberid = fiberid, path = path1)
#        
#        start = timeit.default_timer()
#        # do the fitting
#        q.Fit(name = None,nsmooth = 1, and_or_mask = True, deredden = True, reject_badpix = True, wave_range = None,\
#              wave_mask =None, decomposition_host = False, Mi = None, npca_gal = 5, npca_qso = 20, \
#              Fe_uv_op = True, poly = True, BC = False, rej_abs = True, initial_guess = None, MC = True, \
#              n_trails = 25, linefit = True, tie_lambda = True, tie_width = True, tie_flux_1 = True, tie_flux_2 = True,\
#              save_result = True, plot_fig = True,save_fig = True, plot_line_name = True, plot_legend = True, \
#              dustmap_path = path4, save_fig_path = path3, save_fits_path = path2,save_fits_name = None)
#        
#        end = timeit.default_timer()
#        print 'Fitting {0} finished in : {1} s ' .format(PMF,str(np.round(end-start)))
#        # grey shade on the top is the continuum windiows used to fit.
#    except:
#        unr=open('UnabletodoPyQSOFit.txt','a')
#        print>>unr,'Unable2fit-PyQSOFit: {0}\t{1}\t{2}\t{3}'.format(snfile['name'][i],snfile['plate'][i],snfile['mjd'][i],snfile['fiber'][i])
#        unr.close()
