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
from astropy import units as u
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
from astroquery.irsa import Irsa

def queryWISE(ra,dec):
    wise = Irsa.query_region(SkyCoord(ra,dec, unit=(u.deg,u.deg),frame='fk5'),catalog='allwise_p3as_psd', radius='0d0m2s',selcols='ra,dec,w1mpro,w1sigmpro,w2mpro,w2sigmpro,w3mpro,w3sigmpro,w4mpro,w4sigmpro')
    #print wise['ra'].data[0]
    return (wise['ra'].data[0],wise['dec'].data[0],wise['w1mpro'].data[0],wise['w1sigmpro'].data[0],wise['w2mpro'].data[0],wise['w2sigmpro'].data[0],wise['w3mpro'].data[0],wise['w3sigmpro'].data[0],wise['w4mpro'].data[0],wise['w4sigmpro'].data[0])



def convertWISEmag2Flux(mag,dmag,band):
    if band == 'W1':
        zero_point =  309.540
    elif band == 'W2':
        zero_point =  171.787
    elif band =='W3' :
        zero_point =  31.674
    else :
        zero_point =  8.363
    if mag: 
        fnu = zero_point*10**(-mag/2.5)
    else:
        fnu = -9999.0
    if dmag:
        dfnu = fnu*dmag/1.086
    else:
        dfnu = -9999.0
    #return fnu*1000,dfnu*1000 # for CIGALE
    return fnu,dfnu


def querySDSS(ra,dec):

    sdss = SDSS.query_crossid(coords.SkyCoord(ra,dec,unit=(u.deg,u.deg),frame='fk5'))
    if sdss:
        uprime_mag,uprime_magerr = sdss['psfMag_u'].item(),sdss['psfMagerr_u'].item()
        gprime_mag,gprime_magerr = sdss['psfMag_g'].item(),sdss['psfMagerr_g'].item()
        rprime_mag,rprime_magerr = sdss['psfMag_r'].item(),sdss['psfMagerr_r'].item()
        iprime_mag,iprime_magerr = sdss['psfMag_i'].item(),sdss['psfMagerr_i'].item()
        zprime_mag,zprime_magerr = sdss['psfMag_z'].item(),sdss['psfMagerr_z'].item()
    else:
        uprime_mag,uprime_magerr = -9999,-9999
        gprime_mag,uprime_magerr = -9999,-9999
        rprime_mag,uprime_magerr = -9999,-9999
        iprime_mag,uprime_magerr = -9999,-9999
        zprime_mag,uprime_magerr = -9999,-9999
    psfmag = (uprime_mag,gprime_mag,rprime_mag,iprime_mag,zprime_mag)
    psfmag_err = (uprime_magerr,gprime_magerr,rprime_magerr,iprime_magerr,zprime_magerr)
    return psfmag,psfmag_err

def convertSDSSMag2Flux(mag,dmag):
    if mag > 0:
        fnu =  3.631*10**(-6)* 10**(-0.4*(mag - 22.5))

        dfnu = fnu*dmag/1.086
    else:
        fnu = -9999.0
        dfnu = -9999.0

    return fnu,dfnu

def convertNanoMaggy2Flux(nmagy):
    if nmagy >0:
        #return 3.631*10**(-6)* nmagy*1000 # for CIGALE
        return 3.631*10**(-6)* nmagy
    else:
        return -9999.0

def convertNanoMaggy2Flux_err(nmagy):
    if nmagy >0:
        #return 3.631*10**(-6)* 1.0/np.sqrt(nmagy)*1000 # for CIGALE
        return 3.631*10**(-6)* 1.0/np.sqrt(nmagy) 
    else:
        return -9999.0


def query2mass(ra,dec):
    tmass = Irsa.query_region(SkyCoord(df['RA'][i],df['DEC'][i], unit=(u.deg,u.deg),frame='fk5'),catalog='fp_psc', radius='0d0m2s',selcols='ra,dec,j_m,j_cmsig,h_m,h_cmsig,k_m,k_cmsig')
    if tmass:
        return tmass['ra'].data[0],tmass['dec'].data[0],tmass['j_m'].data[0],tmass['j_cmsig'].data[0],tmass['h_m'].data[0],tmass['h_cmsig'].data[0],tmass['k_m'].data[0],tmass['k_cmsig'].data[0]
    else :
        return ra,dec,-9999.0,-9999.0,-9999.0,-9999.0,-9999.0,-9999.0


def photometry_2mass_mag_to_jy(magnitude,magerr, band):

    """
    This function ...
    :param magnitude:
    :param band:
    :return:

    """
    # 2MASS F_0 (in Jy)
    f_0_2mass = {"2MASS.J": 1594.0, "2MASS.H": 1024.0, "2MASS.K": 666.7}
    #print magnitude,magerr,band 
    if magnitude>0:
        fnu = f_0_2mass["2MASS."+band] * 10.**(-magnitude/2.5)
    else : 
        fnu = -9999.0
    if magerr>0:
        dfnu =  fnu*magerr/1.086
    else :
        dfnu = -9999.0
    #return fnu*1000,dfnu*1000 #for CIGALE
    return fnu,dfnu


def prepare_sed_sdss():
    FUV_wl = 1525
    NUV_wl = 2200
    sdss_u_wl  = 3543
    sdss_g_wl  = 4770
    sdss_r_wl  = 6231
    sdss_i_wl  = 7625
    sdss_z_wl  = 9134
    W1_wl = 33530
    W2_wl = 46030
    W3_wl = 115610
    W4_wl = 220880
    tmass_J_wl = 12350
    tmass_H_wl = 16620
    tmass_K_wl = 21590

    inp=open('inputcatalog_CLQ_nwise_mJy_CIGALE.txt','w')
    #print>>inp,'{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t'.format('#id','redshift', )
    print>>inp,'{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format('#id','redshift','WISE1','WISE1_err','WISE2','WISE2_err','WISE3','WISE3_err','WISE4','WISE4_err','u_prime','u_prime_err','g_prime','g_prime_err','r_prime','r_prime_err','i_prime','i_prime_err','z_prime','z_prime_err','J_2mass','J_2mass_err','H_2mass','H_2mass_err','Ks_2mass','Ks_2mass_err','NUV','NUV_err','FUV','FUV_err')
    
    dr16_clq = Table.read('changing_look_quasar_list_xmatch_dr16.fits',format='fits')
    for i in range(len(dr16_clq)):

        # WISE FLUX
        wise_ra,wise_dec,w1mpro,w1sigmpro,w2mpro,w2sigmpro,w3mpro,w3sigmpro,w4mpro,w4sigmpro =queryWISE(dr16_clq['ra_1'][i],dr16_clq['dec_1'][i])
        #w1mpro,w1sigmpro,w2mpro,w2sigmpro,w3mpro,w3sigmpro,w4mpro,w4sigmpro = dr16_clq['W1_MAG'][i],dr16_clq['W1_MAG_ERR'][i],dr16_clq['W2_MAG'][i],dr16_clq['W2_MAG_ERR'][i],dr16_clq['W3_MAG'][i],dr16_clq['W3_MAG_ERR'][i],dr16_clq['W4_MAG'][i],dr16_clq['W4_MAG_ERR'][i]
        w1flux,w1fluxwerr = convertWISEmag2Flux(w1mpro,w1sigmpro,'W1')
        w2flux,w2fluxwerr = convertWISEmag2Flux(w2mpro,w2sigmpro,'W2')
        w3flux,w3fluxwerr = convertWISEmag2Flux(w3mpro,w3sigmpro,'W3')
        w4flux,w4fluxwerr = convertWISEmag2Flux(w4mpro,w4sigmpro,'W4')
        
        #SDSS FLUX
        sdss_u_flux,sdss_u_flux_err =  convertNanoMaggy2Flux(dr16_clq['PSFFLUX'][i,0]),convertNanoMaggy2Flux(1.0/np.sqrt(dr16_clq['PSFFLUX_IVAR'][i,0]))
        sdss_g_flux,sdss_g_flux_err =  convertNanoMaggy2Flux(dr16_clq['PSFFLUX'][i,1]),convertNanoMaggy2Flux(1.0/np.sqrt(dr16_clq['PSFFLUX_IVAR'][i,1]))
        sdss_r_flux,sdss_r_flux_err =  convertNanoMaggy2Flux(dr16_clq['PSFFLUX'][i,2]),convertNanoMaggy2Flux(1.0/np.sqrt(dr16_clq['PSFFLUX_IVAR'][i,2]))
        sdss_i_flux,sdss_i_flux_err =  convertNanoMaggy2Flux(dr16_clq['PSFFLUX'][i,3]),convertNanoMaggy2Flux(1.0/np.sqrt(dr16_clq['PSFFLUX_IVAR'][i,3]))
        sdss_z_flux,sdss_z_flux_err =  convertNanoMaggy2Flux(dr16_clq['PSFFLUX'][i,4]),convertNanoMaggy2Flux(1.0/np.sqrt(dr16_clq['PSFFLUX_IVAR'][i,4]))

        #2MASS FLUX
        TMASS_J_flux,TMASS_J_flux_err = photometry_2mass_mag_to_jy(float(dr16_clq['JMAG'][i]),float(dr16_clq['JMAG_ERR'][i]),'J')
        TMASS_H_flux,TMASS_H_flux_err = photometry_2mass_mag_to_jy(float(dr16_clq['HMAG'][i]),float(dr16_clq['HMAG_ERR'][i]),'H')
        TMASS_K_flux,TMASS_K_flux_err = photometry_2mass_mag_to_jy(float(dr16_clq['KMAG'][i]),float(dr16_clq['KMAG_ERR'][i]),'K')

        #GALEX FLUX
        if dr16_clq['NUV_IVAR'][i] < 0 :
            nuv_temp_ivar = 0
        else :
            nuv_temp_ivar = dr16_clq['NUV_IVAR'][i]
        if dr16_clq['FUV_IVAR'][i] < 0 :
            fuv_temp_ivar = 0
        else:
            fuv_temp_ivar = dr16_clq['FUV_IVAR'][i]


        GALEX_NUV_flux, GALEX_NUV_flux_err = convertNanoMaggy2Flux(dr16_clq['NUV'][i]),convertNanoMaggy2Flux_err(dr16_clq['NUV_IVAR'][i])
        GALEX_FUV_flux, GALEX_FUV_flux_err = convertNanoMaggy2Flux(dr16_clq['FUV'][i]),convertNanoMaggy2Flux_err(dr16_clq['FUV_IVAR'][i])
        print  GALEX_NUV_flux, GALEX_NUV_flux_err 
        print>>inp,'{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t'.format(dr16_clq['SDSS_NAME'][i][0:8],dr16_clq['Z'][i], 
           W1_wl, w1flux,w1fluxwerr, 
           W2_wl, w2flux,w2fluxwerr, 
           W3_wl, w3flux,w3fluxwerr, 
           W4_wl, w4flux,w4fluxwerr, 
           sdss_u_wl,sdss_u_flux,sdss_u_flux_err,
           sdss_g_wl,sdss_g_flux,sdss_g_flux_err,
           sdss_r_wl,sdss_r_flux,sdss_r_flux_err,
           sdss_i_wl,sdss_i_flux,sdss_i_flux_err,
           sdss_z_wl,sdss_z_flux,sdss_z_flux_err,
           tmass_J_wl, TMASS_J_flux,TMASS_J_flux_err, 
           tmass_H_wl, TMASS_H_flux,TMASS_H_flux_err, 
           tmass_K_wl, TMASS_K_flux,TMASS_K_flux_err, 
           NUV_wl,  GALEX_NUV_flux, GALEX_NUV_flux_err, 
           FUV_wl,  GALEX_FUV_flux, GALEX_FUV_flux_err )


## for CIGALE
#        print>>inp,'{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(dr16_clq['SDSS_NAME'][i],dr16_clq['Z'][i], 
#           w1flux,w1fluxwerr, 
#           w2flux,w2fluxwerr, 
#           w3flux,w3fluxwerr, 
#           w4flux,w4fluxwerr, 
#           sdss_u_flux,sdss_u_flux_err,
#           sdss_g_flux,sdss_g_flux_err,
#           sdss_r_flux,sdss_r_flux_err,
#           sdss_i_flux,sdss_i_flux_err,
#           sdss_z_flux,sdss_z_flux_err,
#           TMASS_J_flux,TMASS_J_flux_err, 
#           TMASS_H_flux,TMASS_H_flux_err, 
#           TMASS_K_flux,TMASS_K_flux_err, 
#           GALEX_NUV_flux, GALEX_NUV_flux_err, 
#           GALEX_FUV_flux, GALEX_FUV_flux_err )


    inp.close()
    

def main():     

    print '-------'
    prepare_sed_sdss()


if __name__== "__main__":
    main()
