
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

from scipy.ndimage.filters import gaussian_filter1d,uniform_filter1d

params = {
   'axes.labelsize': 18,
   'axes.linewidth': 1.5,
   #'text.fontsize': 18,
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



def plot_spectra(i,j,data1,plates,maxflux,minflux,fig,ax):
    if ((i==0) & (j== 1)):
        ax.plot(10**data1['loglam']/(1.0+zvi),gaussian_filter1d(data1['flux']/1.625,2),color=plt.cm.RdYlBu(j*300),alpha=0.95,label='')
        ax.set(ylabel='Flux',xlim=(1400,2900),ylim=(0,30))
        ax.annotate('FeLoBAL changes to LoBAL', xy=(2350, 13), xytext=(2000, 2), fontsize=18,color='blue',
            arrowprops=dict(facecolor='red', shrink=0.05) )
    elif ((i==1) & (j== 1)):
        ax.plot(10**data1['loglam']/(1.0+zvi),gaussian_filter1d(data1['flux']*1.2,5),color=plt.cm.RdYlBu(j*300),alpha=0.95,label='')
        ax.set(ylabel='Flux',xlim=(1300,2500),ylim=(0,max(minflux)+1.0*max(maxflux)))
    elif ((i==1) & (j== 0)):
        ax.plot(10**data1['loglam']/(1.0+zvi),gaussian_filter1d(data1['flux']*1,5),color=plt.cm.RdYlBu(j*300),alpha=0.95,label='')
        ax.set(ylabel='Flux',xlim=(1300,2500),ylim=(0,max(minflux)+1.0*max(maxflux)))
        ax.annotate('Redshifted BAL', xy=(1570, 1.5), xytext=(1600, 0.4), fontsize=18,color='blue',
            arrowprops=dict(facecolor='red', shrink=0.05) )
        ax.annotate('', xy=(1420, 1.85), xytext=(1600, 0.4), fontsize=18,color='red',
            arrowprops=dict(facecolor='red', shrink=0.05) )

    elif ((i==2) & (j== 0)):
        ax.plot(10**data1['loglam']/(1.0+zvi),gaussian_filter1d(data1['flux']/1.15,1),color=plt.cm.RdYlBu(j*300),alpha=0.95,label='')
        ax.set(ylabel='Flux',xlim=(1300,2000),ylim=(0,25))
        ax.annotate('Many narrow trough system', xy=(1390, 2.5), xytext=(1450, 1), fontsize=18,color='blue',
            arrowprops=dict(facecolor='red', shrink=0.05) )


    else:
        ax.plot(10**data1['loglam']/(1.0+zvi),gaussian_filter1d(data1['flux'],1),color=plt.cm.RdYlBu(j*300),alpha=0.95,label='')
    #ax.plot(10**data1['loglam']/(1.0+zvi),1.0/np.sqrt(data1['ivar']),color=plt.cm.RdYlBu(j*300),alpha=0.1)
 
        #ax.set(ylabel='Flux',ylim=(0,max(minflux)+2.0*max(maxflux)))
 
    if i == 2:
        if (len(maxflux) > 0):
            ax.set(xlabel='Rest wavelength ($\AA$)',ylabel='Flux')
    xlim,ylim=ax.get_xlim(),ax.get_ylim()
    ax.tick_params(direction='in')
    #string1 = 'SDSS J{0}\tZ\_VI: {1:4.4f}\tN$\_{{spec}}$: {2}'.format(fsm['SDSS_NAME'][xx[0]],zvi,len(plates))
    #string1 = 'RA: {0:5.4f}\t DEC: {1:5.4f}  \tZ\_VI: {2:4.4f}\tN$\_{{spec}}$: {3}'.format(fsm['ra'][xx[0]],fsm['dec'][xx[0]],zvi,len(plates))
    #print string1,xlim,ylim
    #ax.text(xlim[0]+0.05*(xlim[1] - xlim[0]),ylim[1]-0.05*(ylim[1] - ylim[0]), string1,fontsize=18)
    if j == 1:
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
    
        #ax.legend(loc=1)
        #fig.savefig(pp,format='pdf')
        #sdlfj=raw_input()




text_font = {'fontname':'Times New Roman', 'size':'14'}
#pp = PdfPages('Fullanalysis_BAL_plot_spectra.pdf') 
specdir = 'BALQSOs_Spectra'
linelist = np.genfromtxt('/Users/vivekm/Proposals/linelist_speccy.txt',usecols=(0,1,2),dtype=('|S10',float,'|S5'),names=True)

data = np.load('MasterList_Plate-MJD-Fiber.npz')
#fsm = fits.open('Fullsample_3028_BAL_sources_crossmatch-DR12Q.fits')[1].data
fsm = np.genfromtxt('Master_full_sample_redshift.txt',names=['name','ra','dec','Z_VI','tag'],dtype=('|S25',float,float,float,'|S25'))

plotclass=['Target1622','Target0200','Target0429']

fig,(ax0,ax1,ax2)=plt.subplots(3,1,figsize=(12,8))

for i in range(len(plotclass)):
    x=np.where(data['name']==plotclass[i])[0][0]
    plates = data['plate'][x]
    mjds = data['mjd'][x]
    fibers = data['fiber'][x]
    print data['name'][x]
    print plates
    xx=np.where(fsm['name'] == data['name'][x])[0]
    if len(xx)>0:
        zvi = fsm['Z_VI'][xx[0]]
        print xx,data['name'][i],fsm['name'][xx]
        minflux =[] ;maxflux=[]
        #fig,ax=plt.subplots(figsize=(15,8))
        for j in range(2):
            if plates[j] >=10000:
                PMF = 'spec-{0:05d}-{1:05d}-{2:04d}.fits'.format(plates[j],mjds[j],fibers[j])
            else:
                PMF = 'spec-{0:04d}-{1:05d}-{2:04d}.fits'.format(plates[j],mjds[j],fibers[j])
            print PMF
            if os.path.isfile(os.path.join(specdir,PMF)):
                data1 = fits.open(os.path.join(specdir,PMF))[1].data
                print data1
                minflux.append(np.median((data1['flux']*(data1['and_mask'] == 0)).copy()))
                maxflux.append(np.std((data1['flux']*(data1['and_mask'] == 0)).copy()))
                if i==0:
                    plot_spectra(i,j,data1,plates,maxflux,minflux,fig,ax0)
                if i==1:
                    plot_spectra(i,j,data1,plates,maxflux,minflux,fig,ax1)
                if i==2:
                    plot_spectra(i,j,data1,plates,maxflux,minflux,fig,ax2)
                        
        #ax.legend(loc=1)
        fig.tight_layout()
        fig.savefig('Exampleplots_unusualBALs.jpg',format='jpg')

plt.show()


