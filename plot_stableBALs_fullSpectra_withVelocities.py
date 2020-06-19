import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
import scipy.stats as stats
import seaborn as sns
from matplotlib.ticker import NullFormatter
import os

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


n_v = 1240.81
si_iv = 1399.8
c_iv = 1549.48
al_iii = 1857.4
mg_ii = 2799.117
p_v = 1122.9925#1128.008#1117.977 #


df = np.genfromtxt('stable_dummy_3epoch_info.dat',names=['ra','dec','z','pmf1','pmf2','pmf3'],dtype=(float,float,float,'|S25','|S25','|S25'))
lb_list = np.genfromtxt('name_pmf_match.txt',names=['sdssname','redshift','pmf'],dtype=('|S55',float,'|S25'))
pf = np.genfromtxt('Stable_BAL_troughs_trough_info.dat',names=('index','filename','vmin','vmax','nt','ntroughs','paddmin','paddmax','ew','ewerr'),dtype=(int,'|S45',float,float,float,float,float,float,float,float))

pp = PdfPages('StableBALs_other_ion_troughs_PV.pdf')

for i in range(len(df)):
#for i in range(2):
    file1 = '../Continuumfits_fullspectra/Normspec_'+df['pmf1'][i]+'_powerlaw.txt'
    file2 = '../Continuumfits_fullspectra/Normspec_'+df['pmf2'][i]+'_powerlaw.txt'
    file3 = '../Continuumfits_fullspectra/Normspec_'+df['pmf3'][i]+'_powerlaw.txt'
    fsp1 = np.genfromtxt(file1, names = ['wave','nflux','nerr','cont','flux','err','kate_cont','kate_cont_err'], \
            dtype = (float,float,float,float,float,float,float,float),\
            usecols = [0,1,2,3,4,5,6,7])
    fsp2 = np.genfromtxt(file2, names = ['wave','nflux','nerr','cont','flux','err','kate_cont','kate_cont_err'], \
            dtype = (float,float,float,float,float,float,float,float),\
            usecols = [0,1,2,3,4,5,6,7])
    fsp3 = np.genfromtxt(file3, names = ['wave','nflux','nerr','cont','flux','err','kate_cont','kate_cont_err'], \
            dtype = (float,float,float,float,float,float,float,float),\
            usecols = [0,1,2,3,4,5,6,7])
    os.system('cp {} {}'.format(file1,'PowLaw_contNorm_stableBALs/.'))    
    os.system('cp {} {}'.format(file2,'PowLaw_contNorm_stableBALs/.'))    
    os.system('cp {} {}'.format(file3,'PowLaw_contNorm_stableBALs/.'))    
    legend1 = str(file1.split('/')[2].split('_')[1])
    legend2 = str(file2.split('/')[2].split('_')[1])
    legend3 = str(file3.split('/')[2].split('_')[1])
    
    fig=plt.figure( figsize=(20, 10))
    
    nullfmt = NullFormatter()         # no labels
    
    lm = np.where(lb_list['pmf'] == legend1)[0]
    sdss_name = '{0}\nZ = {1:3.4f}'.format('SDSS J'+str(lb_list['sdssname'][lm[0]]), df['z'][i])

    dim_fs1 = [0.05, 0.8, 0.6, 0.1]
    dim_fs2 = [0.05, 0.55, 0.6, 0.25]
    dim_fs3 = [0.05, 0.3, 0.6, 0.25]
    dim_fs4 = [0.05, 0.05, 0.6, 0.25]
    
    
    dim_vel1 = [0.7, 0.725, 0.25, 0.225]
    dim_vel2 = [0.7, 0.5, 0.25, 0.225]
    dim_vel3 = [0.7, 0.275, 0.25, 0.225]
    dim_vel4 = [0.7, 0.05, 0.25, 0.225]
    
    ax_fs1 = plt.axes(dim_fs1)
    ax_fs2 = plt.axes(dim_fs2)
    ax_fs3 = plt.axes(dim_fs3)
    ax_fs4 = plt.axes(dim_fs4)


    ax_vel1 = plt.axes(dim_vel1)
    ax_vel2 = plt.axes(dim_vel2)
    ax_vel3 = plt.axes(dim_vel3)
    ax_vel4 = plt.axes(dim_vel4)
    
    wavemin,wavemax= min(min(fsp1['wave']),min(fsp2['wave']),min(fsp3['wave']))-50, max(max(fsp1['wave']),max(fsp2['wave']),max(fsp2['wave']) )+50
    
    ax_fs2.plot(fsp1['wave'],fsp1['flux'],color='black',label=legend1,alpha=0.5)
    ax_fs2.plot(fsp1['wave'],fsp1['cont'],color='orange',alpha=0.8,ls='--',lw=2)
    ax_fs2.legend(loc=1)
    ax_fs2.xaxis.set_major_formatter(nullfmt)
    ax_fs2.set_ylabel('Flux')   
    ax_fs2.set_xlim(wavemin,wavemax)
    

    ax_fs3.plot(fsp2['wave'],fsp2['flux'],color='red',label=legend2,alpha=0.5)
    ax_fs3.plot(fsp2['wave'],fsp2['cont'],color='orange',alpha=0.8,ls='--',lw=2)
    ax_fs3.legend(loc=1)
    ax_fs3.xaxis.set_major_formatter(nullfmt)
    ax_fs3.set_ylabel('Flux')    
    ax_fs3.set_xlim(wavemin,wavemax)
    
    ax_fs4.plot(fsp3['wave'],fsp3['flux'],color='blue',label=legend3,alpha=0.5)
    ax_fs4.plot(fsp3['wave'],fsp3['cont'],color='orange',alpha=0.8,ls='--',lw=2)
    ax_fs4.legend(loc=1)
    ax_fs4.set_xlabel('Rest-frame Wavelengths $\AA$') 
    ax_fs4.set_ylabel('Flux')    
    ax_fs4.set_xlim(wavemin,wavemax)
    
    ax_fs1.axis('off')
    ax_fs1.text(0.1,0.65,sdss_name,fontsize =25)


    ax_vel1.plot(velocity(fsp1['wave'],si_iv,0),fsp1['nflux'],color='black',alpha=0.5)
    ax_vel1.plot(velocity(fsp2['wave'],si_iv,0),fsp2['nflux'],color='red',alpha=0.5)
    ax_vel1.plot(velocity(fsp3['wave'],si_iv,0),fsp3['nflux'],color='blue',alpha=0.5)
    ax_vel1.set_xlim(-25000,3000)    
    ax_vel1.set_ylabel('Norm Flux')
    ax_vel1.axvline(0.0,ls='--',lw=2,color='black')
    ax_vel1.xaxis.set_major_formatter(nullfmt)
    ax_vel1.text(-24000,-0.03,'Si IV',fontsize=18)
    ax_vel1.axvspan(pf['vmax'][i],pf['vmin'][i], alpha=0.05, color='red')
    vel1ylim = plt.ylim()

    ax_vel2.plot(velocity(fsp1['wave'],c_iv,0),fsp1['nflux'],color='black',alpha=0.5)
    ax_vel2.plot(velocity(fsp2['wave'],c_iv,0),fsp2['nflux'],color='red',alpha=0.5)
    ax_vel2.plot(velocity(fsp3['wave'],c_iv,0),fsp3['nflux'],color='blue',alpha=0.5)
    ax_vel2.set_xlim(-25000,3000)    
    ax_vel2.set_ylabel('Norm Flux')    
    ax_vel2.axvline(0.0,ls='--',lw=2,color='black')
    ax_vel2.xaxis.set_major_formatter(nullfmt)
    ax_vel2.text(-24000,-0.03,'C IV',fontsize=18)
    ax_vel2.axvspan(pf['vmax'][i],pf['vmin'][i], alpha=0.05, color='red')
    vel2ylim = plt.ylim()

    ax_vel3.plot(velocity(fsp1['wave'],al_iii,0),fsp1['nflux'],color='black',alpha=0.5)
    ax_vel3.plot(velocity(fsp2['wave'],al_iii,0),fsp2['nflux'],color='red',alpha=0.5)
    ax_vel3.plot(velocity(fsp3['wave'],al_iii,0),fsp3['nflux'],color='blue',alpha=0.5)
    ax_vel3.set_xlim(-25000,3000)    
    ax_vel3.set_ylabel('Norm Flux')    
    ax_vel3.axvline(0.0,ls='--',lw=2,color='black')
    ax_vel3.xaxis.set_major_formatter(nullfmt)
    #ax_vel3.set_title('Al III')
    ax_vel3.text(-24000,-0.03,'Al III',fontsize=18)
    ax_vel3.axvspan(pf['vmax'][i],pf['vmin'][i], alpha=0.05, color='red')
    vel3ylim = plt.ylim()

#    ax_vel4.plot(velocity(fsp1['wave'],mg_ii,0),fsp1['nflux'],color='black',alpha=0.5)
#    ax_vel4.plot(velocity(fsp2['wave'],mg_ii,0),fsp2['nflux'],color='red',alpha=0.5)
#    ax_vel4.plot(velocity(fsp3['wave'],mg_ii,0),fsp3['nflux'],color='blue',alpha=0.5)
#    ax_vel4.set_xlim(-25000,3000)    
#    ax_vel4.set_xlabel('Velocity')    
#    ax_vel4.set_ylabel('Norm Flux')    
#    ax_vel4.axvline(0.0,ls='--',lw=2,color='black')
#    ax_vel4.text(-24000,-0.03,'Mg II',fontsize=18)
#    ax_vel4.axvspan(pf['vmax'][i],pf['vmin'][i], alpha=0.05, color='red')
#    vel4ylim = plt.ylim()

    ax_vel4.plot(velocity(fsp1['wave'],p_v,0),fsp1['nflux'],color='black',alpha=0.5)
    ax_vel4.plot(velocity(fsp2['wave'],p_v,0),fsp2['nflux'],color='red',alpha=0.5)
    ax_vel4.plot(velocity(fsp3['wave'],p_v,0),fsp3['nflux'],color='blue',alpha=0.5)
    ax_vel4.set_xlim(-25000,3000)    
    ax_vel4.set_xlabel('Velocity')    
    ax_vel4.set_ylabel('Norm Flux')    
    ax_vel4.axvline(0.0,ls='--',lw=2,color='black')
    ax_vel4.text(-24000,-0.03,'P V',fontsize=18)
    ax_vel4.axvspan(pf['vmax'][i],pf['vmin'][i], alpha=0.05, color='red')
    vel4ylim = plt.ylim()

    
    yllim,yulim = min([vel1ylim[0],vel2ylim[0],vel3ylim[0],vel4ylim[0] ] ), max([vel1ylim[1],vel2ylim[1],vel3ylim[1],vel4ylim[1] ] ) 

    ax_vel4.set_ylim(-0.5,1.85)    
    ax_vel3.set_ylim(-0.5,1.85)    
    ax_vel2.set_ylim(-0.5,1.85)    
    ax_vel1.set_ylim(-0.5,1.85)   
    fig.tight_layout()
    fig.savefig(pp,format='pdf')
    #plt.show()

pp.close()
