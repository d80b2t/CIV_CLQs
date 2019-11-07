'''

Plotting the observed CIV emission lines from the z~2 CLQs and
comparing those to the e.g. QSFIT model parameters

'''

import math
import numpy as np
import scipy.stats as stats

import matplotlib.lines    as mlines
import matplotlib.pyplot   as plt
import matplotlib.gridspec as gridspec

from matplotlib          import colors as mcolors
from matplotlib.gridspec import GridSpec
from matplotlib.patches  import Rectangle

from astropy.io          import ascii
from astropy.io          import fits
from astropy.convolution import convolve, Box1DKernel


##
## R E A D I N G   I N   T H E    D A T A 
##
path = '../../data/QSFIT/'


##  J 1 2 0 5 +  3 4 2 2 
sdssname   = 'spec-2089-53498-0427.fits'
J12_first  = fits.open(path+sdssname)
J12_second = ascii.read(path+'J1205p3422_58538.dat')
J12_third  = ascii.read(path+'J1205p3422_58693.dat')

J12_redshift    = 2.071
J12_first_flux  =       J12_first[1].data.flux
J12_first_wave  = (10**(J12_first[1].data.loglam)) / (1 + J12_redshift)
J12_second_wave = J12_second['wavelength'] / (1 + J12_redshift)
J12_second_flux = J12_second['flux'] 
J12_third_wave  = J12_third['wavelength']  / (1 + J12_redshift)
J12_third_flux  = J12_third['flux']    


##  J 1 6 3 8 +  2 8 2 7
sdssname   = 'spec-2948-54553-0614.fits'
J16_first  = fits.open(path+sdssname)
bossname   = 'spec-5201-55832-0178_v5_10_0.fits'
J16_second = fits.open(path+sdssname)
J16_third  = ascii.read(path+'J1638p2827_58583.dat')

J16_redshift    = 2.182
J16_first_flux  =        J16_first[1].data.flux
J16_first_wave  = (10**( J16_first[1].data.loglam))  / (1 + J16_redshift)
J16_second_flux =       J16_second[1].data.flux
J16_second_wave = (10**(J16_second[1].data.loglam)) / (1 + J16_redshift)
J16_third_wave  = J16_third['wavelength']  / (1 + J16_redshift)
J16_third_flux  = J16_third['flux']    


##  J 2 2 2 8  +  2 2 0 1 
sdssname   = 'spec-6118-56189-0720.fits'
J22_first  = fits.open(path+sdssname)
bossname   = 'spec-7582-56960-0790.fits'
J22_second = fits.open(path+sdssname)
J22_third  = ascii.read(path+'J2228p2201_58693.dat')

J22_redshift    = 2.222
J22_first_flux  =        J22_first[1].data.flux
J22_first_wave  = (10**( J22_first[1].data.loglam))  / (1 + J22_redshift)
J22_second_flux =       J22_second[1].data.flux
J22_second_wave = (10**(J22_second[1].data.loglam)) / (1 + J22_redshift)
J22_third_wave  = J22_third['wavelength']  / (1 + J22_redshift)
J22_third_flux  = J22_third['flux']    


##  R E A D I N G    I N    T H E   L I N E   F I T S    F R O M    Q S F I T
infile = 'QSFIT_CIV_line_params.dat'
qsfit = ascii.read(path+infile)

## c in km/s; nominal wavelength in Ang
qsfit_fwhm_wave  = (qsfit['FWHM']  /3e5) * 1548.
qsfit_sigma      = qsfit_fwhm_wave / 2.*np.sqrt(2*np.log(2)) 


##  S E T T I N G   U P   T H E   P L O T
##
## Adjusting the Whitespace for the plots
left   = 0.06   # the left side of the subplots of the figure
right  = 0.97   # the right side of the subplots of the figure
bottom = 0.16   # the bottom of the subplots of the figure
top    = 0.96   # the top of the subplots of the figure
wspace = 0.26   # the amount of width reserved for blank space between subplots
hspace = 0.26   # the amount of height reserved for white space between subplots

plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
 
## Some NPR defaults
ls              = 'solid'
ms              = 1.
ms_big          = ms*6.
ms_large        = ms*8.
alpha           = 1.0
fontsize        = 16
labelsize       = fontsize
tickwidth       = 2.0
linewidth       = 1.4
tickwidth       = 2.0
ticklength      = 6.0
ticklabelsize   = labelsize
majorticklength = 12
minorticklength = 6
plt.rcParams['text.usetex'] = True

## 3 by 3  multi-panels
fig, ((ax1,ax2,ax3), (ax4,ax5,ax6), (ax7,ax8,ax9)) = plt.subplots(nrows=3, ncols=3, figsize=(12.0, 12.0))
##fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(12.0, 12.0))


## 1549 +/- 81::
xmin = 1468; xmax = 1630.    ## 1468, 1620 to match VdB01, Fig 7;

ax1.set_xlim([xmin, xmax])
ax2.set_xlim([xmin, xmax])
ax3.set_xlim([xmin, xmax])
ax4.set_xlim([xmin, xmax])
ax5.set_xlim([xmin, xmax])
ax6.set_xlim([xmin, xmax])
ax7.set_xlim([xmin, xmax])
ax8.set_xlim([xmin, xmax])
ax9.set_xlim([xmin, xmax])

ymin_J12 = -9.5; ymax_J12 = 64.   
ax1.set_ylim([ymin_J12, ymax_J12])
ax2.set_ylim([ymin_J12, ymax_J12])
ax3.set_ylim([ymin_J12, ymax_J12])

ymin_J16 = -2.9; ymax_J16 = 21.9   
ax4.set_ylim([ymin_J16, ymax_J16])
ax5.set_ylim([ymin_J16, ymax_J16])
ax6.set_ylim([ymin_J16, ymax_J16])

ymin_J22 = -4.9; ymax_J22 = 12.9 
ax7.set_ylim([ymin_J22, ymax_J22])
ax8.set_ylim([ymin_J22, ymax_J22])
ax9.set_ylim([ymin_J22, ymax_J22])

ymin = -0.1; ymax = 1.1
#ax1.set_ylim([ymin, ymax]) ; #ax2.set_ylim([ymin, ymax]) ; #ax3.set_ylim([ymin, ymax])
#ax4.set_ylim([ymin, ymax]) ; #ax5.set_ylim([ymin, ymax]) ; #ax6.set_ylim([ymin, ymax])
#ax7.set_ylim([ymin, ymax]) ; #ax8.set_ylim([ymin, ymax]) ; #ax9.set_ylim([ymin, ymax])

ii=0
## ACTUALLY PLOTTING   J12
##
## trying to pick out the flux peak of the CIV line
norm_civ     = J12_first_flux[800:1200].max()
cont_offsett = J12_first_flux[1150]
ax1.plot(J12_first_wave,  J12_first_flux,  '-b', lw=linewidth, label='J1205+3422 (53498)')
#ax1.plot(J12_first_wave,  ((J12_first_flux-cont_offsett) /norm_civ),  '-b', lw=linewidth, label='MJD 53498')

#mu = 1549.
#linefit_x =  np.linspace(mu - 3*qsfit_sigma[ii], mu + 3*qsfit_sigma[ii], 100)  
#linefit   =  stats.norm.pdf(linefit_x, mu, qsfit_sigma[ii])
#ax1.plot(linefit_x, (linefit*norm_civ))


ax2.plot(J12_second_wave, J12_second_flux, '-r', lw=linewidth, label='J1205+3422 (58538')
ax3.plot(J12_third_wave,  J12_third_flux,  '-k', lw=linewidth, label='J1205+3422 (58693')

## ACTUALLY PLOTTING   J16
ax4.plot(J16_first_wave,  J16_first_flux,  '-b', lw=linewidth, label='J1638+2827 (54553)')
ax5.plot(J16_second_wave, J16_second_flux, '-r', lw=linewidth, label='J1638+2827 (55832)')
ax6.plot(J16_third_wave,  J16_third_flux,  '-k', lw=linewidth, label='J1638+2827 (58583)')

## ACTUALLY PLOTTING   J22
ax7.plot(J22_first_wave,  J22_first_flux,  '-b', lw=linewidth, label='J2228+2201 (56189)')
ax8.plot(J22_second_wave, J22_second_flux, '-r', lw=linewidth, label='J2228+2201 (56960)')
ax9.plot(J22_third_wave,  J22_third_flux,  '-k', lw=linewidth, label='J2228+2201 (58693)')



ax8.set_xlabel(r'Restframe wavelength / ${ \rm \AA}$', fontsize=fontsize*1.4)
ax4.set_ylabel(r'F$_{\lambda}$ / 10$^{-17}$ erg cm$^{-2}$ s$^{-1}$ ${ \rm \AA }^{-1}$', fontsize=fontsize*1.4)


#handles=[neo_w1, neo_w2, crts, ztf_g, ztg_r]
leg = ax1.legend(loc='upper right', fontsize=fontsize/1.25, frameon=True, framealpha=1.0, fancybox=True)
leg = ax2.legend(loc='upper right', fontsize=fontsize/1.25, frameon=True, framealpha=1.0, fancybox=True)
leg = ax3.legend(loc='upper right', fontsize=fontsize/1.25, frameon=True, framealpha=1.0, fancybox=True)
leg = ax4.legend(loc='upper right', fontsize=fontsize/1.25, frameon=True, framealpha=1.0, fancybox=True)
leg = ax5.legend(loc='upper right', fontsize=fontsize/1.25, frameon=True, framealpha=1.0, fancybox=True)
leg = ax6.legend(loc='upper right', fontsize=fontsize/1.25, frameon=True, framealpha=1.0, fancybox=True)
leg = ax7.legend(loc='upper right', fontsize=fontsize/1.25, frameon=True, framealpha=1.0, fancybox=True)
leg = ax8.legend(loc='upper right', fontsize=fontsize/1.25, frameon=True, framealpha=1.0, fancybox=True)
leg = ax9.legend(loc='upper right', fontsize=fontsize/1.25, frameon=True, framealpha=1.0, fancybox=True)


plt.savefig('CIV_fits_temp.png', format='png')
plt.close(fig)
