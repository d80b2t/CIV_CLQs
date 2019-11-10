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
J16_second = fits.open(path+bossname)
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
J22_second = fits.open(path+bossname)
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
linewidth       = 1.0 ## 1.4
beefup          = 2.2
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

ymin_J16 = -2.9; ymax_J16 = 31.9     ## 21.9 in the 'main' figure; larger here to accomodate legend
ax4.set_ylim([ymin_J16,  ymax_J16/1.2])
ax5.set_ylim([ymin_J16,  ymax_J16/2.4])
ax6.set_ylim([ymin_J16,  ymax_J16])

ymin_J22 = -4.9; ymax_J22 = 24.9     ## 12.9 in the 'main' figure; larger here to accomodate legend
ax7.set_ylim([ymin_J22/2.2, ymax_J22/2.8])
ax8.set_ylim([ymin_J22/2.2, ymax_J22])
ax9.set_ylim([ymin_J22,     ymax_J22])

ymin = -0.1; ymax = 1.1
#ax1.set_ylim([ymin, ymax]) ; #ax2.set_ylim([ymin, ymax]) ; #ax3.set_ylim([ymin, ymax])
#ax4.set_ylim([ymin, ymax]) ; #ax5.set_ylim([ymin, ymax]) ; #ax6.set_ylim([ymin, ymax])
#ax7.set_ylim([ymin, ymax]) ; #ax8.set_ylim([ymin, ymax]) ; #ax9.set_ylim([ymin, ymax])

## numerical constants :
c        = 3.e5                              ## SoL in  km/s
conve    = 2. * np.sqrt(2 * np.log(2))       ## FWHM = 2.sqrt(2.ln(2)) \sigma ~  2.355 \sigma
sqrt_2pi = np.sqrt(2*np.pi)                  ## 2.50663
center     = 1549.48                         ## line center for CIV, Angstrom

print()
print("qsfit['FILENAME'][ii], '     ',  qsfit['LUM'][ii], qsfit['FWHM'][ii], qsfit['VOFF'][ii]")
for ii in range(len(qsfit)):
    print(ii, qsfit['FILENAME'][ii], '       ', qsfit['NCOMP'][ii],  qsfit['LUM'][ii], qsfit['FWHM'][ii], qsfit['VOFF'][ii])
    if ii == 0:
        #j=0; k=0
        norm       = qsfit['LUM'][ii]                ## 10^42 erg s^-1
        fwhm       = qsfit['FWHM'][ii]               ## km/s
        v_off      = qsfit['VOFF'][ii]               ## km/s
        cont_off   = 27.                             ## a sophisticated guess at the continnuum level.
        
        x0         = center - ((v_off / c) * center) ##  Angstrom (actual center of the line profile)
        sigma      = (fwhm  / c) * center / conve    ##;  Angstrom (width of the line profile)
        peak_value = norm / sqrt_2pi  / sigma        ##  10^42 erg s^-1 A^-1 (value at the peak)
    
        line_profile = (peak_value * (np.exp( -((J12_first_wave - x0) / sigma)**2 / 2.)))+cont_off

        ax1.plot(J12_first_wave,  J12_first_flux, '-b', lw=linewidth)
        ax1.plot(J12_first_wave,  line_profile,   '-b', lw=beefup)

        norm_str = str(np.around(norm, decimals=3))
        fwhm_str = str(np.around(fwhm, decimals=3))
        voff_str = str(np.around(v_off, decimals=3))
        obj_lab  = mlines.Line2D([], [], label=r'J1205+3422 (53498)',      color='b')
        norm_lab = mlines.Line2D([], [], label=r'Lum.='+norm_str,          color='b', linestyle='None')
        fwhm_lab = mlines.Line2D([], [], label=r'FWHM='+fwhm_str,          color='b', linestyle='None')
        voff_lab = mlines.Line2D([], [], label=r'V$_{\rm off}$='+voff_str, color='b', linestyle='None')
        
        handles  = [obj_lab, norm_lab, fwhm_lab, voff_lab]
        ax1.legend(loc='lower right', handles=handles,
        fontsize=fontsize/1.4, frameon=True, framealpha=1.0,
        fancybox=True)
        
    if ii == 1:
        norm       = qsfit['LUM'][ii]                ## 10^42 erg s^-1
        fwhm       = qsfit['FWHM'][ii]               ## km/s
        v_off      = qsfit['VOFF'][ii]               ## km/s
        cont_off   = 2.5                              ## a sophisticated guess at the continnuum level.
        
        x0           = center - (v_off / c) * center   ##  Angstrom (actual center of the line profile)
        sigma        = (fwhm  / c) * center / conve    ##;  Angstrom (width of the line profile)
        peak_value   = norm / sqrt_2pi  / sigma        ##  10^42 erg s^-1 A^-1 (value at the peak)
        line_profile = (peak_value * (np.exp( -((J12_second_wave - x0) / sigma)**2 / 2.)))+cont_off

        ax2.plot(J12_second_wave, J12_second_flux, '-r', lw=linewidth)
        ax2.plot(J12_second_wave, line_profile,    '-r', lw=beefup)

        norm_str = str(np.around(norm, decimals=3))
        fwhm_str = str(np.around(fwhm, decimals=3))
        voff_lab = str(np.around(v_off, decimals=3))
        obj_lab  = mlines.Line2D([], [], label=r'J1205+3422 (58538)',      color='r')
        norm_lab = mlines.Line2D([], [], label=r'Lum.='+norm_str,          color='r', linestyle='None')
        fwhm_lab = mlines.Line2D([], [], label=r'FWHM='+fwhm_str,          color='r', linestyle='None')
        voff_lab = mlines.Line2D([], [], label=r'V$_{\rm off}$='+voff_lab, color='r', linestyle='None')

        handles  = [obj_lab, norm_lab, fwhm_lab, voff_lab]
        ax2.legend(loc='upper right', handles=handles,
        fontsize=fontsize/1.4, frameon=True, framealpha=1.0,
        fancybox=True)

    if ii == 2:
        norm       = qsfit['LUM'][ii]                ## 10^42 erg s^-1
        fwhm       = qsfit['FWHM'][ii]               ## km/s
        v_off      = qsfit['VOFF'][ii]               ## km/s
        cont_off   = 2.5                             ## a sophisticated guess at the continnuum level.
        
        x0           = center - (v_off / c) * center   ##  Angstrom (actual center of the line profile)
        sigma        = (fwhm  / c) * center / conve    ##;  Angstrom (width of the line profile)
        peak_value   = norm / sqrt_2pi  / sigma        ##  10^42 erg s^-1 A^-1 (value at the peak)
        line_profile = (peak_value * (np.exp( -((J12_third_wave - x0) / sigma)**2 / 2.)))+cont_off
        
        ax3.plot(J12_third_wave, J12_third_flux, '-k', lw=linewidth)
        ax3.plot(J12_third_wave, line_profile,   '-k', lw=beefup)

        norm_str = str(np.around(norm, decimals=3))
        fwhm_str = str(np.around(fwhm, decimals=3))
        voff_lab = str(np.around(v_off, decimals=3))
        obj_lab  = mlines.Line2D([], [], label=r'J1205+3422 (58693)',      color='k')
        norm_lab = mlines.Line2D([], [], label=r'Lum.='+norm_str,          color='k', linestyle='None')
        fwhm_lab = mlines.Line2D([], [], label=r'FWHM='+fwhm_str,          color='k', linestyle='None')
        voff_lab = mlines.Line2D([], [], label=r'V$_{\rm off}$='+voff_lab, color='k', linestyle='None')

        handles  = [obj_lab, norm_lab, fwhm_lab, voff_lab]
        ax3.legend(loc='upper right', handles=handles,
        fontsize=fontsize/1.4, frameon=True, framealpha=1.0,
        fancybox=True)
        print()
        
## ACTUALLY PLOTTING   J16
    if ii == 3:
        norm       = qsfit['LUM'][ii]                ## 10^42 erg s^-1
        fwhm       = qsfit['FWHM'][ii]               ## km/s
        v_off      = qsfit['VOFF'][ii]               ## km/s
        cont_off   = 2.5                              ## a sophisticated guess at the continnuum level.
        
        x0         = center - (v_off / c) * center   ##  Angstrom (actual center of the line profile)
        sigma      = (fwhm  / c) * center / conve    ##;  Angstrom (width of the line profile)
        peak_value = norm / sqrt_2pi  / sigma        ##  10^42 erg s^-1 A^-1 (value at the peak)
    
        line_profile = (peak_value * (np.exp( -((J16_first_wave - x0) / sigma)**2 / 2.)))+cont_off
        ax4.plot(J16_first_wave, J16_first_flux, '-b', lw=linewidth, label='J1638+2827 (54553)')
        ax4.plot(J16_first_wave, line_profile,   '-b', lw=beefup)

        norm_str = str(np.around(norm, decimals=3))
        fwhm_str = str(np.around(fwhm, decimals=3))
        voff_lab = str(np.around(v_off, decimals=3))
        obj_lab  = mlines.Line2D([], [], label=r'J1638+2827 (54553)',      color='b')
        norm_lab = mlines.Line2D([], [], label=r'Lum.='+norm_str,          color='b', linestyle='None')
        fwhm_lab = mlines.Line2D([], [], label=r'FWHM='+fwhm_str,          color='b', linestyle='None')
        voff_lab = mlines.Line2D([], [], label=r'V$_{\rm off}$='+voff_lab, color='b', linestyle='None')

        handles  = [obj_lab, norm_lab, fwhm_lab, voff_lab]
        ax4.legend(loc='upper right', handles=handles,
        fontsize=fontsize/1.4, frameon=True, framealpha=1.0,
        fancybox=True)
        
    if ii == 4:
        norm       = qsfit['LUM'][ii]                ## 10^42 erg s^-1
        fwhm       = qsfit['FWHM'][ii]               ## km/s
        v_off      = qsfit['VOFF'][ii]               ## km/s
        cont_off   = 1.5                              ## a sophisticated guess at the continnuum level.
        
        x0         = center - (v_off / c) * center   ##  Angstrom (actual center of the line profile)
        sigma      = (fwhm  / c) * center / conve    ##;  Angstrom (width of the line profile)
        peak_value = norm / sqrt_2pi  / sigma        ##  10^42 erg s^-1 A^-1 (value at the peak)
    
        line_profile = (peak_value * (np.exp( -((J16_second_wave - x0) / sigma)**2 / 2.)))+cont_off
        ax5.plot(J16_second_wave, J16_second_flux, '-r', lw=linewidth)
        ax5.plot(J16_second_wave, line_profile,    '-r', lw=beefup)

        norm_str = str(np.around(norm, decimals=3))
        fwhm_str = str(np.around(fwhm, decimals=3))
        voff_lab = str(np.around(v_off, decimals=3))
        obj_lab  = mlines.Line2D([], [], label=r'J1638+2827 (55832)',      color='r')
        norm_lab = mlines.Line2D([], [], label=r'Lum.='+norm_str,          color='r', linestyle='None')
        fwhm_lab = mlines.Line2D([], [], label=r'FWHM='+fwhm_str,          color='r', linestyle='None')
        voff_lab = mlines.Line2D([], [], label=r'V$_{\rm off}$='+voff_lab, color='r', linestyle='None')

        handles  = [obj_lab, norm_lab, fwhm_lab, voff_lab]
        ax5.legend(loc='upper right', handles=handles,
        fontsize=fontsize/1.4, frameon=True, framealpha=1.0,
        fancybox=True)
       
    if ii == 5:
        norm       = qsfit['LUM'][ii]                ## 10^42 erg s^-1
        fwhm       = qsfit['FWHM'][ii]               ## km/s
        v_off      = qsfit['VOFF'][ii]               ## km/s
        cont_off   = 5.                              ## a sophisticated guess at the continnuum level.
        
        x0         = center - (v_off / c) * center   ##  Angstrom (actual center of the line profile)
        sigma      = (fwhm  / c) * center / conve    ##;  Angstrom (width of the line profile)
        peak_value = norm / sqrt_2pi  / sigma        ##  10^42 erg s^-1 A^-1 (value at the peak)
        line_profile = (peak_value * (np.exp( -((J16_third_wave - x0) / sigma)**2 / 2.)))+cont_off
            
        ax6.plot(J16_third_wave,  J16_third_flux, '-k', lw=linewidth)
        ax6.plot(J16_third_wave,  line_profile,   '-k', lw=beefup)
        
        norm_str = str(np.around(norm, decimals=3))
        fwhm_str = str(np.around(fwhm, decimals=3))
        voff_lab = str(np.around(v_off, decimals=3))
        obj_lab  = mlines.Line2D([], [], label=r'J1638+2827 (58583)',      color='k')
        norm_lab = mlines.Line2D([], [], label=r'Lum.='+norm_str,          color='k', linestyle='None')
        fwhm_lab = mlines.Line2D([], [], label=r'FWHM='+fwhm_str,          color='k', linestyle='None')
        voff_lab = mlines.Line2D([], [], label=r'V$_{\rm off}$='+voff_lab, color='k', linestyle='None')

        handles  = [obj_lab, norm_lab, fwhm_lab, voff_lab]
        ax6.legend(loc='upper right', handles=handles,
        fontsize=fontsize/1.4, frameon=True, framealpha=1.0,
        fancybox=True)
        print()
         
## ACTUALLY PLOTTING   J22
    if ii == 6:
        norm       = qsfit['LUM'][ii]                ## 10^42 erg s^-1
        fwhm       = qsfit['FWHM'][ii]               ## km/s
        v_off      = qsfit['VOFF'][ii]               ## km/s
        cont_off   = 0.5                              ## a sophisticated guess at the continnuum level.
        
        x0           = center - (v_off / c) * center   ##  Angstrom (actual center of the line profile)
        sigma        = (fwhm  / c) * center / conve    ##;  Angstrom (width of the line profile)
        peak_value   = norm / sqrt_2pi  / sigma        ##  10^42 erg s^-1 A^-1 (value at the peak)
        line_profile = (peak_value * (np.exp( -((J22_first_wave - x0) / sigma)**2 / 2.)))+cont_off

        ax7.plot(J22_first_wave,  J22_first_flux, '-b', lw=linewidth)
        ax7.plot(J22_first_wave,  line_profile,   '-b', lw=beefup)

        norm_str = str(np.around(norm, decimals=3))
        fwhm_str = str(np.around(fwhm, decimals=3))
        voff_lab = str(np.around(v_off, decimals=3))
        obj_lab  = mlines.Line2D([], [], label=r'J2228+2201 (56189)',      color='b')
        norm_lab = mlines.Line2D([], [], label=r'Lum.='+norm_str,          color='b', linestyle='None')
        fwhm_lab = mlines.Line2D([], [], label=r'FWHM='+fwhm_str,          color='b', linestyle='None')
        voff_lab = mlines.Line2D([], [], label=r'V$_{\rm off}$='+voff_lab, color='b', linestyle='None')

        handles  = [obj_lab, norm_lab, fwhm_lab, voff_lab]
        ax7.legend(loc='upper right', handles=handles,
        fontsize=fontsize/1.4, frameon=True, framealpha=1.0,
        fancybox=True)

        
    if ii == 7:
        norm       = qsfit['LUM'][ii]                ## 10^42 erg s^-1
        fwhm       = qsfit['FWHM'][ii]               ## km/s
        v_off      = qsfit['VOFF'][ii]               ## km/s
        cont_off   = 4.5                             ## a sophisticated guess at the continnuum level.
        
        x0           = center - (v_off / c) * center   ##  Angstrom (actual center of the line profile)
        sigma        = (fwhm  / c) * center / conve    ##;  Angstrom (width of the line profile)
        peak_value   = norm / sqrt_2pi  / sigma        ##  10^42 erg s^-1 A^-1 (value at the peak)
        line_profile = (peak_value * (np.exp( -((J22_second_wave - x0) / sigma)**2 / 2.)))+cont_off

        ax8.plot(J22_second_wave, J22_second_flux, '-r', lw=linewidth)
        ax8.plot(J22_second_wave, line_profile,    '-r', lw=beefup)

        norm_str = str(np.around(norm, decimals=3))
        fwhm_str = str(np.around(fwhm, decimals=3))
        voff_lab = str(np.around(v_off, decimals=3))
        obj_lab  = mlines.Line2D([], [], label=r'J2228+2201 (56960)',      color='r')
        norm_lab = mlines.Line2D([], [], label=r'Lum.='+norm_str,          color='r', linestyle='None')
        fwhm_lab = mlines.Line2D([], [], label=r'FWHM='+fwhm_str,          color='r', linestyle='None')
        voff_lab = mlines.Line2D([], [], label=r'V$_{\rm off}$='+voff_lab, color='r', linestyle='None')

        handles  = [obj_lab, norm_lab, fwhm_lab, voff_lab]
        ax8.legend(loc='upper right', handles=handles,
        fontsize=fontsize/1.4, frameon=True, framealpha=1.0,
        fancybox=True)

        
    if ii == 8:
        norm       = qsfit['LUM'][ii]                ## 10^42 erg s^-1
        fwhm       = qsfit['FWHM'][ii]               ## km/s
        v_off      = qsfit['VOFF'][ii]               ## km/s
        cont_off   = 1.5                              ## a sophisticated guess at the continnuum level.
        
        x0           = center - (v_off / c) * center   ##  Angstrom (actual center of the line profile)
        sigma        = (fwhm  / c) * center / conve    ##;  Angstrom (width of the line profile)
        peak_value   = norm / sqrt_2pi  / sigma        ##  10^42 erg s^-1 A^-1 (value at the peak)
        line_profile = (peak_value * (np.exp( -((J22_third_wave - x0) / sigma)**2 / 2.)))+cont_off

        ax9.plot(J22_third_wave,  J22_third_flux, '-k', lw=linewidth)
        ax9.plot(J22_third_wave,  line_profile,   '-k', lw=beefup)

        norm_str = str(np.around(norm, decimals=3))
        fwhm_str = str(np.around(fwhm, decimals=3))
        voff_lab = str(np.around(v_off, decimals=3))
        obj_lab  = mlines.Line2D([], [], label=r'J2228+2201 (58693)',      color='k')
        norm_lab = mlines.Line2D([], [], label=r'Lum.='+norm_str,          color='k', linestyle='None')
        fwhm_lab = mlines.Line2D([], [], label=r'FWHM='+fwhm_str,          color='k', linestyle='None')
        voff_lab = mlines.Line2D([], [], label=r'V$_{\rm off}$='+voff_lab, color='k', linestyle='None')

        handles  = [obj_lab, norm_lab, fwhm_lab, voff_lab]
        ax9.legend(loc='upper right', handles=handles,
        fontsize=fontsize/1.4, frameon=True, framealpha=1.0,
        fancybox=True)

        
ax8.set_xlabel(r'Restframe wavelength / ${ \rm \AA}$', fontsize=fontsize*1.6)
ax4.set_ylabel(r'F$_{\lambda}$ / 10$^{-17}$ erg cm$^{-2}$ s$^{-1}$ ${ \rm \AA }^{-1}$', fontsize=fontsize*1.6)

#handles=[neo_w1, neo_w2, crts, ztf_g, ztg_r]

#leg = ax2.legend(loc='upper right', fontsize=fontsize/1.25, frameon=True, framealpha=1.0, fancybox=True)
#leg = ax3.legend(loc='upper right', fontsize=fontsize/1.25, frameon=True, framealpha=1.0, fancybox=True)
#leg = ax4.legend(loc='upper right', fontsize=fontsize/1.25, frameon=True, framealpha=1.0, fancybox=True)
#leg = ax5.legend(loc='upper right', fontsize=fontsize/1.25, frameon=True, framealpha=1.0, fancybox=True)
#leg = ax6.legend(loc='upper right', fontsize=fontsize/1.25, frameon=True, framealpha=1.0, fancybox=True)
#leg = ax7.legend(loc='upper right', fontsize=fontsize/1.25, frameon=True, framealpha=1.0, fancybox=True)
#leg = ax8.legend(loc='upper right', fontsize=fontsize/1.25, frameon=True, framealpha=1.0, fancybox=True)
#leg = ax9.legend(loc='upper right', fontsize=fontsize/1.25, frameon=True, framealpha=1.0, fancybox=True)


plt.savefig('CIV_fits_temp.png', format='png')
plt.close(fig)
