'''
      J 1 2 0 5 + 3 4 2 2  

12h05m44.7s +34d22m52.4s
181.436164  +34.381229
'''

import numpy as np
import seaborn as sns

import matplotlib as mpl
import matplotlib.lines    as mlines
import matplotlib.pyplot   as plt
import matplotlib.gridspec as gridspec

from matplotlib          import colors as mcolors
from matplotlib.gridspec import GridSpec
from matplotlib.patches  import Rectangle

from astropy.io          import ascii
from astropy.io          import fits
from astropy.convolution import convolve, Box1DKernel

import pylustrator
#pylustrator.start()

##
## R E A D I N G   I N   T H E    D A T A 
##
path = '/cos_pc19a_npr/programs/quasars/CIV_CLQs/data/J1205+3422/'

##    S P E C T R A 
infile   = 'DBSP_J1205p3422_b_58538.dat'
DBSP_b1  = ascii.read(path+infile)
infile   = 'DBSP_J1205p3422_r_58538.dat'
DBSP_r1  = ascii.read(path+infile)
infile   = 'DBSP_J1205p3422_b_58694.dat'
DBSP_b2  = ascii.read(path+infile)
infile   = 'DBSP_J1205p3422_r_58694.dat'
DBSP_r2  = ascii.read(path+infile)

## The Blue and Red arms overlap in wavelength; Blue has better SNR than RED
## in overlapping regions
DBSP_r1 = DBSP_r1[np.where(DBSP_r1['wavelength'] > DBSP_b1['wavelength'].max())]
DBSP_r2 = DBSP_r2[np.where(DBSP_r2['wavelength'] > DBSP_b2['wavelength'].max())]

DBSP_b1_wavelength = DBSP_b1['wavelength']
DBSP_r1_wavelength = DBSP_r1['wavelength']
DBSP_b2_wavelength = DBSP_b2['wavelength']
DBSP_r2_wavelength = DBSP_r2['wavelength']

## Need to convert BUNIT   = 'erg/cm2/s/Hz'  to  'erg/cm2/s/A'
## http://www.stsci.edu/~strolger/docs/UNITS.txt
## 
## [Y erg/cm^2/s/Hz]            = 1000 * [X W/m^2/Hz]
## [Y erg/cm^2/s/A]             = 2.99792458E+21 * [X1 W/m^2/Hz] / [X2 A]^2
## / 1e-17 to put into 'regular' SDSS flux units 
DBSP_b1_flux   = np.array((2.99792458E+21 * (DBSP_b1['flux_density'] / 1000.)) / (DBSP_b1['wavelength']**2) / 1e-17)
DBSP_r1_flux   = np.array((2.99792458E+21 * (DBSP_r1['flux_density'] / 1000.)) / (DBSP_r1['wavelength']**2) / 1e-17)
DBSP_b2_flux   = np.array((2.99792458E+21 * (DBSP_b2['flux_density'] / 1000.)) / (DBSP_b2['wavelength']**2) / 1e-17)
DBSP_r2_flux   = np.array((2.99792458E+21 * (DBSP_r2['flux_density'] / 1000.)) / (DBSP_r2['wavelength']**2) / 1e-17)

##    https://dr15.sdss.org/optical/spectrum/view?plateid=2089&mjd=53498&fiberid=427&run2d=26        
sdssname        = 'spec-2089-53498-0427.fits'
sdss_data       = fits.open(path+sdssname)
sdss_spectrum   = sdss_data[1].data
sdss_flux       = sdss_spectrum.flux
sdss_loglam     = sdss_spectrum.loglam
sdss_wavelength = 10**(sdss_loglam)
sdss_mjd        = sdss_data[2].data['MJD'][0]

## Get the redshift right!
redshift = 2.068
print()
print('Redshift of J1205+3422 is ', redshift)
print()

## Emission Line list
linelist_file = 'emission_lines.dat'
linelist = ascii.read(path+linelist_file)


##  S E T T I N G    U P   T H E   P L O T !!!
##  
## https://jakevdp.github.io/PythonDataScienceHandbook/04.08-multiple-subplots.html
fig = plt.figure(figsize=(8, 4.3))
gs = GridSpec(2, 3, wspace=0.3, hspace=0.2)

## Use ``python slicing syntax'', 
##   gs[which row, which columnn] 
##   with e.g. gs[r:, c] then letting you span the whole row,
##
## https://stackoverflow.com/questions/49323348/matplotlib-gridspec-how-to-specify-the-location-by-numbers
## https://images.app.goo.gl/VQ4ELKbWTdgJNLu39
##
## gs[0 , :]    i.e. row 0,  all columns
## gs[1 , :-1]  i.e. row 1, all columns except last
## gs[1:, -1]   i.e. row 1 untill last row, last column
## gs[-1, 0]    i.e. last row, 0th colum

ax = fig.add_subplot(111)
#ax1 = fig.add_subplot(gs[0:, 0:2]) 
ax2 = fig.add_subplot(gs[0, :])
ax3 = fig.add_subplot(gs[-1, 0])
ax4 = fig.add_subplot(gs[-1, 1])
ax5 = fig.add_subplot(gs[-1, -1])


## Adjusting the Whitespace for the plots
left   = 0.10   # the left side of the subplots of the figure
right  = 0.98   # the right side of the subplots of the figure
bottom = 0.14   # the bottom of the subplots of the figure
top    = 0.98   # the top of the subplots of the figure
wspace = 0.26   # the amount of width reserved for blank space between subplots
hspace = 0.26   # the amount of height reserved for white space between subplots

plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)

 
## Some NPR defaults
lw              = 1.0
ls              = 'solid'
ms              = 1.
ms_large        = ms*8.
alpha           = 1.0
fontsize        = 16
labelsize       = fontsize
tickwidth       = 2.0
linewidth       = 2.4
tickwidth       = 2.0
ticklength      = 6.0
ticklabelsize   = labelsize
majorticklength = 12
minorticklength = 6

##
##  P L O T T I N G    T H E    S P E C T R A 
##
## Turn off axis lines and ticks of the big subplot
ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['right'].set_color('none')
ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)

## full, main spectral plot
xmin = 990;  xmax = 3690    ##used to be 910,3590 
ymin = -9.95; ymax = 99.   
ax2.set_xlim([xmin,xmax])
ax2.set_ylim([ymin, ymax])

## LyA 
xmin_ax3 = 1150; xmax_ax3 = 1285    ## to match VbB01, Fig 7; had been 1080, 1380 previously
ymin_ax3 = -9.5; ymax_ax3 = 92.
ax3.set_xlim([xmin_ax3, xmax_ax3])
ax3.set_ylim([ymin_ax3, ymax_ax3])

## CIV Line
xmin_ax4 = 1468; xmax_ax4 = 1620    ## to match VbB01, Fig 7; had been 1495, 1585 previously
ymin_ax4 = -9.5; ymax_ax4 = 64.   
#ax3.set_xlim([xmin,xmax])
#ax3.set_ylim([ymin, ymax])
ax4.set_xlim([xmin_ax4, xmax_ax4])
ax4.set_ylim([ymin_ax4, ymax_ax4])

## CIII]
#xmin = 1828; xmax = 1985           ## to match VbB01, Fig 7
#ymin = -6.95; ymax = 35.   
#ax4.set_xlim([xmin,xmax])
#ax4.set_ylim([ymin, ymax])

## MgII
xmin_ax5 =  2680; xmax_ax5 = 2910   ## to match VbB01, Fig 7; had been 2659, 2949 previously
ymin_ax5 = -4.95; ymax_ax5 = 24.   
ax5.set_xlim([xmin_ax5, xmax_ax5])
ax5.set_ylim([ymin_ax5, ymax_ax5])


## Smoothing:: https://joseph-long.com/writing/AstroPy-boxcar/
## (with timesteps = np.arange(N)*dt  in above eg.)

DBSP_b1_flux_smoothed = convolve(DBSP_b1_flux, Box1DKernel(5))
DBSP_r1_flux_smoothed = convolve(DBSP_r1_flux, Box1DKernel(5))
DBSP_b2_flux_smoothed = convolve(DBSP_b2_flux, Box1DKernel(5))
DBSP_r2_flux_smoothed = convolve(DBSP_r2_flux, Box1DKernel(5))


## The Spectra curves on the RHS
ax2.plot(   sdss_wavelength/(1+redshift), sdss_flux,           '-b', lw=linewidth/4.0,  label='SDSS MJD 53498')
ax2.plot(DBSP_b1_wavelength/(1+redshift), DBSP_b1_flux_smoothed, '-r', lw=linewidth/4.0,  label='DBSP MJD 58538')
ax2.plot(DBSP_r1_wavelength/(1+redshift), DBSP_r1_flux_smoothed, '-r', lw=linewidth/4.0)
ax2.plot(DBSP_b2_wavelength/(1+redshift), DBSP_b2_flux_smoothed, '-k', lw=linewidth/4.0)
ax2.plot(DBSP_r2_wavelength/(1+redshift), DBSP_r2_flux_smoothed, '-k', lw=linewidth/4.0, label='DBSP MJD 58693')
for ll in range(len(linelist)):
    ax2.axvline(x=linelist['Wavelength'][ll], color='gray', linestyle='--', linewidth=linewidth/3.4)
    label = linelist['LineName'][ll]
    if (np.mod(ll,2) == 0): xylabel = ((linelist['Wavelength'][ll]), 89)
    if (np.mod(ll,2) == 1): xylabel = ((linelist['Wavelength'][ll]), 79)        
#    xylabel = ((linelist['Wavelength'][ll]), 89)
    print(ll, np.mod(ll,2))
    ax2.annotate(label, xy=xylabel, ha='center', va='center', rotation=0, fontsize=fontsize/1.6 )


## Just plotting the spectra again, the axes limits take care of the "zoom in"
ax3.plot(   sdss_wavelength/(1+redshift), sdss_flux,    '-b', lw=linewidth/4.0)
ax3.plot(DBSP_b1_wavelength/(1+redshift), DBSP_b1_flux, '-r', lw=linewidth/4.0)
ax3.plot(DBSP_b2_wavelength/(1+redshift), DBSP_b2_flux, '-k', lw=linewidth/4.0)
#ax3.plot(DBSP_b1_wavelength/(1+redshift), DBSP_b1_flux_smoothed, '-r', lw=linewidth/4.0)
#ax3.plot(DBSP_b2_wavelength/(1+redshift), DBSP_b2_flux_smoothed, '-k', lw=linewidth/4.0)
for ll in range(len(linelist)):
    ax3.axvline(x=linelist['Wavelength'][ll], color='gray', linestyle='--', linewidth=linewidth/3.4)
    label = linelist['LineName'][ll]
    xylabel = ((linelist['Wavelength'][ll]), 78)
    ax3.annotate(label, xy=xylabel, ha='center', va='center', rotation=0, fontsize=fontsize/1.6 )

ax4.plot(   sdss_wavelength/(1+redshift), sdss_flux,    '-b', lw=linewidth/4.0)
ax4.plot(DBSP_b1_wavelength/(1+redshift), DBSP_b1_flux, '-r', lw=linewidth/4.0)
ax4.plot(DBSP_b2_wavelength/(1+redshift), DBSP_b2_flux, '-k', lw=linewidth/4.0)
ax4.plot(DBSP_b1_wavelength/(1+redshift), DBSP_b1_flux, '-r', lw=linewidth/4.0)
ax4.plot(DBSP_b2_wavelength/(1+redshift), DBSP_b2_flux_smoothed, '-k', lw=linewidth/4.0)
#ax4.plot(DBSP_r1_wavelength/(1+redshift), DBSP_r1_flux, '-r', lw=linewidth/4.0)
#ax4.plot(DBSP_r2_wavelength/(1+redshift), DBSP_r2_flux_smoothed, '-k', lw=linewidth/4.0)
for ll in range(len(linelist)):
    ax4.axvline(x=linelist['Wavelength'][ll], color='gray', linestyle='--', linewidth=linewidth/3.4)
    label = linelist['LineName'][ll]
    xylabel = ((linelist['Wavelength'][ll]), ymax_ax4/1.1)
    ax4.annotate(label, xy=xylabel, ha='center', va='center', rotation=0, fontsize=fontsize/1.6 )

ax5.plot(   sdss_wavelength/(1+redshift),    sdss_flux,    '-b', lw=linewidth/4.0)
ax5.plot(DBSP_r1_wavelength/(1+redshift), DBSP_r1_flux, '-r', lw=linewidth/4.0)
ax5.plot(DBSP_r2_wavelength/(1+redshift), DBSP_r2_flux, '-k', lw=linewidth/4.0)
#ax5.plot(DBSP_r1_wavelength/(1+redshift), DBSP_r1_flux_smoothed, '-r', lw=linewidth/4.0)
#ax5.plot(DBSP_r2_wavelength/(1+redshift), DBSP_r2_flux_smoothed, '-k', lw=linewidth/4.0)
for ll in range(len(linelist)):
    ax5.axvline(x=linelist['Wavelength'][ll], color='gray', linestyle='--', linewidth=linewidth/3.4)
    label = linelist['LineName'][ll]
    xylabel = ((linelist['Wavelength'][ll]), ymax_ax5/1.15)
    ax5.annotate(label, xy=xylabel, ha='center', va='center', rotation=0, fontsize=fontsize/1.6 )

   
##  A X E S   L A B E L S
plt.rcParams['text.usetex'] = True

#ax2.set_xlabel(r'Restframe wavelength', fontsize=fontsize)
ax.set_ylabel(r'F$_{\lambda}$ / 10$^{-17}$ erg cm$^{-2}$ s$^{-1}$ ${ \rm \AA }^{-1}$', fontsize=fontsize)
#ax2.set_ylabel(r'F$_{\lambda}$ / 10$^{-17}$ erg cm$^{-2}$ s$^{-1}$ ${ \rm \AA }^{-1}$', fontsize=fontsize/1.2)
ax4.set_xlabel(r'Restframe wavelength', fontsize=fontsize)

## Axes style
ax2.tick_params(labelsize=labelsize/1.15, direction='in', length=ticklength, width=tickwidth)
ax3.tick_params(labelsize=labelsize/1.6,  direction='in', length=ticklength, width=tickwidth)
ax4.tick_params(labelsize=labelsize/1.6,  direction='in', length=ticklength, width=tickwidth)
ax5.tick_params(labelsize=labelsize/1.6,  direction='in', length=ticklength, width=tickwidth)

##
##  L E G E N D S
##
leg = ax2.legend(loc='lower right', fontsize=fontsize/1.2,
                     frameon=True, framealpha=1.0, fancybox=True)
for line in leg.get_lines(): line.set_linewidth(2.2)

## Quasar Lable Name
ax2.text(3000, ymax/1.4, r'{\bf J1205+3422}',  fontsize=fontsize*1.2)

plt.savefig('J1205+3422_landscape_spectra_temp.png', format='png')
plt.close(fig)

# show the plot in a pylustrator window
#plt.show()
