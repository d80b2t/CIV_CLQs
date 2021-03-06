'''
12h05m44.7s +34d22m52.4s
'''

import numpy as np
import seaborn as sns
import matplotlib.lines    as mlines
import matplotlib.pyplot   as plt
import matplotlib.gridspec as gridspec

from matplotlib          import colors as mcolors
from matplotlib.gridspec import GridSpec
from matplotlib.patches  import Rectangle


from astropy.io          import ascii
from astropy.io          import fits
from astropy.convolution import convolve, Box1DKernel

#import pylustrator
#pylustrator.start()

##
## R E A D I N G   I N   T H E    D A T A 
##
path = '/cos_pc19a_npr/programs/quasars/CIV_CLQs/data/J1638+2827/'

##  Light Curve
infile = 'NEOWISER-L1b_J1638+2827.dat'
NEOWISER = ascii.read(path+infile)
NEOWISER_W1_AB = NEOWISER['w1mpro'] + 2.673
NEOWISER_W2_AB = NEOWISER['w2mpro'] + 3.313
infile = 'NEOWISER-L1b_J1638+2827_averaged.dat'
NEOWISER_aver = ascii.read(path+infile)
NEOWISER_aver_W1_AB = NEOWISER_aver['w1mpro_wgt'] + 2.673
NEOWISER_aver_W2_AB = NEOWISER_aver['w2mpro_wgt'] + 3.313

##  Spectra
infile   = 'LRIS_J1638p2827_b_58583.dat'
LRIS_b  = ascii.read(path+infile)
infile   = 'LRIS_J1638p2827_r_58583.dat'
LRIS_r  = ascii.read(path+infile)

LRIS_b_wavelength = LRIS_b['wavelength']
LRIS_r_wavelength = LRIS_r['wavelength']

## Need to convert BUNIT   = 'erg/cm2/s/Hz'  to  'erg/cm2/s/A'
## http://www.stsci.edu/~strolger/docs/UNITS.txt
## 
## [Y erg/cm^2/s/Hz]            = 1000 * [X W/m^2/Hz]
## [Y erg/cm^2/s/A]             = 2.99792458E+21 * [X1 W/m^2/Hz] / [X2 A]^2
## / 1e-17 to put into 'regular' SDSS flux units 
LRIS_b_flux   = (2.99792458E+21 * (LRIS_b['flux_density'] / 1000.)) / (LRIS_b['wavelength']**2)  / 1e-17
LRIS_r_flux   = (2.99792458E+21 * (LRIS_r['flux_density'] / 1000.)) / (LRIS_r['wavelength']**2)  / 1e-17

##    SDSS
sdssname        = 'spec-2948-54553-0614.fits'
sdss_data       = fits.open(path+sdssname)
sdss_spectrum   = sdss_data[1].data
sdss_flux       = sdss_spectrum.flux
sdss_loglam     = sdss_spectrum.loglam
sdss_wavelength = 10**(sdss_loglam)
sdss_mjd        = sdss_data[2].data['MJD'][0]

bossname        = 'spec-5201-55832-0178_v5_10_0.fits'
boss_data       = fits.open(path+bossname)
boss_spectrum   = boss_data[1].data
boss_flux       = boss_spectrum.flux
boss_loglam     = boss_spectrum.loglam
boss_wavelength = 10**(boss_loglam)
boss_mjd        = boss_data[2].data['MJD'][0]
    
## Get the redshift right!
redshift = 2.182
print()
print('Redshift of J1638+2827 is ', redshift)
print()

## Emission Line list
linelist_file = 'emission_lines.dat'
linelist = ascii.read(path+linelist_file)


## https://jakevdp.github.io/PythonDataScienceHandbook/04.08-multiple-subplots.html
fig = plt.figure(figsize=(12, 5))
gs = GridSpec(2, 5, wspace=0.6, hspace=0.3)

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

ax1 = fig.add_subplot(gs[0:, 0:2]) 
ax2 = fig.add_subplot(gs[0, 2:5])
ax3 = fig.add_subplot(gs[-1, 2])
ax4 = fig.add_subplot(gs[-1, 3])
ax5 = fig.add_subplot(gs[-1, -1])


## Adjusting the Whitespace for the plots
left   = 0.06   # the left side of the subplots of the figure
right  = 0.97   # the right side of the subplots of the figure
bottom = 0.16   # the bottom of the subplots of the figure
top    = 0.96   # the top of the subplots of the figure
wspace = 0.26   # the amount of width reserved for blank space between subplots
hspace = 0.26   # the amount of height reserved for white space between subplots

plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)

 
## Some NPR defaults
ms              = 10.
ms_big          = ms*6.
ms_large        = ms*8.
ls              = 'solid'
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


##  T H E    L I G H T     C U R V E S 
xmin = 54500; xmax = 59000
ymin = 19.85; ymax = 15.10   
ax1.set_xlim([xmin,xmax])
ax1.set_ylim([ymin, ymax])

## NEOWISER W1/2 (AB)
ax1.scatter(NEOWISER['mjd'],  NEOWISER_W1_AB, color='k', alpha=alpha, s=ms*1.8)
ax1.scatter(NEOWISER['mjd'],  NEOWISER_W1_AB, color='r', alpha=alpha, s=ms, label='NEOWISE W1')
ax1.scatter(NEOWISER['mjd'],  NEOWISER_W2_AB, color='k', alpha=alpha, s=ms*1.8)
ax1.scatter(NEOWISER['mjd'],  NEOWISER_W2_AB, color='c', alpha=alpha, s=ms, label='NEOWISE W2')

ax1.scatter(NEOWISER_aver['mean_mjd'],  NEOWISER_aver_W1_AB, color='k', alpha=alpha, s=ms_big*1.8)
ax1.scatter(NEOWISER_aver['mean_mjd'],  NEOWISER_aver_W1_AB, color='r', alpha=alpha, s=ms_big)
ax1.scatter(NEOWISER_aver['mean_mjd'],  NEOWISER_aver_W2_AB, color='k', alpha=alpha, s=ms_big*1.8)
ax1.scatter(NEOWISER_aver['mean_mjd'],  NEOWISER_aver_W2_AB, color='c', alpha=alpha, s=ms_big)



xmin = 910;  xmax = 3590
ymin = -9.95; ymax = 48.   
ax2.set_xlim([xmin, xmax])
ax2.set_ylim([ymin, ymax])

## LyA 
xmin = 1095; xmax = 1330
ymin = -9.5; ymax = 48.
ax3.set_xlim([xmin,xmax])
ax3.set_ylim([ymin, ymax])

## CIV Line
xmin = 1480; xmax = 1620
ymin = -4.9; ymax = 14.9   
#ax3.set_xlim([xmin,xmax])
#ax3.set_ylim([ymin, ymax])
ax4.set_xlim([xmin,xmax])
ax4.set_ylim([ymin, ymax])

## CIII]
xmin = 1821; xmax = 1948
ymin = -6.95; ymax = 35.   
#ax4.set_xlim([xmin,xmax])
#ax4.set_ylim([ymin, ymax])

## MgII
xmin = 2700; xmax = 2900
ymin = -4.95; ymax = 10.   
ax5.set_xlim([xmin,xmax])
ax5.set_ylim([ymin, ymax])


## The Spectra curves on the RHS
ax2.plot(   sdss_wavelength/(1+redshift), sdss_flux,         '-b', lw=linewidth/4.0, label='SDSS MJD 54553')
ax2.plot(   boss_wavelength/(1+redshift), boss_flux,         '-r', lw=linewidth/4.0, label='BOSS MJD 55832')
ax2.plot(LRIS_b_wavelength/(1+redshift), LRIS_b_flux,        '-k', lw=linewidth/4.0, label='LRIS MJD 58538')
ax2.plot(LRIS_r_wavelength/(1+redshift), LRIS_r_flux,        '-k', lw=linewidth/4.0)
## Smoothing:: https://joseph-long.com/writing/AstroPy-boxcar/
#LRIS_b_flux_smoothed = convolve(LRIS_b_flux, Box1DKernel(5))
#LRIS_r_flux_smoothed = convolve(LRIS_r_flux, Box1DKernel(5))
#ax2.plot(LRIS_b_wavelength/(1+redshift), LRIS_b_flux_smoothed, '-k', lw=linewidth/4.0)
#ax2.plot(LRIS_r_wavelength/(1+redshift), LRIS_r_flux_smoothed, '-k', lw=linewidth/4.0,     label='LRIS MJD 58538')
## Putting on the Emission Line names and vertical lines
## e.g. https://scipython.com/blog/rotating-text-onto-a-line-in-matplotlib/
for ll in range(len(linelist)):
    ax2.axvline(x=linelist['Wavelength'][ll], color='gray', linestyle='--', linewidth=linewidth/3.4)
    label = linelist['LineName'][ll]
    xylabel = ((linelist['Wavelength'][ll]), 48)
    ax2.annotate(label, xy=xylabel, ha='center', va='center', rotation=90, fontsize=fontsize/1.6 )


## Just plotting the spectra again, the axes limits take care of the "zoom in"
ax3.plot(   sdss_wavelength/(1+redshift), sdss_flux,    '-b', lw=linewidth/4.0)
ax3.plot(   boss_wavelength/(1+redshift), boss_flux,    '-r', lw=linewidth/4.0)
ax3.plot(LRIS_b_wavelength/(1+redshift), LRIS_b_flux, '-k', lw=linewidth/4.0)
#ax3.plot(LRIS_b_wavelength/(1+redshift), LRIS_b_flux_smoothed, '-k', lw=linewidth/4.0)
for ll in range(len(linelist)):
    ax3.axvline(x=linelist['Wavelength'][ll], color='gray', linestyle='--', linewidth=linewidth/3.4)
    label = linelist['LineName'][ll]
    xylabel = ((linelist['Wavelength'][ll]), 22)
    ax3.annotate(label, xy=xylabel, ha='center', va='center', rotation=90, fontsize=fontsize/1.6 )

    
ax4.plot(   sdss_wavelength/(1+redshift), sdss_flux,  '-b', lw=linewidth/4.0)
ax4.plot(   boss_wavelength/(1+redshift), boss_flux,  '-r', lw=linewidth/4.0)
ax4.plot(LRIS_b_wavelength/(1+redshift), LRIS_b_flux, '-k', lw=linewidth/4.0)
ax4.plot(LRIS_r_wavelength/(1+redshift), LRIS_r_flux, '-k', lw=linewidth/4.0)
#ax4.plot(LRIS_b_wavelength/(1+redshift), LRIS_b_flux_smoothed, '-k', lw=linewidth/4.0)
#ax4.plot(LRIS_r_wavelength/(1+redshift), LRIS_r_flux_smoothed, '-k', lw=linewidth/4.0)
for ll in range(len(linelist)):
    ax4.axvline(x=linelist['Wavelength'][ll], color='gray', linestyle='--', linewidth=linewidth/3.4)
    label = linelist['LineName'][ll]
    xylabel = ((linelist['Wavelength'][ll]), 12)
    ax4.annotate(label, xy=xylabel, ha='center', va='center', rotation=90, fontsize=fontsize/1.6 )


ax5.plot(   sdss_wavelength/(1+redshift),    sdss_flux,    '-b', lw=linewidth/4.0)
ax5.plot(   boss_wavelength/(1+redshift),    boss_flux,    '-r', lw=linewidth/4.0)
ax5.plot(LRIS_r_wavelength/(1+redshift), LRIS_r_flux, '-k', lw=linewidth/4.0)
#ax5.plot(LRIS_r_wavelength/(1+redshift), LRIS_r_flux_smoothed, '-k', lw=linewidth/4.0)
for ll in range(len(linelist)):
    ax5.axvline(x=linelist['Wavelength'][ll], color='gray', linestyle='--', linewidth=linewidth/3.4)
    label = linelist['LineName'][ll]
    xylabel = ((linelist['Wavelength'][ll]), 12)
    ax5.annotate(label, xy=xylabel, ha='center', va='center', rotation=90, fontsize=fontsize/1.6 )


## Axes labels
plt.rcParams['text.usetex'] = True

ax1.set_xlabel(r'MJD',          fontsize=fontsize/1.2)
ax1.set_ylabel(r'AB Magnitude', fontsize=fontsize/1.2)
ax2.set_xlabel(r'Restframe wavelength', fontsize=fontsize)
ax2.set_ylabel(r'F$_{\lambda}$ / 10$^{-17}$ erg cm$^{-2}$ s$^{-1}$ ${ \rm \AA }^{-1}$', fontsize=fontsize/1.4)
ax4.set_xlabel(r'Restframe wavelength', fontsize=fontsize)

## Axes style
ax2.tick_params(labelsize=labelsize/1.15, direction='in', length=ticklength,   width=tickwidth)
ax3.tick_params(labelsize=labelsize/1.6, direction='in', length=ticklength,   width=tickwidth)
ax4.tick_params(labelsize=labelsize/1.6, direction='in', length=ticklength,   width=tickwidth)
ax5.tick_params(labelsize=labelsize/1.6, direction='in', length=ticklength,   width=tickwidth)



#% start: automatic generated code from pylustrator
plt.figure(1).ax_dict = {ax.get_label(): ax for ax in plt.figure(1).axes}
import matplotlib as mpl
plt.figure(1).axes[0].get_xaxis().get_label().set_fontsize(fontsize)
plt.figure(1).axes[0].get_yaxis().get_label().set_fontsize(fontsize)
plt.figure(1).axes[0].set_position([0.060000, 0.136250, 0.345, 0.823750])
plt.figure(1).axes[0].xaxis.labelpad = 10.0


plt.figure(1).axes[1].get_xaxis().get_label().set_fontsize(fontsize/1.2)
plt.figure(1).axes[1].get_yaxis().get_label().set_fontsize(fontsize)

x_nudge = 0.012
plt.figure(1).axes[1].get_yaxis().get_label().set_rotation(90.0)
plt.figure(1).axes[1].set_position([0.464541+x_nudge, 0.560924, 0.499459, 0.399076])
plt.figure(1).axes[1].xaxis.labelpad = 323.200000

plt.figure(1).axes[1].yaxis.labelpad = 4.720000
plt.figure(1).axes[1].yaxis.set_label_coords(-0.045, -0.1)

plt.figure(1).axes[2].get_xaxis().get_label().set_text("")
plt.figure(1).axes[2].set_position([0.464541+x_nudge, 0.135000, 0.146836, 0.372826])
plt.figure(1).axes[3].set_position([0.643228+x_nudge, 0.135000, 0.146836, 0.372826])
plt.figure(1).axes[3].xaxis.labelpad = 10.0
plt.figure(1).axes[4].set_position([0.817164+x_nudge, 0.135000, 0.146836, 0.372826])
#% end: automatic generated code from pylustrator
#plt.show()

## LIGHT CURVE Legend
neo_w1 = mlines.Line2D([], [], label='NEOWISE W1', color='red',  marker="o", markeredgecolor='k', markeredgewidth=2.0, markersize=12, linestyle='None')
neo_w2 = mlines.Line2D([], [], label='NEOWISE W2', color='cyan', marker="o", markeredgecolor='k', markeredgewidth=2.0, markersize=12, linestyle='None')
leg = ax1.legend(loc='upper left', fontsize=fontsize/1.2, handles=[neo_w1, neo_w2])

## SPECTRA Legend
leg = ax2.legend(loc='upper right', fontsize=fontsize/1.2)
for line in leg.get_lines(): line.set_linewidth(2.2)

plt.savefig('J1638+2827_temp.png', format='png')
plt.close(fig)
