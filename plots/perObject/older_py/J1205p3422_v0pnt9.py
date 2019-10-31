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
path = '/cos_pc19a_npr/programs/quasars/CIV_CLQs/data/J1205+3422/'

##  Light Curve
infile = 'NEOWISER-L1b_J1205p3422.dat'
NEOWISER = ascii.read(path+infile)
NEOWISER_W1_AB = NEOWISER['w1mpro'] + 2.673
NEOWISER_W2_AB = NEOWISER['w2mpro'] + 3.313
infile = 'NEOWISER-L1b_J1205p3422_averaged.dat'
NEOWISER_aver = ascii.read(path+infile)
NEOWISER_aver_W1_AB = NEOWISER_aver['w1mpro_wgt'] + 2.673
NEOWISER_aver_W2_AB = NEOWISER_aver['w2mpro_wgt'] + 3.313

##  Spectra
infile   = 'DBSP_J1205p3422_b_58538.dat'
DBSP_b1  = ascii.read(path+infile)
infile   = 'DBSP_J1205p3422_r_58538.dat'
DBSP_r1  = ascii.read(path+infile)
infile   = 'DBSP_J1205p3422_b_58693.dat'
DBSP_b2  = ascii.read(path+infile)
infile   = 'DBSP_J1205p3422_r_58693.dat'
DBSP_r2  = ascii.read(path+infile)

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
DBSP_b1_flux   = (2.99792458E+21 * (DBSP_b1['flux_density'] / 1000.)) / (DBSP_b1['wavelength']**2)  / 1e-17
DBSP_r1_flux   = (2.99792458E+21 * (DBSP_r1['flux_density'] / 1000.)) / (DBSP_r1['wavelength']**2)  / 1e-17
DBSP_b2_flux   = (2.99792458E+21 * (DBSP_b2['flux_density'] / 1000.)) / (DBSP_b2['wavelength']**2)  / 1e-17
DBSP_r2_flux   = (2.99792458E+21 * (DBSP_r2['flux_density'] / 1000.)) / (DBSP_r2['wavelength']**2)  / 1e-17


##    https://dr15.sdss.org/optical/spectrum/view?plateid=2089&mjd=53498&fiberid=427&run2d=26        

sdssname        = 'spec-2089-53498-0427.fits'
sdss_data       = fits.open(path+sdssname)
sdss_spectrum   = sdss_data[1].data
sdss_flux       = sdss_spectrum.flux
sdss_loglam     = sdss_spectrum.loglam
sdss_wavelength = 10**(sdss_loglam)
sdss_mjd        = sdss_data[2].data['MJD'][0]


    
redshift = 2.068


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
ymin = -9.95; ymax = 80.   
ax2.set_xlim([xmin,xmax])
ax2.set_ylim([ymin, ymax])

## LyA 
xmin = 1080; xmax = 1380
ymin = -9.5; ymax = 83.
ax3.set_xlim([xmin,xmax])
ax3.set_ylim([ymin, ymax])

## CIV Line
xmin = 1500; xmax = 1580
ymin = -9.5; ymax = 49.   
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
ymin = -4.95; ymax = 20.   
ax5.set_xlim([xmin,xmax])
ax5.set_ylim([ymin, ymax])


## The Spectra curves on the RHS
ax2.plot(   sdss_wavelength/(1+redshift), sdss_flux,          '-b', lw=linewidth/4.0, label='SDSS MJD 53498')
ax2.plot(DBSP_b1_wavelength/(1+redshift), DBSP_b1_flux,        '-r', lw=linewidth/4.0,     label='DBSP MJD 58538')
ax2.plot(DBSP_r1_wavelength/(1+redshift), DBSP_r1_flux,        '-r', lw=linewidth/4.0)
## Smoothing:: https://joseph-long.com/writing/AstroPy-boxcar/
DBSP_b2_flux_smoothed = convolve(DBSP_b2_flux, Box1DKernel(5))
DBSP_r2_flux_smoothed = convolve(DBSP_r2_flux, Box1DKernel(5))
ax2.plot(DBSP_b2_wavelength/(1+redshift), DBSP_b2_flux_smoothed, '-k', lw=linewidth/4.0)
ax2.plot(DBSP_r2_wavelength/(1+redshift), DBSP_r2_flux_smoothed, '-k', lw=linewidth/4.0,     label='DBSP MJD 58693')

## Just plotting the spectra again, the axes limits take care of the "zoom in"
ax3.plot(   sdss_wavelength/(1+redshift), sdss_flux,    '-b', lw=linewidth/4.0)
ax3.plot(DBSP_b1_wavelength/(1+redshift), DBSP_b1_flux, '-r', lw=linewidth/4.0)
ax3.plot(DBSP_b2_wavelength/(1+redshift), DBSP_b2_flux_smoothed, '-k', lw=linewidth/4.0)

ax4.plot(   sdss_wavelength/(1+redshift), sdss_flux,    '-b', lw=linewidth/4.0)
ax4.plot(DBSP_b1_wavelength/(1+redshift), DBSP_b1_flux, '-r', lw=linewidth/4.0)
ax4.plot(DBSP_b2_wavelength/(1+redshift), DBSP_b2_flux_smoothed, '-k', lw=linewidth/4.0)
#ax4.plot(DBSP_r1_wavelength/(1+redshift), DBSP_r1_flux, '-r', lw=linewidth/4.0)
#ax4.plot(DBSP_r2_wavelength/(1+redshift), DBSP_r2_flux_smoothed, '-k', lw=linewidth/4.0)

ax5.plot(   sdss_wavelength/(1+redshift),    sdss_flux,    '-b', lw=linewidth/4.0)
ax5.plot(DBSP_r1_wavelength/(1+redshift), DBSP_r1_flux, '-r', lw=linewidth/4.0)
ax5.plot(DBSP_r2_wavelength/(1+redshift), DBSP_r2_flux_smoothed, '-k', lw=linewidth/4.0)


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

neo_w1_under = mlines.Line2D([], [], color='k',    marker='o', linestyle='None', markersize=ms_big, label='NEOWISE W1')
neo_w2_under = mlines.Line2D([], [], color='k' ,   marker='o', linestyle='None', markersize=ms_big,     label='NEOWISE W2')
neo_w1_over = mlines.Line2D([], [], color='red',  marker='o', linestyle='None', markersize=ms_big, label='')
neo_w2_over = mlines.Line2D([], [], color='cyan', marker='o', linestyle='None', markersize=ms_big,     label='')


neo_w1 = mlines.Line2D([], [], label='NEOWISE W1', color='red',  marker="o", markeredgecolor='k', markeredgewidth=2.0, markersize=12, linestyle='None')
neo_w2 = mlines.Line2D([], [], label='NEOWISE W2', color='cyan', marker="o", markeredgecolor='k', markeredgewidth=2.0, markersize=12, linestyle='None')

#leg = ax1.legend(loc='upper left', fontsize=fontsize/1.2, markerscale=(1/7.6),  handles=[neo_w1_under, neo_w2_under])
#leg = ax1.legend(loc='upper left', fontsize=fontsize/1.2, markerscale=(1/7.8),  handles=[neo_w1_over, neo_w2_over])
leg = ax1.legend(loc='upper left', fontsize=fontsize/1.2, handles=[neo_w1, neo_w2])
#leg = ax1.legend(loc='upper left', fontsize=fontsize/1.2, markerscale=2.8)


#for legend_handle in legend.legendHandles:
#    legend_handle._legmarker.set_markersize(9)
#for line in leg.get_lines(): line.set_markersize(120)
#leg.legendHandles[0]._sizes = [30]


leg = ax2.legend(loc='upper right', fontsize=fontsize/1.2)
for line in leg.get_lines(): line.set_linewidth(2.2)


plt.savefig('J1205+3422_temp.png', format='png')
plt.close(fig)
