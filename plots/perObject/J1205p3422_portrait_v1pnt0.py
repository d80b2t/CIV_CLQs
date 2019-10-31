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

##    L I G H T    C U R V E
## CRTS
infile = 'J120544.67+342252.4_CRTS.dat'
CRTS   = ascii.read(path+infile)
## CRTS MJD offset
CRTS_MJD_offset = 53500
CRTS_MJD = CRTS['mjd'] + CRTS_MJD_offset

infile = 'J120544.67+342252_ztf_g.dat'
ZTF_g  = ascii.read(path+infile)

infile = 'J120544.67+342252_ztf_r.dat'
ZTF_r  = ascii.read(path+infile)


infile = 'NEOWISER-L1b_J1205p3422.dat'
NEOWISER = ascii.read(path+infile)
NEOWISER_W1_AB = NEOWISER['w1mpro'] + 2.673
NEOWISER_W2_AB = NEOWISER['w2mpro'] + 3.313
infile = 'NEOWISER-L1b_J1205p3422_averaged.dat'
NEOWISER_aver = ascii.read(path+infile)
NEOWISER_aver_W1_AB = NEOWISER_aver['w1mpro_wgt'] + 2.673
NEOWISER_aver_W2_AB = NEOWISER_aver['w2mpro_wgt'] + 3.313

##    S P E C T R A 
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
ax2 = fig.add_subplot(gs[0,  2:5])
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


##  T H E    L I G H T     C U R V E S
xmin = 53400; xmax = 59000
ymin = 20.85; ymax = 16.10   
ax1.set_xlim([xmin,xmax])
ax1.set_ylim([ymin, ymax])

## Plotting the SPECTRA as vertical lines
lw = 2.0
ax1.axvline(x=53498, linewidth=lw, linestyle='dotted', color='b')
ax1.axvline(x=58538, linewidth=lw, linestyle='dashed', color='r')
ax1.axvline(x=58693, linewidth=lw, linestyle='dashed', color='k')

## CRTS data
lw = 1
ax1.scatter( CRTS_MJD, CRTS['magnitude'], color='k',       alpha=alpha, s=ms*1.8)
ax1.scatter( CRTS_MJD, CRTS['magnitude'], color='dimgray', alpha=alpha, s=ms, label='CRTS')
ax1.errorbar(CRTS_MJD, CRTS['magnitude'], color='k', yerr=CRTS['magnitude_err'], fmt='o', linewidth=lw, ms=ms)
## ZTF data
ax1.scatter( ZTF_g['mjd'], ZTF_g['mag'], color='k',         alpha=alpha, s=ms*1.8)
ax1.scatter( ZTF_g['mjd'], ZTF_g['mag'], color='olivedrab', alpha=alpha, s=ms, label='ZTF g-band')
ax1.errorbar(ZTF_g['mjd'], ZTF_g['mag'], color='olivedrab', yerr=ZTF_g['magerr'], fmt='o', linewidth=lw, ms=ms)
ax1.scatter( ZTF_r['mjd'], ZTF_r['mag'], color='k',         alpha=alpha, s=ms*1.8)
ax1.scatter( ZTF_r['mjd'], ZTF_r['mag'], color='tomato',    alpha=alpha, s=ms, label='ZTF r-band')
ax1.errorbar(ZTF_r['mjd'], ZTF_r['mag'], color='tomato',    yerr=ZTF_r['magerr'], fmt='o', linewidth=lw, ms=ms)

#olivedrab


## NEOWISER W1/2 (AB)
## indigo and brown were used in Ross et al. (2018)
ms              = 10.
ms_big          = ms*6.
ax1.scatter(NEOWISER['mjd'],  NEOWISER_W1_AB, color='k', alpha=alpha, s=ms*1.8)
ax1.scatter(NEOWISER['mjd'],  NEOWISER_W1_AB, color='r', alpha=alpha, s=ms, label='NEOWISE W1')
ax1.scatter(NEOWISER['mjd'],  NEOWISER_W2_AB, color='k', alpha=alpha, s=ms*1.8)
ax1.scatter(NEOWISER['mjd'],  NEOWISER_W2_AB, color='c', alpha=alpha, s=ms, label='NEOWISE W2')

ax1.scatter(NEOWISER_aver['mean_mjd'],  NEOWISER_aver_W1_AB, color='k', alpha=alpha, s=ms_big*1.8)
ax1.scatter(NEOWISER_aver['mean_mjd'],  NEOWISER_aver_W1_AB, color='r', alpha=alpha, s=ms_big)
ax1.scatter(NEOWISER_aver['mean_mjd'],  NEOWISER_aver_W2_AB, color='k', alpha=alpha, s=ms_big*1.8)
ax1.scatter(NEOWISER_aver['mean_mjd'],  NEOWISER_aver_W2_AB, color='c', alpha=alpha, s=ms_big)


## full, main spectral plot
xmin = 910;  xmax = 3590
ymin = -9.95; ymax = 85.   
ax2.set_xlim([xmin,xmax])
ax2.set_ylim([ymin, ymax])

## LyA 
xmin = 1080; xmax = 1380
ymin = -9.5; ymax = 88.
ax3.set_xlim([xmin,xmax])
ax3.set_ylim([ymin, ymax])

## CIV Line
xmin = 1495; xmax = 1585
ymin = -9.5; ymax =54.   
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
xmin = 2659; xmax = 2949
ymin = -4.95; ymax = 27.   
ax5.set_xlim([xmin,xmax])
ax5.set_ylim([ymin, ymax])


## The Spectra curves on the RHS
ax2.plot(   sdss_wavelength/(1+redshift), sdss_flux,           '-b', lw=linewidth/4.0, label='SDSS MJD 53498')
ax2.plot(DBSP_b1_wavelength/(1+redshift), DBSP_b1_flux,        '-r', lw=linewidth/4.0,     label='DBSP MJD 58538')
ax2.plot(DBSP_r1_wavelength/(1+redshift), DBSP_r1_flux,        '-r', lw=linewidth/4.0)
## Smoothing:: https://joseph-long.com/writing/AstroPy-boxcar/
DBSP_b2_flux_smoothed = convolve(DBSP_b2_flux, Box1DKernel(5))
DBSP_r2_flux_smoothed = convolve(DBSP_r2_flux, Box1DKernel(5))
ax2.plot(DBSP_b2_wavelength/(1+redshift), DBSP_b2_flux_smoothed, '-k', lw=linewidth/4.0)
ax2.plot(DBSP_r2_wavelength/(1+redshift), DBSP_r2_flux_smoothed, '-k', lw=linewidth/4.0,     label='DBSP MJD 58693')
for ll in range(len(linelist)):
    ax2.axvline(x=linelist['Wavelength'][ll], color='gray', linestyle='--', linewidth=linewidth/3.4)
    label = linelist['LineName'][ll]
    xylabel = ((linelist['Wavelength'][ll]), 83)
    ax2.annotate(label, xy=xylabel, ha='center', va='center', rotation=90, fontsize=fontsize/1.6 )


## Just plotting the spectra again, the axes limits take care of the "zoom in"
ax3.plot(   sdss_wavelength/(1+redshift), sdss_flux,    '-b', lw=linewidth/4.0)
ax3.plot(DBSP_b1_wavelength/(1+redshift), DBSP_b1_flux, '-r', lw=linewidth/4.0)
ax3.plot(DBSP_b2_wavelength/(1+redshift), DBSP_b2_flux_smoothed, '-k', lw=linewidth/4.0)
for ll in range(len(linelist)):
    ax3.axvline(x=linelist['Wavelength'][ll], color='gray', linestyle='--', linewidth=linewidth/3.4)
    label = linelist['LineName'][ll]
    xylabel = ((linelist['Wavelength'][ll]), 78)
    ax3.annotate(label, xy=xylabel, ha='center', va='center', rotation=90, fontsize=fontsize/1.6 )


ax4.plot(   sdss_wavelength/(1+redshift), sdss_flux,    '-b', lw=linewidth/4.0)
ax4.plot(DBSP_b1_wavelength/(1+redshift), DBSP_b1_flux, '-r', lw=linewidth/4.0)
ax4.plot(DBSP_b2_wavelength/(1+redshift), DBSP_b2_flux_smoothed, '-k', lw=linewidth/4.0)
#ax4.plot(DBSP_r1_wavelength/(1+redshift), DBSP_r1_flux, '-r', lw=linewidth/4.0)
#ax4.plot(DBSP_r2_wavelength/(1+redshift), DBSP_r2_flux_smoothed, '-k', lw=linewidth/4.0)
for ll in range(len(linelist)):
    ax4.axvline(x=linelist['Wavelength'][ll], color='gray', linestyle='--', linewidth=linewidth/3.4)
    label = linelist['LineName'][ll]
    xylabel = ((linelist['Wavelength'][ll]), 45)
    ax4.annotate(label, xy=xylabel, ha='center', va='center', rotation=90, fontsize=fontsize/1.6 )

ax5.plot(   sdss_wavelength/(1+redshift),    sdss_flux,    '-b', lw=linewidth/4.0)
ax5.plot(DBSP_r1_wavelength/(1+redshift), DBSP_r1_flux, '-r', lw=linewidth/4.0)
ax5.plot(DBSP_r2_wavelength/(1+redshift), DBSP_r2_flux_smoothed, '-k', lw=linewidth/4.0)
for ll in range(len(linelist)):
    ax5.axvline(x=linelist['Wavelength'][ll], color='gray', linestyle='--', linewidth=linewidth/3.4)
    label = linelist['LineName'][ll]
    xylabel = ((linelist['Wavelength'][ll]), 24)
    ax5.annotate(label, xy=xylabel, ha='center', va='center', rotation=90, fontsize=fontsize/1.6 )




## Making all the nice plots boxes the right dimensions and placings
#% start: automatic generated code from pylustrator
plt.figure(1).ax_dict = {ax.get_label(): ax for ax in plt.figure(1).axes}
import matplotlib as mpl
plt.figure(1).set_size_inches(25.470000/2.54, 20.680000/2.54, forward=True)
plt.figure(1).axes[0].get_xaxis().get_label().set_fontsize(fontsize)
plt.figure(1).axes[0].get_yaxis().get_label().set_fontsize(fontsize)
plt.figure(1).axes[0].set_position([0.060000, 0.688887, 0.898036, 0.268065])
plt.figure(1).axes[0].xaxis.labelpad = 10.0
plt.figure(1).axes[0].xaxis.labelpad = -4.080000
#plt.figure(1).axes[0].get_legend()._set_loc((0.055442, 0.054485))
plt.figure(1).axes[1].get_xaxis().get_label().set_fontsize(fontsize/1.2)
plt.figure(1).axes[1].get_yaxis().get_label().set_rotation(90.0)
plt.figure(1).axes[1].set_position([0.060000, 0.391842, 0.898036, 0.248409])
plt.figure(1).axes[1].set_xlim(910.0, 3490.0)
plt.figure(1).axes[1].set_ylim(-9.95, 95.0)
plt.figure(1).axes[1].xaxis.labelpad = 323.200000
plt.figure(1).axes[1].yaxis.labelpad = 4.720000
plt.figure(1).axes[1].yaxis.labelpad = -4.730058
plt.figure(1).axes[1].yaxis.set_label_coords(-0.045, -0.1)
#plt.figure(1).axes[1].get_legend()._set_loc((0.803069, -0.638658))
plt.figure(1).axes[2].get_xaxis().get_label().set_text("")
plt.figure(1).axes[2].set_position([0.060000, 0.117801, 0.205718, 0.225406])
plt.figure(1).axes[2].set_ylim(-9.5, 92.0)
plt.figure(1).axes[3].set_position([0.305213, 0.117820, 0.205887, 0.225386])
plt.figure(1).axes[3].set_ylim(-5.5, 59.0)
plt.figure(1).axes[3].xaxis.labelpad = 10.0
plt.figure(1).axes[4].set_position([0.550626, 0.117801, 0.205887, 0.225386])
#% end: automatic generated code from pylustrator

## Axes labels
plt.rcParams['text.usetex'] = True

ax1.set_xlabel(r'MJD',          fontsize=fontsize/1.2)
ax1.set_ylabel(r'AB Magnitude', fontsize=fontsize/1.2)
ax2.set_xlabel(r'Restframe wavelength', fontsize=fontsize)
ax2.set_ylabel(r'F$_{\lambda}$ / 10$^{-17}$ erg cm$^{-2}$ s$^{-1}$ ${ \rm \AA }^{-1}$', fontsize=fontsize/1.4)
ax4.set_xlabel(r'Restframe wavelength', fontsize=fontsize)

## Axes style
ax2.tick_params(labelsize=labelsize/1.15, direction='in', length=ticklength, width=tickwidth)
ax3.tick_params(labelsize=labelsize/1.6,  direction='in', length=ticklength, width=tickwidth)
ax4.tick_params(labelsize=labelsize/1.6,  direction='in', length=ticklength, width=tickwidth)
ax5.tick_params(labelsize=labelsize/1.6,  direction='in', length=ticklength, width=tickwidth)



##
##  L E G E N D S
##
neo_w1 = mlines.Line2D([], [], label='NEOWISE W1', color='red',
                       marker="o", markeredgecolor='k', markeredgewidth=2.0, markersize=12, linestyle='None')
neo_w2 = mlines.Line2D([], [], label='NEOWISE W2', color='cyan',
                       marker="o", markeredgecolor='k', markeredgewidth=2.0, markersize=12, linestyle='None')
crts   = mlines.Line2D([], [], label='CRTS',       color='k',
                       marker="o", markeredgecolor='k', markeredgewidth=2.0, markersize=5,  linestyle='None')
ztf_g  = mlines.Line2D([], [], label='ZTF g-band', color='olivedrab',
                       marker="o", markeredgecolor='k', markeredgewidth=1.4, markersize=7,  linestyle='None')
ztg_r  = mlines.Line2D([], [], label='ZTF r-band', color='tomato',
                       marker="o", markeredgecolor='k', markeredgewidth=1.4, markersize=7,  linestyle='None')

handles=[neo_w1, neo_w2, crts, ztf_g, ztg_r]
leg = ax1.legend(loc='lower left', fontsize=fontsize/1.4, handles=handles)


leg = ax2.legend(loc='upper right', fontsize=fontsize/1.2)
for line in leg.get_lines(): line.set_linewidth(2.2)

#plt.show()
plt.savefig('J1205+3422_portrait_temp.png', format='png')
plt.close(fig)
