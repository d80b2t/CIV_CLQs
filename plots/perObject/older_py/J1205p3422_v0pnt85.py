
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from matplotlib          import colors as mcolors
from matplotlib.gridspec import GridSpec
from matplotlib.patches  import Rectangle

from astropy.io import ascii
from astropy.io import fits

import pylustrator
pylustrator.start()


## Reading in the data
path = '/cos_pc19a_npr/programs/quasars/CIV_CLQs/data/J1205+3422/'
infile = 'crts1205p3422_b.dat'
crts_b = ascii.read(path+infile)
infile = 'crts1205p3422_r.dat'
crts_r = ascii.read(path+infile)

crts_b_wavelength = crts_b['wavelength']
crts_r_wavelength = crts_r['wavelength']

## Need to convert BUNIT   = 'erg/cm2/s/Hz'  to  'erg/cm2/s/A'
## http://www.stsci.edu/~strolger/docs/UNITS.txt
## 
## [Y erg/cm^2/s/Hz]            = 1000 * [X W/m^2/Hz]
## [Y erg/cm^2/s/A]             = 2.99792458E+21 * [X1 W/m^2/Hz] / [X2 A]^2
## / 1e-17 to put into 'regular' SDSS flux units 
crts_b_flux   = (2.99792458E+21* (crts_b['flux_density'] / 1000.))/(crts_b['wavelength']**2)  / 1e-17
crts_r_flux   = (2.99792458E+21* (crts_r['flux_density'] / 1000.))/(crts_r['wavelength']**2)  / 1e-17

    
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
left   = 0.08   # the left side of the subplots of the figure
right  = 0.96   # the right side of the subplots of the figure
bottom = 0.16   # the bottom of the subplots of the figure
top    = 0.96   # the top of the subplots of the figure
wspace = 0.26   # the amount of width reserved for blank space between subplots
hspace = 0.26   # the amount of height reserved for white space between subplots

plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)

 
## Some NPR defaults
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


## The Light curves on the LHS

## Axes limits
#xmin =  4.90
#xmax =  7.65
#ymin =  0.0
#ymax = 38.0

## Light Curve
xmin = 51000; xmax = 59000
ymin = 19.85; ymax = 15.10   
ax1.set_xlim([xmin,xmax])
ax1.set_ylim([ymin, ymax])

xmin = 1100;  xmax = 3590
ymin = -9.95; ymax = 80.   
ax2.set_xlim([xmin,xmax])
ax2.set_ylim([ymin, ymax])

## CIV Line
xmin = 1500; xmax = 1580
ymin = -9.95; ymax = 20.   
ax3.set_xlim([xmin,xmax])
ax3.set_ylim([ymin, ymax])

## CIII]
xmin = 1890; xmax = 1920
ax4.set_xlim([xmin,xmax])
ax4.set_ylim([ymin, ymax])

## MgII
xmin = 2700; xmax = 2900
ax5.set_xlim([xmin,xmax])
ax5.set_ylim([ymin, ymax])


## The Spectra curves on the RHS
ax2.plot(crts_b_wavelength/(1+redshift), crts_b_flux,  '-k', lw=linewidth/4.0)
ax2.plot(crts_r_wavelength/(1+redshift), crts_r_flux,  '-k', lw=linewidth/4.0)

ax3.plot(crts_b_wavelength/(1+redshift), crts_b_flux,  '-k', lw=linewidth/4.0)
ax3.plot(crts_r_wavelength/(1+redshift), crts_r_flux,  '-k', lw=linewidth/4.0)
ax4.plot(crts_b_wavelength/(1+redshift), crts_b_flux,  '-k', lw=linewidth/4.0)
ax4.plot(crts_r_wavelength/(1+redshift), crts_r_flux,  '-k', lw=linewidth/4.0)
ax5.plot(crts_b_wavelength/(1+redshift), crts_b_flux,  '-k', lw=linewidth/4.0)
ax5.plot(crts_r_wavelength/(1+redshift), crts_r_flux,  '-k', lw=linewidth/4.0)


## Axes labels
ax1.set_xlabel(r'MJD', fontsize=fontsize)
ax1.set_ylabel(r'AB Magnitude', fontsize=fontsize)

ax2.set_xlabel(r'Restframe wavelength', fontsize=fontsize)
ax2.set_ylabel(r'F$_{\lambda}$ (10$^{-17}$ erg cm$^{-2}$ s$^{-1}$ Ang$^{-1}$', fontsize=fontsize/2.)

ax4.set_xlabel(r'Restframe wavelength', fontsize=fontsize)
#ax3.set_ylabel(r'F_{\lambda} (10^{-17} erg cm^{-2} s^{-1} Ang^{-1}', fontsize=fontsize)


## Axes style
ax3.tick_params(labelsize=labelsize/2., direction='in', length=ticklength,   width=tickwidth)
ax4.tick_params(axis='both', which='major', labelsize=labelsize/2., direction='in', length=ticklength,   width=tickwidth)
ax5.tick_params(axis='both', which='major', labelsize=labelsize/2., direction='in', length=ticklength,   width=tickwidth)



#% start: automatic generated code from pylustrator
plt.figure(1).ax_dict = {ax.get_label(): ax for ax in plt.figure(1).axes}
import matplotlib as mpl
plt.figure(1).axes[0].get_xaxis().get_label().set_fontsize(28)
plt.figure(1).axes[0].get_yaxis().get_label().set_fontsize(28)
plt.figure(1).axes[0].set_position([0.040000, 0.136250, 0.357689, 0.823750])
plt.figure(1).axes[0].xaxis.labelpad = 16.720000
plt.figure(1).axes[1].get_xaxis().get_label().set_fontsize(24.000000)
plt.figure(1).axes[1].get_yaxis().get_label().set_rotation(90.0)
plt.figure(1).axes[1].set_position([0.460541, 0.560924, 0.499459, 0.399076])
plt.figure(1).axes[1].xaxis.labelpad = 262.160000
plt.figure(1).axes[1].xaxis.labelpad = 323.200000
plt.figure(1).axes[1].yaxis.labelpad = 4.720000
plt.figure(1).axes[2].get_xaxis().get_label().set_text("")
plt.figure(1).axes[2].set_position([0.460541, 0.135000, 0.146836, 0.372826])
plt.figure(1).axes[3].set_position([0.639228, 0.135000, 0.146836, 0.372826])
plt.figure(1).axes[3].xaxis.labelpad = -3.280000
plt.figure(1).axes[4].set_position([0.813164, 0.135000, 0.146836, 0.372826])
#% end: automatic generated code from pylustrator
plt.show()

plt.savefig('J1205+3422_temp.png', format='png')
plt.close(fig)
