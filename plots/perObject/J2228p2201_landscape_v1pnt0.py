'''
      J 2 2 2 8  +  2 2 0 1
   
22h28m18.7s +22d01m02.9s
337.078194  +22.017478 
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
path = '/cos_pc19a_npr/programs/quasars/CIV_CLQs/data/J2228+2201/'

##    L I G H T    C U R V E
## CRTS
infile = 'J222818.77+220102.6_CRTS.dat'
CRTS   = ascii.read(path+infile)
## CRTS MJD offset
CRTS_MJD_offset = 53500
CRTS_MJD = CRTS['mjd'] + CRTS_MJD_offset
## ZTF
infile = 'J222818+220102_ztf_g.dat'
ZTF_g  = ascii.read(path+infile)
infile = 'J222818+220102_ztf_r.dat'
ZTF_r  = ascii.read(path+infile)

##  Pan-STARRS 
infile    = 'PanSTARRS_DR2_detections.dat'
PanSTARRS = ascii.read(path+infile)
##  g=1, r=2, i=3, z=4, y=5 
PS_gband   = PanSTARRS[np.where(PanSTARRS['filterID'] == 1)]
PS_rband   = PanSTARRS[np.where(PanSTARRS['filterID'] == 2)]
PS_gPSFmag = -2.5*(np.log10(PS_gband['psfFlux']/3631.))
PS_rPSFmag = -2.5*(np.log10(PS_rband['psfFlux']/3631.))

##  ALLWISE
infile  = 'ALLWISE_catalog.dat'
ALLWISE = ascii.read(path+infile)
ALLWISE_W1_AB = ALLWISE['w1mpro'] + 2.673
ALLWISE_W2_AB = ALLWISE['w2mpro'] + 3.313

## NEOWISE-R
infile = 'NEOWISER-L1b_J2228+2201.dat'
NEOWISER = ascii.read(path+infile)
NEOWISER_W1_AB = NEOWISER['w1mpro'] + 2.673
NEOWISER_W2_AB = NEOWISER['w2mpro'] + 3.313
infile = 'NEOWISER-L1b_J2228+2201_averaged.dat'
NEOWISER_aver = ascii.read(path+infile)
NEOWISER_aver_W1_AB = NEOWISER_aver['w1mpro_wgt'] + 2.673
NEOWISER_aver_W2_AB = NEOWISER_aver['w2mpro_wgt'] + 3.313


##    S P E C T R A 
##  It's actaully a BOSS spectrum, but let's just call it SDSS for
##  ease of variable...
sdssname        = 'spec-6118-56189-0720.fits'
sdss_data       = fits.open(path+sdssname)
sdss_spectrum   = sdss_data[1].data
sdss_flux       = sdss_spectrum.flux
sdss_loglam     = sdss_spectrum.loglam
sdss_wavelength = 10**(sdss_loglam)
sdss_mjd        = sdss_data[2].data['MJD'][0]

bossname        = 'spec-7582-56960-0790.fits'
boss_data       = fits.open(path+bossname)
boss_spectrum   = boss_data[1].data
boss_flux       = boss_spectrum.flux
boss_loglam     = boss_spectrum.loglam
boss_wavelength = 10**(boss_loglam)
boss_mjd        = boss_data[2].data['MJD'][0]

##  DBSP spectrum
infile   = 'DBSP_J2228p2201_b_58693.dat'
DBSP_b  = ascii.read(path+infile)
infile   = 'DBSP_J2228p2201_r_58693.dat'
DBSP_r  = ascii.read(path+infile)

## The Blue and Red arms overlap in wavelength;
## Blue has better SNR than RED in overlapping regions
DBSP_r = DBSP_r[np.where(DBSP_r['wavelength'] > DBSP_b['wavelength'].max())]

DBSP_b_wavelength = DBSP_b['wavelength']
DBSP_r_wavelength = DBSP_r['wavelength']

##  Need to convert BUNIT   = 'erg/cm2/s/Hz'  to  'erg/cm2/s/A'
##  http://www.stsci.edu/~strolger/docs/UNITS.txt
## 
##  [Y erg/cm^2/s/Hz]            = 1000 * [X W/m^2/Hz]
##  [Y erg/cm^2/s/A]             = 2.99792458E+21 * [X1 W/m^2/Hz] / [X2 A]^2
##  / 1e-17 to put into 'regular' SDSS flux units 
DBSP_b_flux   = (2.99792458E+21 * (DBSP_b['flux_density'] / 1000.)) / (DBSP_b['wavelength']**2)  / 1e-17
DBSP_r_flux   = (2.99792458E+21 * (DBSP_r['flux_density'] / 1000.)) / (DBSP_r['wavelength']**2)  / 1e-17


## Get the redshift right!
redshift = 2.222
print()
print('Redshift of J2228+2201 is ', redshift)
print()

## Emission Line list
linelist_file = 'emission_lines.dat'
linelist = ascii.read(path+linelist_file)



## https://jakevdp.github.io/PythonDataScienceHandbook/04.08-multiple-subplots.html
fig = plt.figure(figsize=(12, 4.3))
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
lw              = 1.0
ls              = 'solid'
ms              = 1.
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
plt.rcParams['text.usetex'] = True


##  T H E    L I G H T     C U R V E S 
xmin = 53400; xmax = 59000
ymin = 21.20; ymax = 16.20   
ax1.set_xlim([xmin,xmax])
ax1.set_ylim([ymin, ymax])

## Plotting the SPECTRA as vertical lines
lw = 2.0
ax1.axvline(x=56189, linewidth=lw, linestyle='dotted', color='b')
ax1.axvline(x=56960, linewidth=lw, linestyle='dashed', color='r')
ax1.axvline(x=58693, linewidth=lw, linestyle='dashed', color='k')

## CRTS data
lw = 1
ax1.scatter( CRTS_MJD, CRTS['magnitude'], color='k',       alpha=alpha, s=ms*1.8, zorder=0)
ax1.scatter( CRTS_MJD, CRTS['magnitude'], color='dimgray', alpha=alpha, s=ms,     zorder=0, label='CRTS')
ax1.errorbar(CRTS_MJD, CRTS['magnitude'], color='k', yerr=CRTS['magnitude_err'],  zorder=0, fmt='o', linewidth=lw, ms=ms)
## ZTF data
ax1.scatter( ZTF_g['mjd'], ZTF_g['mag'], color='k',         alpha=alpha, s=ms*1.8)
ax1.scatter( ZTF_g['mjd'], ZTF_g['mag'], color='olivedrab', alpha=alpha, s=ms, label='ZTF g-band')
ax1.errorbar(ZTF_g['mjd'], ZTF_g['mag'], color='olivedrab', yerr=ZTF_g['magerr'], fmt='o', linewidth=lw, ms=ms)
ax1.scatter( ZTF_r['mjd'], ZTF_r['mag'], color='k',         alpha=alpha, s=ms*1.8)
ax1.scatter( ZTF_r['mjd'], ZTF_r['mag'], color='tomato',    alpha=alpha, s=ms, label='ZTF r-band')
ax1.errorbar(ZTF_r['mjd'], ZTF_r['mag'], color='tomato',    yerr=ZTF_r['magerr'], fmt='o', linewidth=lw, ms=ms)

## PanSTARRS data
ms     =  8.
ms_big = 36.
ax1.scatter(PS_gband['obsTime'], PS_gPSFmag, color='k',        alpha=alpha, s=ms_big)
ax1.scatter(PS_gband['obsTime'], PS_gPSFmag, color='lime',     alpha=alpha, s=ms, label='Pan-STARRS g-band')
ax1.scatter(PS_rband['obsTime'], PS_rPSFmag, color='k',        alpha=alpha, s=ms_big)
ax1.scatter(PS_rband['obsTime'], PS_rPSFmag, color='deeppink', alpha=alpha, s=ms, label='Pan-STARRS r-band')

## NEOWISER W1/2 (AB)
ms     = 10.
ms_big = ms*6.
ax1.scatter(ALLWISE['w1mjdmean'], ALLWISE_W1_AB, color='k', alpha=alpha, s=ms_big*1.8)
ax1.scatter(ALLWISE['w1mjdmean'], ALLWISE_W1_AB, color='r', alpha=alpha, s=ms_big, label='ALLWISE W1')
ax1.scatter(ALLWISE['w2mjdmean'], ALLWISE_W2_AB, color='k', alpha=alpha, s=ms_big*1.8)
ax1.scatter(ALLWISE['w2mjdmean'], ALLWISE_W2_AB, color='c', alpha=alpha, s=ms_big, label='ALLWISE W2')

ax1.scatter(NEOWISER['mjd'],           NEOWISER_W1_AB, color='k', alpha=alpha, s=ms*1.8)
ax1.scatter(NEOWISER['mjd'],           NEOWISER_W1_AB, color='r', alpha=alpha, s=ms, label='NEOWISE W1')
ax1.scatter(NEOWISER['mjd'],           NEOWISER_W2_AB, color='k', alpha=alpha, s=ms*1.8)
ax1.scatter(NEOWISER['mjd'],           NEOWISER_W2_AB, color='c', alpha=alpha, s=ms, label='NEOWISE W2')
ax1.scatter(NEOWISER_aver['mean_mjd'], NEOWISER_aver_W1_AB, color='k', alpha=alpha, s=ms_big*1.8)
ax1.scatter(NEOWISER_aver['mean_mjd'], NEOWISER_aver_W1_AB, color='r', alpha=alpha, s=ms_big)
ax1.scatter(NEOWISER_aver['mean_mjd'], NEOWISER_aver_W2_AB, color='k', alpha=alpha, s=ms_big*1.8)
ax1.scatter(NEOWISER_aver['mean_mjd'], NEOWISER_aver_W2_AB, color='c', alpha=alpha, s=ms_big)


##  P L O T T I N G    T H E      S P E C T R A  
## full, main spectral plot
xmin_ax2 =  910; xmax_ax2 = 3350
ymin_ax2 = -4.9; ymax_ax2 = 32.
#ymin_ax2 = -1.9; ymax_ax2 = 3.2   ## 32 is 
ax2.set_xlim([xmin_ax2, xmax_ax2])
ax2.set_ylim([ymin_ax2, ymax_ax2])

## LyA 
xmin_ax3 = 1150; xmax_ax3 = 1285    ## to match VdB01, Fig 7; had been 1080, 1380 previously
ymin_ax3 = -4.9; ymax_ax3 = ymax_ax2
ax3.set_xlim([xmin_ax3, xmax_ax3])
ax3.set_ylim([ymin_ax3, ymax_ax3])

## CIV 
xmin_ax4 = 1468; xmax_ax4 = 1620    ## to match VdB01, Fig 7; had been 1495, 1585 previously
ymin_ax4 =  -1.9; ymax_ax4 = 12.9   
ax4.set_xlim([xmin_ax4, xmax_ax4])
ax4.set_ylim([ymin_ax4, ymax_ax4])

## MgII
xmin_ax5 =  2680; xmax_ax5 = 2910   ## to match VdB01, Fig 7; had been 2659, 2949 previously
ymin_ax5 =  -1.9; ymax_ax5 = 8.   
ax5.set_xlim([xmin_ax5, xmax_ax5])
ax5.set_ylim([ymin_ax5, ymax_ax5])


## The Spectra curves on the RHS
ax2.plot(   sdss_wavelength/(1+redshift), sdss_flux,         '-b', lw=linewidth/4.0, label=r'\ BOSS MJD 56189')
ax2.plot(   boss_wavelength/(1+redshift), boss_flux,         '-r', lw=linewidth/4.0, label='eBOSS MJD 56960')
ax2.plot(DBSP_b_wavelength/(1+redshift), DBSP_b_flux,        '-k', lw=linewidth/4.0, label='DBSP MJD 58693')
ax2.plot(DBSP_r_wavelength/(1+redshift), DBSP_r_flux,        '-k', lw=linewidth/4.0)
## Smoothing:: https://joseph-long.com/writing/AstroPy-boxcar/
#DBSP_b_flux_smoothed = convolve(DBSP_b_flux, Box1DKernel(5))
#DBSP_r_flux_smoothed = convolve(DBSP_r_flux, Box1DKernel(5))
#ax2.plot(DBSP_b_wavelength/(1+redshift), DBSP_b_flux_smoothed, '-k', lw=linewidth/4.0)
#ax2.plot(DBSP_r_wavelength/(1+redshift), DBSP_r_flux_smoothed, '-k', lw=linewidth/4.0,     label='DBSP MJD 58538')

## Putting on the Emission Line names and vertical lines
## e.g. https://scipython.com/blog/rotating-text-onto-a-line-in-matplotlib/
for ll in range(len(linelist)):
    ax2.axvline(x=linelist['Wavelength'][ll], color='gray', linestyle='--', linewidth=linewidth/3.4)
    label = linelist['LineName'][ll]
    xylabel = ((linelist['Wavelength'][ll]), ymax_ax2/1.1)
    ax2.annotate(label, xy=xylabel, ha='center', va='center', rotation=90, fontsize=fontsize/1.6 )


## Just plotting the spectra again, the axes limits take care of the "zoom in"
ax3.plot(   sdss_wavelength/(1+redshift), sdss_flux,    '-b', lw=linewidth/4.0)
ax3.plot(   boss_wavelength/(1+redshift), boss_flux,    '-r', lw=linewidth/4.0)
ax3.plot(DBSP_b_wavelength/(1+redshift), DBSP_b_flux, '-k', lw=linewidth/4.0)
#ax3.plot(DBSP_b_wavelength/(1+redshift), DBSP_b_flux_smoothed, '-k', lw=linewidth/4.0)
for ll in range(len(linelist)):
    ax3.axvline(x=linelist['Wavelength'][ll], color='gray', linestyle='--', linewidth=linewidth/3.4)
    label = linelist['LineName'][ll]
    xylabel = ((linelist['Wavelength'][ll]), ymax_ax3/1.1)
    ax3.annotate(label, xy=xylabel, ha='center', va='center', rotation=90, fontsize=fontsize/1.6 )


ax4.plot(   sdss_wavelength/(1+redshift), sdss_flux,  '-b', lw=linewidth/4.0)
ax4.plot(   boss_wavelength/(1+redshift), boss_flux,  '-r', lw=linewidth/4.0)
ax4.plot(DBSP_b_wavelength/(1+redshift), DBSP_b_flux, '-k', lw=linewidth/4.0)
#ax4.plot(DBSP_r_wavelength/(1+redshift), DBSP_r_flux, '-k', lw=linewidth/4.0)
#ax4.plot(DBSP_b_wavelength/(1+redshift), DBSP_b_flux_smoothed, '-k', lw=linewidth/4.0)
#ax4.plot(DBSP_r_wavelength/(1+redshift), DBSP_r_flux_smoothed, '-k', lw=linewidth/4.0)
for ll in range(len(linelist)):
    ax4.axvline(x=linelist['Wavelength'][ll], color='gray', linestyle='--', linewidth=linewidth/3.4)
    label = linelist['LineName'][ll]
    xylabel = ((linelist['Wavelength'][ll]), ymax_ax4/1.1)
    ax4.annotate(label, xy=xylabel, ha='center', va='center', rotation=90, fontsize=fontsize/1.6 )


ax5.plot(   sdss_wavelength/(1+redshift),    sdss_flux,    '-b', lw=linewidth/4.0)
ax5.plot(   boss_wavelength/(1+redshift),    boss_flux,    '-r', lw=linewidth/4.0)
ax5.plot(DBSP_r_wavelength/(1+redshift), DBSP_r_flux, '-k', lw=linewidth/4.0)
#ax5.plot(DBSP_r_wavelength/(1+redshift), DBSP_r_flux_smoothed, '-k', lw=linewidth/4.0)
for ll in range(len(linelist)):
    ax5.axvline(x=linelist['Wavelength'][ll], color='gray', linestyle='--', linewidth=linewidth/3.4)
    label = linelist['LineName'][ll]
    xylabel = ((linelist['Wavelength'][ll]), ymax_ax5/1.15)
    ax5.annotate(label, xy=xylabel, ha='center', va='center', rotation=90, fontsize=fontsize/1.6 )

## Making all the nice plots boxes the right dimensions and placings
#% start: automatic generated code from pylustrator
plt.figure(1).ax_dict = {ax.get_label(): ax for ax in plt.figure(1).axes}
import matplotlib as mpl
plt.figure(1).axes[0].get_xaxis().get_label().set_fontsize(fontsize)
plt.figure(1).axes[0].get_yaxis().get_label().set_fontsize(fontsize)
plt.figure(1).axes[0].set_position([0.060000, 0.136250, 0.345, 0.823750])
plt.figure(1).axes[0].xaxis.labelpad = 10.0

plt.figure(1).axes[1].get_xaxis().get_label().set_fontsize(fontsize/1.2)
plt.figure(1).axes[1].get_yaxis().get_label().set_fontsize(fontsize)
plt.figure(1).axes[1].get_yaxis().get_label().set_rotation(90.0)
x_nudge = 0.012
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


##  A X E S   L A B E L S

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


##
##  L E G E N D S
##
##  LIGHT CURVE Legend
neo_w1 = mlines.Line2D([], [], label='NEOWISE W1', color='red',
                       marker="o", markeredgecolor='k', markeredgewidth=2.0, markersize=12, linestyle='None')
neo_w2 = mlines.Line2D([], [], label='NEOWISE W2', color='cyan',
                       marker="o", markeredgecolor='k', markeredgewidth=2.0, markersize=12, linestyle='None')
crts   = mlines.Line2D([], [], label='CRTS',       color='k',
                       marker="o", markeredgecolor='k', markeredgewidth=2.0, markersize=5,  linestyle='None')
ps_g   = mlines.Line2D([], [], label=r'Pan-STARRS $g$-band', color='lime',
                       marker="o", markeredgecolor='k', markeredgewidth=1.4, markersize=7,  linestyle='None')
ps_r   = mlines.Line2D([], [], label=r'Pan-STARRS $r$-band', color='deeppink',
                       marker="o", markeredgecolor='k', markeredgewidth=1.4, markersize=7,  linestyle='None')
ztf_g  = mlines.Line2D([], [], label=r'ZTF $g$-band', color='olivedrab',
                       marker="o", markeredgecolor='k', markeredgewidth=1.4, markersize=7,  linestyle='None')
ztg_r  = mlines.Line2D([], [], label=r'ZTF $r$-band', color='tomato',
                       marker="o", markeredgecolor='k', markeredgewidth=1.4, markersize=7,  linestyle='None')

handles=[neo_w1, neo_w2, crts,ps_g, ps_r, ztf_g, ztg_r]
leg = ax1.legend(loc='upper left',
                     fontsize=fontsize/1.6,
                     handles=handles,
                    frameon=True, framealpha=1.0, fancybox=True)

## SPECTRA Legend
leg = ax2.legend(loc='upper right', fontsize=fontsize/1.2,
                 frameon=True, framealpha=1.0, fancybox=True)

for line in leg.get_lines(): line.set_linewidth(2.2)

plt.savefig('J2228+2201_landscape_temp.png', format='png')
plt.close(fig)
