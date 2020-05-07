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

## Get the redshift right!
redshift = 2.222
print()
print('Redshift of J2228+2201 is ', redshift)
print()


##
##  S E T T I N G    U P   T H E   P L O T !!!
##  
fig, ax1 = plt.subplots(figsize=(12, 4.3))

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
ls              = 'solid'
alpha           = 1.0
fontsize        = 20
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
ms = 4.
ax1.scatter( CRTS_MJD, CRTS['magnitude'], color='k',       alpha=alpha, s=ms*1.8, zorder=0)
ax1.scatter( CRTS_MJD, CRTS['magnitude'], color='dimgray', alpha=alpha, s=ms,     zorder=0, label='CRTS')
ax1.errorbar(CRTS_MJD, CRTS['magnitude'], color='k', yerr=CRTS['magnitude_err'],  zorder=0, fmt='o', linewidth=lw, ms=ms)

## ZTF data
ms_big          = 28.
ms              = 6.
ax1.errorbar(ZTF_g['mjd'], ZTF_g['mag'], color='olivedrab', yerr=ZTF_g['magerr'], fmt='o', linewidth=lw, ms=ms/2)
ax1.scatter( ZTF_g['mjd'], ZTF_g['mag'], color='k',         alpha=alpha, s=ms_big)
ax1.scatter( ZTF_g['mjd'], ZTF_g['mag'], color='olivedrab', alpha=alpha, s=ms, label='ZTF g-band')
ax1.errorbar(ZTF_r['mjd'], ZTF_r['mag'], color='tomato',    yerr=ZTF_r['magerr'], fmt='o', linewidth=lw, ms=ms/2)
ax1.scatter( ZTF_r['mjd'], ZTF_r['mag'], color='k',         alpha=alpha, s=ms_big)
ax1.scatter( ZTF_r['mjd'], ZTF_r['mag'], color='tomato',    alpha=alpha, s=ms, label='ZTF r-band')

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


##  A X E S   L A B E L S
plt.rcParams['text.usetex'] = True

ax1.set_xlabel(r'MJD',          fontsize=fontsize)
ax1.set_ylabel(r'AB Magnitude', fontsize=fontsize)
ax1.tick_params(labelsize=labelsize, direction='in', length=ticklength, width=tickwidth)

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

## Quasar Lable Name
plt.text(54700, 16.8, r'{\bf J2228+2201}',  fontsize=fontsize*1.2)

plt.savefig('J2228+2201_landscape_LC_temp.png', format='png')
plt.close(fig)
