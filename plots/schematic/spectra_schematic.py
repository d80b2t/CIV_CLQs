'''
A plot with an e.g. the UV/blue/red quasar spectrum (straight from
Vanden Berk et al. (2001) with the given emission lines, and their
Ionization Energies.
'''

import numpy as np
from astropy.io          import ascii

import matplotlib.lines    as mlines
import matplotlib.pyplot   as plt
import matplotlib.gridspec as gridspec

from matplotlib          import colors as mcolors
from matplotlib.gridspec import GridSpec
from matplotlib.patches  import Rectangle



path = '/cos_pc19a_npr/programs/quasars/CIV_CLQs/plots/spectra_schematic/'

##
## Vanden_Berk composite spectrum
##
spectrum_file = 'Vanden_Berk01_Table1.dat'
spectrum = ascii.read(path+spectrum_file)

wave    = spectrum['Wave']
flux    = spectrum['FluxD']
fluxerr = spectrum['e_FluxD']

##
## Emission Line list
## "EMISSION-LINE VELOCITY SHIFTS RELATIVE TO [O III] 5007"
##
linelist_file = 'Vanden_Berk01_Table4.dat'
linelist = ascii.read(path+linelist_file)


##  S E T T I N G    U P   T H E    P L O T 
fig, ax = plt.subplots(figsize=(14, 8), dpi=200, facecolor='w', edgecolor='k')

##  Dark background :-) 
plt.style.use('dark_background')

##  Get things looking "tex-y"
plt.rcParams['text.usetex'] = True

## Adjusting the Whitespace for the plots
left   = 0.10   # the left side of the subplots of the figure
right  = 0.98   # the right side of the subplots of the figure
bottom = 0.18   # the bottom of the subplots of the figure
top    = 0.96   # the top of the subplots of the figure
wspace = 0.26   # the amount of width reserved for blank space between subplots
hspace = 0.06   # the amount of height reserved for white space between subplots
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
 
## Some NPR defaults
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


##  P L O T I N G   T H I N  GS!!
flux_norm = flux / flux.max()

#ax.plot(wave, flux,      alpha=1.0, linewidth=linewidth)
ax.plot(wave, flux_norm, alpha=1.0, linewidth=linewidth)


## Axes limits
xmin =    700.0   # Angstroms  
xmax =   8000.0   # Angstrongs
ymin =      0.0
ymax =      1.2   # Using flux_norm
ax.set_xlim([xmin, xmax])
ax.set_ylim([ymin, ymax])

xmin_log =    750.0   # Angstroms  
xmax_log =   8000.0   # Angstrongs
ymin_log =      0.03   # Angstroms  
ymax_log =      3.4   # Angstrongs

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim([xmin_log, xmax_log])
ax.set_ylim([ymin_log, ymax_log])


## Axes labels
ax.set_xlabel(r'Rest Wavelength / ${ \rm \AA }$',               fontsize=fontsize*1.4)
#ax.set_ylabel(r'Flux Density,  $f_{\lambda}$ / 10$^{-17}$ erg/s/cm$^2$/${ \rm \AA }$', fontsize=fontsize)
ax.set_ylabel(r'Flux Density,  $f_{\lambda}$ (arbitray units)', fontsize=fontsize*1.4)


for ii in range(len(linelist)):
    ax.axvline(x=linelist['lambda_lab'][ii], color='gray', linestyle='--', linewidth=linewidth/3.4)
    label       =    linelist['ion'][ii]
    label_nrg   =    str(linelist['ionization_energy'][ii])
    xylabel     =   (( linelist['lambda_lab'][ii]),  (((ymax+1.25)/1.15)-(0.03*linelist['label_offset'][ii])))
    xylabel_nrg =   (( linelist['lambda_lab'][ii]),  (((ymax+0.90)/1.15)-(0.03*linelist['label_offset'][ii])))
    print(label, xylabel)
    ax.annotate(label,     xy=xylabel,     ha='center', va='center', fontsize=fontsize/1.4)
    ax.annotate(label_nrg, xy=xylabel_nrg, ha='center', va='center', fontsize=fontsize/1.4)

## Axes style
ax.tick_params(axis='both', which='major', labelsize=labelsize, top=True, right=True, direction='in', length=ticklength,   width=tickwidth)
ax.tick_params(axis='both', which='minor', labelsize=labelsize, top=True, right=True, direction='in', length=ticklength/2, width=tickwidth)


plt.savefig('Quasar_spoectrum_lines_and_IP_temp.png', format='png')
plt.close(fig)
