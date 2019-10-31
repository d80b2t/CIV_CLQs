'''
"Baby" python code that will reproduce, e.g. Figure 2 of Wilhite et
al. 2006, ApJ, 641:78
'''

import numpy   as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.lines  as mlines

from matplotlib          import colors as mcolors
from matplotlib.patches  import Rectangle
from matplotlib.gridspec import GridSpec

from astropy.io import ascii
from astropy.io import fits

## Reading in the data from Wilhite et al., 2006, ApJ, 641 78
Wilhite_path = '/cos_pc19a_npr/data/SDSS/Wilhite_CIV/'
Wilhite_file = 'Wilhite_2006_combined_Tables.dat'
Wilhite_input = Wilhite_path+Wilhite_file

data = ascii.read(Wilhite_input) 

## Set up the Line flux ratio variable...
flux_ratio = data['HiSNRlineFlux']/data['LoSNRlineFlux']
## and the Line Width ratio variable...
LineWidth_ratio =  data['HiSNRlineWidth']/data['LoSNRlineWidth']


## J1638+2827
J1638_path  = '/cos_pc19a_npr/programs/quasars/CIV_CLQs/data/CLQ_line_measurements/'
J1638_file  = 'LineStrength_ratios.dat'
J1638_input = J1638_path+J1638_file
J1638_data  = ascii.read(J1638_input)


## Setting up the plot
fig, ax = plt.subplots(figsize=(7,5), dpi=80, facecolor='w', edgecolor='k')

## Adjusting the Whitespace for the plots
left   = 0.14   # the left side of the subplots of the figure
right  = 0.96   # the right side of the subplots of the figure
bottom = 0.14   # the bottom of the subplots of the figure
top    = 0.96   # the top of the subplots of the figure
wspace = 0.26   # the amount of width reserved for blank space between subplots
hspace = 0.06   # the amount of height reserved for white space between subplots
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
 
## Some NPR defaults
ms              = 10.
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

## Some general links/tips on stylings::
##
## https://matplotlib.org/3.1.0/tutorials/introductory/customizing.html
## https://stackoverflow.com/questions/25451294/best-way-to-display-seaborn-matplotlib-plots-with-a-dark-ipython-notebook-profil
## https://smashingtheory.blogspot.com/search/label/matplotlib
## https://earthobservatory.nasa.gov/blogs/elegantfigures/2013/08/05/subtleties-of-color-part-1-of-6/

## ax.errorbar(CRTS['MJD'],   CRTS['mag'],   yerr=CRTS['mag_err'],   fmt='o', linewidth=lw, ms=ms, color='w')
ax.plot(flux_ratio,                 LineWidth_ratio ,              '+', label='Wilhite 2006')
ax.plot(J1638_data['fluxRatio'][2], J1638_data['sigmaRatio'][2]  , 's', label='J1638+2827', color='k', markersize=ms*1.0,)
ax.plot(J1638_data['fluxRatio'][2], J1638_data['sigmaRatio'][2]  , 's', label='J1638+2827', color='lime')

plt.hlines(1.0, 0, 10, linestyles='dashed',   linewidth=linewidth/4.)
plt.vlines(1.0, 0, 10, linestyles='dashed',   linewidth=linewidth/6.)

##  A X E S    L I M I T S 
xmin = 0.20             ; xmax = 3.8              ; ymin = 0.60                  ; ymax = 1.4    ## Wilhite et al. 2006 limits
#xmin = flux_ratio.min() ; xmax = flux_ratio.max() ; ymin = LineWidth_ratio.min() ; ymax = LineWidth_ratio.max() ## 
ax.set_xlim([xmin,xmax])
ax.set_ylim([ymin, ymax])

##  A X E S    L A B E L S 
ax.set_xlabel(r'CIV line flux  ratio', fontsize=fontsize)
ax.set_ylabel(r'CIV line width ratio', fontsize=fontsize)

##  A X E S    S T Y L E 
ax.tick_params(axis='both', which='major', labelsize=labelsize, top=True, right=True, direction='in', length=ticklength,   width=tickwidth)
ax.tick_params(axis='both', which='minor', labelsize=labelsize, top=True, right=True, direction='in', length=ticklength/2, width=tickwidth)


##  L E G E N D
sdss  = mlines.Line2D([], [], label='Wilhite et al. (2006)', color='C0',
                       marker="+", markeredgecolor='C0', markeredgewidth=1.4, markersize=7,  linestyle='None')
j1638  = mlines.Line2D([], [], label='J1638+2827',            color='lime',
                       marker="s", markeredgecolor='k', markeredgewidth=1.4, markersize=7,  linestyle='None')

handles=[sdss, j1638]
leg = ax.legend(loc='upper right', fontsize=fontsize/1.25, handles=handles, 
                frameon=True, framealpha=1.0, fancybox=True)



plt.savefig('Wilhite_2006_Fig2_redux_temp.png', format='png')
plt.close(fig)
