
import numpy   as np
import seaborn as sns
import matplotlib.pyplot as plt

from matplotlib.patches import Rectangle
from matplotlib         import colors as mcolors

from astropy.io import ascii
from astropy.io import fits


path = '/cos_pc19a_npr/programs/quasars/CIV_CLQs/data/Lags'

## Reading in the data
data_NGC5548 = ascii.read(path+infile)


data_Mrk110 = ascii.read(path+infile)




## Setting up the plot
fig, ax = plt.subplots(figsize=(5, 3), dpi=80, facecolor='w', edgecolor='k')

## Adjusting the Whitespace for the plots
left   = 0.18   # the left side of the subplots of the figure
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


## Axes limits
xmin =  4.90
xmax =  7.65
ymin =  0.0
ymax = 38.0
ax.set_xlim([xmin,xmax])
ax.set_ylim([ymin, ymax])

## Axes labels
ax.set_xlabel(r'$z$, redshift')
ax.set_ylabel(r'No. of objects')

## Axes style
ax.tick_params(axis='both', which='major', labelsize=labelsize, top=True, right=True, direction='in', length=ticklength,   width=tickwidth)
ax.tick_params(axis='both', which='minor', labelsize=labelsize, top=True, right=True, direction='in', length=ticklength/2, width=tickwidth)


plt.savefig('matplotlib_template_temp.png', format='png')
plt.close(fig)
