
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from matplotlib         import colors as mcolors
from matplotlib.ticker  import MaxNLocator
from matplotlib.colors  import BoundaryNorm
from matplotlib.patches import Rectangle


from astropy.io import ascii
from astropy.io import fits

## https://docs.scipy.org/doc/scipy/reference/stats.html
## 
from scipy.stats import kde


## Reading in the data
path      = '/cos_pc19a_npr/data/ERQs/Hamann2017_CIVcatalog/'
filename  = 'C4N5REWs_DR12v11_MNRAS.fits'
infile    = path+filename
data_full = fits.open(infile)
tdata     = data_full[1].data

data = tdata[np.where(tdata['rew'] < 10000.)]
## Setting up variable names
REW  = data['rew']        
FWHM = data['fwhm']


data = np.vstack((REW, FWHM)).T                                                                                                                       
x, y = data.T


## Setting up the plot
#fig, ax = plt.subplots(figsize=(5, 3), dpi=80, facecolor='w', edgecolor='k')
## Create a figure with 6 plot areas
fig, axes = plt.subplots(ncols=6, nrows=1, figsize=(21, 5))


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

 
# Everything sarts with a Scatterplot
axes[0].set_title('Scatterplot')
axes[0].plot(x, y, 'ko')
# As you can see there is a lot of overplottin here!
 
# Thus we can cut the plotting window in several hexbins
nbins = 40
axes[1].set_title('Hexbin')
axes[1].hexbin(x, y, gridsize=nbins, mincnt=10, cmap=plt.cm.BuGn_r)
 
# 2D Histogram
axes[2].set_title('2D Histogram')
axes[2].hist2d(x, y, bins=nbins, cmap=plt.cm.BuGn_r)

rew_min = 2.0
# Evaluate a gaussian kde on a regular grid of nbins x nbins over data extents
rew_min  =    2.0
rew_max  =  2000.
fwhm_min =   100.
fwhm_max = 15000.

k = kde.gaussian_kde(data.T)
#xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
xi, yi = np.mgrid[rew_min:rew_max:nbins*1j, fwhm_min:fwhm_max:nbins*1j]
zi     = k(np.vstack([xi.flatten(), yi.flatten()]))

cmap = plt.get_cmap('BuGn_r')
levels = MaxNLocator(nbins=nbins).tick_values(1, 100)
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

# plot a density
axes[3].set_title('Calculate Gaussian KDE')
axes[3].pcolormesh(xi, yi, zi.reshape(xi.shape), cmap=cmap)
 
# add shading
axes[4].set_title('2D Density with shading')
axes[4].pcolormesh(xi, yi, zi.reshape(xi.shape), shading='gouraud', cmap=cmap, norm=norm)
 
# contour
axes[5].set_title('Contour')
axes[5].pcolormesh(xi, yi, zi.reshape(xi.shape), shading='gouraud', cmap=cmap)
axes[5].contour(xi, yi, zi.reshape(xi.shape) )


## Axes limits
xmin =     1.
xmax =  2000.
ymin =     0.0
ymax = 15000.
axes[0].set_xlim([xmin,xmax])
axes[1].set_xlim([xmin,xmax])
axes[2].set_xlim([xmin,xmax])
axes[3].set_xlim([xmin,xmax])
axes[4].set_xlim([xmin,xmax])
axes[5].set_xlim([xmin,xmax])
#ax.set_ylim([ymin, ymax])
axes[0].set_ylim([ymin,ymax])
axes[1].set_ylim([ymin,ymax])
axes[2].set_ylim([ymin,ymax])
axes[3].set_ylim([ymin,ymax])
axes[4].set_ylim([ymin,ymax])
axes[5].set_ylim([ymin,ymax])

## Axes labels
#ax.set_xlabel(r'$z$, redshift')
#$ax.set_ylabel(r'No. of objects')
axes[0].set_xscale('log')
axes[1].set_xscale('log')
axes[2].set_xscale('log')
axes[3].set_xscale('log')
axes[4].set_xscale('log')
axes[5].set_xscale('log')

## Axes style
#ax.tick_params(axis='both', which='major', labelsize=labelsize, top=True, right=True, direction='in', length=ticklength,   width=tickwidth)
#ax.tick_params(axis='both', which='minor', labelsize=labelsize, top=True, right=True, direction='in', length=ticklength/2, width=tickwidth)


plt.savefig('CIV_CLQs_contours_temp.png', format='png')
plt.close(fig)
