
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns


from matplotlib.patches import Rectangle
from matplotlib         import colors as mcolors

from astropy.io import ascii
from astropy.io import fits
#from astroquery.sdss import SDSS



def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        print(len(x.ndim))
        #raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        print(len(x.ndim))
        #raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        print(len(x.ndim))
        #raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y







print()
whichQuasar = input("Which quasar??    (a) J1205+3422;  (b) J1638+2827;  (c) J2228+2201 ")
print()


datapath = '/cos_pc19a_npr/programs/quasars/CIV_CLQs/data/'

##
## Data model: spec
## https://data.sdss.org/datamodel/files/BOSS_SPECTRO_REDUX/RUN2D/spectra/PLATE4/spec.html
##
if (whichQuasar == 'a'):
    print(' J1205+3422 ' )
    name = 'J1205+3422'
    initial_spectrum = fits.open(datapath+name+'/spec-2089-53498-0427.fits')
    #secondEpoch_spectrum = ''
    #thirdEpoch_spectrum = ''
    #final_spectrum = ''
    
elif (whichQuasar == 'b'):
    print(' J1638+2827 ' )
    name = 'J1638+2827'
    initial_spectrum = fits.open(datapath+name+'/spec-2948-54553-0614.fits')
    second_spectrum = fits.open(datapath+name+'/spec-5201-55832-0178_v5_10_0.fits')
else :
    print(' J2228+2201 ' )
    name = 'J2228+2201/'
    initial_spectrum = fits.open(datapath+name+'spec-2089-53498-0427.fits')

## Getting the redshift...
zstruc = initial_spectrum[2].data
redshift = zstruc.field('Z')

initial_spectrum_data       = initial_spectrum[1].data
initial_spectrum_flux       = initial_spectrum_data.flux
initial_spectrum_loglam     = initial_spectrum_data.loglam
initial_spectrum_wavelength = 10**(initial_spectrum_loglam)

second_spectrum_data        = second_spectrum[1].data
second_spectrum_flux        = second_spectrum_data.flux
second_spectrum_loglam      = second_spectrum_data.loglam
second_spectrum_wavelength  = 10**(second_spectrum_loglam)

## In [12]: initial_spectrum_loglam[0]
## Out[12]: 3.5824
##
## In [13]: second_spectrum_loglam[311]
## Out[13]: 3.5824

## In [16]: initial_spectrum_loglam[-1]
## Out[16]: 3.9642
##
## In [23]: second_spectrum_loglam[4129]
## Out[23]: 3.9642

second_spectrum_trimmed            = second_spectrum[311:4130]
second_spectrum_trimmed_flux       = second_spectrum_data[311:4130].flux
second_spectrum_trimmed_loglam     = second_spectrum_data[311:4130].loglam
second_spectrum_trimmed_wavelength = 10**(second_spectrum_trimmed_loglam)

fluxRatio = (initial_spectrum_flux/second_spectrum_trimmed_flux)**(-1)


## Setting up the plot
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True,
                               gridspec_kw={'height_ratios': [2, 1]}, 
                               figsize=(12, 8), dpi=80, facecolor='w', edgecolor='k')

## Adjusting the Whitespace for the plots
left   = 0.12   # the left side of the subplots of the figure
right  = 0.96   # the right side of the subplots of the figure
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


##
##   O B S E R V E D
##
# Plot the resulting spectrum
#fig, ax = plt.subplots(figsize=(5, 3.75))
ax1.plot(initial_spectrum_wavelength, initial_spectrum_flux, '-k', lw=linewidth/3.0)
ax1.plot(second_spectrum_wavelength,   second_spectrum_flux, '-r', lw=linewidth/3.0)
if (whichQuasar == 'b'):
    #ax2.plot(second_spectrum_trimmed_wavelength, fluxRatio, '-k', lw=linewidth/3.0)
    ax2.plot(smooth(second_spectrum_trimmed_wavelength,window_len=21), smooth(fluxRatio, window_len=21,window='hanning'), color='tab:blue')
    

## Axes limits
xmin =  3000.  ## Angstroms
xmax = 10000.  ## Angstroms
ymin =   -5.0    ## f_lambda  /  10^-17  erg/s/cm^2/Ang  
ymax =   second_spectrum_flux.max()*1.4   ## f_lambda  /  10^-17  erg/s/cm^2/Ang

ax1.set_xlim([xmin,xmax])
ax1.set_ylim([ymin, ymax])

ymin =    0.0    ## f_lambda  /  10^-17  erg/s/cm^2/Ang  
ymax =    2.0

ax2.set_xlim([xmin,xmax])
ax2.set_ylim([ymin, ymax])


fig.subplots_adjust(hspace=0)

## Axes labels
matplotlib.rcParams['text.usetex'] = True
ax2.set_xlabel(r'Wavelength / ${ \rm \AA }$')
#ax2.set_xlabel(r'Wavelength $\lambda {(\rm \AA)}$')
ax1.set_ylabel(r'Flux, $f_{\lambda}$ / 10$^{-17}$ erg/s/cm$^2$/${ \rm \AA }$')
## Title label
#ax.set_title(name)

## Axes style
ax1.tick_params(axis='both', which='major', labelsize=labelsize, top=True, right=True, direction='in', length=ticklength,   width=tickwidth)
ax1.tick_params(axis='both', which='minor', labelsize=labelsize, top=True, right=True, direction='in', length=ticklength/2, width=tickwidth)
ax2.tick_params(axis='both', which='major', labelsize=labelsize, top=True, right=True, direction='in', length=ticklength,   width=tickwidth)
ax2.tick_params(axis='both', which='minor', labelsize=labelsize, top=True, right=True, direction='in', length=ticklength/2, width=tickwidth)


plt.savefig('plotSpectra_temp.png', format='png')
plt.close(fig)

''''
##
##   R E S T   -   F R A M E
##
# Plot the resulting spectrum
fig, ax = plt.subplots(figsize=(5, 3.75))
ax.plot((initial_spectrum_wavelength/(1.+redshift)), initial_spectrum_flux, '-k', lw=linewidth/2.0)

## Axes limits
xmin =  3000.  ## Angstroms
xmax = 10000.  ## Angstroms
ymin =   -10.0    ## f_lambda  /  10^-17  erg/s/cm^2/Ang  
ymax =   second_spectrum_flux.max()*1.4   ## f_lambda  /  10^-17  erg/s/cm^2/Ang

ax.set_xlim([xmin,xmax])
ax.set_ylim([ymin, ymax])

## Axes labels
ax.set_xlabel(r'Wavelegnth $\lambda {(\rm \AA)}$')
ax.set_ylabel(r'Flux, f_{\lambda} / 10^{-17} erg/s/cm$^2$/Ang')
## Title label
#ax.set_title(name)

## Axes style
ax.tick_params(axis='both', which='major', labelsize=labelsize, top=True, right=True, direction='in', length=ticklength,   width=tickwidth)
ax.tick_params(axis='both', which='minor', labelsize=labelsize, top=True, right=True, direction='in', length=ticklength/2, width=tickwidth)

plt.savefig('plotSpectraRest_temp.png', format='png')
plt.close(fig)
'''





