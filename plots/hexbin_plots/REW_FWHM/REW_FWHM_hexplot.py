'''
A ``hexbin'' plot for the SDSS-III BOSS data from the Hamann et
al. (2017) catalog and analysis. This code is for rest equivalent
widths (REWs), and Full width at half maximum (FWHM). 
'''

import numpy as np 
from astropy.io import fits
from astropy.io import ascii

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.lines  as mlines

##
## READING IN THE DATA
##
##  F U L L    B O S S    D R 1 2    catalog
path      = '/cos_pc19a_npr/data/ERQs/Hamann2017_CIVcatalog/'
filename  = 'C4N5REWs_DR12v11_MNRAS.fits'
infile    = path+filename
data_full = fits.open(infile)
tdata     = data_full[1].data
## knocking out a couple of objects with bad REW values 
data = tdata[np.where( (tdata['rew'] >0.) & (tdata['rew'] < 10000.)  )]

## Setting up variable names
REW  = data['rew']        
FWHM = data['fwhm']

## Putting te Equiv. Widths into log and
## the FWHMs into '000s of kms^-1.
log_REW = np.log10(REW)
FWHM    = FWHM/1000.


##  T H E    C I V    C L Q    S A M P L E 
path      = '/cos_pc19a_npr/programs/quasars/CIV_CLQs/data/CLQ_line_measurements/'
filename  = 'REW_FWHM_CIV.dat'
CLQs      = ascii.read(path+filename)

CLQ_FWHM    = CLQs['FWHM']    / 1000.
CLQ_errFWHM = CLQs['errFWHM'] / 1000.
CLQ_REW     = CLQs['REW']
CLQ_errREW  = CLQs['errREW'] 

CLQ_logREW = np.log10(CLQs['REW'])

## have to put the x-errors into log-space
#log_REW_hi = np.log10(CLQs['REW'] np.log10(CLQs['errREW']))
#log_REW_lo = np.log10(CLQs['REW']-CLQs['errREW'])
## and then have them as
#CLQ_asym_xerr = [log_REW_lo, log_REW_hi] 
log_errREW = np.log10(CLQs['errREW'])



##   S E T T I N G   U P   T H E    P L O T  
matplotlib.rc('text', usetex=True)
fig, ax = plt.subplots(figsize=(5, 5), dpi=80, facecolor='w', edgecolor='k')

## Adjusting the Whitespace for the plots
left   = 0.16   # the left side of the subplots of the figure
right  = 0.96   # the right side of the subplots of the figure
bottom = 0.14   # the bottom of the subplots of the figure
top    = 0.96   # the top of the subplots of the figure
wspace = 0.26   # the amount of width reserved for blank space between subplots
hspace = 0.06   # the amount of height reserved for white space between subplots
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
 
## Some NPR defaults
s               = 0.5
lw              = 1.0
alpha           = 0.4
fontsize        = 16
labelsize       = fontsize
tickwidth       = 2.0
linewidth       = 2.4
tickwidth       = 2.0
ticklength      = 6.0
ticklabelsize   = labelsize
majorticklength = 12
minorticklength = 6

## Furthe defaults and  hexbin params
mincnt    = 3.     ##   3 works well for log(REW)
gridsize  = 200    ## 200 works well for log(REW)
gridsize  = 320    
color_map = plt.cm.Spectral_r


##  P L O T T I N G    T H E    D A T A
##  LOG IN THE REW:
#ax.scatter(log_REW, FWHM, s=s, alpha=alpha)
#hb = plt.hexbin(log_REW, FWHM, bins='log',gridsize=gridsize, cmap=color_map, mincnt=mincnt)
ax.scatter(REW, FWHM, s=s, alpha=alpha)
hb = plt.hexbin(REW, FWHM, bins='log', gridsize=gridsize, cmap=color_map, mincnt=mincnt)


for ii in range(len(CLQs)):
    ms      =  12.
    ms_back =   7.
    l_back  = 2
    alpha  =  1.0
    
    if str(CLQs['name'][ii]) == 'J1205+3422':
        if (CLQs['mjd'][ii] ==  53498):
            ax.errorbar(CLQ_REW[ii], CLQ_FWHM[ii], color='fuchsia', xerr=CLQ_errREW[ii], yerr=CLQ_errFWHM[ii], fmt='o', linewidth=lw)
            ax.scatter( CLQ_REW[ii], CLQ_FWHM[ii], color='k',       alpha=alpha, marker='o', s=ms*ms_back)
            ax.scatter( CLQ_REW[ii], CLQ_FWHM[ii], color='fuchsia', alpha=alpha, marker='o', s=ms)
            #ax.scatter(CLQ_logREW[ii],  CLQ_FWHM[ii], color='k',       alpha=alpha, marker='o', s=s*2.2)
            #ax.scatter(CLQ_logREW[ii],  CLQ_FWHM[ii], color='fuchsia', alpha=alpha, marker='o', s=s)
            #ax.errorbar(CLQ_logREW[ii], CLQ_FWHM[ii], color='fuchsia', xerr=log_errREW[ii], yerr=errFWHM[ii], fmt='o', linewidth=lw)
            J12_53498  = mlines.Line2D([], [], label='J1205+3422 (53498)', color='fuchsia',
                       marker="o", markeredgecolor='k', markeredgewidth=1.4, markersize=7,  linestyle='None')
        if (CLQs['mjd'][ii] ==  58693):
            ax.errorbar(CLQ_REW[ii], CLQ_FWHM[ii], color='k',       xerr=CLQ_errREW[ii], yerr=CLQ_errFWHM[ii], fmt='s', linewidth=lw*l_back)
            ax.errorbar(CLQ_REW[ii], CLQ_FWHM[ii], color='fuchsia', xerr=CLQ_errREW[ii], yerr=CLQ_errFWHM[ii], fmt='s', linewidth=lw)
            ax.scatter( CLQ_REW[ii], CLQ_FWHM[ii], color='k',       alpha=alpha, marker='s', s=ms*ms_back)
            ax.scatter( CLQ_REW[ii], CLQ_FWHM[ii], color='fuchsia', alpha=alpha, marker='s', s=ms)
            J12_58693 = mlines.Line2D([], [], label='J1205+3422 (58693)', color='fuchsia',
                       marker="s", markeredgecolor='k', markeredgewidth=1.4, markersize=7,  linestyle='None')
            
    if str(CLQs['name'][ii]) == 'J1638+2827':
        if (CLQs['mjd'][ii] ==  54553):
            ax.errorbar(CLQ_REW[ii], CLQ_FWHM[ii], color='lime', xerr=CLQ_errREW[ii], yerr=CLQ_errFWHM[ii], fmt='o', linewidth=lw)
            ax.scatter( CLQ_REW[ii], CLQ_FWHM[ii], color='k',    alpha=alpha, marker='o', s=ms*ms_back)
            ax.scatter( CLQ_REW[ii], CLQ_FWHM[ii], color='lime', alpha=alpha, marker='o', s=ms)

            J16_54553 = mlines.Line2D([], [], label='J1638+2827 (54553)', color='lime',
                       marker="o", markeredgecolor='k', markeredgewidth=1.4, markersize=7,  linestyle='None')
        if (CLQs['mjd'][ii] ==  55832): 
            ax.scatter(CLQ_REW[ii],  CLQ_FWHM[ii], color='k',    alpha=alpha, marker='s', s=ms*ms_back)
            ax.scatter(CLQ_REW[ii],  CLQ_FWHM[ii], color='lime', alpha=alpha, marker='s', s=ms)
            ax.errorbar(CLQ_REW[ii], CLQ_FWHM[ii], color='lime', xerr=log_errREW[ii], yerr=CLQ_errFWHM[ii], fmt='s', linewidth=lw)
            J16_55832 = mlines.Line2D([], [], label='J1638+2827 (55832)', color='lime',
                       marker="s", markeredgecolor='k', markeredgewidth=1.4, markersize=7,  linestyle='None')
        if (CLQs['mjd'][ii] ==  58583): 
            ax.scatter(CLQ_REW[ii],  CLQ_FWHM[ii], color='k',    alpha=alpha, marker='D', s=ms*ms_back)
            ax.scatter(CLQ_REW[ii],  CLQ_FWHM[ii], color='lime', alpha=alpha, marker='D', s=ms)
            ax.errorbar(CLQ_REW[ii], CLQ_FWHM[ii], color='lime', xerr=log_errREW[ii], yerr=CLQ_errFWHM[ii], fmt='D', linewidth=lw)
            J16_58583 = mlines.Line2D([], [], label='J1638+2827 (58583)', color='lime',
                       marker="D", markeredgecolor='k', markeredgewidth=1.4, markersize=7,  linestyle='None')
            
    if str(CLQs['name'][ii]) == 'J2228+2201':
        if (CLQs['mjd'][ii] ==  56189):
            ax.scatter(CLQ_REW[ii],  CLQ_FWHM[ii], color='k',    alpha=alpha, marker='o', s=ms*ms_back)
            ax.scatter(CLQ_REW[ii],  CLQ_FWHM[ii], color='cyan', alpha=alpha, marker='o', s=ms)
            ax.errorbar(CLQ_REW[ii], CLQ_FWHM[ii], color='cyan', xerr=log_errREW[ii], yerr=CLQ_errFWHM[ii], fmt='o', linewidth=lw)
            J22_56189  = mlines.Line2D([], [], label='J2228+2201 (56189)', color='cyan',
                        marker="o", markeredgecolor='k', markeredgewidth=1.4, markersize=7,  linestyle='None')
        if (CLQs['mjd'][ii] ==  56960):
            ax.scatter(CLQ_REW[ii],  CLQ_FWHM[ii], color='k',    alpha=alpha, marker='s', s=ms*ms_back)
            ax.scatter(CLQ_REW[ii],  CLQ_FWHM[ii], color='cyan', alpha=alpha, marker='s', s=ms)
            ax.errorbar(CLQ_REW[ii], CLQ_FWHM[ii], color='cyan', xerr=log_errREW[ii], yerr=CLQ_errFWHM[ii], fmt='s', linewidth=lw)
            J22_56960  = mlines.Line2D([], [], label='J2228+2201 (56960)', color='cyan',
                       marker="s", markeredgecolor='k', markeredgewidth=1.4, markersize=7,  linestyle='None')
        if (CLQs['mjd'][ii] ==  58693):
            ax.scatter(CLQ_REW[ii],  CLQ_FWHM[ii], color='k',    alpha=alpha, marker='D', s=ms*ms_back) 
            ax.scatter(CLQ_REW[ii],  CLQ_FWHM[ii], color='cyan', alpha=alpha, marker='D', s=ms)
            ax.errorbar(CLQ_REW[ii], CLQ_FWHM[ii], color='cyan', xerr=log_errREW[ii], yerr=CLQ_errFWHM[ii], fmt='D', linewidth=lw)
            J22_58693  = mlines.Line2D([], [], label='J2228+2201 (58693)', color='cyan',
                        marker="D", markeredgecolor='k', markeredgewidth=1.4, markersize=7,  linestyle='None')

## AXES LIMITS 
xmin =   -20.  ## -0.75 when REW is log
xmax =  375.  ## 3.15 when REW is log
ymin = -1.0
ymax = 19.7  ## when in '000 km s^-1; 22.7 good when legend in top left corner and x is in log

ax.axis([xmin, xmax, ymin, ymax])

## AXES LABELS
#ax.set_xlabel(r'log$_{10}$(REW)',         fontsize=fontsize)
#ax.set_xlabel(r'rest EW / ${ \rm \AA }$', fontsize=fontsize)
ax.set_xlabel(r'Equiv. Width / ${ \rm \AA }$', fontsize=fontsize)
ax.set_ylabel(r"FWHM / '000 km s$^{-1}$",      fontsize=fontsize)

## AXES TICK FORMAT
ax.tick_params(axis='both', which='major', labelsize=labelsize)
ax.tick_params(axis='both', which='minor', labelsize=labelsize)

##
##  L E G E N D S
boss  = mlines.Line2D([], [], label='BOSS quasars', color='skyblue',
                        marker=".", markeredgecolor='k', markeredgewidth=1.1,
                          markersize=7,  linestyle='None')

handles=[boss, J12_53498, J12_58693,
               J16_54553, J16_55832, J16_58583, 
               J22_56189, J22_56960, J22_58693]
leg = ax.legend(loc='upper right',
                fontsize=fontsize/1.6, handles=handles, 
                frameon=True, framealpha=1.0, fancybox=True)


## SAVING THE FIGURE
plt.savefig('CIV_CLQs_REWvsFWHM_temp.png', format='png')
plt.close()
