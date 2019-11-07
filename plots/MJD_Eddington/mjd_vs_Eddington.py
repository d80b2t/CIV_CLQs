'''

Plot of MJD of Quasar observation, and Eddington ratio, 
e.g. from Shen et al. (2011) or 

'''
import math
import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.lines  as mlines

from matplotlib import colors as mcolors
from astropy.io import ascii
from astropy.io import fits

##
##  S h e n   et al.  (2 0 1 1)   D R 7 Q  
##
path = '../../data/QSO_CIV_catalogs/'
filename ='Shen_dr7_bh_May_2010.fits'
infile = path+filename

## Okay, what we really want to do ;-)
#data_table = pyfits.getdata(infile)
data_in    = fits.open(infile)
data_table = data_in[1].data

## Quick check on the format/dimensions of the FITS table file...
print(type(data_table), '\n')
print('The number of rows of is.... ',     data_table.shape,  '\n')
print('The number of columns is...  ', len(data_table.names), '\n\n')


## Making some senible selections and cuts 
mjd_full    = data_table.field('MJD')
logtau_full = data_table.field('LOGEDD_RATIO')

dr7q        = data_table[np.where(data_table['LOGEDD_RATIO'] > -9.0)]
dr7q_mjd    = dr7q['MJD']
dr7q_logtau = dr7q['LOGEDD_RATIO']
dr7q_tau    = (10**dr7q['LOGEDD_RATIO'])*100

    
## The data for our CIV CLQs
path   = '../../data/CLQ_line_measurements/'
infile = 'MJD_Eddington.dat'
CLQs   = ascii.read(path+infile)


## Setting up the plot
fig, ax = plt.subplots(figsize=(10.0, 6.0), dpi=80, facecolor='w', edgecolor='k')   # was 14.0, 8.0

## Adjusting the Whitespace for the plots
left   = 0.08   # the left side of the subplots of the figure
right  = 0.96   # the right side of the subplots of the figure
bottom = 0.12   # the bottom of the subplots of the figure
top    = 0.96   # the top of the subplots of the figure
wspace = 0.22   # the amount of width reserved for blank space between subplots
hspace = 0.06   # the amount of height reserved for white space between subplots
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)

## Some NPR defaults
lw              = 5.0
ls              = 'solid'
ms              = 80.
ms_large        = ms*3.
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


# COLOR MAP
color_map = plt.cm.get_cmap('Greys')

## For the hexbinning
gridsize  = 30  #use log bins
mincnt    = .1

##  T H E     D R 7 Q     q u a s a r s 
ax.hexbin(dr7q_mjd, dr7q_logtau,         bins='log',  gridsize=gridsize, cmap=color_map, mincnt=mincnt)




##  T H E     C I V     C L Q s
for ii in range(len(CLQs)):
    
    if str(CLQs['Object'][ii]) == 'J1205+3422':
        if (CLQs['MJD'][ii] ==  53498):
            ax.scatter(CLQs['MJD'][ii], CLQs['eta_MgII'][ii], color='k',       alpha=alpha, marker='o', s=ms_large*1.6)
            ax.scatter(CLQs['MJD'][ii], CLQs['eta_MgII'][ii], color='fuchsia', alpha=alpha, marker='o', s=ms_large)
            J12_53498  = mlines.Line2D([], [], label='J1205+3422 (53498)', color='fuchsia',
                       marker="o", markeredgecolor='k', markeredgewidth=1.4, markersize=7,  linestyle='None')
            
        if (CLQs['MJD'][ii] ==  58693):
            ax.scatter(CLQs['MJD'][ii], CLQs['eta_MgII'][ii], color='k',       alpha=alpha, marker='s', s=ms_large*1.6)
            ax.scatter(CLQs['MJD'][ii], CLQs['eta_MgII'][ii], color='fuchsia', alpha=alpha, marker='s', s=ms_large)
            J12_58693 = mlines.Line2D([], [], label='J1205+3422 (58693)', color='fuchsia',
                       marker="s", markeredgecolor='k', markeredgewidth=1.4, markersize=7,  linestyle='None')
            
    if str(CLQs['Object'][ii]) == 'J1638+2827':
        if (CLQs['MJD'][ii] ==  54553):
            ax.scatter(CLQs['MJD'][ii], CLQs['eta_MgII'][ii], color='k',    alpha=alpha, marker='o', s=ms_large*1.6)
            ax.scatter(CLQs['MJD'][ii], CLQs['eta_MgII'][ii], color='lime', alpha=alpha, marker='o', s=ms_large)

            J16_54553 = mlines.Line2D([], [], label='J1638+2827 (54553)', color='lime',
                       marker="o", markeredgecolor='k', markeredgewidth=1.4, markersize=7,  linestyle='None')
        if (CLQs['MJD'][ii] ==  55832): 
            ax.scatter(CLQs['MJD'][ii], CLQs['eta_MgII'][ii], color='k',    alpha=alpha, marker='s', s=ms_large*1.6)
            ax.scatter(CLQs['MJD'][ii], CLQs['eta_MgII'][ii], color='lime', alpha=alpha, marker='s', s=ms_large)
            J16_55832 = mlines.Line2D([], [], label='J1638+2827 (55832)', color='lime',
                       marker="s", markeredgecolor='k', markeredgewidth=1.4, markersize=7,  linestyle='None')
        if (CLQs['MJD'][ii] ==  58583): 
            ax.scatter(CLQs['MJD'][ii], CLQs['eta_MgII'][ii], color='k',    alpha=alpha, marker='D', s=ms_large*1.6)
            ax.scatter(CLQs['MJD'][ii], CLQs['eta_MgII'][ii], color='lime', alpha=alpha, marker='D', s=ms_large)
            J16_58583 = mlines.Line2D([], [], label='J1638+2827 (58583)', color='lime',
                       marker="D", markeredgecolor='k', markeredgewidth=1.4, markersize=7,  linestyle='None')
            
    if str(CLQs['Object'][ii]) == 'J2228+2201':
        if (CLQs['MJD'][ii] ==  56189):
            ax.scatter(CLQs['MJD'][ii], CLQs['eta_MgII'][ii], color='k',    alpha=alpha, marker='o', s=ms_large*1.6)
            ax.scatter(CLQs['MJD'][ii], CLQs['eta_MgII'][ii], color='cyan', alpha=alpha, marker='o', s=ms_large)
            J22_56189  = mlines.Line2D([], [], label='J2228+2201 (56189)', color='cyan',
                        marker="o", markeredgecolor='k', markeredgewidth=1.4, markersize=7,  linestyle='None')
        if (CLQs['MJD'][ii] ==  56960):
            ax.scatter(CLQs['MJD'][ii], CLQs['eta_MgII'][ii], color='k',    alpha=alpha, marker='s', s=ms_large*1.6)
            ax.scatter(CLQs['MJD'][ii], CLQs['eta_MgII'][ii], color='cyan', alpha=alpha, marker='s', s=ms_large)
            J22_56960  = mlines.Line2D([], [], label='J2228+2201 (56960)', color='cyan',
                       marker="s", markeredgecolor='k', markeredgewidth=1.4, markersize=7,  linestyle='None')
        if (CLQs['MJD'][ii] ==  58693):
            ax.scatter(CLQs['MJD'][ii], CLQs['eta_MgII'][ii], color='k',    alpha=alpha, marker='D', s=ms_large*1.6) 
            ax.scatter(CLQs['MJD'][ii], CLQs['eta_MgII'][ii], color='cyan', alpha=alpha, marker='D', s=ms_large)
            J22_58693  = mlines.Line2D([], [], label='J2228+2201 (58693)', color='cyan',
                        marker="D", markeredgecolor='k', markeredgewidth=1.4, markersize=7,  linestyle='None')

    if str(CLQs['Object'][ii]) == 'M87':
            ax.scatter(CLQs['MJD'][ii], CLQs['eta_MgII'][ii], color='k',          alpha=alpha, marker='H', s=ms_large*1.6) 
            ax.scatter(CLQs['MJD'][ii], CLQs['eta_MgII'][ii], color='darkorange', alpha=alpha, marker='H', s=ms_large)
            M87_leg  = mlines.Line2D([], [], label='M87 (57854)', color='darkorange',
                        marker="H", markeredgecolor='k', markeredgewidth=1.4, markersize=7,  linestyle='None')

#ax.scatter(CLQs['MJD'], CLQs['eta_MgII'], marker='o', s=ms_large*1.6, color='black')
#ax.scatter(CLQs['MJD'], CLQs['eta_MgII'], marker='o', s=ms_large,     color='dodgerblue')

#ax.scatter(risers_data['mjd'], np.log10(risers_data['EDD_RATIO']/100),  marker='o', s=ms_large*1.6, color='black')
#ax.scatter(risers_data['mjd'], np.log10(risers_data['EDD_RATIO']/100),  marker='o', s=ms_large, color='red')



## Tidy up the figure
xmin     = 51400      ## 51100
xmax     = 58664.+360
ymin     =   -5.1     ## -3.4 
ymax     =   1.4      ##  1.2
ymin_log =   0.1      ##  0.1
ymax_log = 100.0    

ax.set_xlim((xmin, xmax))
ax.set_ylim((ymin, ymax))
#ax.set_yscale('log')
#ax.set_ylim((ymin_log, ymax_log))

ax.tick_params('x', direction='in')
ax.tick_params('y', direction='in')
#ax.minorticks_on('x', direction='in')
#ay.minorticks_on()
ax.tick_params('x', direction='in', which='major', bottom='True', top='True', left='True', right='True', labelsize=fontsize/1.2)
ax.tick_params('x', direction='in', which='minor', bottom='True', top='True', left='True', right='True', labelsize=fontsize/1.2)
ax.tick_params('y', direction='in', which='both',  bottom='True', top='True', left='True', right='True', labelsize=fontsize)
ax.minorticks_on() 

## ALLWISE   timespan
NodaDone         = 0.02
NodaDone_range   = 1.5
log_NodaDone_min = np.log10(NodaDone / NodaDone_range)
log_NodaDone_max = np.log10(NodaDone * NodaDone_range)

ax.axhspan(log_NodaDone_min, log_NodaDone_max, alpha=0.6, color='red')
ax.text( log_NodaDone_min,   0.8, 'Noda-Done',   style='italic',    fontsize=fontsize/1.2, rotation=270)

ax.set_xlabel('MJD',                       fontsize=fontsize)
ax.set_ylabel(r'log$_{10}$ Eddington Ratio', fontsize=fontsize)



##
##  L E G E N D S
##
boss  = mlines.Line2D([], [], label='DR7 quasars', color='grey',
                          marker=".", markeredgecolor='k', markeredgewidth=1.2,
                          markersize=8,  linestyle='None')

handles=[boss, J12_53498, J12_58693,
               J16_54553, J16_55832, J16_58583, 
               J22_56189, J22_56960, J22_58693, M87_leg]
leg = ax.legend(loc='lower left',
                fontsize=fontsize/1.4, handles=handles, 
                frameon=True, framealpha=1.0, fancybox=True, ncol=2)





##plt.show()
plt.savefig('MJD_vs_Eddington_temp.png',format='png')
plt.close(fig)
