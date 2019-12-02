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

## Reading in the data
##     C I V     C  L Q s
path     = '../../data/CLQ_line_measurements/'
infile   = 'MJD_Eddington.dat'
CLQs     = ascii.read(path+infile)
CLQs_MJD = CLQs['MJD']
CLQs_MBH = CLQs['MBH']
CLQs_eta = CLQs['eta']

##
##  D R 7 Q    S h e n   et al.  (2011)     
##
path       = '../../data/QSO_CIV_catalogs/'
filename   = 'Shen_dr7_bh_May_2010.fits'
dr7q_in    = fits.open(path+filename)
dr7q_table = dr7q_in[1].data

## Quick check on the format/dimensions of the FITS table file...
print(type(dr7q_table), '\n')
print('The number of rows of is.... ',     dr7q_table.shape,  '\n')
print('The number of columns is...  ', len(dr7q_table.names), '\n\n')

## Making some senible selections and cuts 
mjd_full    = dr7q_table.field('MJD')
logtau_full = dr7q_table.field('LOGEDD_RATIO')
dr7q        = dr7q_table[np.where(dr7q_table['LOGEDD_RATIO'] > -9.0)]
dr7q_mjd    = dr7q['MJD']
dr7q_logtau = dr7q['LOGEDD_RATIO']
dr7q_tau    = (10**dr7q['LOGEDD_RATIO'])*100

##
##  D R 1 2 Q    Kozlowski (2017)  and  Hamann et al. (2017)
##
filename   = 'BOSS_Ham17Koz17_DR12Q.fits'
dr12q_in    = fits.open(path+filename)
dr12q_table = dr12q_in[1].data

## Quick check on the format/dimensions of the FITS table file...
print(type(dr12q_table), '\n')
print('The number of rows of is.... ',     dr12q_table.shape,  '\n')
print('The number of columns is...  ', len(dr12q_table.names), '\n\n')

## Making some senible selections and cuts 
dr12q        = dr12q_table[np.where(dr12q_table['nEdd'] > -9.0)]
dr12q_mjd    = dr12q['MJD_x']
dr12q_etaEdd = dr12q['nEdd']

    


## Setting up the plot
fig, ax = plt.subplots(figsize=(10.0, 6.0), dpi=80, facecolor='w', edgecolor='k')   # was 14.0, 8.0

## Adjusting the Whitespace for the plots
left   = 0.10   # the left side of the subplots of the figure
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
ax.hexbin( dr7q_mjd,  dr7q_logtau, bins='log',  gridsize=gridsize, cmap=color_map, mincnt=mincnt)
ax.hexbin(dr12q_mjd, dr12q_etaEdd, bins='log',  gridsize=gridsize, cmap=color_map, mincnt=mincnt)


##  T H E     C I V     C L Q s
for ii in range(len(CLQs)):
    
    if str(CLQs['Object'][ii]) == 'J1205+3422':
        if (CLQs_MJD[ii] ==  53498):
            ax.scatter(CLQs_MJD[ii], CLQs_eta[ii], color='k',       alpha=alpha, marker='o', s=ms_large*1.6)
            ax.scatter(CLQs_MJD[ii], CLQs_eta[ii], color='fuchsia', alpha=alpha, marker='o', s=ms_large)
            J12_53498  = mlines.Line2D([], [], label='J1205+3422 (53498)', color='fuchsia',
                       marker="o", markeredgecolor='k', markeredgewidth=1.4, markersize=7,  linestyle='None')
        if (CLQs_MJD[ii] ==  58538):
            ax.scatter(CLQs_MJD[ii], CLQs_eta[ii], color='k',       alpha=alpha, marker='s', s=ms_large*1.6)
            ax.scatter(CLQs_MJD[ii], CLQs_eta[ii], color='fuchsia', alpha=alpha, marker='s', s=ms_large)
            J12_58538  = mlines.Line2D([], [], label='J1205+3422 (58538)', color='fuchsia',
                       marker="s", markeredgecolor='k', markeredgewidth=1.4, markersize=7,  linestyle='None')
        if (CLQs_MJD[ii] ==  58693):
            ax.scatter(CLQs_MJD[ii], CLQs_eta[ii], color='k',       alpha=alpha, marker='D', s=ms_large*1.6, zorder=10)
            ax.scatter(CLQs_MJD[ii], CLQs_eta[ii], color='fuchsia', alpha=alpha, marker='D', s=ms_large, zorder=10)
            J12_58693  = mlines.Line2D([], [], label='J1205+3422 (58693)', color='fuchsia',
                       marker="D", markeredgecolor='k', markeredgewidth=1.4, markersize=7,  linestyle='None')
            
    if str(CLQs['Object'][ii]) == 'J1638+2827':
        if (CLQs_MJD[ii] ==  54553):
            ax.scatter(CLQs_MJD[ii], CLQs_eta[ii], color='k',    alpha=alpha, marker='o', s=ms_large*1.6)
            ax.scatter(CLQs_MJD[ii], CLQs_eta[ii], color='lime', alpha=alpha, marker='o', s=ms_large)
            J16_54553 = mlines.Line2D([], [], label='J1638+2827 (54553)', color='lime',
                       marker="o", markeredgecolor='k', markeredgewidth=1.4, markersize=7,  linestyle='None')
        if (CLQs_MJD[ii] ==  55832): 
            ax.scatter(CLQs_MJD[ii], CLQs_eta[ii], color='k',    alpha=alpha, marker='s', s=ms_large*1.6)
            ax.scatter(CLQs_MJD[ii], CLQs_eta[ii], color='lime', alpha=alpha, marker='s', s=ms_large)
            J16_55832 = mlines.Line2D([], [], label='J1638+2827 (55832)', color='lime',
                       marker="s", markeredgecolor='k', markeredgewidth=1.4, markersize=7,  linestyle='None')
        if (CLQs_MJD[ii] ==  58583): 
            ax.scatter(CLQs_MJD[ii], CLQs_eta[ii], color='k',    alpha=alpha, marker='D', s=ms_large*1.6)
            ax.scatter(CLQs_MJD[ii], CLQs_eta[ii], color='lime', alpha=alpha, marker='D', s=ms_large)
            J16_58583 = mlines.Line2D([], [], label='J1638+2827 (58583)', color='lime',
                       marker="D", markeredgecolor='k', markeredgewidth=1.4, markersize=7,  linestyle='None')
            
    if str(CLQs['Object'][ii]) == 'J2228+2201':
        if (CLQs_MJD[ii] ==  56189):
            ax.scatter(CLQs_MJD[ii], CLQs_eta[ii], color='k',    alpha=alpha, marker='o', s=ms_large*1.6)
            ax.scatter(CLQs_MJD[ii], CLQs_eta[ii], color='cyan', alpha=alpha, marker='o', s=ms_large)
            J22_56189  = mlines.Line2D([], [], label='J2228+2201 (56189)', color='cyan',
                        marker="o", markeredgecolor='k', markeredgewidth=1.4, markersize=7,  linestyle='None')
        if (CLQs_MJD[ii] ==  56960):
            ax.scatter(CLQs_MJD[ii], CLQs_eta[ii], color='k',    alpha=alpha, marker='s', s=ms_large*1.6)
            ax.scatter(CLQs_MJD[ii], CLQs_eta[ii], color='cyan', alpha=alpha, marker='s', s=ms_large)
            J22_56960  = mlines.Line2D([], [], label='J2228+2201 (56960)', color='cyan',
                       marker="s", markeredgecolor='k', markeredgewidth=1.4, markersize=7,  linestyle='None')
        if (CLQs_MJD[ii] ==  58693):
            ax.scatter(CLQs_MJD[ii], CLQs_eta[ii], color='k',    alpha=alpha, marker='D', s=ms_large*1.6) 
            ax.scatter(CLQs_MJD[ii], CLQs_eta[ii], color='cyan', alpha=alpha, marker='D', s=ms_large)
            J22_58693  = mlines.Line2D([], [], label='J2228+2201 (58693)', color='cyan',
                        marker="D", markeredgecolor='k', markeredgewidth=1.4, markersize=7,  linestyle='None')

## Tidy up the figure
xmin     = 51400      ## 51100
xmax     = 58664.+360
ymin     =   -2.8     ##  used to be -5.1
ymax     =   1.4      ##  1.2
ymin_log =   0.1      ##  0.1
ymax_log = 100.0    

ax.set_xlim((xmin, xmax))
ax.set_ylim((ymin, ymax))

ax.tick_params('x', direction='in')
ax.tick_params('y', direction='in')
ax.tick_params('x', direction='in', which='major', bottom='True', top='True', left='True', right='True', labelsize=fontsize/1.2)
ax.tick_params('x', direction='in', which='minor', bottom='True', top='True', left='True', right='True', labelsize=fontsize/1.2)
ax.tick_params('y', direction='in', which='both',  bottom='True', top='True', left='True', right='True', labelsize=fontsize)
ax.minorticks_on() 

## Eddington  ratio ranges
text_min = 51800
NodaDone         = 0.02
NodaDone_range   = 1.35
log_NodaDone_min = np.log10(NodaDone / NodaDone_range)
log_NodaDone_max = np.log10(NodaDone * NodaDone_range)
ax.axhspan(log_NodaDone_min, log_NodaDone_max, alpha=0.6, color='red', zorder=1)
#ax.text(log_NodaDone_min,   0.8, 'Noda-Done',   style='italic',                   fontsize=fontsize/1.2, rotation=270)


## From Ruan (2019a) 
log_HighSoft   =  -1  
log_LowHard_one = -1.8
LowHard_one_range = 0.1
log_LowHard_min = log_LowHard_one + LowHard_one_range
log_LowHard_max = log_LowHard_one - LowHard_one_range
log_LowHard_two = -2

#ax.hlines(  log_HighSoft,    xmin, xmax, color='c', linestyles='--', linewidth=linewidth/2.6)
##ax.axhspan(log_LowHard_min, log_LowHard_max, alpha=0.6, color='c')
#ax.hlines(  log_LowHard_two, xmin, xmax, color='c', linestyles='--', linewidth=linewidth/2.6)
ax.hlines(  log_HighSoft,    xmin, xmax, color='darkturquoise', linestyles='--', linewidth=linewidth/2.6)
##ax.axhspan(log_LowHard_min, log_LowHard_max, alpha=0.6, color='c')
ax.hlines(  log_LowHard_two, xmin, xmax, color='darkcyan', linestyles='--', linewidth=linewidth/2.6)

ax.arrow(51700,   -0.90, 0.0,  0.4, width=15., head_width=80., head_length=0.15, color='darkturquoise' )
ax.text(text_min, -0.95, '"High/Soft"', style='italic', weight='bold', fontsize=fontsize, color='darkturquoise')
ax.text(text_min, -1.75, 'Noda-Done',   style='italic', weight='bold', fontsize=fontsize)
ax.text(text_min, -2.18, '"Low/Hard"',  style='italic', weight='bold', fontsize=fontsize, color='darkcyan')
ax.arrow(51700,   -2.05, 0.0, -0.4, width=20., head_width=80., head_length=0.15, color='darkcyan' )

ax.set_xlabel('MJD',                       fontsize=fontsize)
ax.set_ylabel(r'log$_{10}$ Eddington Ratio', fontsize=fontsize)


##
##  L E G E N D S
##
boss  = mlines.Line2D([], [], label='SDSS/BOSS quasars', color='grey',
                          marker=".", markeredgecolor='k', markeredgewidth=1.2,
                          markersize=8,  linestyle='None')

handles=[J12_53498, J12_58538, J12_58693,
         J16_54553, J16_55832, J16_58583, 
         J22_56189, J22_56960, J22_58693]

leg = ax.legend(loc='upper right',
                fontsize=fontsize/1.6, handles=handles, 
                frameon=True, framealpha=1.0, fancybox=True, ncol=3)
#handles=[boss]
#leg = ax.legend(loc='upper left',
#                fontsize=fontsize/1.6, handles=handles, 
#                frameon=False, framealpha=1.0, fancybox=True, ncol=1)

##plt.show()
plt.savefig('MJD_vs_Eddington_temp.png',format='png')
plt.close(fig)

