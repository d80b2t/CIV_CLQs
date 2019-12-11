'''
A ``hexbin'' plot for the SDSS-III BOSS data from the Hamann et
al. (2017) catalog and analysis.  This initially concentrated on
'just' f1450 converted into Lum_1450.  However, matching with the
Kozłowski (2017) catalog means we can do the Baldwin plot with LBol
directly.
'''

import numpy as np
from scipy             import stats
from astropy.io        import fits
from astropy.io        import ascii
from astropy.cosmology import FlatLambdaCDM

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.mlab   as mlab
import matplotlib.lines  as mlines

## Setting up the cosmology
cosmo = FlatLambdaCDM(H0=71, Om0=0.27, Tcmb0=2.725)                                                                                                                          
## Reading in the data
##  
##  T H E    C L Q s
##
path      = '../../../data/CLQ_line_measurements/'
filename  = 'QSFIT_CIV_line_params.dat'
CLQs      = ascii.read(path+filename)

CLQs_logREW     = np.log10(CLQs['EW'])
CLQs_logREWerr  = np.log10(CLQs['EW_err'])
CLQs_MBH        = CLQs['MBH']


##
##    Q S F I T
##
## CONT*__LUM::  1450, 2245, 3000, 4210 or 5100:
## AGN continuum (Section 2.2) νLν luminosities and their uncertainties, in units of 1042 erg s−1
## 
path        = '../../../data/QSO_CIV_catalogs/'
filename    = 'qsfit_1.2.4.fits'
infile      = path+filename
sdata_full  = fits.open(infile)
sdata       = sdata_full[1].data

data_qsofit = sdata[np.where((sdata['BR_CIV_1549__LUM']  > 0.) &
                             (sdata['BR_CIV_1549__FWHM'] > 0.) &
                             (sdata['BR_CIV_1549__EW']   > 0.))] 

EW_QSFit        =           data_qsofit['BR_CIV_1549__EW']
log_EW_QSFit    = np.log10(data_qsofit['BR_CIV_1549__EW'])
## Virial BH estimate:: "just"  (nu.L_nu)^0.5  *  FWHM^2 
logBH_virial    = np.log10((data_qsofit['BR_CIV_1549__LUM']**.5)*(data_qsofit['BR_CIV_1549__FWHM']**2))  


##
##  H A M A N N et al. (2017)   and   K O Z L O W S K I   (2017)  merged BOSS DR12 catalog
##
path      = '../../../data/QSO_CIV_catalogs/'
filename  = 'BOSS_Ham17Koz17_DR12Q.fits'
infile    = path+filename
data_full = fits.open(infile)
tdata     = data_full[1].data

## Somewhat amazingly, there's only object with rew=0 (`next' minimum
## object has rew=0.5006 and there's one object with rew=2710846 (next
## maximum object has rew=3417.)
data = tdata[np.where( (tdata['rew'] > 0.) & (tdata['rew'] < 10000.))]

log_tREW        = np.log10(data['rew'])
log_tMBH        = data['MBH_CVI']
eta_Edd         = data['nEdd']


## 
##   S E T T I N G   U P   T H E    P L O T
##
matplotlib.rc('text', usetex=True)
fig, ax = plt.subplots(figsize=(5, 5), dpi=80, facecolor='w', edgecolor='k')  ## used to be 7,5

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
mincnt    = 3.
gridsize  = 200
color_map = plt.cm.Spectral_r

##
## Setting things all up... 
## Comment in/out which values you’d like
##
## 20K objects::
x = logBH_virial
y = log_EW_QSFit
## 200K objects::
#x = log_tMBH
#y = log_tREW

##
##  P L O T T I N G    T H E    D A T A
##
ax.scatter(     x, y, s=s, alpha=alpha)
hb = plt.hexbin(x, y, bins='log', gridsize=gridsize, cmap=color_map, mincnt=mincnt)

## docs.scipy.org/doc/scipy/reference/generated/scipy.stats.linregress.html
## Calculate a linear least-squares regression for two sets of measurements.
slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
print()
print("Performed linear least-squares regression...")
print("The slope     is: %f" % slope)
print("The intercept is: %f" % intercept)
print("R-squared       : %f" % r_value**2)
print("p-value         : %f" % p_value)
print("Standard errro  : %f" % std_err)
print()

#plt.plot(x, (slope*x + intercept), 'r')
## for the label...
slope_str     = str(np.around(slope, decimals=3))
intercept_str = str(np.around(intercept, decimals=3))
R_sq          = str(np.around((r_value**2), decimals=3))               
        

## AXES LIMITS 
xmin =   7.30   ## Log MBH
xmax =  10.88   ## Log MBH
ymin =  0.4
ymax =  2.95    
ax.axis([xmin, xmax, ymin, ymax])

## AXES LABELS
ax.set_xlabel(r"log$_{10}$(M$_{\rm BH}$/M$_{\odot}$)",         fontsize=fontsize)
ax.set_ylabel(r"log$_{10}$(REW CIV / {\rm \AA })", fontsize=fontsize)
## AXES TICK FORMAT
ax.tick_params(axis='both', which='major', labelsize=labelsize)
ax.tick_params(axis='both', which='minor', labelsize=labelsize)

##
##
##
for ii in range(len(CLQs)):
    ms         =  32.  # 16 is decent
    ms_back    =  72.  # 16:45 is decent
    l_back     =   2
    alpha      =  1.0
    markersize = 6.0
    mrkedwd    = 1.1 ## markeredgewidth
    
    if str(CLQs['Object'][ii]) == 'J1205p3422':
        if (CLQs['MJD'][ii] ==  53498):
            ax.scatter(CLQs_MBH[ii], CLQs_logREW[ii], color='k',       alpha=alpha, marker='o', s=ms_back)
            ax.scatter(CLQs_MBH[ii], CLQs_logREW[ii], color='fuchsia', alpha=alpha, marker='o', s=ms)
            J12_53498 = mlines.Line2D([], [], label=r'J1205+3422 (53498)', color='fuchsia',
                        marker="o", markeredgecolor='k', markeredgewidth=mrkedwd, markersize=markersize,  linestyle='None')
        if (CLQs['MJD'][ii] ==  58538):
            ax.scatter (CLQs_MBH[ii], CLQs_logREW[ii], color='k',       alpha=alpha, marker='s', s=ms_back)
            ax.scatter (CLQs_MBH[ii], CLQs_logREW[ii], color='fuchsia', alpha=alpha, marker='s', s=ms)
            J12_58538 = mlines.Line2D([], [], label=r'J1205+3422 (58538)$^{*}$', color='fuchsia',
                        marker="s", markeredgecolor='k', markeredgewidth=mrkedwd, markersize=markersize,  linestyle='None')
        if (CLQs['MJD'][ii] ==  58693):
            ax.scatter (CLQs_MBH[ii], CLQs_logREW[ii], color='k',       alpha=alpha, marker='D', s=ms_back)
            ax.scatter (CLQs_MBH[ii], CLQs_logREW[ii], color='fuchsia', alpha=alpha, marker='D', s=ms)
            J12_58693 = mlines.Line2D([], [], label=r'J1205+3422 (58693)', color='fuchsia',
                        marker="D", markeredgecolor='k', markeredgewidth=mrkedwd, markersize=markersize,  linestyle='None')
            
    if str(CLQs['Object'][ii]) == 'J1638p2827':
        if (CLQs['MJD'][ii] == 54553):
            ax.scatter(CLQs_MBH[ii], CLQs_logREW[ii], color='k',    alpha=alpha, marker='o', s=ms_back)
            ax.scatter(CLQs_MBH[ii], CLQs_logREW[ii], color='lime', alpha=alpha, marker='o', s=ms)
            J16_54553 = mlines.Line2D([], [], label='J1638+2827 (54553)', color='lime',
                        marker="o", markeredgecolor='k', markeredgewidth=mrkedwd, markersize=markersize,  linestyle='None')
        if (CLQs['MJD'][ii] == 55832):
            ax.scatter(CLQs_MBH[ii], CLQs_logREW[ii], color='k',    alpha=alpha, marker='s', s=ms_back)
            ax.scatter(CLQs_MBH[ii], CLQs_logREW[ii], color='lime', alpha=alpha, marker='s', s=ms)
            J16_55832 = mlines.Line2D([], [], label='J1638+2827 (55832)', color='lime',
                        marker="s", markeredgecolor='k', markeredgewidth=mrkedwd, markersize=markersize,  linestyle='None')
        if (CLQs['MJD'][ii] == 58583):
            ax.scatter(CLQs_MBH[ii], CLQs_logREW[ii], color='k',    alpha=alpha, marker='D', s=ms_back)
            ax.scatter(CLQs_MBH[ii], CLQs_logREW[ii], color='lime', alpha=alpha, marker='D', s=ms)
            J16_58583 = mlines.Line2D([], [], label='J1638+2827 (58583)', color='lime',
                        marker="D", markeredgecolor='k', markeredgewidth=mrkedwd, markersize=markersize,  linestyle='None')
            
    if str(CLQs['Object'][ii]) == 'J2228p2201':
        if (CLQs['MJD'][ii] == 56189):
            ax.scatter(CLQs_MBH[ii], CLQs_logREW[ii], color='k',   alpha=alpha, marker='o', s=ms_back)
            ax.scatter(CLQs_MBH[ii], CLQs_logREW[ii], color='cyan',alpha=alpha, marker='o', s=ms)
            J22_56189 = mlines.Line2D([], [], label=r'J2228+2201 (56189)$^{*}$', color='cyan',
                        marker="o", markeredgecolor='k', markeredgewidth=mrkedwd, markersize=markersize,  linestyle='None')
        if (CLQs['MJD'][ii] == 56960):
            ax.scatter(CLQs_MBH[ii], CLQs_logREW[ii], color='k',       alpha=alpha, marker='s', s=ms_back)
            ax.scatter(CLQs_MBH[ii], CLQs_logREW[ii], color='cyan', alpha=alpha, marker='s', s=ms)
            J22_56960 = mlines.Line2D([], [], label=r'J2228+2201 (56960)', color='cyan',
                        marker="s", markeredgecolor='k', markeredgewidth=mrkedwd, markersize=markersize,  linestyle='None')
        if (CLQs['MJD'][ii] == 58693):
            ax.scatter(CLQs_MBH[ii], CLQs_logREW[ii], color='k',       alpha=alpha, marker='D', s=ms_back)
            ax.scatter(CLQs_MBH[ii], CLQs_logREW[ii], color='cyan', alpha=alpha, marker='D', s=ms)
            J22_58693 = mlines.Line2D([], [], label=r'J2228+2201 (58693)', color='cyan',
                        marker="D", markeredgecolor='k', markeredgewidth=mrkedwd, markersize=markersize,  linestyle='None')

##
##    L E G E N D S
##  
#boss = mlines.Line2D([], [], label='163,437 BOSS quasars', color='skyblue',
#boss = mlines.Line2D([], [], label=r' 20,374 {\sc qsfit} quasars',
#boss = mlines.Line2D([], [], label=r'{\sc qsfit} quasars',
boss  = mlines.Line2D([], [], label='SDSS quasars',                        
                         color='skyblue', marker=".", markeredgecolor='k', markeredgewidth=1.1,
                         markersize=7, linestyle='None')

#bestfit_line  = mlines.Line2D([], [], label=r'$\beta = $'+slope_str, color='r')
bestfit_line  = mlines.Line2D([], [], label=r'$\beta$ = -0.1997', color='r')

#handles=[boss, bestfit_line,
handles=[boss, J12_53498, J12_58538, J12_58693, 
               J16_54553, J16_55832, J16_58583,
               J22_56189, J22_56960, J22_58693]
#             bestfit_line]
    
leg = ax.legend(loc='upper right',
                fontsize=fontsize/2.00, handles=handles,
                ncol=2,
                frameon=True, framealpha=1.0, fancybox=True)

#ax.text(0.31, 0.91, r'$\beta=%.4f$' % (slope, ), fontsize=14, color='r',
#            verticalalignment='bottom', horizontalalignment='right',
#            transform=ax.transAxes,
#            bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})
#            bbox={'facecolor':'white', 'alpha':0.5, 'boxstyle':'round'})

## SAVING THE FIGURE
plt.savefig('CIV_CLQs_MBHvsREW_temp.png', format='png')
plt.close()


