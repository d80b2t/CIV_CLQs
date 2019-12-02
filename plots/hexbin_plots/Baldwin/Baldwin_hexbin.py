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
path      = '../../../data/CLQ_line_measurements/'
filename  = 'QSFIT_CIV_line_params.dat'
CLQs      = ascii.read(path+filename)

CLQs_FWHM       = CLQs['FWHM']    / 1000.
CLQs_FWHMerr    = CLQs['FWHM_err'] / 1000.
CLQs_logREW     = np.log10(CLQs['EW'])
CLQs_logREWerr  = np.log10(CLQs['EW_err'])
CLQs_logContLum = np.log10(CLQs['ContLum']*1e42)
CLQs_logLineLum = np.log10(CLQs['LineLum_err']*1e42)
CLQs_logContLum = CLQs_logContLum - 42.   ## To ``normalize'' to 10^42 erg s^-1
CLQs_logLineLum = CLQs_logLineLum - 42.   ## To ``normalize'' to 10^42 erg s^-1



##  Hamann et al. (2017) and Kozłowski (2017) merged BOSS DR12 catalog
path      = '../../../data/QSO_CIV_catalogs/'
filename  = 'BOSS_Ham17Koz17_DR12Q.fits'
infile    = path+filename
data_full = fits.open(infile)
tdata     = data_full[1].data

## Knocking out a couple of objects with bad REW values;
## No BALs
## qflag; quality flag for the C IV fit: 0 = no problems
## 2.0 < ze < 3.4  to match "Full Sample" in Hamann++17
## and kick out the LBol
data = tdata[np.where( (tdata['rew'] > 0.     )     &
                       (tdata['rew'] < 10000. )     & 
                       (tdata['bal_flag_vi'] == 0 ) &
                       (tdata['qflag']       == 0 ) & 
                       (tdata['z_dr12'] > 2.0)      &
                       (tdata['z_dr12'] < 3.4)      & 
                       (tdata['Lbol']   > 0.0))]

## Setting up variable names
f1450_Ham17     = data['f1450']
alpha_civ_Ham17 = data['alpha_civ']
REW_Ham17       = data['rew']        
FWHM_Ham17      = data['fwhm']
z_dr12          = data['z_dr12']
log_LBol_Koz17  = data['LBol']
log_L1350_Koz17 = data['L1350']


##    Q S F I T 
## CONT*__LUM::  1450, 2245, 3000, 4210 or 5100:
## AGN continuum (Section 2.2) νLν luminosities and their uncertainties, in units of 1042 erg s−1
## 
path        = '/cos_pc19a_npr/programs/quasars/CIV_CLQs/data/QSO_CIV_catalogs/'
filename    = 'qsfit_1.2.4.fits'
infile      = path+filename
sdata_full  = fits.open(infile)
sdata       = sdata_full[1].data

data_qsofit = sdata[np.where( (sdata['BR_CIV_1549__LUM'] > 0.   )     &
                             (sdata['BR_CIV_1549__EW']   > 0.   ))    ] 

log_L1450_QSFit = np.log10((data_qsofit['BR_CIV_1549__LUM'] * 1e42))
log_LCont_QSFit = np.log10((data_qsofit['CONT1__LUM']       * 1e42))
Cont1__WAVE     =           data_qsofit['CONT1__WAVE']
EW_QSFit        =           data_qsofit['BR_CIV_1549__EW']
log_EW_QSFit    = np.log10( data_qsofit['BR_CIV_1549__EW'])



## Putting te Equiv. Widths into log and
## the FWHMs into '000s of kms^-1.
log_REW_Ham17   = np.log10(REW_Ham17)
FWHM_kilo_Ham17 = FWHM_Ham17/1000.

## The f1450 values are fluxes in the uncorrected BOSS spectrum at 1450 Å rest (10−17 ergs s−1 cm−2 Å−1 )
## and straigt from Hamannn, Zakamska, Ross et al. (2017; Appendix A) but I've never been able to reproduce
## e.g. Figure 3 from that/our paper... (!!!) 
f_lam   = f1450_Ham17*((1450/1450)**alpha_civ_Ham17)     

## L = 4 * pi * (D_L^2) * F
##   where L is in W and F is in W/m^2
##   3.086e+22 meters in a megaparsec
## Set-up the Luminosity Distance 
DL      = cosmo.luminosity_distance(z_dr12)             
DL_inM  = DL.value * 3.08567758128e22
DL_incm = DL.value * 3.08567758128e24

## DL in cm; f in 10−17 ergs s−1 cm−2 Å−1
Lum1450     = (4 * np.pi * (DL_incm*DL_incm) * 1e-17 * f_lam) * 1450.
log_Lum1450 = np.log10(Lum1450)


 
##   S E T T I N G   U P   T H E    P L O T  
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
## Setting things all up... :-)
##
##  x - a x i s
##
## L_1450 values taking f1450 from Hamann et al.
## L_1350          from Kozłowski (2017)
## L_Bol           from Kozłowski (2017)
## log_L1450_QSFit from QSFit (v1.2.4; Calderone et al. 2017) 
#
#x = log_LBol_Koz17
#x = log_L1350_Koz17
#x = log_Lum1450
#x = log_L1450_QSFit
x = log_LCont_QSFit - 42.  ## -42.0 for more direct comparison to Fig 15 of Calderone+2017.

##
##  y  - a x i s 
##
#y = log_REW
y = log_EW_QSFit 


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
## "Fitting", by eye, the the CLQs points
plt.plot(x, ((-.38*x) + (intercept+1.1)), 'r', linestyle=':', zorder=0)

## for the label...
slope_str     = str(np.around(slope, decimals=3))
intercept_str = str(np.around(intercept, decimals=3))
R_sq          = str(np.around((r_value**2), decimals=3))               
        

## AXES LIMITS 
xmin =  2.30   
xmax =  5.55    ## log_LCont_QSFit.max() ## log_LBol.max()   ## can be e.g. 48.4 for log_LBol
ymin =  0.4
ymax =  2.95    ## when in '000 km s^-1;  3.4 when the Legend handles has the # and slope line
ax.axis([xmin, xmax, ymin, ymax])

## AXES LABELS
#ax.set_xlabel(r"log$_{10}[\lambda$L$_{\lambda}$(1450{\rm \AA })] ergs/s ",  fontsize=fontsize)
ax.set_xlabel(r"log$_{10}[\nu    $L$_{\nu}    $(1450{\rm \AA })] / 10$^{42}$ erg s$^{-1}$",  fontsize=fontsize)
#ax.set_xlabel(r"log$_{10}$(L$_{\rm 1350}$) [ergs/s] ",  fontsize=fontsize)
#ax.set_xlabel(r"log$_{10}$(L$_{\rm Bol}$) [ergs/s] ",  fontsize=fontsize)
ax.set_ylabel(r"log$_{10}$(REW CIV / {\rm \AA })",     fontsize=fontsize)

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
            #ax.errorbar(CLQs_logContLum[ii], CLQs_logREW[ii], color='k',      xerr=CLQs_logREWerr[ii], yerr=CLQs_FWHMerr[ii], fmt='o', linewidth=lw*l_back)
            #ax.errorbar(CLQs_logContLum[ii], CLQs_logREW[ii], color='fuchsia',xerr=CLQs_logREWerr[ii], yerr=CLQs_FWHMerr[ii], fmt='o', linewidth=lw)
            ax.scatter(CLQs_logContLum[ii], CLQs_logREW[ii], color='k',       alpha=alpha, marker='o', s=ms_back)
            ax.scatter(CLQs_logContLum[ii], CLQs_logREW[ii], color='fuchsia', alpha=alpha, marker='o', s=ms)
            J12_53498 = mlines.Line2D([], [], label='J1205+3422 (53498)', color='fuchsia',
                        marker="o", markeredgecolor='k', markeredgewidth=mrkedwd, markersize=markersize,  linestyle='None')
        if (CLQs['MJD'][ii] ==  58538):
            #ax.errorbar(CLQs_logContLum[ii], CLQs_logREW[ii], color='k',      xerr=CLQs_logREWerr[ii], yerr=CLQs_FWHMerr[ii], fmt='s', linewidth=lw*l_back)
            #ax.errorbar(CLQs_logContLum[ii], CLQs_logREW[ii], color='fuchsia',xerr=CLQs_logREWerr[ii], yerr=CLQs_FWHMerr[ii], fmt='s', linewidth=lw)
            ax.scatter (CLQs_logContLum[ii], CLQs_logREW[ii], color='k',       alpha=alpha, marker='s', s=ms_back)
            ax.scatter (CLQs_logContLum[ii], CLQs_logREW[ii], color='fuchsia', alpha=alpha, marker='s', s=ms)
            J12_58538 = mlines.Line2D([], [], label='J1205+3422 (58538)', color='fuchsia',
                        marker="s", markeredgecolor='k', markeredgewidth=mrkedwd, markersize=markersize,  linestyle='None')
        if (CLQs['MJD'][ii] ==  58693):
            #ax.errorbar(CLQs_logContLum[ii], CLQs_logREW[ii], color='k',      xerr=CLQs_logREWerr[ii], yerr=CLQs_FWHMerr[ii], fmt='D', linewidth=lw*l_back)
            #ax.errorbar(CLQs_logContLum[ii], CLQs_logREW[ii], color='fuchsia',xerr=CLQs_logREWerr[ii], yerr=CLQs_FWHMerr[ii], fmt='D', linewidth=lw)
            ax.scatter (CLQs_logContLum[ii], CLQs_logREW[ii], color='k',       alpha=alpha, marker='D', s=ms_back)
            ax.scatter (CLQs_logContLum[ii], CLQs_logREW[ii], color='fuchsia', alpha=alpha, marker='D', s=ms)
            J12_58693 = mlines.Line2D([], [], label='J1205+3422 (58693)', color='fuchsia',
                        marker="D", markeredgecolor='k', markeredgewidth=mrkedwd, markersize=markersize,  linestyle='None')
            
    if str(CLQs['Object'][ii]) == 'J1638p2827':
        if (CLQs['MJD'][ii] == 54553):
            ax.scatter(CLQs_logContLum[ii], CLQs_logREW[ii], color='k',    alpha=alpha, marker='o', s=ms_back)
            ax.scatter(CLQs_logContLum[ii], CLQs_logREW[ii], color='lime', alpha=alpha, marker='o', s=ms)
            J16_54553 = mlines.Line2D([], [], label='J1638+2827 (54553)', color='lime',
                        marker="o", markeredgecolor='k', markeredgewidth=mrkedwd, markersize=markersize,  linestyle='None')
        if (CLQs['MJD'][ii] == 55832):
            ax.scatter(CLQs_logContLum[ii], CLQs_logREW[ii], color='k',    alpha=alpha, marker='s', s=ms_back)
            ax.scatter(CLQs_logContLum[ii], CLQs_logREW[ii], color='lime', alpha=alpha, marker='s', s=ms)
            J16_55832 = mlines.Line2D([], [], label='J1638+2827 (55832)', color='lime',
                        marker="s", markeredgecolor='k', markeredgewidth=mrkedwd, markersize=markersize,  linestyle='None')
        if (CLQs['MJD'][ii] == 58583):
            ax.scatter(CLQs_logContLum[ii], CLQs_logREW[ii], color='k',    alpha=alpha, marker='D', s=ms_back)
            ax.scatter(CLQs_logContLum[ii], CLQs_logREW[ii], color='lime', alpha=alpha, marker='D', s=ms)
            J16_58583 = mlines.Line2D([], [], label='J1638+2827 (58583)', color='lime',
                        marker="D", markeredgecolor='k', markeredgewidth=mrkedwd, markersize=markersize,  linestyle='None')
            
    if str(CLQs['Object'][ii]) == 'J2228p2201':
        if (CLQs['MJD'][ii] == 56189):
            ax.scatter(CLQs_logContLum[ii], CLQs_logREW[ii], color='k',   alpha=alpha, marker='o', s=ms_back)
            ax.scatter(CLQs_logContLum[ii], CLQs_logREW[ii], color='cyan',alpha=alpha, marker='o', s=ms)
            J22_56189 = mlines.Line2D([], [], label='J2228+2201 (56189)', color='cyan',
                        marker="o", markeredgecolor='k', markeredgewidth=mrkedwd, markersize=markersize,  linestyle='None')
        if (CLQs['MJD'][ii] == 56960):
            ax.scatter(CLQs_logContLum[ii], CLQs_logREW[ii], color='k',       alpha=alpha, marker='s', s=ms_back)
            ax.scatter(CLQs_logContLum[ii], CLQs_logREW[ii], color='cyan', alpha=alpha, marker='s', s=ms)
            J22_56960 = mlines.Line2D([], [], label='J1205+3422 (56960)', color='cyan',
                        marker="s", markeredgecolor='k', markeredgewidth=mrkedwd, markersize=markersize,  linestyle='None')
        if (CLQs['MJD'][ii] == 58693):
            ax.scatter(CLQs_logContLum[ii], CLQs_logREW[ii], color='k',       alpha=alpha, marker='D', s=ms_back)
            ax.scatter(CLQs_logContLum[ii], CLQs_logREW[ii], color='cyan', alpha=alpha, marker='D', s=ms)
            J22_58693 = mlines.Line2D([], [], label='J1205+3422 (58693)', color='cyan',
                        marker="D", markeredgecolor='k', markeredgewidth=mrkedwd, markersize=markersize,  linestyle='None')

##
##    L E G E N D S
##  
#boss  = mlines.Line2D([], [], label='163,437 BOSS quasars', color='skyblue',
#boss = mlines.Line2D([], [], label=r' 20,374 {\sc qsfit} quasars',
#boss = mlines.Line2D([], [], label=r'{\sc qsfit} quasars',
boss = mlines.Line2D([], [], label='SDSS quasars',                        
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
plt.savefig('CIV_CLQs_Baldwin_temp.png', format='png')
plt.close()





##
## Just keen to see what the distribution of CONT1__WAVE -- from QSFit
## and Calderone et al. (2017) -- is.  This is the value of where the
## continuum luminosities and slopes havw been computed. There are
## five such quantities whose values are equally spaced in the
## logarithmic wavelength range [λmin , λmax ];
## 
##   S E T T I N G   U P   T H E    P L O T
##
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
mincnt    = 3.
gridsize  = 200
color_map = plt.cm.Spectral_r

# the histogram of the data
n, bins, patches = plt.hist(Cont1__WAVE, 100, facecolor='green', alpha=0.75)
l = plt.plot(bins)

plt.xlabel('CONT1\_\_WAVELENGTH')
plt.ylabel('No. of QSFit Quasars')
#plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=100,\ \sigma=15$')
plt.axis([bins.min()-10, bins.max()+10, 0,  n.max()*1.1])
plt.grid(True, linestyle='--', zorder=0)

plt.savefig('CIV_CLQs_Baldwin_WAVE_temp.png', format='png')
plt.close()
