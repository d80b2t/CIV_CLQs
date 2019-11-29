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
##  T H E   C L Qs   
##
path      = '../../../data/CLQ_line_measurements/'
#filename  = 'Baldwin_CIV.dat'
filename  = 'QSFIT_CIV_line_params.dat'
CLQs      = ascii.read(path+filename)

CLQs_FWHM       = CLQs['FWHM']     / 1000.
CLQs_errFWHM    = CLQs['FWHM_err'] / 1000.
CLQs_logREW     = np.log10(CLQs['EW'])
CLQs_logerrREW  = np.log10(CLQs['EW_err'])
CLQs_loglineLum = np.log10(CLQs['LineLum']*1e42)
CLQs_logContLum = np.log10(CLQs['ContLum']*1e42)
CLQs_MBH        = CLQs['MBH']
CLQs_MBHerr     = CLQs['MBH_err']
CLQs_eta        = CLQs['eta']  
#CLQs_MBH        = CLQs['MBH_CIV']  ## MBH_CIV from Kozłowski (2017)
#CLQs_eta        = CLQs['eta_CIV']  ## MBH_CIV from Kozłowski (2017)
CLQs_loglineLum = CLQs_loglineLum - 42.  ## To ``normalize'' to 10^42 erg s^-1
CLQs_logContLum = CLQs_logContLum - 42.  ## To ``normalize'' to 10^42 erg s^-1

##
##     Q S F I T 
## CONT*__LUM::  1450, 2245, 3000, 4210 or 5100:
## AGN continuum (Section 2.2) νLν luminosities and their uncertainties, in units of 1042 erg s−1
##
path      = '../../../data/QSO_CIV_catalogs/'
filename    = 'qsfit_1.2.4.fits'
infile      = path+filename
sdata_full  = fits.open(infile)
sdata       = sdata_full[1].data

data_qsofit = sdata[np.where((sdata['BR_CIV_1549__LUM'] > 0. ) &
                             (sdata['BR_CIV_1549__EW']  > 0. ) )] 

log_L1450_QSFit = np.log10((data_qsofit['BR_CIV_1549__LUM'] * 1e42))
log_LCont_QSFit = np.log10((data_qsofit['CONT1__LUM']       * 1e42))
Cont1__WAVE     =           data_qsofit['CONT1__WAVE']
EW_QSFit        =           data_qsofit['BR_CIV_1549__EW']
log_EW_QSFit    = np.log10( data_qsofit['BR_CIV_1549__EW'])
logBH_virial    = np.log10((data_qsofit['BR_CIV_1549__LUM']**.5)*(data_qsofit['BR_CIV_1549__FWHM']**2))  

##
##  Q S O    C A T A L O Gs
## Hamann et al. (2017) and Kozłowski (2017) merged BOSS DR12 catalog
filename  = 'BOSS_Ham17Koz17_DR12Q.fits'
infile    = path+filename
data_full = fits.open(infile)
tdata     = data_full[1].data

## Knocking out a couple of objects with bad REW values;
## No BALs; qflag; quality flag for the C IV fit: 0 = no problems
## 2.0 < ze < 3.4  to match "Full Sample" in Hamann++17
## and kick out the LBol < 0.
data = tdata[np.where( (tdata['rew'] > 0.     )     &
                       (tdata['rew'] < 10000. )     &   ## there's one point at like 2,700,000! 
                       (tdata['bal_flag_vi'] == 0 ) &
                       (tdata['qflag']       == 0 ) & 
                       (tdata['z_dr12']    > 2.0 )  &
                       (tdata['z_dr12']    < 3.4 )  & 
                       (tdata['Lbol']      > 0.0 )  &
                       (tdata['MBH_MgII']  > 4.0 ))]

## Setting up variable names
f1450_Ham17     = data['f1450']
alpha_civ_Ham17 = data['alpha_civ']
REW_Ham17       = data['rew']        
FWHM_Ham17      = data['fwhm']
z_dr12          = data['z_dr12']
log_LBol_Koz17  = data['LBol']
log_L1350_Koz17 = data['L1350']
log_MBH         = data['MBH_MgII']
log_REW         = np.log10(data['rew'])
log_REW_Ham17   = np.log10(REW_Ham17)    ## same thing, diff name!
FWHM_kilo_Ham17 = FWHM_Ham17/1000.
#log_tMBH        = tdata['MBH_MgII']
log_tMBH        = tdata['MBH_CVI']
eta_Edd         = tdata['nEdd']
log_tREW        = np.log10(tdata['rew'])


 
##   S E T T I N G   U P   T H E    P L O T  
matplotlib.rc('text', usetex=True)
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(14, 4.5), sharey=True, dpi=80, facecolor='w', edgecolor='k')  ## used to be 7,5

## Adjusting the Whitespace for the plots
left   = 0.06   # the left side of the subplots of the figure
right  = 0.96   # the right side of the subplots of the figure
bottom = 0.14   # the bottom of the subplots of the figure
top    = 0.96   # the top of the subplots of the figure
wspace = 0.26   # the amount of width reserved for blank space between subplots
hspace = 0.06   # the amount of height reserved for white space between subplots
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)

fig.subplots_adjust(wspace=0)

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
majorticklength = 6
minorticklength = 3


##
##   X - a x i s
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
x_one = log_LCont_QSFit - 42.  ## -42.0 for more direct comparison to Fig 15 of Calderone+2017.

##
##  y  - a x i s 
##
#y = log_REW
y_one = log_EW_QSFit 


##
##  P L O T T I N G    T H E    D A T A
##
## Further defaults and  hexbin params
color_map = plt.cm.Spectral_r

mincnt    = 3.
gridsize  = 200
ax1.scatter(     x_one, y_one, s=s, alpha=alpha)
hb1 = ax1.hexbin(x_one, y_one, bins='log', gridsize=gridsize, cmap=color_map, mincnt=mincnt)

mincnt    = 5.
gridsize  = 320
color_map = plt.cm.Spectral_r
ax2.scatter(     log_tMBH, log_tREW, s=s, alpha=alpha)
hb2 = ax2.hexbin(log_tMBH, log_tREW, bins='log', gridsize=gridsize, cmap=color_map, mincnt=mincnt)
#ax2.scatter(     logBH_virial, y_one, s=s, alpha=alpha)
#hb2 = ax2.hexbin(logBH_virial, y_one, bins='log', gridsize=gridsize, cmap=color_map, mincnt=mincnt)
ax3.scatter(     eta_Edd,  log_tREW, s=s, alpha=alpha)
hb3 = ax3.hexbin(eta_Edd,  log_tREW, bins='log', gridsize=gridsize, cmap=color_map, mincnt=mincnt)


## docs.scipy.org/doc/scipy/reference/generated/scipy.stats.linregress.html
## Calculate a linear least-squares regression for two sets of measurements.
slope, intercept, r_value, p_value, std_err = stats.linregress(x_one, y_one)
print()
print("Performed linear least-squares regression...")
print("The slope     is: %f" % slope)
print("The intercept is: %f" % intercept)
print("R-squared       : %f" % r_value**2)
print("p-value         : %f" % p_value)
print("Standard errro  : %f" % std_err)
print()
plt.plot(x_one, (slope*x_one + intercept), 'r')
#plt.plot(x, (slope*x), color='mediumturquoise', linestyle='--')
## for the label...
slope_str     = str(np.around(slope, decimals=3))
intercept_str = str(np.around(intercept, decimals=3))
R_sq          = str(np.around((r_value**2), decimals=3))               
        

##  A X E S    L I M I T S 
##
## Y-axis the same across all three plots
y_min =  0.4
y_max =  2.999    ## when in '000 km s^-1;  3.4 when the Legend handles has the # and slope line

##   EW vs. Cont. Lumin
x1_min =  2.30   
x1_max =  5.55    ## log_LCont_QSFit.max() ## log_LBol.max()   ## can be e.g. 48.4 for log_LBol
##   EW vs. Log MBH
x2_min =   7.30   
x2_max =  10.88
##   EW vs. eta  (Lbol/LEdd)
x3_min =  -3.2   
x3_max =   1.7
ax1.axis([x1_min, x1_max, y_min, y_max])
ax2.axis([x2_min, x2_max, y_min, y_max])
ax3.axis([x3_min, x3_max, y_min, y_max])

no_MBH_points = len(log_tMBH[np.where( (log_tMBH > x2_min) & (log_tMBH < x2_max) & (log_tREW > y_min) * (log_tREW < y_max))])
                        


## AXES LABELS
ax1.set_xlabel(r"log$_{10}[\nu    $L$_{\nu}    $(1450{\rm \AA })] / 10$^{42}$ erg s$^{-1}$",  fontsize=fontsize)
ax2.set_xlabel(r"log$_{10}$[M$_{\rm BH}$]",                           fontsize=fontsize)
ax3.set_xlabel(r"$\eta$ = log$_{10}$[L$_{\rm Bol}$ / L$_{\rm Edd}$]", fontsize=fontsize)
ax1.set_ylabel(r"log$_{10}$(REW CIV / {\rm \AA })",                  fontsize=fontsize)

## AXES TICK FORMAT
ax1.tick_params(axis='both', which='major', labelsize=labelsize, right=True, length=majorticklength)
ax1.tick_params(axis='both', which='minor', labelsize=labelsize)
ax2.tick_params(axis='both', which='major', labelsize=labelsize, right=True, length=majorticklength)
ax2.tick_params(axis='both', which='minor', labelsize=labelsize)
ax3.tick_params(axis='both', which='major', labelsize=labelsize, right=True, length=majorticklength)
ax3.tick_params(axis='both', which='minor', labelsize=labelsize)

## Plotting each spectral epoch on its own (there's a nicer way to
## code this, but I don't know how!)
for ii in range(len(CLQs)):
    ms         =  16.
    ms_back    =  45.
    l_back     =   2
    alpha      =  1.0
    markersize = 6.0
    
    if str(CLQs['Object'][ii]) == 'J1205p3422':
        if (CLQs['MJD'][ii] ==  53498):
            ax1.scatter(CLQs_logContLum[ii], CLQs_logREW[ii], color='k',       alpha=alpha, marker='o', s=ms_back)
            ax1.scatter(CLQs_logContLum[ii], CLQs_logREW[ii], color='fuchsia', alpha=alpha, marker='o', s=ms)
            ax2.scatter(CLQs_MBH[ii],        CLQs_logREW[ii], color='k',       alpha=alpha, marker='o', s=ms_back)
            ax2.scatter(CLQs_MBH[ii],        CLQs_logREW[ii], color='fuchsia', alpha=alpha, marker='o', s=ms)
            ax3.scatter(CLQs_eta[ii],        CLQs_logREW[ii], color='k',       alpha=alpha, marker='o', s=ms_back)
            ax3.scatter(CLQs_eta[ii],        CLQs_logREW[ii], color='fuchsia', alpha=alpha, marker='o', s=ms)
            J12_53498 = mlines.Line2D([], [], label='J1205+3422 (53498)', color='fuchsia',
                        marker="o", markeredgecolor='k', markeredgewidth=1.4, markersize=markersize,  linestyle='None')

        if (CLQs['MJD'][ii] ==  58538):
            ax1.scatter(CLQs_logContLum[ii], CLQs_logREW[ii], color='k',       alpha=alpha, marker='s', s=ms_back)
            ax1.scatter(CLQs_logContLum[ii], CLQs_logREW[ii], color='fuchsia', alpha=alpha, marker='s', s=ms)
            ax2.scatter(CLQs_MBH[ii],        CLQs_logREW[ii], color='k',       alpha=alpha, marker='s', s=ms_back)
            ax2.scatter(CLQs_MBH[ii],        CLQs_logREW[ii], color='fuchsia', alpha=alpha, marker='s', s=ms)
            ax3.scatter(CLQs_eta[ii],        CLQs_logREW[ii], color='k',       alpha=alpha, marker='s', s=ms_back)
            ax3.scatter(CLQs_eta[ii],        CLQs_logREW[ii], color='fuchsia', alpha=alpha, marker='s', s=ms)
            J12_58538 = mlines.Line2D([], [], label='J1205+3422 (58538)', color='fuchsia',
                        marker="s", markeredgecolor='k', markeredgewidth=1.4, markersize=markersize,  linestyle='None')

        if (CLQs['MJD'][ii] ==  58693):
            ax1.scatter(CLQs_logContLum[ii], CLQs_logREW[ii], color='k',       alpha=alpha, marker='D', s=ms_back)
            ax1.scatter(CLQs_logContLum[ii], CLQs_logREW[ii], color='fuchsia', alpha=alpha, marker='D', s=ms)
            ax2.scatter(CLQs_MBH[ii],        CLQs_logREW[ii], color='k',       alpha=alpha, marker='D', s=ms_back)
            ax2.scatter(CLQs_MBH[ii],        CLQs_logREW[ii], color='fuchsia', alpha=alpha, marker='D', s=ms)
            ax3.scatter(CLQs_eta[ii],        CLQs_logREW[ii], color='k',       alpha=alpha, marker='D', s=ms_back)
            ax3.scatter(CLQs_eta[ii],        CLQs_logREW[ii], color='fuchsia', alpha=alpha, marker='D', s=ms)
            J12_58693 = mlines.Line2D([], [], label='J1205+3422 (58693)', color='fuchsia',
                        marker="D", markeredgecolor='k', markeredgewidth=1.4, markersize=markersize,  linestyle='None')
            

    if str(CLQs['Object'][ii]) == 'J1638p2827':
        if (CLQs['MJD'][ii] ==  54553):
            ax1.scatter(CLQs_logContLum[ii], CLQs_logREW[ii], color='k',    alpha=alpha, marker='o', s=ms_back)
            ax1.scatter(CLQs_logContLum[ii], CLQs_logREW[ii], color='lime', alpha=alpha, marker='o', s=ms)
            ax2.scatter(CLQs_MBH[ii],        CLQs_logREW[ii], color='k',    alpha=alpha, marker='o', s=ms_back)
            ax2.scatter(CLQs_MBH[ii],        CLQs_logREW[ii], color='lime', alpha=alpha, marker='o', s=ms)
            ax3.scatter(CLQs_eta[ii],        CLQs_logREW[ii], color='k',    alpha=alpha, marker='o', s=ms_back)
            ax3.scatter(CLQs_eta[ii],        CLQs_logREW[ii], color='lime', alpha=alpha, marker='o', s=ms)
            J16_54553 = mlines.Line2D([], [], label='J1638+2827 (54553)', color='lime',
                        marker="o", markeredgecolor='k', markeredgewidth=1.4, markersize=markersize,  linestyle='None')
            
        if (CLQs['MJD'][ii] ==  55832):
            ax1.scatter(CLQs_logContLum[ii], CLQs_logREW[ii], color='k',    alpha=alpha, marker='s', s=ms_back)
            ax1.scatter(CLQs_logContLum[ii], CLQs_logREW[ii], color='lime', alpha=alpha, marker='s', s=ms)
            ax2.scatter(CLQs_MBH[ii],        CLQs_logREW[ii], color='k',    alpha=alpha, marker='s', s=ms_back)
            ax2.scatter(CLQs_MBH[ii],        CLQs_logREW[ii], color='lime', alpha=alpha, marker='s', s=ms)
            ax3.scatter(CLQs_eta[ii],        CLQs_logREW[ii], color='k',    alpha=alpha, marker='s', s=ms_back)
            ax3.scatter(CLQs_eta[ii],        CLQs_logREW[ii], color='lime', alpha=alpha, marker='s', s=ms)
            J16_55832 = mlines.Line2D([], [], label='J1638+2827 (55832)', color='lime',
                        marker="s", markeredgecolor='k', markeredgewidth=1.4, markersize=markersize,  linestyle='None')
            
        if (CLQs['MJD'][ii] ==  58583): 
            ax1.scatter(CLQs_logContLum[ii], CLQs_logREW[ii], color='k',    alpha=alpha, marker='D', s=ms_back)
            ax1.scatter(CLQs_logContLum[ii], CLQs_logREW[ii], color='lime', alpha=alpha, marker='D', s=ms)
            ax2.scatter(CLQs_MBH[ii],        CLQs_logREW[ii], color='k',    alpha=alpha, marker='D', s=ms_back)
            ax2.scatter(CLQs_MBH[ii],        CLQs_logREW[ii], color='lime', alpha=alpha, marker='D', s=ms)
            ax3.scatter(CLQs_eta[ii],        CLQs_logREW[ii], color='k',    alpha=alpha, marker='D', s=ms_back)
            ax3.scatter(CLQs_eta[ii],        CLQs_logREW[ii], color='lime', alpha=alpha, marker='D', s=ms)
            J16_58583 = mlines.Line2D([], [], label='J1638+2827 (58583)', color='lime',
                        marker="D", markeredgecolor='k', markeredgewidth=1.4, markersize=markersize,  linestyle='None')
            
    if str(CLQs['Object'][ii]) == 'J2228p2201':
        if (CLQs['MJD'][ii] ==  56189):
            ax1.scatter(CLQs_logContLum[ii], CLQs_logREW[ii], color='k',    alpha=alpha, marker='o', s=ms_back)
            ax1.scatter(CLQs_logContLum[ii], CLQs_logREW[ii], color='cyan', alpha=alpha, marker='o', s=ms)
            ax2.scatter(CLQs_MBH[ii],        CLQs_logREW[ii], color='k',    alpha=alpha, marker='o', s=ms_back)
            ax2.scatter(CLQs_MBH[ii],        CLQs_logREW[ii], color='cyan', alpha=alpha, marker='o', s=ms)
            ax3.scatter(CLQs_eta[ii],        CLQs_logREW[ii], color='k',    alpha=alpha, marker='o', s=ms_back)
            ax3.scatter(CLQs_eta[ii],        CLQs_logREW[ii], color='cyan', alpha=alpha, marker='o', s=ms)
            J22_56189 = mlines.Line2D([], [], label='J2228+2201 (56189)', color='cyan',
                        marker="o", markeredgecolor='k', markeredgewidth=1.4, markersize=markersize,  linestyle='None')
            
        if (CLQs['MJD'][ii] ==  56960):
            ax1.scatter(CLQs_logContLum[ii], CLQs_logREW[ii], color='k',    alpha=alpha, marker='s', s=ms_back)
            ax1.scatter(CLQs_logContLum[ii], CLQs_logREW[ii], color='cyan', alpha=alpha, marker='s', s=ms)
            ax2.scatter(CLQs_MBH[ii],        CLQs_logREW[ii], color='k',    alpha=alpha, marker='s', s=ms_back)
            ax2.scatter(CLQs_MBH[ii],        CLQs_logREW[ii], color='cyan', alpha=alpha, marker='s', s=ms)
            ax3.scatter(CLQs_eta[ii],        CLQs_logREW[ii], color='k',    alpha=alpha, marker='s', s=ms_back)
            ax3.scatter(CLQs_eta[ii],        CLQs_logREW[ii], color='cyan', alpha=alpha, marker='s', s=ms)
            J22_56960 = mlines.Line2D([], [], label='J2228+2201 (56960)', color='cyan',
                        marker="s", markeredgecolor='k', markeredgewidth=1.4, markersize=markersize,  linestyle='None')

        if (CLQs['MJD'][ii] ==  58693):
            ax1.scatter(CLQs_logContLum[ii], CLQs_logREW[ii], color='k',    alpha=alpha, marker='D', s=ms_back)
            ax1.scatter(CLQs_logContLum[ii], CLQs_logREW[ii], color='cyan', alpha=alpha, marker='D', s=ms)
            ax2.scatter(CLQs_MBH[ii],        CLQs_logREW[ii], color='k',    alpha=alpha, marker='D', s=ms_back)
            ax2.scatter(CLQs_MBH[ii],        CLQs_logREW[ii], color='cyan', alpha=alpha, marker='D', s=ms)
            ax3.scatter(CLQs_eta[ii],        CLQs_logREW[ii], color='k',    alpha=alpha, marker='D', s=ms_back)
            ax3.scatter(CLQs_eta[ii],        CLQs_logREW[ii], color='cyan', alpha=alpha, marker='D', s=ms)
            J22_58693 = mlines.Line2D([], [], label='J2228+2201 (58693)', color='cyan',
                        marker="D", markeredgecolor='k', markeredgewidth=1.4, markersize=markersize,  linestyle='None')

##
##    L E G E N D S
##  
boss = mlines.Line2D([], [], label='SDSS quasars',                        
                         color='skyblue', marker=".", markeredgecolor='k', markeredgewidth=1.1,
                         markersize=7, linestyle='None')

#bestfit_line  = mlines.Line2D([], [], label=r'$\beta = $'+slope_str, color='r')
bestfit_line  = mlines.Line2D([], [], label=r'$\beta$ = -0.1997', color='r')

handles=[boss, J12_53498, J12_58538, J12_58693,
#handles=[    J12_53498, J12_58538, J12_58693,
             J16_54553, J16_55832, J16_58583,
             J22_56189, J22_56960, J22_58693]
#             bestfit_line]
    
#leg = ax1.legend(loc='upper right',
 #                    fontsize=fontsize/1.85, handles=handles,
  #                   ncol=2,
   #                  frameon=True, framealpha=1.0, fancybox=True)

leg = ax3.legend(loc='upper right',
                     fontsize=fontsize/1.85, handles=handles,
                     ncol=1,
                     frameon=True, framealpha=1.0, fancybox=True)

## SAVING THE FIGURE
plt.savefig('CIV_CLQs_Baldwin_triptych_temp.png', format='png')
plt.close()



