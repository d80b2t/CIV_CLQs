'''
A ``hexbin'' plot for the SDSS-III BOSS data from the Hamann et al. (2017) catalog and analysis. 
This initially concentrated on 'just' f1450 converted into Lum_1450.
However, matching with the Kozłowski (2017) catalog means we can do the Baldwin plot with LBol directly. 
'''


import numpy as np
from scipy             import stats
from astropy.io        import fits
from astropy.io        import ascii
from astropy.cosmology import FlatLambdaCDM

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.lines  as mlines

## Setting up the cosmology
cosmo = FlatLambdaCDM(H0=71, Om0=0.27, Tcmb0=2.725)                                                                                                                         

## Reading in the data
path      = '/cos_pc19a_npr/programs/quasars/CIV_CLQs/data/CLQ_line_measurements/'
filename  = 'REW_FWHM.dat'
CLQs      = ascii.read(path+filename)

CLQ_FWHM   = CLQs['FWHM']    / 1000.
errFWHM    = CLQs['errFWHM'] / 1000.
CLQ_logREW = np.log10(CLQs['REW'])
## have to put the x-errors into log-space
#log_REW_hi = np.log10(CLQs['REW'] np.log10(CLQs['errREW']))
#log_REW_lo = np.log10(CLQs['REW']-CLQs['errREW'])
## and then have them as
#CLQ_asym_xerr = [log_REW_lo, log_REW_hi] 
log_errREW = np.log10(CLQs['errREW'])


##  Hamann et al. (2017) and Kozłowski (2017) merged BOSS DR12 catalog
path      = '/cos_pc19a_npr/programs/quasars/CIV_CLQs/data/QSO_CIV_catalogs/'
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
f1450     = data['f1450']
alpha_civ = data['alpha_civ']
REW       = data['rew']        
FWHM      = data['fwhm']
redshift  = data['z_dr12']
log_LBol  = data['LBol']
log_L1350 = data['L1350']


## Putting te Equiv. Widths into log and
## the FWHMs into '000s of kms^-1.
log_REW = np.log10(REW)
FWHM    = FWHM/1000.

## f1450 = flux in the uncorrected BOSS spectrum at 1450 Å rest (10−17 ergs s−1 cm−2 Å−1 )
f_lam   = f1450*((1450/1450)**alpha_civ)     

## L = 4 * pi * (D_L^2) * F
##   where L is in W and F is in W/m^2
##   3.086e+22 meters in a megaparsec
## Set-up the Luminosity Distance 
DL      = cosmo.luminosity_distance(redshift)             
DL_inM  = DL.value * 3.08567758128e22
DL_incm = DL.value * 3.08567758128e24

## DL in cm; f in 10−17 ergs s−1 cm−2 Å−1
Lum1450     = (4 * np.pi * (DL_incm*DL_incm) * 1e-17 * f_lam) * 1450.
log_Lum1450 = np.log10(Lum1450)


 
##   S E T T I N G   U P   T H E    P L O T  
matplotlib.rc('text', usetex=True)
fig, ax = plt.subplots(figsize=(7, 5), dpi=80, facecolor='w', edgecolor='k')

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


## Setting things all up... :-)
## L_1450 values taking f1450 from Hamann et al.
## L_1350 from Kozłowski (2017)
## L_Bol from Kozłowski (2017)
x = log_LBol
#x = log_L1350
#x = log_Lum1450


##
##  P L O T T I N G    T H E    D A T A
##
ax.scatter(     x, log_REW, s=s, alpha=alpha)
hb = plt.hexbin(x, log_REW, bins='log', gridsize=gridsize, cmap=color_map, mincnt=mincnt)

y = log_REW

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
plt.plot(x, (slope*x + intercept), 'r')
#plt.plot(x, (slope*x), color='mediumturquoise', linestyle='--')
## for the label...
slope_str     = str(np.around(slope, decimals=3))
intercept_str = str(np.around(intercept, decimals=3))
R_sq          = str(np.around((r_value**2), decimals=3))               

        

## AXES LIMITS 
xmin =  44.8              ##  - ( log_LBol.max()-log_L1350.max())
xmax =  log_LBol.max()    ## can be e.g. 48.4 for log_LBol
ymin =  0.5
ymax =  2.8  ## when in '000 km s^-1
ax.axis([xmin, xmax, ymin, ymax])

## AXES LABELS
#ax.set_xlabel(r"log$_{10}[\lambda$L$_{\lambda}$(1450{\rm \AA })] ergs/s ",  fontsize=fontsize)
#ax.set_xlabel(r"log$_{10}$(L$_{\rm 1350}$) [ergs/s] ",  fontsize=fontsize)
ax.set_xlabel(r"log$_{10}$(L$_{\rm Bol}$) [ergs/s] ",  fontsize=fontsize)
ax.set_ylabel(r"log$_{10}$(REW CIV / {\rm \AA })",     fontsize=fontsize)


## AXES TICK FORMAT
ax.tick_params(axis='both', which='major', labelsize=labelsize)
ax.tick_params(axis='both', which='minor', labelsize=labelsize)

##
##  L E G E N D S
#boss  = mlines.Line2D([], [], label='163,450 BOSS quasars', color='skyblue',
boss  = mlines.Line2D([], [], label='163,437 BOSS quasars', color='skyblue',
                       marker=".", markeredgecolor='k', markeredgewidth=1.4, markersize=7,  linestyle='None')

J12  = mlines.Line2D([], [], label='J1205+3422', color='fuchsia',
                       marker="o", markeredgecolor='k', markeredgewidth=1.4, markersize=7,  linestyle='None')
J16  = mlines.Line2D([], [], label='J1638+2827', color='lime',
                       marker="s", markeredgecolor='k', markeredgewidth=1.4, markersize=7,  linestyle='None')
J22  = mlines.Line2D([], [], label='J2228+2201', color='cyan',
                       marker="D", markeredgecolor='k', markeredgewidth=1.4, markersize=7,  linestyle='None')
## 
bestfit_line  = mlines.Line2D([], [], label=r'$\beta = $'+slope_str, color='r')

#handles=[boss, J12, J16, J22]
handles=[boss, bestfit_line ]
leg = ax.legend(loc='upper right', fontsize=fontsize/1.25, handles=handles, 
                frameon=True, framealpha=1.0, fancybox=True)


## SAVING THE FIGURE
plt.savefig('CIV_CLQs_Baldwin_temp.png', format='png')
plt.close()
