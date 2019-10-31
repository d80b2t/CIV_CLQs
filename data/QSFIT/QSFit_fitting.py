'''
The beginning of some code that takes the outputs from QSFit
e.g. the QSFIT.log file, and then overplots the model fit 
on top of the data spectrum. 

This would be very very much in the same vein as Figure 6 of
Calderone et al., 2017, MNRAS, 472, 4051 

Useful links/URLs::
  https://stackoverflow.com/questions/10582795/finding-the-full-width-half-maximum-of-a-peak
'''


import numpy as np
import pylab as pl
import scipy.optimize as opt

from scipy.interpolate import UnivariateSpline


spectrum_name = 'spec-6118-56189-0720'

infile_data  = spectrum_name+'.fits'
infile_model = spectrum_name+'_QSFIT.log'


def make_norm_dist(x, mean, sd):
    return 1.0/(sd*np.sqrt(2*np.pi))*np.exp(-(x - mean)**2/(2*sd**2))

x = np.linspace(10, 110, 1000)
green = make_norm_dist(x, 50, 10)
pink = make_norm_dist(x, 60, 10)

blue = green + pink   

# create a spline of x and blue-np.max(blue)/2 
spline = UnivariateSpline(x, blue-np.max(blue)/2, s=0)
r1, r2 = spline.roots() # find the roots

pl.plot(x, blue)
pl.axvspan(r1, r2, facecolor='g', alpha=0.5)
pl.show()


def FWHM(X,Y):
    half_max = max(Y) / 2.
    
    #find when function crosses line half_max (when sign of diff flips)
    #take the 'derivative' of signum(half_max - Y[])
    d = sign(half_max - array(Y[0:-1])) - sign(half_max - array(Y[1:]))
    
    #plot(X[0:len(d)],d) #if you are interested
    #find the left and right most indexes
    left_idx = find(d > 0)[0]
    right_idx = find(d < 0)[-1]
    
    #return the difference (full width)
    return X[right_idx] - X[left_idx]


def gauss(x, p): # p[0]==mean, p[1]==stdev
    return 1.0/(p[1]*np.sqrt(2*np.pi))*np.exp(-(x-p[0])**2/(2*p[1]**2))

# Create some sample data
known_param = np.array([2.0, .7])
xmin,xmax = -1.0, 5.0
N = 1000
X = np.linspace(xmin,xmax,N)
Y = gauss(X, known_param)

# Add some noise
Y += .10*np.random.random(N)

# Renormalize to a proper PDF
Y /= ((xmax-xmin)/N)*Y.sum()

# Fit a guassian
p0 = [0,1] # Inital guess is a normal distribution
errfunc = lambda p, x, y: gauss(x, p) - y # Distance to the target function
p1, success = opt.leastsq(errfunc, p0[:], args=(X, Y))

fit_mu, fit_stdev = p1

FWHM = 2*np.sqrt(2*np.log(2))*fit_stdev
print("FWHM", FWHM)



## for sake of ease, set all the inital compoments to zero
     CONTINUUM = 0   ## :  qsfit_comp_sbpowerlaw            
        GALAXY = 0   ## :  qsfit_comp_galaxytemplate     
        BALMER = 0   ## :  qsfit_comp_balmer             
        IRONUV = 0   ## :  qsfit_comp_ironuv             
       IRONOPT = 0   ## :  qsfit_comp_ironoptical        
        NA_LYA = 0   ## :  qsfit_comp_emline             
        BR_LYA = 0   ## :  qsfit_comp_emline             
    NA_NV_1241 = 0   ## :  qsfit_comp_emline             
    BR_OI_1306 = 0   ## :  qsfit_comp_emline             
   BR_CII_1335 = 0   ## :  qsfit_comp_emline             
  BR_SIIV_1400 = 0   ## :  qsfit_comp_emline             
   BR_CIV_1549 = 0   ## :  qsfit_comp_emline             
  BR_CIII_1909 = 0   ## :  qsfit_comp_emline             
        BR_CII = 0   ## :  qsfit_comp_emline             
  BR_MGII_2798 = 0   ## :  qsfit_comp_emline             
  NA_NEVI_3426 = 0   ## :  qsfit_comp_emline             
   NA_OII_3727 = 0   ## :  qsfit_comp_emline             
 NA_NEIII_3869 = 0   ## :  qsfit_comp_emline             
         BR_HD = 0   ## :  qsfit_comp_emline             
         BR_HG = 0   ## :  qsfit_comp_emline             
         NA_HB = 0   ## :  qsfit_comp_emline             
         BR_HB = 0   ## :  qsfit_comp_emline             
  NA_OIII_4959 = 0   ## :  qsfit_comp_emline             
  NA_OIII_5007 = 0   ## :  qsfit_comp_emline             
   BR_HEI_5876 = 0   ## :  qsfit_comp_emline             
   NA_NII_6549 = 0   ## :  qsfit_comp_emline             
         NA_HA = 0   ## :  qsfit_comp_emline             
         BR_HA = 0   ## :  qsfit_comp_emline             
   NA_NII_6583 = 0   ## :  qsfit_comp_emline             
   NA_SII_6716 = 0   ## :  qsfit_comp_emline             
   NA_SII_6731 = 0   ## :  qsfit_comp_emline             
          UNK1 = 0   ## :  qsfit_comp_emline             
          UNK2 = 0   ## :  qsfit_comp_emline             
          UNK3 = 0   ## :  qsfit_comp_emline             
          UNK4 = 0   ## :  qsfit_comp_emline             
          UNK5 = 0   ## :  qsfit_comp_emline             
          UNK6 = 0   ## :  qsfit_comp_emline             
          UNK7 = 0   ## :  qsfit_comp_emline             
          UNK8 = 0   ## :  qsfit_comp_emline             
          UNK9 = 0   ## :  qsfit_comp_emline             
         UNK10 = 0   ## :  qsfit_comp_emline             
  LINE_HA_BASE = 0   ## :  qsfit_comp_emline             
  LINE_OIII_BW = 0   ## :  qsfit_comp_emline             


  
br_CIV_1549 = 


EXPR_CONTINUUM = continuum
EXPR_GALAXY    = galaxy
EXPR_CONTGALAX = continuum + galaxy
EXPR_BALMER    = balmer
EXPR_IRON      = ironuv + ironopt
EXPR_BROADLINE = line_Ha_base + br_Lya + br_OI_1306 + br_CII_1335 + br_SiIV_1400 + br_CIV_1549 + br_CIII_1909 + br_CII + br_MgII_2798 + br_Hd + br_Hg + br_Hb + br_HeI_5876 + br_Ha
EXPR_NARROWLIN :  na_Lya + na_NV_1241 + na_NeVI_3426 + na_OII_3727 + na_NeIII_3869 + na_Hb + na_OIII_4959 + na_OIII_5007 + na_NII_6549 + na_Ha + na_NII_6583 + na_SII_6716 + na_SII_6731 + line_OIII_bw
EXPR_ABSLINES :  0
EXPR_UNKNOWN :  unk1 + unk2 + unk3 + unk4 + unk5 + unk6 + unk7 + unk8 + unk9 + unk10


FULL_FIT =  CONTINUUM + GALAXY + BALMER + IRONUV + IRONOPT + NA_LYA + BR_LYA + NA_NV_1241 + BR_OI_1306 + BR_CII_1335 + BR_SIIV_1400 + BR_CIV_1549 + BR_CIII_1909 + BR_CII + BR_MGII_2798 + NA_NEVI_3426 + NA_OII_3727 + NA_NEIII_3869 + BR_HD + BR_HG + NA_HB + BR_HB + NA_OIII_4959 + NA_OIII_5007 + BR_HEI_5876 + NA_NII_6549 + NA_HA + BR_HA + NA_NII_6583 + NA_SII_6716 + NA_SII_6731 + UNK1 + UNK2 + UNK3 + UNK4 + UNK5 + UNK6 + UNK7 + UNK8 + UNK9 + UNK10 + LINE_HA_BASE + LINE_OIII_BW

