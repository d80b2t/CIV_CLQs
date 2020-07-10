'''
  BLACK HOLE VIRIAL MASS ESTIMATES

The virial BH mass calibrations and scaling relations used in this
paper [Shen et al. 2011] are from McLure & Dunlop (2004, Hβ and Mg
II), Vestergaard & Peterson (2006, Hβ and C IV), and Vestergaard &
Osmer (2009, Mg II). 
These are for quasars:
 SDSS  J160544.68+342252.4

'''

import numpy as np


a_MD04_MgII = 0.505
b_MD04_MgII = 0.62

a_VO09_MgII = 0.860
b_VO09_MgII = 0.50

a_S10_MgII = 0.740
b_S10_MgII = 0.62

a_VP06_CIV = 0.660
b_VP06_CIV = 0.53

'''
	From Shen et al. (2011) 
	3.8. The Spectral Catalog
	``BC_1350 = 3.81 from the composite spectral energy distribution (SED) in Richards et al. (2006a).''
'''


'''
     J 1 2 0 5 + 3 4 2 2 

Data from Shen et al. (2011) 
Shen_dr7_bh_May_2010.fits
((PLATE == 2089) && (FIBER == 427) && (MJD == 53498))
'''

log_L3000     = 46.368801937169536
log_L3000_err = 0.008489474504106198
log_L1350     = 46.634680779574154
log_L1350err  = 0.004262343930953705

log_LMgII     = 44.45386264728122
log_LMgII_err = 0.045117652479984116
log_LCIV      = 44.896004599940575
log_LCIV_err  = 0.01388080688529314

FWHM_MgII     = 4729.812569514398
FWHM_MgII_err = 313.64594572537635
FWHM_CIV      = 5230.7138671875
FWHM_CIV_err  = 219.283447265625

LOGBH_MGII_MD04 = 9.639619130064474
LOGBH_MGII_VO09 = 9.394088830762948
LOGBH_MGII_S10  = 9.558345063223292
LOGBH_CIV_VP06  = 9.493502740592843

log_LBol      = 47.215605755249776
log_LBol_err  = 0.004262343930953705
log_BHmass    = 9.493502740592843               ##   LOGBH_CIV_VP06
log_LEdd      = np.log10(1.26e38)+log_BHmass    ## ``For pure ionized hydrogen''

log_EddRatio = -0.3782675304606329

## Values from QSFit
## These are the 'precise' values...
log_L1450_J16_1st = 46.61937511301522       ##  log10((41627*1e42))     
log_L1450_J16_2nd = 45.457881896733994		##  log10(( 2870*1e42))   
log_L1450_J16_3rd = 45.95640857119583	    ##  log10(( 9045*1e42))   
log_LCIV_J16_1st  = 44.84757265914211       ##  log10((  704*1e42))     
log_LCIV_J16_2nd  = 44.227886704613674		##  log10(( 169*1e42))   
log_LCIV_J16_3rd  = 44.59714648783373	    ##  log10((395.5*1e42))   
FWHM_CIV_J16_1st  = 4700.                   ## +/- 120 
FWHM_CIV_J16_2nd  = 14900.					## +/- 2800 
FWHM_CIV_J16_3rd  = 6983. 					## +/- 89 

#M_BH = a
#np.log10((3.8503777798695754e+45/1e44)**0.53) + (2*np.log10(4180.7216796875))+0.66

log_MBH_MD04_MgII_J16  = a_MD04_MgII + np.log10(( 10**log_L3000 /1e44)**b_MD04_MgII) + (2*np.log10(FWHM_MgII))
log_MBH_VO09_MgII_J16  = a_VO09_MgII + np.log10(( 10**log_L3000 /1e44)**b_VO09_MgII) + (2*np.log10(FWHM_MgII))
log_MBH_S10_MgII_J16   = a_S10_MgII  + np.log10(( 10**log_L3000 /1e44)**b_S10_MgII)  + (2*np.log10(FWHM_MgII))
log_MBH_VP06_CIV_J16   = a_VP06_CIV  + np.log10(( 10**log_L1350 /1e44)**b_VP06_CIV) + (2*np.log10(FWHM_CIV))


print()
print('  J 1 2 0 5  + 3 4 2 2   ')
print()
print('log_MBH_MD04_MgII_J16::  ', log_MBH_MD04_MgII_J16, '  vs. ', LOGBH_MGII_MD04 ,' in Shen et al. (2011)')
print('log_MBH_VO09_MgII_J16::  ', log_MBH_VO09_MgII_J16, '  vs. ', LOGBH_MGII_VO09 ,' in Shen et al. (2011)')
print('log_MBH_S10_MgII_J16::   ', log_MBH_S10_MgII_J16,  '  vs. ', LOGBH_MGII_S10  ,' in Shen et al. (2011)')
print('log_MBH_VP06_CIV_J16::   ', log_MBH_VP06_CIV_J16,  '  vs. ', LOGBH_CIV_VP06  ,' in Shen et al. (2011)')
print()
print(' From Shen et al. (2011), log(L_CIV)               =', log_LCIV,   '% of L_bol  = ', (10**(log_LCIV -log_LBol))*100.   )
print(' From Shen et al. (2011), log(log_L1350)           =', log_L1350,  '% of L_bol  = ', (10**(log_L1350-log_LBol))*100.   )
print(' From Shen et al. (2011), log(L_bol)               =', log_LBol,   '% of L_bol  = ', (10**(log_LBol-log_LBol))*100.   )
print(' From Shen et al. (2011), ``final'' log(BH_mass)   =', log_BHmass)
print(' ==> ')
print(' L_Edd = (4pi.G.M.m_p.c) / sigma_T = 1.26x10^38(M/Msol) erg/s = ', log_LEdd    )
print()
print(' LEdd / Lbol = ', log_LBol/log_LEdd, ' and log10(L_Edd) = ', np.log10((10**log_LBol)/(10**log_LEdd)), ', vs ', log_EddRatio, 'in Shen et al. (2011)')
print()

log_MBH_QSFit_J16_1st  = a_VP06_CIV  + np.log10(( 10**log_L1450_J16_1st /1e44)**b_VP06_CIV) + (2*np.log10(FWHM_CIV_J16_1st))
log_MBH_QSFit_J16_2nd  = a_VP06_CIV  + np.log10(( 10**log_L1450_J16_2nd /1e44)**b_VP06_CIV) + (2*np.log10(FWHM_CIV_J16_2nd))
log_MBH_QSFit_J16_3rd  = a_VP06_CIV  + np.log10(( 10**log_L1450_J16_3rd /1e44)**b_VP06_CIV) + (2*np.log10(FWHM_CIV_J16_3rd))

print('log_MBH_QSFit_J16 1st epoch (using VP06) ',  log_MBH_QSFit_J16_1st, '   vs.   9.66  currently in Table 3 (',  9.66-log_MBH_QSFit_J16_1st, 'diff)')
print('log_MBH_QSFit_J16 2nd epoch (using VP06) ',  log_MBH_QSFit_J16_2nd, '   vs   10.08  currently in Table 3 (', 10.08-log_MBH_QSFit_J16_2nd, 'diff)')
print('log_MBH_QSFit_J16 3rd epoch (using VP06) ',  log_MBH_QSFit_J16_3rd, '   vs.   9.67  currently in Table 3 (',  9.67-log_MBH_QSFit_J16_3rd, 'diff)')
print()

## Eddington Luminosity. For pure ionized hydrogen, 
## L_Edd = (4pi.G.M.m_p.c) / sigma_T 
##       =  1.26x10^38   (M/Msol) erg/s 
log_LEdd_J16_1st = np.log10(1.26e38)+log_MBH_QSFit_J16_1st 
log_LEdd_J16_2nd = np.log10(1.26e38)+log_MBH_QSFit_J16_2nd
log_LEdd_J16_3rd = np.log10(1.26e38)+log_MBH_QSFit_J16_3rd 

## Bolometric luminosity, LBol, computed from L1350/L1450 using the spectral fits and 
## bolometric corrections:: BC1350 = 3.81 from the composite spectral energy distribution (
## SED) in Richards et al. (2006a). 
##   log10(3.81) = 0.5809249756756193
log_LBol_J16_1st = log_L1450_J16_1st + 0.5809249756756216
log_LBol_J16_2nd = log_L1450_J16_2nd + 0.5809249756756216   
log_LBol_J16_3rd = log_L1450_J16_3rd + 0.5809249756756216   

print()
print('J16 1st epoch ::  log_LBol = ', log_LBol_J16_1st  ,' and LEdd = ', log_LEdd_J16_1st, ' and  eta =log10(L_Bol/L_Edd) = ', np.log10((10**log_LBol_J16_1st)/(10**log_LEdd_J16_1st)))	
print('J16 2nd epoch ::  log_LBol = ', log_LBol_J16_2nd  ,' and LEdd = ', log_LEdd_J16_2nd, ' and  eta =log10(L_Bol/L_Edd) = ', np.log10((10**log_LBol_J16_2nd)/(10**log_LEdd_J16_1st)))	
print('J16 3rd epoch ::  log_LBol = ', log_LBol_J16_3rd  ,' and LEdd = ', log_LEdd_J16_3rd, ' and  eta =log10(L_Bol/L_Edd) = ', np.log10((10**log_LBol_J16_3rd)/(10**log_LEdd_J16_3rd)))	
print()




''' 
    J 1 6 3 8 + 2 8 2 7
'''
## Values from Shen et al. (2011)
## ((PLATE == 2948) && (MJD == 54553) && (FIBER == 614))

log_L3000     = 45.49768513378036
log_L3000_err = 0.0835457146161076
log_L1350     = 45.58550334240817
log_L1350err  = 0.04035432415603424

log_LMgII     = 44.14148506498093
log_LMgII_err = 0.07201572951290913
log_LCIV      = 44.44932090425458
log_LCIV_err  = 0.035716244627341354

FWHM_MgII     = 4757.42541816648
FWHM_MgII_err = 2224.264991733375
FWHM_CIV      = 4180.7216796875
FWHM_CIV_err  = 735.8486328125

LOGBH_MGII_MD04 = 9.161413814175969
LOGBH_MGII_VO09 = 8.963586544105137
LOGBH_MGII_S10  = 9.031423035597872
LOGBH_CIV_VP06  = 8.742819284525387

log_LBol      = 46.16642831808379
log_LBol_err  = 0.04035432415603424
log_BHmass    = 8.742819284525387
log_LEdd      = np.log10(1.26e38)+log_BHmass    ## ``For pure ionized hydrogen''
log_EddRatio = -0.6767615115591639

## Values from Kozłowski, 2017, ApJS, 228, 9 
log_L3000_Koz17    = 45.942   ##  +/- 0.128
log_L1350_Koz17    = 46.172	  ##  +/- 0.089
log_LBol_Koz17     = 46.721	  ##  +/- 0.073
log_MBH_MgII_Koz17 =  9.338   ##  no error reported
log_MBH_CIV_Koz17  =  9.128   ##  no error reported
log_eta_Koz17      = -0.717   ##  no error reported


## Values from QSFit
log_L1450_J16_1st = 45.60852603357719		##  log10((4060   *1e42))   
log_L1450_J16_2nd = 45.40925665203891		##  log10((2566   *1e42))   
log_L1450_J16_3rd = 45.934700401715425		##  log10((8604   *1e42))   
log_LCIV_J16_1st  = 44.84757265914211       ##  log10(( 329   *1e42))     
log_LCIV_J16_2nd  = 44.227886704613674		##  log10((  99.8 *1e42))  
log_LCIV_J16_3rd  = 44.59714648783373	    ##  log10(( 395.5 *1e42))   
FWHM_CIV_J16_1st  = 4630.                    ##  $\pm$  190
FWHM_CIV_J16_2nd  = 4990.                    ##  $\pm$  300
FWHM_CIV_J16_3rd  = 4616.                    ##  $\pm$  72

#M_BH = a
#np.log10((3.8503777798695754e+45/1e44)**0.53) + (2*np.log10(4180.7216796875))+0.66

log_MBH_MD04_MgII_J16 = a_MD04_MgII + np.log10(( 10**log_L3000 /1e44)**b_MD04_MgII) + (2*np.log10(FWHM_MgII))
log_MBH_VO09_MgII_J16 = a_VO09_MgII + np.log10(( 10**log_L3000 /1e44)**b_VO09_MgII) + (2*np.log10(FWHM_MgII))
log_MBH_S10_MgII_J16  = a_S10_MgII  + np.log10(( 10**log_L3000 /1e44)**b_S10_MgII)  + (2*np.log10(FWHM_MgII))
log_MBH_VP06_CIV_J16  = a_VP06_CIV  + np.log10(( 10**log_L1350 /1e44)**b_VP06_CIV)  + (2*np.log10(FWHM_CIV))

print()
print()
print('   J 1 6 3 8 + 2 8 2 7   ')
print()
print('log_MBH_MD04_MgII_J16::  ', log_MBH_MD04_MgII_J16, ' vs. ', LOGBH_MGII_MD04 ,' in Shen et al. (2011)')
print('log_MBH_VO09_MgII_J16::  ', log_MBH_VO09_MgII_J16, ' vs. ', LOGBH_MGII_VO09 ,' in Shen et al. (2011)')
print('log_MBH_S10_MgII_J16 ::  ', log_MBH_S10_MgII_J16,  ' vs. ', LOGBH_MGII_S10  ,' in Shen et al. (2011)')
print('log_MBH_VP06_CIV_J16 ::  ', log_MBH_VP06_CIV_J16,  ' vs. ', LOGBH_CIV_VP06  ,' in Shen et al. (2011)')
print()
print(' From Shen et al. (2011), log(L_CIV)               =', log_LCIV,   '% of L_bol  = ', (10**(log_LCIV -log_LBol))*100.   )
print(' From Shen et al. (2011), log(log_L1350)           =', log_L1350,  '% of L_bol  = ', (10**(log_L1350-log_LBol))*100.   )
print(' From Shen et al. (2011), log(L_bol)               =', log_LBol,   '% of L_bol  = ', (10**(log_LBol-log_LBol))*100.   )
print(' From Shen et al. (2011), ``final'' log(BH_mass)   =', log_BHmass)
print(' ==> ')
print(' L_Edd = (4pi.G.M.m_p.c) / sigma_T = 1.26x10^38(M/Msol) erg/s = ', log_LEdd    )
print()
print(' LEdd / Lbol = ', log_LBol/log_LEdd, ' and log10(L_Edd) = ', np.log10((10**log_LBol)/(10**log_LEdd)), ', vs ', log_EddRatio, 'in Shen et al. (2011)')
print()

print(' From Shen et al. (2011), log(L_CIV)               =', log_LCIV,   '% of L_bol  = ', (10**(log_LCIV -log_LBol))*100.   )
print(' From Shen et al. (2011), log(log_L1350)           =', log_L1350,  '% of L_bol  = ', (10**(log_L1350-log_LBol))*100.   )
print(' From Shen et al. (2011), log(L_bol)               =', log_LBol,   '% of L_bol  = ', (10**(log_LBol-log_LBol))*100.   )
print(' From Shen et al. (2011), ``final'' log(BH_mass)   =', log_BHmass)
print(' ==> ')
print(' L_Edd = (4pi.G.M.m_p.c) / sigma_T = 1.26x10^38(M/Msol) erg/s = ', log_LEdd    )
print()
print(' LEdd / Lbol = ', log_LBol/log_LEdd, ' and log10(L_Edd) = ', np.log10((10**log_LBol)/(10**log_LEdd)), ', vs ', log_EddRatio, 'in Shen et al. (2011)')
print()


log_MBH_QSFit_J16_1st  = a_VP06_CIV  + np.log10(( 10**log_L1450_J16_1st /1e44)**b_VP06_CIV) + (2*np.log10(FWHM_CIV_J16_1st))
log_MBH_QSFit_J16_2nd  = a_VP06_CIV  + np.log10(( 10**log_L1450_J16_2nd /1e44)**b_VP06_CIV) + (2*np.log10(FWHM_CIV_J16_2nd))
log_MBH_QSFit_J16_3rd  = a_VP06_CIV  + np.log10(( 10**log_L1450_J16_3rd /1e44)**b_VP06_CIV) + (2*np.log10(FWHM_CIV_J16_3rd))

print('log_MBH_QSFit_J16_1st epoch',  log_MBH_QSFit_J16_1st, '  (using VP06)  vs.  9.14  in Table 3  (', 9.14-log_MBH_QSFit_J16_1st, 'diff)')
print('log_MBH_QSFit_J16_2nd epoch',  log_MBH_QSFit_J16_2nd, '  (using VP06)  vs.  9.10  in Table 3  (', 9.10-log_MBH_QSFit_J16_2nd, 'diff)')
print('log_MBH_QSFit_J16_3rd epoch',  log_MBH_QSFit_J16_3rd, '  (using VP06)  vs.  9.30  in Table 3  (', 9.30-log_MBH_QSFit_J16_3rd, 'diff)')

## Eddington Luminosity. For pure ionized hydrogen, 
## L_Edd = (4pi.G.M.m_p.c) / sigma_T 
##       =  1.26x10^38   (M/Msol) erg/s 
log_LEdd_J16_1st = np.log10(1.26e38)+log_MBH_QSFit_J16_1st 
log_LEdd_J16_2nd = np.log10(1.26e38)+log_MBH_QSFit_J16_2nd
log_LEdd_J16_3rd = np.log10(1.26e38)+log_MBH_QSFit_J16_3rd 

## Bolometric luminosity, LBol, computed from L1350/L1450 using the spectral fits and 
## bolometric corrections:: BC1350 = 3.81 from the composite spectral energy distribution (
## SED) in Richards et al. (2006a). 
##   log10(3.81) = 0.5809249756756193
log_LBol_J16_1st = log_L1450_J16_1st + 0.5809249756756216
log_LBol_J16_2nd = log_L1450_J16_2nd + 0.5809249756756216   
log_LBol_J16_3rd = log_L1450_J16_3rd + 0.5809249756756216   

print()
print('J16 1st epoch ::  log_LBol = ', log_LBol_J16_1st  ,' and LEdd = ', log_LEdd_J16_1st, ' and  eta =log10(L_Bol/L_Edd) = ', np.log10((10**log_LBol_J16_1st)/(10**log_LEdd_J16_1st)))	
print('J16 2nd epoch ::  log_LBol = ', log_LBol_J16_2nd  ,' and LEdd = ', log_LEdd_J16_2nd, ' and  eta =log10(L_Bol/L_Edd) = ', np.log10((10**log_LBol_J16_2nd)/(10**log_LEdd_J16_1st)))	
print('J16 3rd epoch ::  log_LBol = ', log_LBol_J16_3rd  ,' and LEdd = ', log_LEdd_J16_3rd, ' and  eta =log10(L_Bol/L_Edd) = ', np.log10((10**log_LBol_J16_3rd)/(10**log_LEdd_J16_3rd)))	
print()



''' 
   J 2 2 2 8  +   2 2 0 1 
'''
## Values from Kozłowski, 2017, ApJS, 228, 9 
log_L3000_Koz17    = 45.477   ##  +/- 0.116
log_L1350_Koz17    = 45.678	  ##  +/- 0.094
log_LBol_Koz17     = 46.231	  ##  +/- 0.073
log_MBH_MgII_Koz17 =  9.448   ##  no error reported
log_MBH_CIV_Koz17  =  8.733	  ##  no error reported
log_eta_Koz17      = -1.317	  ##  no error reported


## Values from QSFit
log_L1450_1st = 44.93449845124357
log_L1450_2nd = 45.961421094066445
log_L1450_3rd = 45.44870631990508
FWHM_CIV_1st = 5930.
FWHM_CIV_2nd = 7000.
FWHM_CIV_3rd = 5930.

print()
print()
print()
print('   J 2 2 2 8  + 2 2 0 1   ')
print()
print('log_MBH_MD04_MgII_J16::  ', log_MBH_MD04_MgII_J16, ' vs. ', LOGBH_MGII_MD04 ,' in Shen et al. (2011)')
print('log_MBH_VO09_MgII_J16::  ', log_MBH_VO09_MgII_J16, ' vs. ', LOGBH_MGII_VO09 ,' in Shen et al. (2011)')
print('log_MBH_S10_MgII_J16 ::  ', log_MBH_S10_MgII_J16,  ' vs. ', LOGBH_MGII_S10  ,' in Shen et al. (2011)')
print('log_MBH_VP06_CIV_J16 ::  ', log_MBH_VP06_CIV_J16,  ' vs. ', LOGBH_CIV_VP06  ,' in Shen et al. (2011)')
print()
print()
print(' From Shen et al. (2011), log(L_CIV)               =', log_LCIV,   '% of L_bol  = ', (10**(log_LCIV -log_LBol))*100.   )
print(' From Shen et al. (2011), log(log_L1350)           =', log_L1350,  '% of L_bol  = ', (10**(log_L1350-log_LBol))*100.   )
print(' From Shen et al. (2011), log(L_bol)               =', log_LBol,   '% of L_bol  = ', (10**(log_LBol-log_LBol))*100.   )
print(' From Shen et al. (2011), ``final'' log(BH_mass)   =', log_BHmass)
print(' ==> ')
print(' L_Edd = (4pi.G.M.m_p.c) / sigma_T = 1.26x10^38(M/Msol) erg/s = ', log_LEdd    )
print()
print(' LEdd / Lbol = ', log_LBol/log_LEdd, ' and log10(L_Edd) = ', np.log10((10**log_LBol)/(10**log_LEdd)), ', vs ', log_EddRatio, 'in Shen et al. (2011)')
print()

log_MBH_QSFit_J16_1st  = a_VP06_CIV  + np.log10(( 10**log_L1450_1st /1e44)**b_VP06_CIV) + (2*np.log10(FWHM_CIV_1st))
log_MBH_QSFit_J16_2nd  = a_VP06_CIV  + np.log10(( 10**log_L1450_2nd /1e44)**b_VP06_CIV) + (2*np.log10(FWHM_CIV_2nd))
log_MBH_QSFit_J16_3rd  = a_VP06_CIV  + np.log10(( 10**log_L1450_3rd /1e44)**b_VP06_CIV) + (2*np.log10(FWHM_CIV_3rd))

print
print('log_MBH_QSFit_J16_1st epoch',  log_MBH_QSFit_J16_1st, '  (using VP06)  vs.  9.01  currently in the paper ( ', 9.01-log_MBH_QSFit_J16_1st, 'diff)')
print('log_MBH_QSFit_J16_2nd epoch',  log_MBH_QSFit_J16_2nd, '  (using VP06)  vs.  9.67  currently in the paper ( ', 9.67-log_MBH_QSFit_J16_2nd, 'diff)')
print('log_MBH_QSFit_J16_3rd epoch',  log_MBH_QSFit_J16_3rd, '  (using VP06)  vs.  9.27  currently in the paper ( ', 9.26-log_MBH_QSFit_J16_3rd, 'diff)')
print()

## Eddington Luminosity. For pure ionized hydrogen, 
## L_Edd = (4pi.G.M.m_p.c) / sigma_T 
##       =  1.26x10^38   (M/Msol) erg/s 
log_LEdd_J16_1st = np.log10(1.26e38)+log_MBH_QSFit_J16_1st 
log_LEdd_J16_2nd = np.log10(1.26e38)+log_MBH_QSFit_J16_2nd
log_LEdd_J16_3rd = np.log10(1.26e38)+log_MBH_QSFit_J16_3rd 

## Bolometric luminosity, LBol, computed from L1350/L1450 using the spectral fits and 
## bolometric corrections:: BC1350 = 3.81 from the composite spectral energy distribution (
## SED) in Richards et al. (2006a). 
##   log10(3.81) = 0.5809249756756193
log_LBol_J16_1st = log_L1450_J16_1st + 0.5809249756756216
log_LBol_J16_2nd = log_L1450_J16_2nd + 0.5809249756756216   
log_LBol_J16_3rd = log_L1450_J16_3rd + 0.5809249756756216   

print()
print('J16 1st epoch ::  log_LBol = ', log_LBol_J16_1st  ,' and LEdd = ', log_LEdd_J16_1st, ' and  eta =log10(L_Bol/L_Edd) = ', np.log10((10**log_LBol_J16_1st)/(10**log_LEdd_J16_1st)))	
print('J16 2nd epoch ::  log_LBol = ', log_LBol_J16_2nd  ,' and LEdd = ', log_LEdd_J16_2nd, ' and  eta =log10(L_Bol/L_Edd) = ', np.log10((10**log_LBol_J16_2nd)/(10**log_LEdd_J16_1st)))	
print('J16 3rd epoch ::  log_LBol = ', log_LBol_J16_3rd  ,' and LEdd = ', log_LEdd_J16_3rd, ' and  eta =log10(L_Bol/L_Edd) = ', np.log10((10**log_LBol_J16_3rd)/(10**log_LEdd_J16_3rd)))	
print()

