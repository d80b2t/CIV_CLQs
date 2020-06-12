
import numpy as np

'''
Scaling relations
The virial BH mass calibrations used in this paper [Shen et al. 2011]
are from McLure & Dunlop (2004, Hβ and Mg ii), Vestergaard & Peterson
(2006, Hβ and C iv), and Vestergaard & Osmer (2009, Mg ii). These
calibrations have parameters:
'''

a_MD04_MgII = 0.505
b_MD04_MgII = 0.62

a_VO09_MgII = 0.860
b_VO09_MgII = 0.50

a_S10_MgII = 0.740
b_S10_MgII = 0.62

a_VP06_CIV = 0.660
b_VP06_CIV = 0.53




'''
     J 1 2 0 5 + 3 4 2 2 

Okay, data from Shen et al. (2011) for
((PLATE == 2089) && (FIBER == 427) && (MJD == 53498))
J120544.68+342252.4
'''
logL_MgII     = 44.45386264728122
logL_MgII_err = 0.045117652479984116
logL_CIV      = 44.896004599940575
logL_CIV_err  = 0.01388080688529314

logL_3000     = 46.368801937169536
logL_3000_err = 0.008489474504106198
logL_1350     = 46.634680779574154
logL_1350err  = 0.004262343930953705


FWHM_MgII     = 4729.812569514398
FWHM_MgII_err = 313.64594572537635
FWHM_CIV      = 5230.7138671875
FWHM_CIV_err  = 219.283447265625

## Values from QSFit
logL_1450_1st = 46.619406410886775
logL_1450_2nd = 45.457881896733994
logL_1450_3rd = 45.9566485792052
FWHM_CIV_1st = 4700.
FWHM_CIV_2nd = 14900.
FWHM_CIV_3rd = 6980.

#M_BH = a
#np.log10((3.8503777798695754e+45/1e44)**0.53) + (2*np.log10(4180.7216796875))+0.66

log_MBH_MD04_MgII_J12  = a_MD04_MgII + np.log10(( 10**logL_3000 /1e44)**b_MD04_MgII) + (2*np.log10(FWHM_MgII))
log_MBH_VO09_MgII_J12  = a_VO09_MgII + np.log10(( 10**logL_3000 /1e44)**b_VO09_MgII) + (2*np.log10(FWHM_MgII))
log_MBH_S10_MgII_J12   = a_S10_MgII  + np.log10(( 10**logL_3000 /1e44)**b_S10_MgII)  + (2*np.log10(FWHM_MgII))
log_MBH_VP06_CIV_J12   = a_VP06_CIV  + np.log10(( 10**logL_1350 /1e44)**b_VP06_CIV) + (2*np.log10(FWHM_CIV))


print()
print('  J 1 2 0 5  + 3 4 2 2   ')
print()
print('log_MBH_MD04_MgII_J12::  ', log_MBH_MD04_MgII_J12, ' vs. 9.639619 in Shen et al. (2011)')
print('log_MBH_VO09_MgII_J12::  ', log_MBH_VO09_MgII_J12, ' vs. 9.394088 in Shen et al. (2011)')
print('log_MBH_S10_MgII_J12::   ',  log_MBH_S10_MgII_J12,  ' vs. 9.55834 in Shen et al. (2011)')
print('log_MBH_VP06_CIV_J12::   ',  log_MBH_VP06_CIV_J12,  ' vs. 9.49350 in Shen et al. (2011)')
print()

log_MBH_QSFit_J12_1st  = a_VP06_CIV  + np.log10(( 10**logL_1450_1st /1e44)**b_VP06_CIV) + (2*np.log10(FWHM_CIV_1st))
log_MBH_QSFit_J12_2nd  = a_VP06_CIV  + np.log10(( 10**logL_1450_2nd /1e44)**b_VP06_CIV) + (2*np.log10(FWHM_CIV_2nd))
log_MBH_QSFit_J12_3rd  = a_VP06_CIV  + np.log10(( 10**logL_1450_3rd /1e44)**b_VP06_CIV) + (2*np.log10(FWHM_CIV_3rd))

print
print('log_MBH_QSFit_J12_1st epoch',  log_MBH_QSFit_J12_1st, '  (using VP06)  vs.  9.66  currently in the paper  (', 9.66-log_MBH_QSFit_J12_1st,  'diff)')
print('log_MBH_QSFit_J12_2nd epoch',  log_MBH_QSFit_J12_2nd, '  (using VP06)  vs. 10.08  currently in the paper  (', 10.08-log_MBH_QSFit_J12_2nd, 'diff)')
print('log_MBH_QSFit_J12_3rd epoch',  log_MBH_QSFit_J12_3rd, '  (using VP06)  vs.  9.67  currently in the paper  (',  9.67-log_MBH_QSFit_J12_3rd, 'diff)' )
print()





''' 
  J 1 6 3 8 + 2 8 2 7
'''
logL_3000     = 45.49768513378036
logL_3000_err = 0.0835457146161076
logL_1350     = 45.58550334240817
logL_1350err  = 0.04035432415603424


FWHM_MgII     = 4757.42541816648
FWHM_MgII_err = 2224.264991733375
FWHM_CIV      = 4180.7216796875
FWHM_CIV_err  = 735.8486328125

## Values from QSFit
logL_1450_1st = 45.60852603357719
logL_1450_2nd = 45.40993312333129
logL_1450_3rd = 45.93449845124357
FWHM_CIV_1st = 4630.
FWHM_CIV_2nd = 4990.
FWHM_CIV_3rd = 4620.

#M_BH = a
#np.log10((3.8503777798695754e+45/1e44)**0.53) + (2*np.log10(4180.7216796875))+0.66

log_MBH_MD04_MgII_J12  = a_MD04_MgII + np.log10(( 10**logL_3000 /1e44)**b_MD04_MgII) + (2*np.log10(FWHM_MgII))
log_MBH_VO09_MgII_J12  = a_VO09_MgII + np.log10(( 10**logL_3000 /1e44)**b_VO09_MgII) + (2*np.log10(FWHM_MgII))
log_MBH_S10_MgII_J12   = a_S10_MgII  + np.log10(( 10**logL_3000 /1e44)**b_S10_MgII)  + (2*np.log10(FWHM_MgII))
log_MBH_VP06_CIV_J12   = a_VP06_CIV  + np.log10(( 10**logL_1350 /1e44)**b_VP06_CIV) + (2*np.log10(FWHM_CIV))

print()
print('   J 1 6 3 8 + 2 8 2 7   ')
print()
print('log_MBH_MD04_MgII_J12::  ', log_MBH_MD04_MgII_J12, ' vs. 9.16141 in Shen et al. (2011)')
print('log_MBH_VO09_MgII_J12::  ', log_MBH_VO09_MgII_J12, ' vs. 8.96358 in Shen et al. (2011)')
print('log_MBH_S10_MgII_J12::   ',  log_MBH_S10_MgII_J12,  ' vs. 9.03142 in Shen et al. (2011)')
print('log_MBH_VP06_CIV_J12::   ',  log_MBH_VP06_CIV_J12,  ' vs. 8.74281 in Shen et al. (2011)')
print()

log_MBH_QSFit_J12_1st  = a_VP06_CIV  + np.log10(( 10**logL_1450_1st /1e44)**b_VP06_CIV) + (2*np.log10(FWHM_CIV_1st))
log_MBH_QSFit_J12_2nd  = a_VP06_CIV  + np.log10(( 10**logL_1450_2nd /1e44)**b_VP06_CIV) + (2*np.log10(FWHM_CIV_2nd))
log_MBH_QSFit_J12_3rd  = a_VP06_CIV  + np.log10(( 10**logL_1450_3rd /1e44)**b_VP06_CIV) + (2*np.log10(FWHM_CIV_3rd))

print
print('log_MBH_QSFit_J12_1st epoch',  log_MBH_QSFit_J12_1st, '  (using VP06)  vs.  9.14  currently in the paper  (', 9.14-log_MBH_QSFit_J12_1st, 'diff)')
print('log_MBH_QSFit_J12_2nd epoch',  log_MBH_QSFit_J12_2nd, '  (using VP06)  vs.  9.10  currently in the paper  (', 9.10-log_MBH_QSFit_J12_2nd, 'diff)')
print('log_MBH_QSFit_J12_3rd epoch',  log_MBH_QSFit_J12_3rd, '  (using VP06)  vs.  9.30  currently in the paper  (', 9.30-log_MBH_QSFit_J12_3rd, 'diff)')
print()



''' 
   J 2 2 2 8  +   2 2 0 1 
'''

## Values from QSFit
logL_1450_1st = 44.93449845124357
logL_1450_2nd = 45.961421094066445
logL_1450_3rd = 45.44870631990508
FWHM_CIV_1st = 5930.
FWHM_CIV_2nd = 7000.
FWHM_CIV_3rd = 5930.

print()
print()
print()
print('   J 2 2 2 8  + 2 2 0 1   ')
print()

log_MBH_QSFit_J12_1st  = a_VP06_CIV  + np.log10(( 10**logL_1450_1st /1e44)**b_VP06_CIV) + (2*np.log10(FWHM_CIV_1st))
log_MBH_QSFit_J12_2nd  = a_VP06_CIV  + np.log10(( 10**logL_1450_2nd /1e44)**b_VP06_CIV) + (2*np.log10(FWHM_CIV_2nd))
log_MBH_QSFit_J12_3rd  = a_VP06_CIV  + np.log10(( 10**logL_1450_3rd /1e44)**b_VP06_CIV) + (2*np.log10(FWHM_CIV_3rd))

print
print('log_MBH_QSFit_J12_1st epoch',  log_MBH_QSFit_J12_1st, '  (using VP06)  vs.  9.01  currently in the paper ( ', 9.01-log_MBH_QSFit_J12_1st, 'diff)')
print('log_MBH_QSFit_J12_2nd epoch',  log_MBH_QSFit_J12_2nd, '  (using VP06)  vs.  9.67  currently in the paper ( ', 9.67-log_MBH_QSFit_J12_2nd, 'diff)')
print('log_MBH_QSFit_J12_3rd epoch',  log_MBH_QSFit_J12_3rd, '  (using VP06)  vs.  9.27  currently in the paper ( ', 9.26-log_MBH_QSFit_J12_3rd, 'diff)')
print()
