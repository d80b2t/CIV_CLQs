'''
Calculate the Schwarzschild radius, R_S, that marks the event horizon
of a nonrotating black hole: 
    R_{s} = 2 r_{g} =  2 GM/c^{2}  


'''

import numpy as np

import astropy.units as u     
from astropy import constants as const


m_M87 = 6.5e9   ## from the EHT papers

## From EHT Paper I::
## ``When viewed from infinity, a nonrotating Schwarzschild (1916) black
## hole has a photon capture radius Rc = sqrt(27) rg, where rg = GM/c2 is
## the characteristic lengthscale of a black hole. The photon capture
## radius is larger than the Schwarzschild radius RS that marks the
## event horizon of a nonrotating black hole, RS ≡ 2 rg.

rs =  (const.G * m_M87 * const.M_sun) / (const.c * const.c)
rs_in_cm = rs*1e2


## Eq. 6 from EHT (2019) Paper V:
## LEdd = 4.pi.GM.c.m_p / sigma_T where, sigma_T is the Thomson cross section
#
## Mdot_edd = L_Edd / epsilonn/c^2
##
## emission arises from 5rg 
## mp is the proton mass
## n_e is electron number density

n_e  = 1.0547e11 / (u.meter  * u.meter  * u.meter )   ##  1.05995e11 m^-3

Mdot            =  4. * np.pi * ((5*(rs/2.))**2) * n_e * const.m_p * (const.c / np.sqrt(5))
Mdot_solperyer  = (Mdot.value*3.1556952e7)/2e30
    
## Mdot quoted from Eq. (5) in EHT Paper V
Mdot_quoted = 2.7e-3  ## Msol /yr 

## Equation (6) from EHT Paper V
## epsilon = radidative efficiency 

epsilion = 0.1
Mdot_Edd = (2.2/epsilion) * (m_M87 / 1e9) 


print()
print(' Schwarzschild radius, R_S, of M87 (using EHT data), is :   ',  rs_in_cm / 1e15 , ' x 10^15 cm')
print()
print(' Mdot is: ', Mdot)
print('   which is : ',  Mdot_solperyer,   ' M_sol / yr')
print('   vs.      : ',  Mdot_solperyer/Mdot_quoted)

print
print(' Eddington accretion rate, of M87 is                    :   ',  Mdot_Edd, ' M_sol yr^-1  (assuming  $\epsilon = 0.1$).')
print('   thus Mdot / Mdot_Edd                 is :: ',          (Mdot_solperyer/  Mdot_Edd), ' and ')
print('   eta where eta = log(Mdot / Mdot_Edd) is :: ', np.log10((Mdot_solperyer/  Mdot_Edd)))

print()



''' 
 E [Joules] = M [kg]  . c^2   [m^2 s^-2]
 L [Watts]  = eta  .    Mdot   . c^2
 L          = 0.1   * ((2.7e-3*2e30)/(24*3600*365.25 ))*(9e16)       
'''
