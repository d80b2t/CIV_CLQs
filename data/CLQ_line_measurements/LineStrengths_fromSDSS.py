''' 
Super simple program to take in e.g. (emission, line) fluxes, that
is the Flux density in units of e.g W/m2 or erg/s/cm2 (dimensions of M
T−^-3) and convert, using a cosmological Luminosity Distance to just
e.g. line luminosity
 
https://en.wikipedia.org/wiki/Luminosity
http://www.astro.ucla.edu/~wright/cosmo_02.htm#DL
'''
import numpy as np
from astropy.cosmology import FlatLambdaCDM


## Setting up the cosmology
cosmo = FlatLambdaCDM(H0=71, Om0=0.27, Tcmb0=2.725)

print()
print('Setting up the cosmology')
print(cosmo)
print()

lineFlux = float(input ("Enter line flux, [units of 10-17 erg/cm2/s] :  "))
redshift = float(input ("Enter redshift                              :  "))
print()

## https://docs.astropy.org/en/stable/cosmology/
DLum = cosmo.luminosity_distance(redshift)

## Converting (Mega)parsecs to centimeters
Mpc_to_cm = 3.08567782e18 * 1e6

## Lum = 4.pi.R^2.flux
lineLum = (((DLum.value * Mpc_to_cm)**2) * (lineFlux * 1e-17) * (4*np.pi)) / 1e42
    
print()
print('The Line Luminosity, in units of 10^42 erg/s is:    ', lineLum)
print()
