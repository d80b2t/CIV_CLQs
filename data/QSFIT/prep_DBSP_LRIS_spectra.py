'''
Just a wee script that takes e.g. the DBSP and LRIS data 
from the directories above and converts things into 
columns of wavelength, flux and flux_error with 
the flux units being 10^-17 erg s^-1 cm^-2 Ang^-1. 
https://qsfit.inaf.it/cat_1.30/onlinefit.php
'''

import numpy as np
from astropy.io    import ascii
from astropy.table import Table, vstack

##
## R E A D I N G   I N   T H E    D A T A 
##
path = '/cos_pc19a_npr/programs/quasars/CIV_CLQs/data/'


##      S P E C T R A
##
##  J1205+3422
#object_name = 'J1205+3422/'
#infile = 'DBSP_J1205p3422_b_58538.dat'
#infile = 'DBSP_J1205p3422_b_58693.dat'
#blue   = ascii.read(path+object_name+infile)
#infile = 'DBSP_J1205p3422_r_58538.dat'
#infile = 'DBSP_J1205p3422_r_58693.dat'
#red    = ascii.read(path+object_name+infile)
##
##  J1638+2827
#object_name = 'J1638+2827/'
#infile = 'LRIS_J1638p2827_b_58583.dat'
#blue   = ascii.read(path+object_name+infile)
#infile = 'LRIS_J1638p2827_r_58583.dat'
#red    = ascii.read(path+object_name+infile)
##
##  J2228+3422
object_name = 'J2228+2201/'
infile = 'DBSP_J2228p2201_b_58693.dat'
blue   = ascii.read(path+object_name+infile)
infile = 'DBSP_J2228p2201_r_58693.dat'
red    = ascii.read(path+object_name+infile)


print()
print('Read-in', object_name+infile)
print()


wavelength_blue = blue['wavelength']
wavelength_red  =  red['wavelength']


## Need to convert BUNIT   = 'erg/cm2/s/Hz'  to  'erg/cm2/s/A'
## http://www.stsci.edu/~strolger/docs/UNITS.txt
## 
## [Y erg/cm^2/s/Hz]            = 1000 * [X W/m^2/Hz]
## [Y erg/cm^2/s/A]             = 2.99792458E+21 * [X1 W/m^2/Hz] / [X2 A]^2
## / 1e-17 to put into 'regular' SDSS flux units 
flux_blue = (2.99792458E+21 * (blue['flux_density'] / 1000.)) / (blue['wavelength']**2) / 1e-17
flux_red  = (2.99792458E+21 * ( red['flux_density'] / 1000.)) / ( red['wavelength']**2) / 1e-17


## We don’t have flux uncertainties but the SNR
## is typically between 10 - 30.
flux_error_blue = flux_blue / 10.
flux_error_red  =  flux_red / 10.


## basically use all the Blue arm info;
## cut out Red arm info that have Blue arm coverage
good_wave = np.where(wavelength_red > wavelength_blue.max())

wavelength = vstack([wavelength_blue, wavelength_red[good_wave]])
flux       = vstack([flux_blue,       flux_red[good_wave]])             
err_flux   = vstack([flux_error_blue, flux_error_red[good_wave]])

#output = np.array(Table([wavelength,flux,err_flux], names=['wavelength', 'flux', 'err_flux']))
output = Table( [wavelength,flux,err_flux], names=['wavelength', 'flux', 'err_flux'])


ascii.write(output, 'spec_output_temp.dat',  overwrite=True) # names=['wavelength', 'flux', 'err_flux'])

