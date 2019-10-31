#!/usr/bin/env python

from astropy.io import ascii
import matplotlib.pyplot as plt 

data_first  = ascii.read("J1205p3422_58538.dat")
data_second = ascii.read("J1205p3422_58693.dat")

plt.plot((data_first['wavelength']/(1+2.071)),  data_first['flux'], color='r')
plt.plot((data_second['wavelength']/(1+2.071)), data_second['flux'],  color='k')
plt.ylim((-10,50))

#plt.savefig("test.png")
plt.ylim((-1,50))
plt.show()  



