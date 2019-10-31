
import numpy as np 
from astropy.io import fits
from astropy.io import ascii

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.lines  as mlines


##  Hamann et al. (2017) ERQ BOSS DR12 catalog
path       = '/cos_pc19a_npr/data/ERQs/Hamann2017_CIVcatalog/'
filename   = 'C4N5REWs_DR12v11_MNRAS.fits'
infile     = path+filename
data_full  = fits.open(infile)
Ham17_full = data_full[1].data

## knocking out a couple of objects with bad REW values 
#Ham17 = Ham17_full[np.where( (Ham17['rew'] >0.) & (Ham17['rew'] < 10000.)  )]
Ham17 = Ham17_full


##  Kozłowski_2017_ApJS_228_9.  BOSS DR12 "Value Added" catalog
path       = '/cos_pc19a_npr/programs/quasars/CIV_CLQs/data/QSO_CIV_catalogs/'
filename   = 'SDSS-DR12Q-BH_extra.fits'
infile     = path+filename
data_full  = fits.open(infile)
Koz17_full = data_full[1].data

Koz17      = Koz17_full



file1 = open("temp.txt","w+") 

    
#for ii in range(len(Koz17)):
for ii in range(100):
    print(ii)

    name        = Ham17[np.where((Koz17['SDSS_NAME'][ii]  == Ham17['SDSS_NAME']))]['SDSS_NAME']
    REW         = Ham17[np.where((Koz17['SDSS_NAME'][ii]  == Ham17['SDSS_NAME']))]['rew']
    bal_flag_vi = Ham17[np.where((Koz17['SDSS_NAME'][ii]  == Ham17['SDSS_NAME']))]['bal_flag_vi']
    f1450       = Ham17[np.where((Koz17['SDSS_NAME'][ii]  == Ham17['SDSS_NAME']))]['f1450']
    
    if (len(name) > 0 and  bal_flag_vi <1) :
        #print(ii, Koz17['SDSS_NAME'][ii], name, REW, f1450)
        file1.write(str(ii)+str(Koz17['SDSS_NAME'][ii])+str(we_name)+str(we_rew))
        file1.write(" {} {} {}  {} {} {}  {} {} {}  \n".format(ii, name, Koz17['RA'][ii], Koz17['DEC'][ii],
                                                    #bal_flag_vi,
                                                    REW, f1450, 
                                                    Koz17['L1350'][ii], Koz17['LBol'][ii], #Koz17['eLBol'][ii],
                                                    Koz17['nEdd'][ii]))
file1.close() 
