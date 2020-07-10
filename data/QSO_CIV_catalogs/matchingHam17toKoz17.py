
import numpy  as np 
import pandas as pd
from   pandas import DataFrame
from astropy.io    import fits
from astropy.io    import ascii
from astropy.table import Table


##  Hamann et al. (2017) ERQ BOSS DR12 catalog
path       = '/cos_pc19a_npr/data/ERQs/Hamann2017_CIVcatalog/'
filename   = 'C4N5REWs_DR12v11_MNRAS.fits'
infile     = path+filename
data_full  = fits.open(infile)
Ham17_full = data_full[1].data

#with fits.open(infile) as data:
#    df_Ham17 = pd.DataFrame(data[0].data)

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


## astropy Table to pandas DataFrame #2804
## https://github.com/astropy/astropy/issues/2804
Ham17_table = Table(Ham17)
Ham17_df    = DataFrame(np.array(Ham17_table))
print('len(Ham17_df)', len(Ham17_df))


Koz17_table = Table(Koz17)

## Have to remove this '5 value' columns in order for the
## DataFrame swith-aroo to work..
Koz17_table.remove_column('PSFFLUX')
Koz17_table.remove_column('IVAR_PSFFLUX')
Koz17_table.remove_column('PSFMAG')
Koz17_table.remove_column('ERR_PSFMAG')


## Some nice DataFrame polish/manipulation
Koz17_df    = DataFrame(np.array(Koz17_table))
print('len(Koz17_df)', len(Koz17_df))

## Columns names are case senstive!
Koz17_df.rename(columns={'SDSS_NAME':'sdss_name'}, inplace=True)              

## Testing on a wee bit of the DataFrame...
mini_merge = pd.merge(Ham17_df[0:100], Koz17_df[0:100], on="sdss_name")

The_DF = pd.merge(Ham17_df, Koz17_df, on="sdss_name")


## Just wanting to write a few things out to a simple text file. 
file1 = open("temp.txt","w+") 
    
#for ii in range(len(Koz17)):     # if the full catalog is wanted!  
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

