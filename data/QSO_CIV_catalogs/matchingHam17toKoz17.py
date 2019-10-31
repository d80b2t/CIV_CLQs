
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

Koz17_df    = DataFrame(np.array(Koz17_table))
print('len(Koz17_df)', len(Koz17_df))

## COlumns names are case senstive!
Koz17_df.rename(columns={'SDSS_NAME':'sdss_name'}, inplace=True)              

## Testing on a wee bit of the DataFrame...
mini_merge = pd.merge(Ham17_df[0:100], Koz17_df[0:100], on="sdss_name")

The_DF = pd.merge(Ham17_df, Koz17_df, on="sdss_name")
