from astropy.io import fits
import numpy as np

hdul = fits.open("Data/PROMISE/Full_Data.fits")
hdul.info()
hdr = hdul[0].header.copy()
print(repr(hdr))
data = hdul[0].data
print(np.shape(data))