from astropy.io import fits
import numpy as np

hdul = fits.open(r"Data\PROMISE\Clouds\Cloud_1_Slice_0.fits")
hdul.info()
hdr = hdul[0].header.copy()
print(repr(hdr))
data = hdul[0].data
print(np.shape(data))