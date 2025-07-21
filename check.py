from astropy.io import fits
import numpy as np

hdul = fits.open("Data/dist_10_v1.fits")
data = hdul[0].data

print(np.shape(data))
print(np.sum(data[0]!=0), np.sum(data[1]!=0))