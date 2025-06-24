import cloud
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

hdul = fits.open("Data/PROMISE/Full_Data.fits")
hdr = hdul[0].header
data = hdul[0].data
distances = np.zeros_like(data)

clouds = cloud.loadclouds("Data/clouds_15kpc.pkl")

keys = []

central_x = hdr["CRPIX1"]
central_y = hdr["CRPIX2"]

for key in clouds:
    keys.append(key)
    if clouds[key].n == 1:
        #print(int(central_x - clouds[key].xc)- int(central_x - clouds[key].xc + clouds[key].delta_x))
        #print(int(central_y - clouds[key].yc)- int(central_y - clouds[key].yc + clouds[key].delta_y))
        #print(np.shape(clouds[key].data))
        mask = clouds[key].data != 0
        mask = mask*clouds[key].mean
        distances[int(central_y - clouds[key].yc): int(central_y - clouds[key].yc + clouds[key].delta_y),int(central_x - clouds[key].xc): int(central_x - clouds[key].xc + clouds[key].delta_x)] = mask
    elif clouds[key].n == 0:
        pass
    elif clouds[key].diff > 0.5:
        mask = clouds[key].data != 0
        mask = mask*clouds[key].mean
        distances[int(central_y - clouds[key].yc): int(central_y - clouds[key].yc + clouds[key].delta_y),int(central_x - clouds[key].xc): int(central_x - clouds[key].xc + clouds[key].delta_x)] = mask

hdu_new = fits.PrimaryHDU(distances, hdr)
hdul_new = fits.HDUList([hdu_new])
hdu_new.writeto(fr"Data/dist_test_15_diff.fits", overwrite = True)