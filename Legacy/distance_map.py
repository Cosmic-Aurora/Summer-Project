import numpy as np
from astropy.io import fits
import skimage as sk

data = np.loadtxt(r"Data/Distances_HiGAL/Useful.dat")
hdul = fits.open("Data/PROMISE/Masked_Data.fits")
mask = hdul[0].data
hdr = hdul[0].header

l0 = hdr["CRVAL1"]
dl = hdr["CD1_1"]
x0 = hdr["CRPIX1"]

b0 = hdr["CRVAL2"]
db = hdr["CD2_2"]
y0 = hdr["CRPIX2"]

l = data[:,2]
b = data[:,3]

x = x0 + (l-l0)/dl - 1
y = y0 + (b-b0)/db - 1
print(len(x),max(y))
dist = np.zeros_like(mask)
for d in range(len(data[:,0])):
    dist[sk.draw.disk((int(y[d]),int(x[d])),30, shape = (7000, 120000))] = data[d,4]
hdu_new = fits.PrimaryHDU(dist, hdr)
hdul_new = fits.HDUList([hdu_new])
hdu_new.writeto(fr"Data/distances.fits", overwrite = True)