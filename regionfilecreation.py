import numpy as np
from astropy.io import fits

hdul = fits.open(r"Data/PROMISE/Masked_Data.fits")
hdr = hdul[0].header

l0 = hdr["CRVAL1"]
x0 = hdr["CRPIX1"]
dl = hdr["CD1_1"]

b0 = hdr["CRVAL2"]
y0 = hdr["CRPIX2"]
db = hdr["CD2_2"]


data = np.loadtxt(r"Data/Distances_HiGAL/Useful.dat")
lista = []
for i in range(len(data[:,0])):
    new_element = f"circle({x0+(data[i,2]-l0)/dl},{y0+(data[i,3]-b0)/db},30)"
    lista.append(new_element)



sparalista = np.array(lista, dtype = str)
np.savetxt("Data/region_image.reg", sparalista, fmt = "%s")