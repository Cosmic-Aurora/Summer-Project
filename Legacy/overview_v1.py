import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

hdul = fits.open("Data/dist_10_v1.fits")
hdr = hdul[0].header
data = hdul[0].data[0]
confidence = hdul[0].data[1]

hdul = fits.open("Data/PROMISE/Masked_Data.fits")
magnitude = hdul[0].data

conf = 1

mask = confidence >= conf
distances = data*mask

scaling_factors = np.array([20,50,100,500])
topdown = []

for s in scaling_factors:
    topdown.append(np.zeros((10*s, int(6.5*s))))

dl = hdr["CD1_1"]
db = hdr["CD2_2"]

cl = hdr["CRVAL1"]
cb = hdr["CRVAL2"]

cx = hdr["CRPIX1"]
cy = hdr["CRPIX2"]


for i in range(len(distances)):
    for j,p in enumerate(distances[i]):
        if mask[i,j]:
            l = (cl + (j-cx)*dl)*np.pi/180
            b = (cb + (i-cy)*db)*np.pi/180
            for k,s in enumerate(scaling_factors):
                flattened_distance = p*np.cos(b)*s
                y = int(flattened_distance*np.cos(l))
                x = int(6.5*s-1-flattened_distance*np.sin(l))
                topdown[k][y,x] += magnitude[i,j]
    print(i)

for i, t in enumerate(topdown):
    plt.figure(figsize = (16,8))
    plt.imshow(t, origin="lower", extent=(-6.5,0,0,10))
    plt.savefig(f"Data/Figures/milky_way_{10*scaling_factors[i]}_px_masked.png", dpi = 1000, bbox_inches="tight")
    plt.show()