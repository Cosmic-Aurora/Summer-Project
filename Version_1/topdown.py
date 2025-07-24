import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

hdul = fits.open("Product/distance_mask.fits") # Open distance and confidence data
hdr = hdul[0].header
data = hdul[0].data[0]
confidence = hdul[0].data[1]

hdul = fits.open("Data/Full_Data.fits") # Open full PROMISE data file
magnitude = hdul[0].data

conf = 1 # Lowest confidence level distance to be included

mask = confidence >= conf # Removes all distance data with too low confidence
distances = data*mask

scaling_factors = np.array([20,50,100,500]) # Resolution of inmages to be created (value x will result in a resolution of 10x*6.5x)
topdown = []

for s in scaling_factors:
    topdown.append(np.zeros((10*s, int(6.5*s))))

dl = hdr["CD1_1"] # Extract values from header
db = hdr["CD2_2"]

cl = hdr["CRVAL1"]
cb = hdr["CRVAL2"]

cx = hdr["CRPIX1"]
cy = hdr["CRPIX2"]


for i in range(len(distances)): # Iterates through each pixel in the distance map
    for j,p in enumerate(distances[i]):
        if mask[i,j]:
            l = (cl + (j-cx)*dl)*np.pi/180 # Converts pixel indices to galactic coordinates
            b = (cb + (i-cy)*db)*np.pi/180
            for k,s in enumerate(scaling_factors): # Adds the magnitude to the correct pixel in the top-down image for each scaling factor
                flattened_distance = p*np.cos(b)*s
                y = int(flattened_distance*np.cos(l))
                x = int(6.5*s-1-flattened_distance*np.sin(l))
                topdown[k][y,x] += magnitude[i,j]

for i, t in enumerate(topdown): # Creates and saves an image for each scaling factor
    plt.figure(figsize = (16,8))
    plt.imshow(np.log(t+1), origin="lower", extent=(-6.5,0,0,10))
    plt.savefig(f"Product/Images/milky_way_{10*scaling_factors[i]}_px.png", dpi = 1000, bbox_inches="tight")