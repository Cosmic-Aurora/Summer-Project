import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

milky_way_map = np.loadtxt("Data/finalmap.dat",usecols = (0,1,2,3)) # Load in distance data files
useless = np.where(np.any((milky_way_map[:,2] <=0., milky_way_map[:,2] >= 10000.,milky_way_map[:,0] < 0, milky_way_map[:,0] > 40, np.abs(milky_way_map[:,1]) > 1.16), axis = 0))[0] # Find useless distance row indices (No measurement or outside of PROMISE range)
useful_map = np.delete(milky_way_map, useless, axis = 0) # Remove useless rows

hdul = fits.open("Product/distance_mask.fits") # Open distance and confidence data
hdr = hdul[0].header
data = hdul[0].data[0]
confidence = hdul[0].data[1]

hdul = fits.open("Data/Full_Data.fits") # Open full PROMISE data file
magnitude = hdul[0].data

conf = 3 # Lowest confidence level distance to be included

mask = confidence >= conf # Removes all distance data with too low confidence
distances = data*mask

scaling_factors = np.array([2,4,6,8]) # Resolution of inmages to row[1]e created (value x will result in a resolution of 10x*6.5x)
topdown = []
map = []

for s in scaling_factors:
    topdown.append(np.zeros((10*s, int(6.5*s))))
    map.append(np.zeros((10*s, int(6.5*s))))

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

for row in useful_map:
    for k,s in enumerate(scaling_factors): # Adds the magnitude to the correct pixel in the top-down image for each scaling factor
        flattened_distance = row[2]*np.cos(row[1]*np.pi/180)*s/1000
        y = int(flattened_distance*np.cos(row[0]*np.pi/180))
        x = int(int(6.5*s)-1-flattened_distance*np.sin(row[0]*np.pi/180))
        map[k][y,x] += row[3]

for i, t in enumerate(topdown): # Creates and saves an image for each scaling factor
    fig, ax = plt.subplots(1,3,figsize = (16,16))
    ax[0].imshow(np.log(t+1), origin="lower", extent=(-6.5,0,0,10))
    ax[1].imshow(map[i], origin = "lower", extent = (-6.5,0,0,10), cmap = "autumn_r")
    ax[1].imshow(np.log(t+1), origin="lower", extent=(-6.5,0,0,10), alpha = 0.5)
    ax[2].imshow(map[i], origin = "lower", extent = (-6.5,0,0,10), cmap = "autumn_r")
    plt.savefig(f"Product/Images/milky_way_{10*scaling_factors[i]}_px_comparison.png", dpi = 1000, bbox_inches="tight")