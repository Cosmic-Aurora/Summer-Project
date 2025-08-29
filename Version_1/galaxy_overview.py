import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy import constants

compare = True

try:
    milky_way_map = np.loadtxt("Data/finalmap.dat",usecols = (0,1,2,3)) # Load in distance data files
    useless = np.where(np.any((milky_way_map[:,2] <=0., milky_way_map[:,2] >= 10000.,milky_way_map[:,0] < 0, milky_way_map[:,0] > 40, np.abs(milky_way_map[:,1]) > 1.16), axis = 0))[0] # Find useless distance row indices (No measurement or outside of PROMISE range)
    useful_map = np.delete(milky_way_map, useless, axis = 0) # Remove useless rows
except:
    compare = False

hdul = fits.open("Product/distance_mask.fits") # Open distance and confidence data
hdr = hdul[0].header
data = hdul[0].data[0]
confidence = hdul[0].data[1]

hdul = fits.open("Data/Masked_Data.fits") # Open full PROMISE data file
magnitude = hdul[0].data

conf = 3 # Lowest confidence level distance to be included

mask = confidence >= conf # Removes all distance data with too low confidence
distances = data*mask

scaling_factors = np.array([6,10,20,40]) # Modify for different resolutions

comp_scaling_factor = 8 # Resolution of inmages to row[1]e created (value x will result in a resolution of 10x*6.5x)

topdown = []
map = np.zeros((10*comp_scaling_factor, int(6.5*comp_scaling_factor)))

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
                pixel_area = 4*(p*1000)**2*np.abs(np.tan(dl*np.pi/360)*np.tan(db*np.pi/360))
                pixel_mass = pixel_area*magnitude[i,j]*((3.086e18)**2)*9.4e20*2*constants.m_p*1.4/(2e30)
                column_density = pixel_mass*s**2/1e6 # Calculates the column density that the pixel adds in solar masses/pc^2
                topdown[k][y,x] += column_density
if compare:
    for row in useful_map:
        flattened_distance = row[2]*np.cos(row[1]*np.pi/180)*comp_scaling_factor/1000
        y = int(flattened_distance*np.cos(row[0]*np.pi/180))
        x = int(int(6.5*comp_scaling_factor)-1-flattened_distance*np.sin(row[0]*np.pi/180))
        map[y,x] += row[3]




for i, t in enumerate(topdown): # Creates and saves an two images for each scaling factor, one with the comparison and one without
    plt.figure(figsize = (16,10))
    plt.imshow(np.log(t), origin="lower", extent=(-6.5,0,0,10))
    plt.colorbar(label = r"Logarithmized column density [M$_{\odot}$/pc$^2$]")
    plt.xlabel("kpc")
    plt.ylabel("kpc")
    plt.savefig(f"Analysis/Images/milky_way_{10*scaling_factors[i]}_px.png", dpi = 1000, bbox_inches="tight")
    
    if compare:
        plt.figure(figsize=(16,10))
        plt.imshow(np.log(t), origin="lower", extent=(-6.5,0,0,10))
        plt.colorbar(label = r"Logarithmized column density [M$_{\odot}$/pc$^2$]")
        plt.xlabel("kpc")
        plt.ylabel("kpc")
        plt.contour(map, origin = "lower", extent = (-6.5,0,0,10), cmap = "autumn_r")
        plt.savefig(f"Analysis/Images/milky_way_{10*scaling_factors[i]}_px_comp.png", dpi = 1000, bbox_inches="tight")
