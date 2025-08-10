from astropy.io import fits
import numpy as np
from scipy import ndimage
import cloud

distances = np.loadtxt("Data/distances.dat",usecols = (3,4,8)) # Load in distance data files
useless = np.where(np.any((distances[:,2] <=0, distances[:,2] > 10, distances[:,0] < 0, distances[:,0] > 40, np.abs(distances[:,1]) > 1.16), axis = 0))[0] # Find useless distance row indices (No measurement or outside of PROMISE range)
useful_distances = np.delete(distances, useless, axis = 0) # Remove useless rows

hdul = fits.open("Data/Masked_Data.fits") # Load in masked PROMISE data from FITS file
data = hdul[0].data # Extract the masked data
hdr = hdul[0].header
hdul.close()

sections_vertical = [0]
for x in range(0, np.shape(data)[1]-1): # Finds places to slice the total data file
    if sum(data[:,x]) == 0 and sections_vertical[-1] + 8000 < x:
        sections_vertical.append(x)
sections_vertical.append(120000)

xdists = []
unsorted_clouds = []

for j in range(len(sections_vertical)-1):

    data_slice = data[:,sections_vertical[j]:sections_vertical[j+1]]
    
    if data_slice.dtype.byteorder == ">": # Used to avoid "ValueError: Big-endian buffer not supported on little-endian compiler"
        data_slice = data_slice.view(data_slice.dtype.newbyteorder()).byteswap()

    
    labels, num = ndimage.label(data_slice) # num = number of clouds, labels = labeled matrix the same size as data

    for i in range(1,num + 1):
        cloud_data = np.trim_zeros(np.where(labels == i, data_slice, 0)) # Keep data that is in he cloud mask, remove all other
        
        mask = np.where(labels == i, True, False) # Create a mask for the data in the full file in order to extract the position of the lower left corner of the mask
        argmaxy = np.argmax(mask, axis = 0)
        argmaxx = np.argmax(mask, axis = 1)
        deltay = np.min(argmaxy[argmaxy != 0])
        deltax = np.min(argmaxx[argmaxx != 0]) + sections_vertical[j]
        
        xdists.append(deltax) # Appends the cloud and its distance to the y-axis (in pixels)
        unsorted_clouds.append(cloud.Cloud(cloud_data.copy(), useful_distances, hdr, deltax, deltay)) 
        print(f"Done with {i}/{num} in slice {j}.")


sorted_clouds = [c for _,_, c in sorted(zip(xdists,range(len(xdists)),unsorted_clouds))] # Sorts the cloud list based on the xdist list

clouds = {}

for i, c in enumerate(sorted_clouds):
    clouds[f"cloud_{i}"] = c # Appends the clouds to a dictionary in order

cloud.saveclouds(clouds, "Data/clouds.pkl") # Saves the cloud dictionary