import numpy as np
import matplotlib.pyplot as plt

distances = np.loadtxt("Data/finalmap.dat",usecols = (0,1,2,3)) # Load in distance data files
useless = np.where(np.any((distances[:,2] <=0., distances[:,2] >= 10000., distances[:,0] < 0, distances[:,0] > 40, np.abs(distances[:,1]) > 1.16), axis = 0))[0] # Find useless distance row indices (No measurement or outside of PROMISE range)
useful_distances = np.delete(distances, useless, axis = 0) # Remove useless rows


scaling_factors = np.array([2,4,6,8]) # Resolution of inmages to row[1]e created (value x will result in a resolution of 10x*6.5x)
topdown = []

for s in scaling_factors:
    topdown.append(np.zeros((10*s, int(6.5*s))))

i = 0
for row in useful_distances:
    for k,s in enumerate(scaling_factors): # Adds the magnitude to the correct pixel in the top-down image for each scaling factor
        flattened_distance = row[2]*np.cos(row[1]*np.pi/180)*s/1000
        y = int(flattened_distance*np.cos(row[0]*np.pi/180))
        x = int(int(6.5*s)-1-flattened_distance*np.sin(row[0]*np.pi/180))
        topdown[k][y,x] += row[3]
    i += 1

for i, t in enumerate(topdown): # Creates and saves an image for each scaling factor
    plt.figure(figsize = (16,8))
    plt.imshow(t, origin="lower", extent=(-6.5,0,0,10), cmap = "autumn_r")
    plt.savefig(f"Product/Images/milky_way_{10*scaling_factors[i]}_px_reference.png", dpi = 1000, bbox_inches="tight")
    
xlist = []
ylist = []

for row in useful_map:
    flattened_distance = row[2]*np.cos(row[1]*np.pi/180)/1000
    ylist.append(flattened_distance*np.cos(row[0]*np.pi/180))
    xlist.append(int(6.5)-flattened_distance*np.sin(row[0]*np.pi/180))

for i, t in enumerate(topdown): # Creates and saves an image for each scaling factor
    fig, ax = plt.subplots(1,3,figsize = (16,16))
    thing = ax[0].imshow(t, origin="lower", extent=(-6.5,0,0,10))
    ax[1].contour(xlist, ylist, useful_map[:,3], origin = "lower", extent = (-6.5,0,0,10))
    ax[1].imshow(t, origin="lower", extent=(-6.5,0,0,10), alpha = 0.5)
    ax[2].contour(xlist, ylist, useful_map[:,3], origin = "lower", extent = (-6.5,0,0,10))
    fig.colorbar(thing, ax = ax[0])
    plt.savefig(f"Images/milky_way_{10*scaling_factors[i]}_px_comp_lin.png", dpi = 1000, bbox_inches="tight")