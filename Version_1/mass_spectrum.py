import numpy as np
import cloud
import scipy as sc
import skimage as sk
import matplotlib.pyplot as plt

def quadratic(x,A,B,C):
    return A*x**2 + B*x + C

def gaussian(x,mu,sigma,A):
    return A*np.exp(-(x-mu)**2/(2*sigma**2))

clouds = cloud.loadclouds("Data/clouds.pkl") # Loads the clouds
masses = []

for key in clouds: # Loops through all the cloud objects
    cl = clouds[key]
    distances = np.unique(cl.distance_matrix)[1:] # Finds all unique distances in the cloud object

    full_filtered_data = np.zeros_like(cl.data)
    for d in distances:
        sigma = 1e-4/(2*d*np.tan(cl.db*np.pi/360))
        pixel_area = (2*d*np.tan(cl.db*np.pi/360)*1000)**2
        if sigma < 2:
            continue

        mask = (cl.distance_matrix == d)*(cl.confidence_matrix >= 3) # Creates a mask that only includes the distance we are currently looking at and areas with sufficient confidence class
        data = cl.data*mask # Extracts the relevant data from the cloud data
        data_mass = data*pixel_area*((3.086e18)**2)*9.4e20*2*sc.constants.m_p*1.4/(2e30) # Creates a mass matrix where each pixel value is the mass of the pixel in solar masses
        smoothed_data = sc.ndimage.gaussian_filter(data, sigma) # Smoothes the data
        filtered_data = data - smoothed_data # Removes the smoothed data, in effect creating a high-pass filter keeping only structures that are small enough
        filter_mask = filtered_data > 3 # Applies a threshold that removes everything below
        filtered_data_mask = sk.morphology.area_opening(mask*filter_mask, int(sigma**2)) # Removes areas that are smaller then sigma squared
        big_sized = sk.morphology.area_opening(filtered_data_mask, 4*int(sigma**2)) # Removes areas that are bigger then 4 sigma squared
        final_mask = np.logical_xor(filtered_data_mask,big_sized)
        full_filtered_data = np.where(final_mask, data_mass, full_filtered_data) # Fills out the mask with the actual mass data
        
    labels, num = sc.ndimage.label(full_filtered_data) # Labels each core
    for l in range(num-1): # Iterates through each core
        mass_matrix = np.where(labels == l+1, full_filtered_data, 0)
        mass = np.sum(mass_matrix) # Finds the mass of the core
        if mass < 0.5: # Skips the rest of the loop if core mass is below 0.5 solar masses
            continue
        masses.append(mass) # Appends mass to a list
    print(f"Done with cloud {key[6:]} out of {len(clouds)-1}.")


factor = 2
bins = np.arange(factor*np.max(masses))/factor + 1/2 # Calculates the bin ranges for the histogram
plt.figure(figsize = (16,10))
n, b, p = plt.hist(masses, bins) # Creates the histogram
xlist = bins[:int(3.5*factor)] + 1/(2*factor)
ylist = n[:len(xlist)]
param_q, cov = sc.optimize.curve_fit(quadratic, xlist, ylist, p0 = [-200,1200,-1460], bounds = ([-np.inf,0,-np.inf],[0,np.inf,np.inf])) # Makes a quadratic approximation to the earlier part of the histogram
param_g, cov = sc.optimize.curve_fit(gaussian, xlist, ylist, p0 = [2, 0.5, 340], bounds = ([0,-np.inf, 0],[5,np.inf, np.inf])) # Makes a gaussian approximation to the earlier part of the histogram
x = np.linspace(0,np.max(masses),1000)
plt.xlabel(r"Mass [M$_{\odot}$]")
plt.ylabel("Count")
plt.ylim(0, np.max(n)+20)
plt.xlim(0, len(n)/factor + 1)
plt.savefig("Analysis/Images/histogram.png", dpi = 1000, bbox_inches="tight") # Saves an image with just the histogram
plt.plot(x, quadratic(x, param_q[0], param_q[1], param_q[2]), color = "orange", label = "Quadratic approximation") # Plots the quadratic approximation
plt.plot(x, gaussian(x, param_g[0], param_g[1], param_g[2]), color = "green", label = "Gaussian approximation") # Plots the gaussian approximation
plt.legend()
plt.savefig("Analysis/Images/histogram_approx.png", dpi = 1000, bbox_inches="tight") # Saves an image with the approximations as well
