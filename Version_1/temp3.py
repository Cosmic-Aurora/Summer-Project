import numpy as np
import cloud
import scipy as sc
import matplotlib.pyplot as plt
import cmath
import skimage as sk
from methods import *



clouds = cloud.loadclouds("Data/clouds.pkl")

cl = clouds["cloud_0"]
distances = np.unique(cl.distance_matrix)[1:] # Finds all unique distances in the cloud object

full_filtered_data = np.zeros_like(cl.data)
fourier_filtered_data = np.zeros_like(cl.data)

for d in distances:
    sigma = 1e-4/(d*np.tan(cl.db*np.pi/360))
    print(sigma)
    if sigma < 2:
        continue

    mask = (cl.distance_matrix == d)*(cl.confidence_matrix >= 3) # Creates a mask that only includes the distance we are currently looking at and areas with sufficient confidence class
    data = cl.data*mask # Extracts the relevant data from the cloud data
    ft_data = np.fft.fft2(data)
    M,N = np.shape(ft_data)
    filter = gaussian_fourier_filter(M,N,sigma)
    ift_data = np.fft.ifft2(ft_data*filter)
    smoothed_data = sc.ndimage.gaussian_filter(data, sigma)
    filtered_data = data - smoothed_data
    filter_mask = filtered_data>3
    fourier_mask = ift_data>3
    filtered_data = filtered_data*filter_mask
    fourier_data = ift_data*fourier_mask
    full_filtered_data = np.where(mask == True, filtered_data, full_filtered_data)
    fourier_filtered_data = np.where(mask == True, fourier_data, fourier_filtered_data)

fig, ax = plt.subplots(1,3,figsize = (24,10))
img1 = ax[0].imshow(cl.data)
img2 = ax[1].imshow(full_filtered_data)
img3 = ax[2].imshow(fourier_filtered_data)
fig.colorbar(img1, ax = ax[0])
fig.colorbar(img2, ax = ax[1])
fig.colorbar(img3, ax = ax[2])
plt.show()