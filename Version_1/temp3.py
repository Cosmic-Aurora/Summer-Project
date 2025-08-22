import numpy as np
import cloud
import scipy as sc
import matplotlib.pyplot as plt
import cmath
import skimage as sk
from methods import *


def quadratic(x,A,B,C):
    return A*x**2 + B*x + C

clouds = cloud.loadclouds("Data/clouds.pkl")

cl = clouds["cloud_0"]
distances = np.unique(cl.distance_matrix)[1:] # Finds all unique distances in the cloud object

full_filtered_data = np.zeros_like(cl.data)
fourier_filtered_data = np.zeros_like(cl.data)



for d in distances:
    sigma = 1e-4/(2*d*np.tan(cl.db*np.pi/360))
    pixel_area = (2*d*np.tan(cl.db*np.pi/360)*1000)**2
    print(sigma)
    if sigma < 2:
        continue

    mask = (cl.distance_matrix == d)*(cl.confidence_matrix >= 3) # Creates a mask that only includes the distance we are currently looking at and areas with sufficient confidence class
    data = cl.data*mask # Extracts the relevant data from the cloud data
    data_mass = data*pixel_area*((3.086e18)**2)*9.4e20*2*sc.constants.m_p*1.4/(2e30) # Creates a mass matrix where each pixel value is the mass of the pixel in solar masses
    #ft_data = np.fft.fft2(data)
    #M,N = np.shape(ft_data)
    #filter = gaussian_fourier_filter(M,N,sigma)
    #gaussian_matrix = gaussian_kernel(np.zeros_like(cl.data), sigma)
    #plt.imshow(gaussian_matrix)
    #plt.colorbar()
    #plt.show()
    #ift_data = np.fft.ifft2(ft_data*filter)
    smoothed_data = sc.ndimage.gaussian_filter(data, sigma)
    filtered_data = data - smoothed_data
    filter_mask = filtered_data>3
    #fourier_mask = ift_data>3
    filtered_data_mask = sk.morphology.area_opening(mask*filter_mask, int(sigma**2)/2)
    big_sized = sk.morphology.area_opening(filtered_data_mask, 4*int(sigma**2))
    final_mask = np.logical_xor(filtered_data_mask,big_sized)
    #fourier_data = ift_data*fourier_mask
    full_filtered_data = np.where(final_mask, data_mass, full_filtered_data)
    #fourier_filtered_data = np.where(mask == True, fourier_data, fourier_filtered_data)
labels, num = sc.ndimage.label(full_filtered_data)
masses = []
for l in range(num-1):
    mass_matrix = np.where(labels == l+1, full_filtered_data, 0)
    mass = np.sum(mass_matrix)
    if mass < 0.5:
        continue
    masses.append(mass)

factor = 1
bins = np.arange(factor*np.max(masses))/factor + 1/2
n, b, p = plt.hist(masses, bins)
xlist = bins[:-1] + 1/(2*factor)
param, cov = sc.optimize.curve_fit(quadratic, xlist, n)
x = np.linspace(0.5,np.max(masses),1000)
plt.savefig(f"Analysis/Images/histogram_{factor}_1.png")
plt.plot(x, quadratic(x, param[0], param[1], param[2]))
plt.show()


fig, ax = plt.subplots(1,2,figsize = (24,10))
img1 = ax[0].imshow(cl.data)
img2 = ax[1].imshow(full_filtered_data)
#img3 = ax[2].imshow(fourier_filtered_data)
fig.colorbar(img1, ax = ax[0])
fig.colorbar(img2, ax = ax[1])
#fig.colorbar(img3, ax = ax[2])
plt.show()