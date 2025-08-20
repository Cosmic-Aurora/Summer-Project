import numpy as np
import cloud
import scipy as sc

clouds = cloud.loadclouds("Data/clouds.pkl")

for key in clouds:
    cl = clouds[key]
    distances = np.unique(cl.distance_matrix)[1:] # Finds all unique distances in the cloud object
    for d in distances:
        sigma = 1e-4/(d*np.tan(cl.db/2))
        if sigma < 2:
            continue
        mask = (cl.distance_matrix == d)*(cl.confidence_matrix >= 3) # Creates a mask that only includes the distance we are currently looking at and areas with sufficient confidence class
        data = cl.data*mask # Extracts the relevant data from the cloud data
        smoothed_data = sc.ndimage.gaussian_filter(data, sigma)
        filtered_data = data - smoothed_data
