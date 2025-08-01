# ID, l, b (of highest extinction value within the cloud), Area (physical area), distance, flag, mass (sum over pixel magnitudes*area of each pixel (use distance)) (assume square pixels, small enough angels that b and l dont matter), # of points

# Comparing topdown to other images (sent by Jouni)
import cloud
import astrodendro as ad
from methods import *
import scipy as sc
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

clouds = cloud.loadclouds("Data/clouds.pkl")

cl = clouds["cloud_0"]

table = np.zeros(8)
counter = 1

points = cl.distances[:,:2] # Extracts pixel position of distance mearuements within cloud
distances = cl.distances[:,4] # Extracts distance values
data = cl.data # Exctracts data image
mask = data != 0 # Creates a mask of the data
g_data = sc.ndimage.gaussian_filter(data,5)*mask # Applies a gaussian smoothing to the data
d = ad.Dendrogram.compute(g_data, min_delta = 0.1, min_value = 0.1, is_independent = ad.pruning.contains_seeds(points[:,::-1].T.astype(int))) # Computes a dendrogram for the smoothed image

distance_matrix, confidence_matrix = dendro_to_distance(d, mask, points, distances) # Computed the distance and confidence matrix of the cloud (see methods.py for more detail)
id_matrix, table, counter = id_assignment(np.zeros_like(confidence_matrix), distance_matrix, confidence_matrix, table, points, counter, cl)

np.savetxt("Product/temp.txt", table[1:],fmt="%s", delimiter="\t") # Saves the matrix

df = pd.DataFrame(table[1:], columns = ["ID", "l", "b", "Area [pc^2]", "Distance", "Flag", "Mass", "# of distance points"])

max_lengths = df.map(str).map(len).max()
max_lengths = np.max(max_lengths, df.columns.str.len(), axis = 0)

df.columns = df.columns.str.pad(max_lengths, side='right')

df = df.apply(lambda col: col.str.pad(max_lengths[col.name], side='right'))

df.to_csv("Product/temp.txt", index = False, sep = "\t")