import numpy as np
from astropy.io import fits
import cloud
import matplotlib.pyplot as plt
import scipy as sc
import sklearn as sk

clouds = cloud.loadclouds("Data/clouds_10kpc.pkl")

cl = clouds["slice_0_cloud_1.fits"]
# print(cl.n)
# #print(sc.cluster.hierarchy.linkage(cl.distances[:,4]))
# print(sc.cluster.hierarchy.fclusterdata(cl.distances[:,4].reshape(-1,1),0.031, criterion = "distance"))
# print(cl.distances[:,4])

# clusters = sk.cluster.DBSCAN(eps = 0.5, min_samples = 1).fit(cl.distances[:,4].reshape(-1,1))
# print(clusters.labels_)
# print(clusters.core_sample_indices_)

distance_divide = 0.5
print(cl.n)
dist = np.ma.array(cl.distances[:,4], mask = False)
remove = False
for i in range(len(dist)):
    dist.mask[i] = True
    if np.abs(np.max(dist)-np.min(dist)-cl.diff) > 3:
        if remove != False:
            remove = False
            continue
        remove = i
    dist.mask[i] = False

if remove != False:
    cl.remove_row(remove)
print(cl.n)
cl.show_data(plt)