import numpy as np
from astropy.io import fits
import cloud
import matplotlib.pyplot as plt
import scipy as sc
import sklearn as sk

clouds = cloud.loadclouds("Data/clouds_10kpc.pkl")

cl = clouds["slice_2_cloud_64.fits"]
print(cl.n)
#print(sc.cluster.hierarchy.linkage(cl.distances[:,4]))
print(sc.cluster.hierarchy.fclusterdata(cl.distances[:,4].reshape(-1,1),0.031, criterion = "distance"))
print(cl.distances[:,4])

clusters = sk.cluster.DBSCAN(eps = 0.031, min_samples = 2).fit(cl.distances[:,4].reshape(-1,1))
print(clusters.labels_)
print(clusters.core_sample_indices_)


cl.show_data(plt)