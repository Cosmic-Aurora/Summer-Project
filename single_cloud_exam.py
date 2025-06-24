import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import cloud


slice = 5
cloud_index = 280
distance_limit = 15

clouds = cloud.loadclouds(f"Data/clouds_{distance_limit}kpc.pkl")
cl = clouds[f"slice_{slice}_cloud_{cloud_index}.fits"]

print(f"Variance is: {cl.var}")
print(f"Mean is: {cl.mean}")
print(f"Number of datapoints: {cl.n}")
print(f"Absolute difference: {cl.diff}")

#cl.property_list()

#plt.scatter(pruned_data[:,2],pruned_data[:,3], s = 30, color = "red")
plt.figure(figsize = (10,8))
plt.imshow(cl.data, extent = (cl.l0,cl.l1,cl.b1,cl.b0), aspect = "equal", cmap = "YlGn_r")
#plt.show()
# print(np.shape(cloud))
# print(x0,(l-lc)/dl ,x0 + (l-lc)/dl)
# print(y0 + (b-bc)/db)


plt.scatter(cl.distances[:,2],cl.distances[:,3], s = 30, c = cl.distances[:,4], cmap = "cool")
plt.colorbar()
plt.show()
plt.figure(figsize = (10,8))
plt.hist(cl.distances[:,4], bins = np.arange(11))
plt.show()