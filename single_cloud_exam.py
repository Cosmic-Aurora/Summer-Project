import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt



slice = 3
cloud_index = 342


hdul = fits.open(f"Data/PROMISE/Clouds/slice_{slice}_cloud_{cloud_index}.fits")
cloud = hdul[0].data
hdr = hdul[0].header

data = np.loadtxt("Data/Distances_HiGAL/Useful.dat")

lc = hdr["CRVAL1"]
dl = hdr["CD1_1"]
x0 = hdr["CRPIX1"]
lx = hdr["NAXIS1"]

bc = hdr["CRVAL2"]
db = hdr["CD2_2"]
y0 = hdr["CRPIX2"]
by = hdr["NAXIS2"]


l0 = lc - dl*x0
b1 = bc - db*y0

l1 = l0 + lx*dl
b0 = b1 + by*db

remove = np.where(np.any((data[:,2] < min(l0,l1), data[:,2] > max(l0,l1), data[:,3] < min(b0,b1), data[:,3] > max(b0,b1)), axis = 0))
pruned_data = np.delete(data, remove, axis = 0)

l = pruned_data[:,2]
b = pruned_data[:,3]

remove2 = np.where(cloud[(y0 + (b-bc)/db).astype(int), (x0 + (l-lc)/dl).astype(int)] == 0)

data = np.delete(pruned_data, remove2, axis = 0)

print(f"Mean is {np.mean(data[:,4])}")
print(f"Variance is {np.var(data[:,4])}")

#plt.scatter(pruned_data[:,2],pruned_data[:,3], s = 30, color = "red")
plt.figure(figsize = (10,10))
plt.imshow(cloud, extent = (l0,l1,b0,b1), aspect = "equal")
#plt.show()
# print(np.shape(cloud))
# print(x0,(l-lc)/dl ,x0 + (l-lc)/dl)
# print(y0 + (b-bc)/db)


plt.scatter(data[:,2],data[:,3], s = 30, c = data[:,4])
plt.colorbar()
plt.show()
plt.figure(figsize = (10,10))
plt.hist(data[:,4], bins = np.arange(11))
plt.show()