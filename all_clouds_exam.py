import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import os

directory = "Data/PROMISE/Clouds"

slice = 0
cloud_index = 3
hits = []
variances = []
means = []
i = 0
print(os.listdir(directory))
for file in os.listdir(directory):
    #hdul = fits.open(f"Data/PROMISE/Clouds/slice_{slice}_cloud_{cloud_index}.fits")
    hdul = fits.open(os.path.join(directory, file))
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
    hits.append(len(data[:,0]))
    variances.append(np.var(data[:,4]))
    means.append(np.mean(data[:,4]))
    print(i)
    i += 1

#print(hits)
clearzeros = np.where(np.array(hits) == 0)[0]
print(len(clearzeros))
#variances = np.delete(variances, clearzeros)
#hits = np.delete(hits, clearzeros)
#means = np.delete(means, clearzeros)

lookat = np.where(np.all((np.array(variances) > 10,np.array(hits) > 10), axis = 0))[0]
lookat = np.where(np.all((np.array(hits) < 3,np.array(hits) > 0, np.array(variances) > 0.1), axis = 0))[0]
print(len(lookat))
print(np.array(os.listdir(directory))[lookat])
#ax, fig = plt.sub
plt.figure(figsize = (10,10))
plt.hist(hits, bins = np.arange(20)-0.5)
plt.show()

plt.figure(figsize = (10,10))
plt.hist(variances)
plt.show()
plt.figure(figsize = (10,10))
plt.hist(means)
plt.show()
plt.figure(figsize = (10,10))
plt.hist(np.divide(variances,means))
plt.show()