import numpy as np
from astropy.io import fits
import cloud
import matplotlib.pyplot as plt
import scipy as sc

clouds = cloud.loadclouds("Data/clouds_10kpc.pkl")

cl = clouds["slice_2_cloud_64.fits"]
print(sc.cluster.hierarchy.fclusterdata(cl.distances[:,4].reshape(-1,1),0.5))

medians = []
keys = []
fracs = []

for key in clouds:
    if clouds[key].diff > 0.5 and clouds[key].n >=1:
        keys.append(key)
        inside = np.where(np.abs(clouds[key].median - clouds[key].distances[:,4]) < 0.5)[0]
        frac = len(inside)/len(clouds[key].distances[:,4])
        fracs.append(frac)
print(fracs)
y = []
for i in range(101):
    y.append(sum(np.array(fracs) >= i/100)/len(fracs))
print(y)
x = np.linspace(0,1,101)
plt.plot(x,y)
plt.vlines([1/2,1/3,2/3,1/4,3/4], 0, 1, color = "red")
plt.show()

# # arr = np.array([0, 0, 0, 0, 0, 0, 8, 83, 120, 111, 31, 37, 10, 0, 0, 0, 0, 0, 0, 0])

# # sak = np.matrix([[0,1,0,4],[3,0,2,1],[5,6,0,0]])
# # sak = np.pad(sak,2)
# # print(sak)
# # s = sak != 0
# # print(s)
# # print(np.argmax(s[::-1],axis = 0))
# # print(np.argmax(s[::],axis = 1))


# # m = arr!=0
# # print(m)    
# # print(m.argmax())

# # arr = np.trim_zeros(arr)
# # print(arr)
# import matplotlib.pyplot as plt

# hdul = fits.open("Data/PROMISE/Clouds/slice_0_cloud_1.fits")
# data = np.loadtxt("Data/Distances_HiGAL/Useful.dat")
# c1 = cloud.Cloud(hdul, data)
# #c1.show_data(plt)
# #c1.property_list()
# hdul.close()

# hdul = fits.open("Data/PROMISE/Clouds/slice_0_cloud_2.fits")
# data = np.loadtxt("Data/Distances_HiGAL/Useful.dat")
# c2 = cloud.Cloud(hdul, data)
# #c1.show_data(plt)
# #c2.property_list()
# hdul.close()

# clouds = [c1,c2]
# #print(clouds[0])
# cloud.saveclouds(clouds, "Data/temp.pickle")

# loadedclouds = cloud.loadclouds("Data/temp.pickle")
# #print(loadedclouds[0])

# clouds[0].show_data(plt)
# loadedclouds[0].show_data(plt)