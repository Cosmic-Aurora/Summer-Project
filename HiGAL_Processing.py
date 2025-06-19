import numpy as np

#HiGAL columns representation
#0 - Right ascention
#1 - Declination
#2 - Galactic longitude
#3 - Galactic latitude
#4 - Distance (kpc)
#5 - Distance uncertainty (kpc)

data = np.loadtxt(r"Data/Distances_HiGAL/table3.dat",usecols = (3,4,8))

remove = np.where(np.any((data[:,2] == -999., data[:,0] < 0, data[:,0] > 40, np.abs(data[:,1] > 1.16)), axis = 0))[0] #Extracts which data points are useless

useful_data = np.delete(data, remove, axis = 0)

np.savetxt("Data/Distances_HiGAL/promise_slice.dat", useful_data)