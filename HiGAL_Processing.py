import numpy as np

#HiGAL columns representation
#0 - Right ascention
#1 - Declination
#2 - Galactic longitude
#3 - Galactic latitude
#4 - Distance (kpc)
#5 - Distance uncertainty (kpc)

data = np.loadtxt(r"Data/Distances_HiGAL/table3.dat",usecols = (1,2,3,4,8,9))

remove = np.where(np.any((data[:,4] == -999., data[:,2] < 0, data[:,2] > 40, np.abs(data[:,3] > 1.16), data[:,4] > 15), axis = 0))[0] #Extracts which data points are useless

useful_data = np.delete(data, remove, axis = 0)

np.savetxt("Data/Distances_HiGAL/Useful.dat", useful_data)