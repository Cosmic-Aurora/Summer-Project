import numpy as np

#HiGAL columns representation
#0 - Right ascention
#1 - Declination
#2 - Galactic longitude
#3 - Galactic latitude
#4 - Distance (kpc)
#5 - Distance uncertainty (kpc)

data = np.loadtxt(r"Data\Distances HiGAL\table3.dat",usecols = (1,2,3,4,8,9))

remove = np.where(data[:,4] == -999. or data[:,2] < 0 or data[:,2] > 40 or np.abs(data[:,3] > 1.2))[0] #Extracts which data points are useless