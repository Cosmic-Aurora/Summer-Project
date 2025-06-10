import numpy as np
import matplotlib.pyplot as plt

#HiGAL columns representation
#0 - Right ascention
#1 - Inclination
#2 - Galactic longitude
#3 - Galactic latitude
#4 - Distance (kpc)
#5 - Distance uncertainty (kpc)

HiGAL5 = np.loadtxt(r"C:\Users\aurora\Documents\Sommarprojekt 2025\Data\HiGAL\tablee1.dat\table3.dat", usecols = (1,2,3,4,8,9)) #Loads Data
HiGAL3 = np.loadtxt(r"C:\Users\aurora\Documents\Sommarprojekt 2025\Data\HiGAL\tablee1b.dat\table3b.dat", usecols = (1,2,3,4,8,9)) 


remove5 = np.where(HiGAL5[:,4] == -999.)[0] #Extracts which data points are useless
remove3 = np.where(HiGAL3[:,4] == -999.)[0]

HiGALdat5 = np.delete(HiGAL5, remove5, axis = 0) #Removes those data points
HiGALdat3 = np.delete(HiGAL3, remove5, axis = 0)

print(HiGALdat5, np.shape(HiGALdat5))
print(HiGALdat3, np.shape(HiGALdat3))