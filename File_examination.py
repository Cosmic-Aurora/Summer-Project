import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#HiGAL columns representation
#0 - Right ascention
#1 - Declination
#2 - Galactic longitude
#3 - Galactic latitude
#4 - Distance (kpc)
#5 - Distance uncertainty (kpc)

HiGAL5 = np.loadtxt(r"Data\Distances HiGAL\table3.dat", usecols = (1,2,3,4,8,9)) #Loads Data
HiGAL3 = np.loadtxt(r"C:\Users\aurora\Documents\Sommarprojekt 2025\Data\HiGAL\tablee1b.dat\table3b.dat", usecols = (1,2,3,4,8,9)) 
HiGALh = np.loadtxt(r"C:\Users\aurora\Documents\Sommarprojekt 2025\Data\HiGAL\hcatalog.dat\hcatalog.dat", usecols = (4,5,2,3,57))
HiGALl = np.loadtxt(r"C:\Users\aurora\Documents\Sommarprojekt 2025\Data\HiGAL\lcatalog.dat\lcatalog.dat", usecols = (4,5,2,3,57))

remove5 = np.where(HiGAL5[:,4] == -999.)[0] #Extracts which data points are useless
remove3 = np.where(HiGAL3[:,4] == -999.)[0]
removeh = np.where(HiGALh[:,4] == 0.)[0]
removel = np.where(HiGALl[:,4] == 0.)[0]

HiGALdat5 = np.delete(HiGAL5, remove5, axis = 0) #Removes those data points
HiGALdat3 = np.delete(HiGAL3, remove3, axis = 0)
HiGALdath = np.delete(HiGALh, removeh, axis = 0)
HiGALdatl = np.delete(HiGALl, removel, axis = 0)

remove2 = np.where(HiGALdat5[:,4] > 33)
HiGAL52 = np.delete(HiGALdat5, remove2, axis = 0)

print("HiGAL5:",len(remove5),np.shape(HiGAL5), np.shape(HiGALdat5))
print("HiGAL3:",len(remove3),np.shape(HiGAL3), np.shape(HiGALdat3))
print("HiGALh:",len(removeh), np.shape(HiGALh), np.shape(HiGALdath))
print("HiGALl:",len(removel), np.shape(HiGALl), np.shape(HiGALdatl))

data = HiGALdath
def polartocathartic(r,theta):
    return -r*np.sin(theta), r*np.cos(theta)

r = data[:,4]
truer = r*np.cos(data[:,3])
theta = data[:,2]*np.pi/180
x,y = polartocathartic(truer,theta)

length = int(np.max(truer)/2)
liner = np.arange(0,length)
theta1 = np.ones(length)*1.5*np.pi/180
theta2 = np.ones(length)*40*np.pi/180
x1,y1 = polartocathartic(liner, theta1)
x2,y2 = polartocathartic(liner, theta2)


plt.figure(figsize = (16,16))
#plt.axis('square')
plt.scatter(x,y, s = 0.1)
plt.plot(x1,y1, color = "red")
plt.plot(x2,y2, color = "red")
#plt.xlim(-1000,1000)
#plt.ylim(-1000,1000)

plt.grid()
plt.savefig(r"C:\Users\aurora\Downloads\temp.png")
plt.show()
#plt.scatter(r,theta, s = 1)
#plt.show()
