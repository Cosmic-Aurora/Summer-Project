import cloud
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sklearn as sk
from voronoi import voronoi

hdul = fits.open("Data/PROMISE/Full_Data.fits")
hdr = hdul[0].header
data = hdul[0].data
distances = np.zeros_like(data)

clouds = cloud.loadclouds("Data/clouds_10kpc.pkl")

keys = []

central_x = hdr["CRPIX1"]
central_y = hdr["CRPIX2"]
distance_divide = 0.5
count1 = 0
count2 = 0

for key in clouds:
    keys.append(key)
    if clouds[key].n == 1:
        #print(int(central_x - clouds[key].xc)- int(central_x - clouds[key].xc + clouds[key].delta_x))
        #print(int(central_y - clouds[key].yc)- int(central_y - clouds[key].yc + clouds[key].delta_y))
        #print(np.shape(clouds[key].data))
        mask = clouds[key].data != 0
        mask = mask*clouds[key].mean
        count1 += 1
        distances[int(central_y - clouds[key].yc): int(central_y - clouds[key].yc + clouds[key].delta_y),int(central_x - clouds[key].xc): int(central_x - clouds[key].xc + clouds[key].delta_x)] = mask
        continue
    elif clouds[key].n == 0:
        continue
    

    if clouds[key].n >= 10:
        dist = np.ma.array(clouds[key].distances[:,4], mask = False)
        remove = False
        for i in range(len(dist)):
            dist.mask[i] = True
            if np.abs(np.max(dist)-np.min(dist)-clouds[key].diff) > 3:
                if remove != False:
                    remove = False
                    continue
                remove = i
            dist.mask[i] = False

        if remove != False:
            clouds[key].remove_row(remove)


    if clouds[key].diff < distance_divide:
        mask = clouds[key].data != 0
        mask = mask*clouds[key].mean
        count2 += 1
        distances[int(central_y - clouds[key].yc): int(central_y - clouds[key].yc + clouds[key].delta_y),int(central_x - clouds[key].xc): int(central_x - clouds[key].xc + clouds[key].delta_x)] = mask
        continue
    points = clouds[key].distances[:,:2]
    data = clouds[key].data
    vor = voronoi(points, data)



print(count1, count2)
hdu_new = fits.PrimaryHDU(distances, hdr)
hdul_new = fits.HDUList([hdu_new])
hdu_new.writeto(fr"Data/dist_test_10.fits", overwrite = True)