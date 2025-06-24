import numpy as np
from astropy.io import fits
import os
import re
import cloud

def sorted_alphanumeric(data):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(data, key=alphanum_key)

directory = "Data/PROMISE/Clouds"
data = np.loadtxt("Data/Distances_HiGAL/promise_slice.dat")

distance = 5

remove = np.where(data[:,2] > distance)[0]
pruned_data = np.delete(data, remove, axis = 0)

clouds = {}

filepaths = sorted_alphanumeric(os.listdir(directory))
for file in filepaths:
    hdul = fits.open(os.path.join(directory, file))
    c = cloud.Cloud(hdul, pruned_data)
    clouds[file] = c

cloud.saveclouds(clouds, f"Data/clouds_{distance}kpc.pkl")