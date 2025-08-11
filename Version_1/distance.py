import cloud
import numpy as np
from astropy.io import fits
import scipy as sc
import astrodendro as ad
from methods import *
import pandas as pd

hdul = fits.open("Data/Full_Data.fits") # Opens the PROMISE data (this is used purely to determine the size of the matrix that is to be created, as well as to copy the header and extract some values)
hdr = hdul[0].header
full_data = hdul[0].data
full_distances = np.zeros_like(full_data) # Creates three copies with only zeroes from the data, one for distance data and one for confidence
confidence = np.zeros_like(full_data)
identities = np.zeros_like(full_data)

smoothing_level = 5

clouds = cloud.loadclouds("Data/clouds.pkl") # Opens the cloud dictionary

central_x = hdr["CRPIX1"]
central_y = hdr["CRPIX2"]

counter = 1 # Used to assign ID:s to clouds, increases by 1 every time is assigns an ID, thus making sure each cloud gets a different one
table = np.zeros(8)

for key in clouds: # Iterates through each cloud
    cl = clouds[key]
    if cl.n == 0: # Does nothing if the cloud contains no distance measurements
        print(f"Done with cloud {key[6:]} out of {len(clouds)-1}.")
        continue
    
    points = cl.distances[:,:2] # Extracts pixel position of distance mearuements within cloud
    distances = cl.distances[:,4] # Extracts distance values
    data = cl.data # Exctracts data image
    mask = data != 0 # Creates a mask of the data
    g_data = sc.ndimage.gaussian_filter(data,smoothing_level)*mask # Applies a gaussian smoothing to the data
    d = ad.Dendrogram.compute(g_data, min_delta = 0.1, min_value = 0.1, is_independent = ad.pruning.contains_seeds(points[:,::-1].T.astype(int))) # Computes a dendrogram for the smoothed image

    distance_matrix, confidence_matrix = dendro_to_distance(d, mask, points, distances) # Computed the distance and confidence matrix of the cloud (see methods.py for more detail)
    id_matrix, table, counter = id_assignment(np.zeros_like(confidence_matrix), distance_matrix, confidence_matrix, table, points, counter, cl)
    full_distances[int(central_y - cl.yc): int(central_y - cl.yc + cl.delta_y),int(central_x - cl.xc): int(central_x - cl.xc + cl.delta_x)] += distance_matrix # Adds the cloud distance to the full matrix
    confidence[int(central_y - cl.yc): int(central_y - cl.yc + cl.delta_y),int(central_x - cl.xc): int(central_x - cl.xc + cl.delta_x)] += confidence_matrix # Adds the cloud confidence to the full matrix
    identities[int(central_y - cl.yc): int(central_y - cl.yc + cl.delta_y),int(central_x - cl.xc): int(central_x - cl.xc + cl.delta_x)] += id_matrix # Adds the cloud id to the full matrix
    
    print(f"Done with cloud {key[6:]} out of {len(clouds)-1}.")

table = table[1:]
df = pd.DataFrame(table, columns = ["ID", "l [deg]", "b [deg]", "Area [pc^2]", "Distance [kpc]", "Flag", "Mass [M_o]", "# of distance points"])

max_lengths = df.map(str).map(len).max()
df = df.apply(lambda col: col.str.pad(max_lengths[col.name], side='right'))

df.to_csv("Product/ID_catalogue.tsv", index = False, sep = "\t")

hdr.set("NAXIS", 3) # Update headers
hdr.set("NAXIS3", 3, after = "NAXIS2")

data = np.array([full_distances, confidence, identities]).astype(np.float32) # Combine 2 2D matrices to a 3D one

hdu_new = fits.PrimaryHDU(data, hdr) # Save the combines distance and confidence matrix to a FITS file
hdul_new = fits.HDUList([hdu_new])
hdu_new.writeto(fr"Product/distance_mask.fits", overwrite = True)