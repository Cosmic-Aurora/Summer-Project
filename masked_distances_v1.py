def area_spread(labeled_image, mask, step_size = 5): # can be used a svoronoi with step_size = big (image size)
    import numpy as np
    import skimage as sk
    if not np.all(labeled_image == 0):
        while True:
            unmasked_image = sk.segmentation.expand_labels(labeled_image, distance = step_size)
            labeled_image = unmasked_image*mask
            if np.all((labeled_image != 0) == mask):
                return labeled_image
    else:
        raise TypeError("Labeled image must have labels.")

def voronoi(labeled_image, mask):
    return area_spread(labeled_image, mask, np.inf)

def points_to_labels(points, data, min_value = 0):
    labels = np.zeros_like(data).astype(int)
    for i,p in enumerate(points.astype(int)):
        labels[p[1],p[0]] = i + 1 + min_value
    return labels

def dendro_to_distance(d, mask, points, distances, counter = 0):
    import numpy as np
    distance_matrix = np.zeros_like(mask).astype(float)
    confidence_matrix = np.zeros_like(mask).astype(int)
    for s in d.trunk:
        distance_matrix, confidence_matrix, counter = main_loop(d, s, mask, distance_matrix, confidence_matrix, points, distances, counter)
    return distance_matrix, confidence_matrix, counter

def main_loop(d, s, mask, distance_matrix, confidence_matrix, points, distances, counter):
    import numpy as np
    import sklearn as sk
    import skimage as si
    A = 4 # Non-conflicting clouds
    B = 3 # Non-conflicting leaves/branches with conflicts below
    C = 2 # Conflicting leaves
    D = 1 # Conflicting branches (the rest)
    if s.is_branch:
        for c in s.children:
            distance_matrix, confidence_matrix, counter = main_loop(d, c, mask, distance_matrix, confidence_matrix, points, distances, counter)
    ind_x = s.indices()[0]
    ind_y = s.indices()[1]
    structure_mask = np.zeros_like(mask)
    structure_mask[ind_x,ind_y] = 1
    points_within = []
    for i, p in enumerate(points):
        if structure_mask[int(p[1]),int(p[0])]:
            points_within.append(i)
    points_within = np.array(points_within)
    
    if np.abs(np.max(distances[points_within])-np.min(distances[points_within])) < 0.5:
        distance_matrix[ind_x,ind_y] = float(np.mean(distances[points_within]))
        confidence_matrix[ind_x,ind_y] = A
        return distance_matrix, confidence_matrix, counter
    else: 
        if s.is_leaf:
            groups = sk.cluster.DBSCAN(eps = 0.5, min_samples = 1).fit(distances[points_within].reshape(-1,1))
            distance_means = distances[points_within]
            for i in range(np.max(groups.labels_)):    
                distance_means = np.where(groups.labels_ == i, np.mean(distances[np.where(groups.labels_==i)]),distance_means) 
            labels = points_to_labels(points[points_within], structure_mask)
            labeled_image = voronoi(labels, structure_mask)
            for i in range(len(distance_means)):
                distance_matrix = np.where(labeled_image == i+1, distance_means[i], distance_matrix)
            confidence_matrix[ind_x, ind_y] = C
            return distance_matrix, confidence_matrix, counter
        elif s not in [d.structure_at(points[i,::-1].astype(int)) for i in points_within]: # no new points on this branch that arent from an above branch/leaf
            unique_distances = np.unique(distance_matrix)[1:]
            groups = sk.cluster.DBSCAN(eps = 0.5, min_samples = 1).fit(unique_distances.reshape(-1,1))
            distance_means = unique_distances
            for i in range(np.max(groups.labels_)):    
                distance_means = np.where(groups.labels_ == i, np.mean(unique_distances[np.where(groups.labels_==i)]),distance_means) 
            labeled_image = np.zeros_like(mask).astype(int)
            for i in range(len(unique_distances)):
                labeled_image = np.where(distance_matrix == unique_distances[i], i + 1, labeled_image)
            labeled_image = area_spread(labeled_image, structure_mask)
            for i in range(len(distance_means)):
                distance_matrix = np.where(labeled_image == i+1, distance_means[i], distance_matrix)
            confidence_matrix = np.where(confidence_matrix == A, B, confidence_matrix)
            confidence_matrix = np.where(confidence_matrix == 0, D*structure_mask, confidence_matrix)
            return distance_matrix, confidence_matrix, counter
        else:
            counter += 1
            
            new_indices = np.where(np.array([d.structure_at(points[i,::-1].astype(int)) for i in points_within]) == s)[0]
            new_points = points[points_within[new_indices]]
            new_distances = distances[points_within[new_indices]]
            unique_distances = np.unique(distance_matrix)[1:]
            labels = points_to_labels(new_points, structure_mask, len(unique_distances)-len(new_distances))
            labeled_image_new = si.segmentation.expand_labels(labels, 25)*structure_mask
            all_distances = np.append(unique_distances, new_distances)

            labeled_image = np.zeros_like(mask).astype(int)
            for i in range(len(unique_distances)):
                labeled_image = np.where(distance_matrix == unique_distances[i], i + 1, labeled_image)
            labeled_image = np.where(labeled_image == 0, labeled_image_new, labeled_image)
            labeled_image = area_spread(labeled_image, structure_mask)
            for i in range(len(all_distances)):
                distance_matrix = np.where(labeled_image == i+1, all_distances[i], distance_matrix)
                
            groups = sk.cluster.DBSCAN(eps = 0.5, min_samples = 1).fit(all_distances.reshape(-1,1))
            distance_means = all_distances
            for i in range(np.max(groups.labels_)):    
                distance_means = np.where(groups.labels_ == i, np.mean(all_distances[np.where(groups.labels_==i)]),distance_means)
            for i in range(len(unique_distances)):
                distance_matrix = np.where(distance_matrix == all_distances[i], distance_means[i], distance_matrix)
            confidence_matrix = np.where(confidence_matrix == A, B, confidence_matrix)
            confidence_matrix = np.where(confidence_matrix == 0, D*structure_mask, confidence_matrix)
            return distance_matrix, confidence_matrix, counter

import cloud
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sklearn as sk
import scipy as sc
import astrodendro as ad

hdul = fits.open("Data/PROMISE/Full_Data.fits")
hdr = hdul[0].header
full_data = hdul[0].data
full_distances = np.zeros_like(full_data)
confidence = np.zeros_like(full_data)

clouds = cloud.loadclouds("Data/clouds_10kpc.pkl")
counter = 0
keys = []

central_x = hdr["CRPIX1"]
central_y = hdr["CRPIX2"]

for key in clouds:
    cl = clouds[key]
    if cl.n == 0:
        print(f"Done with cloud {key} and counter is at {counter}.")
        continue
    
    points = cl.distances[:,:2]
    distances = cl.distances[:,4]
    data = cl.data
    mask = data != 0
    g_data = sc.ndimage.gaussian_filter(data,5)*mask
    d = ad.Dendrogram.compute(g_data, min_delta = 0.1, min_value = 0.1, is_independent = ad.pruning.contains_seeds(points[:,::-1].T.astype(int)))

    distance_matrix, confidence_matrix, counter = dendro_to_distance(d, mask, points, distances, counter)
    full_distances[int(central_y - clouds[key].yc): int(central_y - clouds[key].yc + clouds[key].delta_y),int(central_x - clouds[key].xc): int(central_x - clouds[key].xc + clouds[key].delta_x)] += distance_matrix
    confidence[int(central_y - clouds[key].yc): int(central_y - clouds[key].yc + clouds[key].delta_y),int(central_x - clouds[key].xc): int(central_x - clouds[key].xc + clouds[key].delta_x)] += confidence_matrix
    
    print(f"Done with cloud {key} and counter is at {counter}.")

hdr.set("NAXIS", 3)
hdr.set("NAXIS3", 2, after = "NAXIS2")

data = np.array([full_distances, confidence]).astype(np.float32)


hdu_new = fits.PrimaryHDU(data, hdr)
hdul_new = fits.HDUList([hdu_new])
hdu_new.writeto(fr"Data/dist_10_v1.fits", overwrite = True)