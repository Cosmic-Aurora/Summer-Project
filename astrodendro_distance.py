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
    for s in d.trunk:
        distance_matrix, counter = main_loop(d, s, mask, distance_matrix, points, distances, counter)
    print(f"Counter: {counter}")
    return distance_matrix

def main_loop(d, s, mask, distance_matrix, points, distances, counter):
    import numpy as np
    import sklearn as sk
    if s.is_branch:
        for c in s.children:
            distance_matrix, counter = main_loop(d, c, mask, distance_matrix, points, distances, counter)
   
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
        return distance_matrix, counter
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
            return distance_matrix, counter
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
            return distance_matrix, counter
        else:
            counter += 1
            structure_mean = np.mean(distances[points_within[np.where(np.array([d.structure_at(points[i,::-1].astype(int)) for i in points_within]) == s)[0]]])
            distance_matrix = np.where(distance_matrix == 0, structure_mean*structure_mask, distance_matrix)
            unique_distances = np.unique(distance_matrix)[1:]
            groups = sk.cluster.DBSCAN(eps = 0.5, min_samples = 1).fit(unique_distances.reshape(-1,1))
            distance_means = unique_distances
            for i in range(np.max(groups.labels_)):    
                distance_means = np.where(groups.labels_ == i, np.mean(unique_distances[np.where(groups.labels_==i)]),distance_means)
            for i in range(len(unique_distances)):
                distance_matrix = np.where(distance_matrix == unique_distances[i], distance_means[i], distance_matrix)
            return distance_matrix, counter

def allstructures(s, matrix, i):
    ind_x = s.indices()[0]
    ind_y = s.indices()[1]
    matrix[ind_x,ind_y] = i + 10
    i += 1
    if s.is_branch:
        for l in s.children:
            matrix, i = allstructures(l, matrix, i)
    return matrix, i






import numpy as np
import cloud
import matplotlib.pyplot as plt
import sklearn as sk
import astrodendro as ad
import scipy as sc

slice = 8
cloud_i = 91
distance_limit = 10

clouds = cloud.loadclouds(f"Data/clouds_{distance_limit}kpc.pkl")
cl = clouds[f"slice_{slice}_cloud_{cloud_i}.fits"]
points = cl.distances[:,:2]
distances = cl.distances[:,4]
data = cl.data
mask = data != 0


g_data = sc.ndimage.gaussian_filter(data,5)
d = ad.Dendrogram.compute(g_data, min_delta = 0.1, min_value = 0.1, is_independent = ad.pruning.contains_seeds(points[:,::-1].T.astype(int)))

distance_matrix = dendro_to_distance(d, mask, points, distances)

v = d.viewer()
v.show()

test_mask = np.zeros_like(data)
i = 0
for l in d:
    test_mask, i =  allstructures(l, test_mask, i)

fig, ax =  plt.subplots(1,2,figsize = (8,8))
ax[0].imshow(data, cmap = "YlGn_r")
ax[0].imshow(distance_matrix, alpha = 0.3)
scatter = ax[0].scatter(points[:,0], points[:,1], c = cl.distances[:,4], s = 5)
fig.colorbar(scatter, ax = ax[0])

ax[1].imshow(data, cmap = "YlGn_r")
ax[1].imshow(test_mask, alpha = 0.3)
ax[1].scatter(points[:,0], points[:,1], c = cl.distances[:,4])

plt.show()