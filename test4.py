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

def points_to_labels(points, data):
    labels = np.zeros_like(data)
    for i,p in enumerate(points.astype(int)):
        labels[p[1],p[0]] = i + 1
    return labels
import numpy as np

test = np.zeros(5,5)
mask = np.ones(5,5)
test[2,2] = 2
test[1,1:3] = 1
test[3, 1:3] = 1
test[2,1] = 1
test[2,3] = 1
print(test)

import numpy as np
import cloud
import matplotlib.pyplot as plt
import sklearn as sk

slice = 0
cloud_i = 35
distance_limit = 10

clouds = cloud.loadclouds(f"Data/clouds_{distance_limit}kpc.pkl")
cl = clouds[f"slice_{slice}_cloud_{cloud_i}.fits"]
points = cl.distances[:,:2]
data = cl.data
mask = data != 0

lables = points_to_labels(points, data)

vor = area_spread(lables, mask)

masked_vor = np.zeros_like(vor)
for i in range(len(points)):
    masked_vor = np.where(vor.astype(int) == i+1, cl.distances[i,4],masked_vor)

fig, ax = plt.subplots(1,2)
imgage = ax[0].imshow(masked_vor, cmap = "YlGn_r")
scatter = ax[0].scatter(points[:,0], points[:,1], c = cl.distances[:,4], cmap = "cool")
fig.colorbar(scatter, ax = ax[1])
fig.colorbar(imgage, ax = ax[0], location = "left")
ax[1].imshow(data, cmap = "YlGn_r")
ax[1].scatter(points[:,0], points[:,1], c = cl.distances[:,4], cmap = "cool", s = 2)
for i in (0,1):
    ax[i].set_xticks([])
    ax[i].set_yticks([])
fig.tight_layout()
plt.show()

groups =  sk.cluster.DBSCAN(eps = 0.5, min_samples = 1).fit(cl.distances[:,4].reshape(-1,1))
data_means = cl.distances[:,4]
for i in range(np.max(groups.labels_)):    
    data_means = np.where(groups.labels_ == i, np.mean(cl.distances[np.where(groups.labels_==i),4]),data_means) 
print(data_means)
masked_vor = np.zeros_like(vor)
for i in range(len(points)):
    masked_vor = np.where(vor.astype(int) == i+1, data_means[i],masked_vor)

fig, ax = plt.subplots(1,2)
imgage = ax[0].imshow(masked_vor, cmap = "YlGn_r")
scatter = ax[0].scatter(points[:,0], points[:,1], c = cl.distances[:,4], cmap = "cool")
fig.colorbar(scatter, ax = ax[1])
fig.colorbar(imgage, ax = ax[0], location = "left")
ax[1].imshow(data, cmap = "YlGn_r")
ax[1].scatter(points[:,0], points[:,1], c = cl.distances[:,4], cmap = "cool", s = 2)
for i in (0,1):
    ax[i].set_xticks([])
    ax[i].set_yticks([])
fig.tight_layout()
plt.show()
