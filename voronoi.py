def voronoi(points, data):
    import numpy as np
    vor = np.zeros_like(data)
    mask = data != 0
    vor = np.array([[np.argmin(np.linalg.norm(points[:, ::-1] - np.array([i, j]), axis=1)) + 1 for j in range(vor.shape[1])]for i in range(vor.shape[0])])
    masked_voronai = vor*mask
    return masked_voronai

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
vor = voronoi(points, data)

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
