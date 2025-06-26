def voronai(points, data):
    import numpy as np
    vor = np.zeros_like(data)
    mask = data != 0
    for i in range(len(vor)):
        for j in range(len(vor.T)):
            vor[i,j] = np.argmin([np.linalg.norm((i,j) - points[k,::-1]) for k in range(0,len(points))]) + 1
    masked_voronai = vor*mask
    return masked_voronai

import numpy as np
import cloud
import matplotlib.pyplot as plt

slice = 2
cloud_i = 64
distance_limit = 10

clouds = cloud.loadclouds(f"Data/clouds_{distance_limit}kpc.pkl")
cl = clouds[f"slice_{slice}_cloud_{cloud_i}.fits"]
points = cl.distances[:,:2]
data = cl.data
vor = voronai(points, data)

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
