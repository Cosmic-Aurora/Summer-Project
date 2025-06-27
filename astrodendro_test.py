
# %%

import cloud
import astrodendro as ad
import numpy as np
import matplotlib.pyplot as plt
import scipy as sc

def iterations(s, matrix, i):
    if s.is_leaf:
        ind_x = s.indices()[0]
        ind_y = s.indices()[1]
        matrix[ind_x,ind_y] = i + 10
        i += 1
        return matrix, i
    else:
        for l in s.children:
            matrix, i = iterations(l, matrix, i)
        return matrix, i

clouds = cloud.loadclouds("Data/clouds_10kpc.pkl")
slice = 2
cloud_i =183
cl = clouds[f"slice_{slice}_cloud_{cloud_i}.fits"]

data = cl.data
g_data = sc.ndimage.gaussian_filter(data,8)

combined_data = np.array([data,g_data])
minima = np.amin(combined_data)
maxima = np.amax(combined_data)
points = cl.distances[:,:2]

fig, ax = plt.subplots(1,2, figsize = (12,8))
imgage = ax[0].imshow(g_data, cmap = "YlGn_r", vmin = minima, vmax = maxima)
scatter = ax[0].scatter(points[:,0], points[:,1], c = cl.distances[:,4], cmap = "cool")
fig.colorbar(scatter, ax = ax[1])
fig.colorbar(imgage, ax = ax[0], location = "left")
smooth_image = ax[1].imshow(data, cmap = "YlGn_r", vmin = minima, vmax = maxima)
ax[1].scatter(points[:,0], points[:,1], c = cl.distances[:,4], cmap = "cool", s = 2)
for i in (0,1):
    ax[i].set_xticks([])
    ax[i].set_yticks([])
fig.tight_layout()
plt.show()

#print(cl.distances[:,:2].T)
d = ad.Dendrogram.compute(g_data, min_delta = 0.1, min_value = 0.1, is_independent = ad.pruning.contains_seeds(points[:,::-1].T.astype(int)))
#d = ad.Dendrogram.compute(g_data, min_delta = 1, min_value = 0.1)
#print(np.array(d.trunk[0].children))
print(d.trunk[0].parent)
v = d.viewer()
v.show()

test_mask = np.zeros_like(data)
i = 0
for l in d:
    test_mask, i =  iterations(l, test_mask, i)

# for i, l in enumerate(d.trunk[0].descendants):
#     if l.is_leaf:
#         print(l)
#         ind_x = l.indices()[0]
#         ind_y = l.indices()[1]
#         test_mask[ind_x,ind_y] = i + 10
plt.imshow(data, cmap = "YlGn_r")
plt.imshow(test_mask, alpha = 0.3)
plt.scatter(points[:,0], points[:,1], c = cl.distances[:,4])
plt.colorbar()
plt.show()