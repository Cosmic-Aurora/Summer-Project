import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import os
import cloud

hits = []
variances = []
means = []
diff = []

clouds = cloud.loadclouds("Data/clouds_10kpc.pkl")
keys = []
print(len(clouds))
for key in clouds:
    keys.append(key)
    c = clouds[key]
    hits.append(c.n)
    variances.append(c.var)
    means.append(c.mean)
    diff.append(c.diff)

keys = np.array(keys)

n_x = 4
n_y = 2
n_tot = n_x*n_y

min_num = 1
max_num = np.inf

min_var = 10
max_var = np.inf

min_mean = -np.inf
max_mean = np.inf

min_diff = 0
max_diff = np.inf

random = True

within = np.where(np.all((np.array(hits) >= min_num, np.array(hits) <= max_num, np.array(variances) >= min_var, np.array(variances) <= max_var, np.array(means) >= min_mean, np.array(means) <= max_mean, np.array(diff) >= min_diff, np.array(diff) <= max_diff), axis = 0))[0]
print(len(within))

if len(within) < n_tot:
    n_tot = len(within)

selected_keys = keys[within[:n_tot]]
if random:
    selected_keys = keys[np.random.choice(within,size = n_tot, replace = False)]


fig, ax = plt.subplots(n_y,n_x, figsize = (20,8))
for i, k in enumerate(selected_keys):
    x = int(i/n_x)
    y = i%n_x
    ax[x,y].imshow(clouds[k].data, cmap = "YlGn_r")
    ax[x,y].set_title(f"{k}, N = {clouds[k].n}")
    ax[x,y].set_xticks([])
    ax[x,y].set_yticks([])
    ax[x,y].set_xlabel(f"mean: {clouds[k].mean}.")
    ax[x,y].set_ylabel(f"var: {clouds[k].var}")
    scatter = ax[x,y].scatter(clouds[k].distances[:,0], clouds[k].distances[:,1], c = clouds[k].distances[:,4], cmap = "cool", s = 5)
    fig.colorbar(scatter,ax = ax[x,y])
    #print(k, clouds[k].diff)
fig.tight_layout()
plt.show()
    









#print(hits)
#clearzeros = np.where(np.array(hits) == 0)[0]
#print(len(clearzeros))
#variances = np.delete(variances, clearzeros)
#hits = np.delete(hits, clearzeros)
#means = np.delete(means, clearzeros)

# lookat = np.where(np.all((np.array(variances) > 10,np.array(hits) > 10), axis = 0))[0]
# lookat = np.where(np.all((np.array(hits) < 3,np.array(hits) > 0, np.array(variances) > 0.1), axis = 0))[0]
# print(len(lookat))
# print(keys[lookat])
# #ax, fig = plt.sub
# plt.figure(figsize = (10,10))
# plt.hist(hits, bins = np.arange(20)-0.5)
# plt.show()

# plt.figure(figsize = (10,10))
# plt.hist(variances)
# plt.show()
# plt.figure(figsize = (10,10))
# plt.hist(means)
# plt.show()
# plt.figure(figsize = (10,10))
# plt.hist(np.divide(variances,means))
# plt.show()