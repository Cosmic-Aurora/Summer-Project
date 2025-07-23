import cloud
import astrodendro as ad
import numpy as np
import matplotlib.pyplot as plt
import scipy as sc

def iterations(l, leaflist):
    if l.is_leaf:
        leaflist.append(l)
        return leaflist
    else:
        for l in l.children:
            leaflist = iterations(l, leaflist)
        return leaflist


clouds = cloud.loadclouds("Data/selected_clouds.pkl")
keys = []

smoothing = 6

total_datapoints = 0

conflicted_leafs = 0
conflicted_leafs_ids = []

empty_leaves = 0

nones = 0
nones_ids = []

trunks = 0
trunks_w_conflicts = 0

leaves = 0
leaves_w_conflicts = 0

for key in clouds:
    keys.append(key)
    print(key)
    cl = clouds[key]
    points = cl.distances[:,:2]
    total_datapoints += len(points)
    data = cl.data
    g_data = sc.ndimage.gaussian_filter(data,smoothing)
    d = ad.Dendrogram.compute(g_data, min_delta = 0.1, min_value = 0.1, is_independent = ad.pruning.contains_seeds(points[:,::-1].T.astype(int)))
    #d = ad.Dendrogram.compute(g_data, min_delta = 0.1, min_value = 0.1)
    trunks += len(d.trunk)
    leaflist = []
    leaves_w_points = []
    for l in d.trunk:
        leaflist = iterations(l, leaflist)
    leaves += len(leaflist)
    usedtrunks = []
    usedleaves = []
    for i0, p0 in enumerate(points):
        #print(d.structure_at(p0[::-1].astype(int)))
        if d.structure_at(p0[::-1].astype(int)) == None:
            nones += 1
            nones_ids.append(key)
        else:
            for i1, p1 in enumerate(points):
                if d.structure_at(p0[::-1].astype(int)) == d.structure_at(p1[::-1].astype(int)) and i0 != i1 and np.abs(cl.distances[i0,4] - cl.distances[i1,4]) > 0.5 and d.structure_at(p0[::-1].astype(int)).is_leaf:
                    conflicted_leafs += 1
                    conflicted_leafs_ids.append(key)
                    if d.structure_at(p0[::-1].astype(int)) not in usedleaves:
                        usedleaves.append(d.structure_at(p0[::-1].astype(int)))
                        leaves_w_conflicts += 1
                if d.structure_at(p1[::-1].astype(int)) == None:
                    pass
                elif d.structure_at(p0[::-1].astype(int)).ancestor == d.structure_at(p1[::-1].astype(int)).ancestor and i0 != i1 and np.abs(cl.distances[i0,4] - cl.distances[i1,4]) > 0.5 and d.structure_at(p0[::-1].astype(int)).ancestor not in usedtrunks:
                    trunks_w_conflicts += 1
                    usedtrunks.append(d.structure_at(p0[::-1].astype(int)).ancestor)
            if d.structure_at(p0[::-1].astype(int)).is_leaf:
                leaves_w_points.append(d.structure_at(p0[::-1].astype(int)))

    empty_leaves += len(list(set(leaflist) - set(leaves_w_points)))

print(conflicted_leafs)
#print(conflicted_leafs_ids)
print(nones)
#print(nones_ids)
print(total_datapoints)
print(empty_leaves)
print(trunks)
print(trunks_w_conflicts)
print(leaves)
print(leaves_w_conflicts)

plt.imshow(data)
plt.scatter(points[:,0], points[:,1], c = np.arange(len(points)))
plt.show()