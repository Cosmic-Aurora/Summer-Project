import cloud
import matplotlib.pyplot as plt
import numpy as np

clouds5 = cloud.loadclouds("Data/clouds_5kpc.pkl")
clouds10 = cloud.loadclouds("Data/clouds_10kpc.pkl")
clouds15 = cloud.loadclouds("Data/clouds_15kpc.pkl")
clouds20 = cloud.loadclouds("Data/clouds_20kpc.pkl")

clouds = [clouds5, clouds10, clouds15, clouds20]

n = [[],[],[],[]]
means = [[],[],[],[]]
vars = [[],[],[],[]]
diff = [[],[],[],[]]

# n5 = []
# n10 = []
# n15 = []
# n20 = []

# means5 = []
# means10 = []
# means15 = []
# means20 = []

# vars5 = []
# vars10 = []
# vars15 = []
# vars20 = []

# diff5 = []
# diff10 = []
# diff15 = []
# diff20 = []

zeros = np.zeros(4)
ones = np.zeros(4)
complication = np.zeros(4)
additional = np.zeros(4)
keys = []
increases = []

for key in clouds5:
    keys.append(key)
    for i in range(4):
        n[i].append(clouds[i][key].n)
        means[i].append(clouds[i][key].mean)
        vars[i].append(clouds[i][key].var)
        diff[i].append(clouds[i][key].diff)
        if clouds[i][key].n == 0:
            zeros[i] += 1
        elif clouds[i][key].n == 1:
            ones[i] +=1
        if clouds[i][key].n != 1 and np.any([clouds[k][key].n == 1 for k in range(0,i)]):
            complication[i] += 1
        if clouds[i][key].n != 0 and np.any([clouds[k][key].n == 0 for k in range(0,i)]):
            additional[i] += 1
        if np.any([clouds[k][key].n < clouds[i][key].n for k in range(0,i)]):
            increases.append(key)
    
    
    # n5.append(clouds5[key].n)
    # means5.append(clouds5[key].mean)
    # vars5.append(clouds5[key].var)
    # diff5.append(clouds5[key].diff)

    # n10.append(clouds10[key].n)
    # means10.append(clouds10[key].mean)
    # vars10.append(clouds10[key].var)
    # diff10.append(clouds10[key].diff)

    # n15.append(clouds15[key].n)
    # means15.append(clouds15[key].mean)
    # vars15.append(clouds15[key].var)
    # diff15.append(clouds15[key].diff)

    # n20.append(clouds20[key].n)
    # means20.append(clouds20[key].mean)
    # vars20.append(clouds20[key].var)
    # diff20.append(clouds20[key].diff)

colors = ["blue", "green", "purple", "red"]
fig, ax = plt.subplots(4,5, figsize = (20,8))
#fig, ax = plt.subplots(4,5, figsize = (30,13))
labels = ["5 kpc", "10 kpc", "15 kpc", "20 kpc"]
for i in range(4):
    ax[i,0].hist(vars[i], bins = np.arange(max(vars[i])+1), color = colors[i])
    ax[i,0].set_title(f"Variance {5*(i+1)} kpc")
    ax[i,1].hist(means[i], bins = np.arange(max(means[i])+1),color = colors[i])
    ax[i,1].set_title(f"Mean {5*(i+1)} kpc")
    ax[i,2].hist(n[i], bins = np.arange(30),color = colors[i])
    ax[i,2].set_title(f"# of datapoints {5*(i+1)} kpc")
    ax[i,3].hist(diff[i], bins = np.arange(max(diff[i])+1),color = colors[i])
    ax[i,3].set_title(f"Absolute difference {5*(i+1)} kpc")
ax[0,4].bar(labels, zeros, color = colors)
ax[1,4].bar(labels,ones, color = colors)
ax[2,4].bar(labels, additional, color = colors)
ax[3,4].bar(labels, complication, color = colors)
fig.tight_layout()
#plt.show()

n_x = 4
n_y = 3
n_tot = n_x*n_y
random = True

selected_keys = increases[:n_y]
if random:
    selected_keys = np.random.choice(increases,size = n_y, replace = False)
selected_keys = np.repeat(selected_keys, n_x)

fig, ax = plt.subplots(n_y,n_x, figsize = (20,8))
for i, k in enumerate(selected_keys):
    x = int(i/n_x)
    y = int(i%n_x)
    ax[x,y].imshow(clouds[x][k].data, cmap = "YlGn_r")
    ax[x,y].set_title(f"{k} {5*(y+1)} kpc, N = {clouds[y][k].n}")
    ax[x,y].set_xticks([])
    ax[x,y].set_yticks([])
    ax[x,y].set_xlabel(f"mean: {clouds[y][k].mean}.")
    ax[x,y].set_ylabel(f"var: {clouds[y][k].var}")
    scatter = ax[x,y].scatter(clouds[y][k].distances[:,0], clouds[y][k].distances[:,1], c = clouds[y][k].distances[:,4], cmap = "cool", s = 5)
    fig.colorbar(scatter,ax = ax[x,y])
fig.tight_layout()
plt.show()
