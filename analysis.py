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

keys = []

for key in clouds5:
    keys.append(key)
    for i in range(4):
        n[i].append(clouds[i][key].n)
        means[i].append(clouds[i][key].mean)
        vars[i].append(clouds[i][key].var)
        diff[i].append(clouds[i][key].diff)
    
    
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


fig, ax = plt.subplots(4,5, figsize = (30,13))

for i in range(4):
    ax[i,0].hist(vars[i], bins = np.arange(max(vars[i])+1))
    ax[i,1].hist(means[i], bins = np.arange(max(means[i])+1))
    ax[i,2].hist(n[i], bins = np.arange(max(n[i])+1))
    ax[i,3].hist(diff[i], bins = np.arange(max(diff[i])+1))

fig.tight_layout()
plt.show()