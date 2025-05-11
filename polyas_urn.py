import numpy as np
import matplotlib.pyplot as plt

rng = np.random.default_rng()
# Draw 10 marbles of 3 color
urn = rng.integers(0,3,10)
_, counts =np.unique(urn, return_counts=True)
stats = [counts]
for i in range(1000):
    size = len(urn)
    index = rng.integers(0,size)
    marble = urn[index]
    added_marbels = marble*np.ones(3)
    urn = np.concatenate((urn, added_marbels))
    _, counts =np.unique(urn, return_counts=True)
    stats.append(counts)

stats = np.array(stats)
greens = stats[:,0]/np.sum(stats,axis=1)
reds = stats[:,1]/np.sum(stats,axis=1)
blues = stats[:,2]/np.sum(stats,axis=1)
y = np.vstack([greens,reds,blues])

fig, ax = plt.subplots()
x = np.arange(len(reds))
ax.stackplot(x,y, colors=["green", "red", "blue"])
ax.set(xlim=(0,len(reds)-1), ylim=(0,1))
plt.show()
