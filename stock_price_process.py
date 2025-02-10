import numpy as np
import matplotlib.pyplot as plt

mu = .15
sigma = .3
days = 252
dt = 1/252
years = 2
N = int(years/dt)

stock = np.ones(N)
# stock = np.ones(time)
# print(type(stock[0]))
for index in range(N-1):
    rand = np.random.normal(0,1)
    stock[index+1] = stock[index]*(1 + mu*dt + rand*np.sqrt(dt)*sigma)


fig, ax = plt.subplots()
ax.plot(np.array([i for i in range(N)]),stock)
plt.show()