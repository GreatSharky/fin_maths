import numpy as np
from scipy import stats

m = 100
n = 252
S0 = 50
K = 50
r = .10
sigma = .15
dt = 1/n
S1_maxes = []
S2_maxes = []
V1 = []
V2 = []
X = 60
W = []
for i in range(m):
    stocks1 = [S0]
    stocks2 = [S0]
    for j in range(n):
        e = np.random.normal(0,1)
        stock = stocks1[j]*np.exp((r-1/2*sigma**2)*dt + sigma*np.sqrt(dt)*e)
        stocks1.append(stock)
        stock = stocks2[j]*np.exp((r-1/2*sigma**2)*dt - sigma*np.sqrt(dt)*e)
        stocks2.append(stock)
    S1_maxes.append(np.max(stocks1))
    S2_maxes.append(np.max(stocks2))
    if S1_maxes[i] > X:
        V1.append(np.exp(-r*n*dt)*np.max([stocks1[-1]-K,0]))
    else:
        V1.append(0)
    if S2_maxes[i] > X:
        V2.append(np.exp(-r*n*dt)*np.max([stocks2[-1]-K,0]))
    else:
        V2.append(0)
    W.append(1/2*(V1[-1]+V2[-1]))

print(np.var(V1)/np.var(W))
