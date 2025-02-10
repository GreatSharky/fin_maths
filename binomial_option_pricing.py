import numpy as np
import math
S0 = 100
U = 1.07
D = 1/U
r = .03
K = 105
T = .5
dt = 1/12
m = int(T/dt)

def binompmf(n,k,p):
    nck = math.factorial(n)/(math.factorial(k)*math.factorial(n-k))
    return nck*p**k*(1-p)**(n-k)

diskont = np.exp(-r*T)
p = (np.exp(r*dt)-D)/(U-D)

payoff = lambda x: np.max([x,0])
sum = 0

for k in range(m+1):
    pay = payoff(K-S0 * U**k * D**(m-k))
    sum += binompmf(m,k,p)*payoff(K-S0 * U**k * D**(m-k))

price = diskont*sum

print(price)
