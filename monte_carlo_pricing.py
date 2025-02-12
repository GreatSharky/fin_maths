import numpy as np
import matplotlib.pyplot as plt

N = int(10E3)

S0 = 50

days = int(.5*252)
T = .5
volatiltiy = .3
expected_growth = .13

dt = 1/126


sims = []
exact = []
for sim in range(N):
    stock_approx = np.ones(days)
    stock_approx[0] = S0
    stock_exact = np.ones(days)
    stock_exact[0] = S0
    for day in range(days-1):
        epsilon = np.random.normal(0,1) 
        stock_approx[day+1] = stock_approx[day]*(1+expected_growth*dt + epsilon*np.sqrt(dt)*volatiltiy)
        stock_exact[day+1] = stock_exact[day]*np.exp((expected_growth - 1/2*volatiltiy**2)*(dt) + volatiltiy*epsilon*np.sqrt(dt))
    sims.append(stock_approx)
    exact.append(stock_exact)

fig, ax = plt.subplots()
x_axis = np.array([x for x in range(days)])
print(exact[-1].shape)
ax.plot(x_axis, exact[-1], x_axis, sims[-1], "--")
plt.show()