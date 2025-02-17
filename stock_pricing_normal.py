import numpy as np
import matplotlib.pyplot as plt

N = int(10E5)

S0 = 50

days = int(.5*252)
T = .5
volatiltiy = .3
expected_growth = .13

dt = T

sims = []
exact = []
for sim in range(N):
    epsilon = np.random.normal(0,1) 
    stock_approx = S0*(1+expected_growth*dt + epsilon*np.sqrt(dt)*volatiltiy)
    stock_exact = S0*np.exp((expected_growth - 1/2*volatiltiy**2)*(dt) + volatiltiy*epsilon*np.sqrt(dt))
    sims.append(stock_approx)
    exact.append(stock_exact)

exact = np.array(exact)
count_sims, bin_sims = np.histogram(sims, 100)
count_exact, bin_exact = np.histogram(exact,100)
plt.stairs(count_sims, bin_sims, fill=True)
plt.stairs(count_exact, bin_exact)
plt.show()

print("Sim mean",np.mean(exact))
true_mean = S0*np.exp(dt*expected_growth)
print("True mean", true_mean)
print("Sim log mean", np.mean(np.log(exact)))
print("True log mean", np.log(S0)+(expected_growth-.5*volatiltiy**2)*dt)
