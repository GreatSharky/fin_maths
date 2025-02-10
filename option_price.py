import numpy as np

K = .9
S0 = 1
U = 1.2
D = 1/U

r = .05
dt = 1/12
p = (np.exp(r*dt)-D)/(U-D)
option_pay = lambda x: np.max([0,x])
C = (p*option_pay(S0*U-K)+(1-p)*option_pay(S0*D-K))* np.exp(-r*dt)
print(C,"\t", p)