import numpy as np

K = 105
s = 100
r = .03
U = 1.07
D = 1/U
months = 6
dt = 1/12
p = (np.exp(r*dt)-D)/(U-D)

stock_tree = [[s]]
for i in range(months):
    next_price = []
    for j in stock_tree[i]:
        next_price.append(U*j)
    next_price.append(stock_tree[i][-1]*D)
    stock_tree.append(next_price)

option_tree = [[] for i in range(len(stock_tree))]
option_pay = lambda x: float(np.max([0,x]))
for i in stock_tree[-1]:
    option_value = option_pay(K-i)
    option_tree[0].append(option_value)

option_price = lambda U, D: float((p*U + (1-p)*D)*np.exp(-r*dt))
for index, iteration in enumerate(option_tree):
    for i, v in enumerate(iteration):
        if i + 1 == len(iteration):
            break
        lower = iteration[i + 1]
        upper = iteration[i]
        val = option_price(upper, lower)
        option_tree[index+1].append(val)
    if index + 1 == len(option_tree):
        break
    


print(option_tree[-1][-1])



# homework European options:
