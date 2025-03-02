import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

def stock_sim(stock, time, mu, sigma):
    epsilon = np.random.normal(0,1)
    return stock*np.exp((mu-1/2*sigma**2)*time + sigma*epsilon*np.sqrt(time)), epsilon

def bs_formula(stock, strike, time, riskless_rate, volatility):
    d1 = (np.log(stock/strike) + (riskless_rate + 1/2*volatility**2)*(time))/(volatility*np.sqrt(time))
    d2 = d1-volatility*np.sqrt(time)
    cdf = sp.stats.norm.cdf
    return stock*cdf(d1) - strike*np.exp(-riskless_rate*time)*cdf(d2)

def a():
    dt = 1/252
    T = 1
    N = int(1/dt)
    S0 = 50
    mu = .15
    sigma = .5
    K = 50
    q = 0
    r = .13
    stocks = [S0]
    epsilons = []
    calls = []
    for i in range(N):
        time = i*dt
        stock, epsilon = stock_sim(stocks[i], dt, mu, sigma)
        stocks.append(stock)
        epsilons.append(epsilon)
        call_price = bs_formula(stock, strike=K, time=1-time, riskless_rate=r, volatility=sigma)
        calls.append(call_price)
    
    fig, ax = plt.subplots()
    ax.plot(np.arange(N), stocks[1:], np.arange(N), calls, np.arange(N), K*np.ones(N))
    plt.show()

a()
