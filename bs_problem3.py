import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

def stock_sim(stock, time, mu, sigma):
    epsilon = np.random.norm(0,1)
    return stock*np.exp((mu-1/2*sigma**2)*time + sigma*epsilon*np.sqrt(time)), epsilon

def bs_formula(stock, strike, time, riskless_rate, volatility):
    d1 = (np.log(stock/strike) + (riskless_rate + 1/2*volatility**2)*(time))/(volatility*np.sqrt(time))
    d2 = d1-volatility*np.sqrt(time)
    cdf = sp.stats.norm.cdf
    return stock*cdf(d1) - strike*np.exp(-riskless_rate*time)*cdf(d2)

