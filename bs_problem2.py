import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

expected_growth = 0.15
volatility = 0.5
dividend_yield = .08
r = .13
dt = 1/252
maturity_time = 5
N = int(maturity_time/dt)

def theta(stock, strike, time):
    cdf = sp.stats.norm.cdf
    pdf = sp.stats.norm.pdf
    d1 = (np.log(stock/strike) +  (r - dividend_yield + 1/2*volatility**2)*(time))/(volatility*np.sqrt(time))
    d2 = (np.log(stock/strike) + (r - dividend_yield-1/2*volatility**2)*(time))
    one = -(stock*pdf(d1)*volatility*np.exp(-dividend_yield*(time)))/(2*np.sqrt(time))
    two = dividend_yield*stock*cdf(d1)*np.exp(-dividend_yield*(time))
    three = -r*strike*np.exp(-r*(time))*cdf(d2)
    return one + two + three

def delta(stock, strike, time):
    d1 = (np.log(stock/strike) +  (r - dividend_yield + 1/2*volatility**2)*(time))/(volatility*np.sqrt(time))
    return np.exp(-dividend_yield*(time))*sp.stats.norm.cdf(d1)

def gamma(stock, strike, time):
    pdf = sp.stats.norm.pdf
    d1 = (np.log(stock/strike) +  (r - dividend_yield + 1/2*volatility**2)*(time))/(volatility*np.sqrt(time))
    return (pdf(d1)*np.exp(-dividend_yield*(time)))/(stock*volatility*np.sqrt(time))

def bs_formula(stock, strike, time):
    d1 = (np.log(stock/strike) + (r + 1/2*volatility**2)*(time))/(volatility*np.sqrt(time))
    d2 = d1-volatility*np.sqrt(time)
    cdf = sp.stats.norm.cdf
    return stock*cdf(d1) - strike*np.exp(-r*time)*cdf(d2)

def call_dynamic(option, stock, strike, time):
    mu_call = 1/option*(theta(stock, strike, time) + expected_growth*stock*delta(stock, strike, time) + 1/2*volatility**2*stock**2*gamma(stock, strike, time))
    sigma_call = 1/option*volatility*stock*delta(stock, strike, time)
    return mu_call, sigma_call

def a():
    stock = 50
    strike = 50
    time = np.arange(1/252,5,1/252)
    call_prices = bs_formula(stock, strike, time)
    mu, sigma = call_dynamic(call_prices, stock, strike, time)
    fig, ax = plt.subplots()
    ax.plot(time, mu, time, sigma)
    plt.show()

def b():
    stock = np.arange(0.01,100,0.1)
    time = 1
    strike = 50
    call_prices = bs_formula(stock, strike, time)
    mu, sigma = call_dynamic(call_prices, stock, strike, time)
    fig, ax = plt.subplots()
    ax.plot(stock[150:], mu[150:],stock[:], sigma[:])
    plt.show()

b()