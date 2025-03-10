import numpy as np
from scipy import stats

def C_geom(S, K, r, sigma, T):
    sigma_g = sigma/np.sqrt(3)
    c = 1/2 * (r-1/2*sigma_g**2)
    d1 = (np.log(S/K)+ (c + 1/2*sigma_g**2)*T)/(sigma_g*np.sqrt(T))
    d2 = d1 - (sigma_g*np.sqrt(T))
    return S*np.exp((c-r)*T)*stats.norm.cdf(d1) - K*np.exp(-r*T)*stats.norm.cdf(d2)

def b(x, y):
    cov = np.cov(x,y)
    var = np.var(x)
    return cov/var

def asian_call(m, n, S0, K, r, sigma, dt):
    ave_stocks = []
    ave_geom_stocks = []
    dis_payoffs = []
    dis_geom_payoffs = []
    T = n*dt
    c_geom = C_geom(S0,K,r,sigma,T)
    for i in range(m):
        stocks = [S0]
        for j in range(n):
            epsilon = np.random.normal(0,1)
            stock = stocks[j]*np.exp((r-1/2*sigma**2)*dt + sigma * np.sqrt(dt)*epsilon)
            stocks.append(stock)
        stocks_ave = np.average(stocks)
        ave_stocks.append(stocks_ave)
        stocks_geom_ave = stats.gmean(stocks)
        ave_geom_stocks.append(stocks_geom_ave)
        dis_payoff = np.exp(-r*n*dt)*np.max([stocks_ave-K, 0])
        dis_payoffs.append(dis_payoff)
        dis_geom_payoff = np.exp(-r*n*dt)*np.max([stocks_geom_ave-K,0])
        dis_geom_payoffs.append(dis_geom_payoff)
    c_avg = np.average(dis_payoffs)-b(dis_geom_payoffs, c_geom)*(np.average(dis_geom_payoffs)-c_geom)
    return ave_stocks, c_avg
