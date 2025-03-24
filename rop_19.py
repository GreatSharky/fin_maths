import numpy as np

def option_payoff(stocks, barrier):
    if np.min(stocks) < barrier:
        strike = np.average(stocks)
        pay = np.max([stocks[-1]-strike,0])
    else:
        pay = 0
    return pay

def asian_call_avgK(m=100000, S0=1, r=.035, sigma=.3, dt=1/252, barrier=.85, T = 1):
    payoffs_B = []
    payoffs_S = []
    n = int(T/dt)
    S0 = float(S0)
    for i in range(m):
        stocks_p = [S0]
        stocks_n = [S0]
        for j in range(n):
            epsilon = np.random.normal(0,1)
            stockp = stocks_p[j]*np.exp((r-1/2*sigma**2)*dt + sigma * np.sqrt(dt)*epsilon)
            stockn = stocks_n[j]*np.exp((r-1/2*sigma**2)*dt - sigma * np.sqrt(dt)*epsilon)
            stocks_p.append(stockp)
            stocks_n.append(stockn)
        pay_p = option_payoff(stocks_p, barrier)
        pay_n = option_payoff(stocks_n, barrier)
        payoffs_B.append(np.exp(-r*T)*pay_p)
        payoffs_B.append(np.exp(-r*T)*pay_n)
        payoffs_S.append(pay_n/stockn)
        payoffs_S.append(pay_p/stockp)
    ave_pay_B = np.average(payoffs_B)
    ave_pay_S = np.average(payoffs_S)
    return ave_pay_B, ave_pay_S

print(asian_call_avgK())