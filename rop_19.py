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
        stocks_pB = [S0]
        stocks_nB = [S0]
        stocks_nS = [S0]
        stocks_pS = [S0]
        for j in range(n):
            epsilon = np.random.normal(0,1)
            stockpB = stocks_pB[j]*np.exp((r-1/2*sigma**2)*dt + sigma * np.sqrt(dt)*epsilon)
            stocknB = stocks_nB[j]*np.exp((r-1/2*sigma**2)*dt - sigma * np.sqrt(dt)*epsilon)
            stocks_pB.append(stockpB)
            stocks_nB.append(stocknB)
            stockpS = stocks_pS[j]*np.exp((r+sigma**2-.5*sigma*2)*dt + sigma*np.sqrt(dt)*epsilon)
            stocknS = stocks_nS[j]*np.exp((r+sigma**2-.5*sigma*2)*dt - sigma*np.sqrt(dt)*epsilon)
            stocks_nS.append(stocknS)
            stocks_pS.append(stockpS)
        pay_p = option_payoff(stocks_pB, barrier)
        pay_n = option_payoff(stocks_nB, barrier)
        payoffs_B.append(np.exp(-r*T)*pay_p)
        payoffs_B.append(np.exp(-r*T)*pay_n)
        payoffs_S.append(option_payoff(stocks_nS, barrier)/stocknS)
        payoffs_S.append(option_payoff(stocks_pS, barrier)/stockpS)
    ave_pay_B = np.average(payoffs_B)
    ave_pay_S = np.average(payoffs_S)
    return ave_pay_B, ave_pay_S

s =asian_call_avgK()
print(s, s[0]-s[1])