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
        epsilon = np.random.normal(0,1,n)
        dynamicspB = np.ones(n)*((r-1/2*sigma**2)*dt + epsilon *sigma*np.sqrt(dt))
        dynamicsnB = np.ones(n)*((r-1/2*sigma**2)*dt - epsilon *sigma*np.sqrt(dt))
        stocks_pB = S0*np.exp(np.cumsum(dynamicspB))
        stocks_nB = S0*np.exp(np.cumsum(dynamicsnB))
        dynamicspS = np.ones(n)*((r+1/2*sigma**2)*dt + epsilon*sigma*np.sqrt(dt))
        dynamicsnS = np.ones(n)*((r+1/2*sigma**2)*dt - epsilon*sigma*np.sqrt(dt))
        stocks_pS = S0*np.exp(np.cumsum(dynamicspS))
        stocks_nS = S0*np.exp(np.cumsum(dynamicsnS))
        pay_p = option_payoff(stocks_pB, barrier)
        pay_n = option_payoff(stocks_nB, barrier)
        payoffs_B.append(np.exp(-r*T)*pay_p)
        payoffs_B.append(np.exp(-r*T)*pay_n)
        payoffs_S.append(option_payoff(stocks_nS, barrier)/stocks_nS[-1])
        payoffs_S.append(option_payoff(stocks_pS, barrier)/stocks_pS[-1])
    ave_pay_B = np.average(payoffs_B)
    ave_pay_S = np.average(payoffs_S)
    return ave_pay_B, ave_pay_S

s =asian_call_avgK()
print(s, s[0]-s[1])