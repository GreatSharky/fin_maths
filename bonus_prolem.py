import numpy as np

def option_payoff(stocks, barrier):
    """Option payoff with barrier to reduce code repetition
    """
    if np.min(stocks) < barrier:
        strike = np.average(stocks)
        pay = np.max([stocks[-1]-strike,0])
    else:
        pay = 0
    return pay

def asian_call_avgK(m=100000, S0=1, r=.035, sigma=.3, dt=1/252, barrier=.85, T = 1):
    """Calculates the price of asian down and in call option with given parametrs as default.
    m int: number of stock sims.
    S0 float: intial stock price.
    r float: risk free interest
    sigma float: volatility
    dt float: differnece in time steps
    barrier float: barrier for the option to become active
    T float: time to maturity in years.

    A good idea would be to pass mu gotten from the martingale process as a parameter to make the function generic.
    Here both numeraires are used and returned
    """
    payoffs_B = [] # Payoffs for B as numeraire
    payoffs_S = [] # Payoofs for stock as numeraire
    n = int(T/dt) # Number of stock evaluations
    S0 = float(S0) # Make sure intial stock price is float 
    for i in range(m):
        # initialize stock
        stocks_pB = [S0]
        stocks_nB = [S0]
        stocks_nS = [S0]
        stocks_pS = [S0]
        # sim stock processes
        for j in range(n):
            epsilon = np.random.normal(0,1)
            # stock process with B
            stockpB = stocks_pB[j]*np.exp((r-1/2*sigma**2)*dt + sigma * np.sqrt(dt)*epsilon)
            stocknB = stocks_nB[j]*np.exp((r-1/2*sigma**2)*dt - sigma * np.sqrt(dt)*epsilon) # antitethic option
            stocks_pB.append(stockpB)
            stocks_nB.append(stocknB)
            # stock process with S
            stockpS = stocks_pS[j]*np.exp((r + .5*sigma**2)*dt + sigma*np.sqrt(dt)*epsilon)
            stocknS = stocks_nS[j]*np.exp((r + .5*sigma**2)*dt - sigma*np.sqrt(dt)*epsilon) # antitethic option
            stocks_nS.append(stocknS)
            stocks_pS.append(stockpS)
        # calculate payoffs
        pay_p = option_payoff(stocks_pB, barrier)
        pay_n = option_payoff(stocks_nB, barrier)
        payoffs_B.append(np.exp(-r*T)*pay_p)
        payoffs_B.append(np.exp(-r*T)*pay_n)
        payoffs_S.append(option_payoff(stocks_nS, barrier)/stocknS)
        payoffs_S.append(option_payoff(stocks_pS, barrier)/stockpS)
    ave_pay_B = np.average(payoffs_B) # price of the option with B as the numeraire
    ave_pay_S = np.average(payoffs_S) # price of the option with S as the numeraire
    return ave_pay_B, ave_pay_S

soltion = asian_call_avgK()
print(soltion, soltion[0]-soltion[1]) 
# Did million sims and the results were: B=0.023733, S=0.023724) Diff=8.602830e-06