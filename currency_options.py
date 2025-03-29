import numpy as np
from scipy import stats

def payoff_down_n_in(stock, strike, barrier):
    if np.min(stock) < barrier:
        return np.max([stock[-1] - strike, 0])
    else:
        return 0
    
def payoff_down_n_out(stock, strike, barrier):
    if np.min(stock) < barrier:
        return 0
    else:
        return np.max([stock[-1] - strike, 0])
    

def problem1(T = 8/12, K = .5, X0 = .52, sigma = .12, rq = .04, rb = .08):
    """
    Calculate the values of an 8 month European call and put options on a currency
    with a strike K of 0.50. The current exchange rate X0 is 0.52, the volatility of the
    exchange rate sigma is 12%, the domestic risk-free rate rq is 4% per annum, and the
    foreign risk-free rate rb is 8% per annum.
    (i) Implement the analytic solution shown at the lecture1
    (ii) Determine the prices
    by using Matlab's blsprice to double check the solution in (i).
    i) H_c(t,x) = x*exp(-rb*dt)*N(d1) - K*exp(-rq*dt)*N(d2)
    F = S-D*K
    C-P = S-DK
    P = -S+DK+C = DK-S+(max(S-DK, 0)) DK > s: C = 0
    P = DK-S
    H_p(t,x) = K*exp(-rq*dt)*N(d1) - x*exp(-rb*dt)*N(d2)
    """
    d1 = (np.log(X0/K) + (rq-rb+1/2*sigma**2)*(T))/(sigma*np.sqrt(T))
    d2 = (np.log(X0/K) + (rq-rb-1/2*sigma**2)*(T))/(sigma*np.sqrt(T))
    call = (X0*np.exp(-rb*T)*stats.norm.cdf(d1)-K*np.exp(-rq*T)*stats.norm.cdf(d2))
    put = (K*np.exp(-rq*T)*stats.norm.cdf(d1)-X0*np.exp(-rb*T)*stats.norm.cdf(d2))
    return call, put

def problem2(m: int):
    """Calculate the price of a down-and-out European call option on non-dividend paying
    stock using Monte Carlo simulation (as you learned in the previous exercises). 
    Current stock price S0 is $50, the strike price K is $50, the risk-free rate r is 5% p.a.
    the volatility sigma is 30%, time to maturity T is one year and the barrier H is $45.
    Check that your code is pricing the call option correctly by calculating the price
    of corresponding down-and-in option and comparing the sum of the prices to the
    price of plain vanilla European call (remember that a plain vanilla European call
    option is worth exactly the same as the sum of corresponding down-and-out and
    down-and-in call options).
    """
    S0 = 50
    K = 50
    r = .05
    sigma = .3
    barrier = 45
    T = 1
    days = 252
    d1 = (np.log(S0/K) + (r+1/2*sigma**2)*(T))/(sigma*np.sqrt(T))
    d2 = (np.log(S0/K) + (r-1/2*sigma**2)*(T))/(sigma*np.sqrt(T))
    vanilla_call = S0*stats.norm.cdf(d1)-np.exp(-r*T)*K*stats.norm.cdf(d2)
    dno_payoffs = []
    dni_payoffs = []
    for i in range(m):
        epsilon = np.random.normal(0,1,days)
        dynamics = np.ones(days)*((r-1/2*sigma**2)*1/days + epsilon*sigma*np.sqrt(1/days))
        dynamicsn = np.ones(days)*((r-1/2*sigma**2)*1/days - epsilon*sigma*np.sqrt(1/days))
        ST = S0*np.exp(np.cumsum(dynamics))
        STn = S0*np.exp(np.cumsum(dynamicsn))
        diskount = np.exp(-r*T)
        pay_dno = diskount*payoff_down_n_out(ST,K, barrier)
        pay_dni = diskount*payoff_down_n_in(ST, K, barrier)
        payn_dno = diskount*payoff_down_n_out(STn, K, barrier)
        payn_dni = diskount*payoff_down_n_in(STn,K, barrier)
        dno_payoffs.append(pay_dno)
        dno_payoffs.append(payn_dno)
        dni_payoffs.append(pay_dni)
        dni_payoffs.append(payn_dni)
    dni_call = np.average(dni_payoffs)
    dno_call = np.average(dno_payoffs)
    return dni_call, dno_call, vanilla_call

def problem3(m: int):
    """Consider the following structured knock-out option. The option has five observation
    dates ti = 1, 2, 3, 4, 5 years after t0. On those dates the payoff is max(Sti - K, 0).
    Moreover, if Sti < 0.7 * St0, the contract is terminated after paying the (possible)
    payoff for that observation date.
    Use Monte Carlo simulation to price the contract. St0 = $100, strike price is $110,
    risk free interest rate 3 % and volatility 30 %. Calculate also the percentages of how
    many times the contract was terminated at each observation date.
    """
    ti = [1,2,3,4,5]
    payoff = lambda s,k: np.max([s-k,0])
    S0 = 100
    K = 110
    barrier = S0*.7
    r = .03
    sigma = .3
    dt = 1
    kills = 0
    vals = 0
    options = []
    for i in range(m):
        pays = 0
        paysn = 0
        stocks = [S0]
        stocksn = [S0]
        for year in range(len(ti)):
            vals += 2
            epsilon = np.random.normal(0,1)
            stok = stocks[year]*np.exp((r-1/2*sigma**2)*dt + sigma*epsilon*np.sqrt(dt))
            stokn = stocksn[year]*np.exp((r-1/2*sigma**2)*dt - epsilon *sigma *np.sqrt(dt))
            call = payoff(stok, K)
            pays += np.exp(-r*year+1)*call
            if stok !=0 and stok < barrier:
                kills += 1
                stok = 0
            calln = payoff(stokn, K)
            paysn += np.exp(-r*year+1)*calln
            if stokn !=0 and stokn < barrier:
                kills += 1
                stokn = 0
            stocks.append(stok)
            stocksn.append(stokn)
        options.append(pays)
        options.append(paysn)
    return np.average(options), kills/vals

print(problem3(int(10000)))
