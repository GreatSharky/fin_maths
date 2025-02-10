import numpy as np
import math
import csv

payoff = lambda x : np.max([x,0])

def read_data(file_name, key):
    with open(file_name, "r") as file:
        reader = csv.DictReader(file)
        data = []
        for row in reader:
            data.append(float(row[key]))

    return data

def binompmf(n,k,p):
    nck = math.factorial(n)/(math.factorial(k)*math.factorial(n-k))
    return nck*p**k*(1-p)**(n-k)

def get_UD(data):
    days = len(data) # 251
    N = days-1  # 250
    dt = 1/days
    data = np.array(data)
    v_hat_d = 1/(N)*np.log(data[-1]/data[0])
    v_hat = v_hat_d*days
    sigma2_hat_d = 1/(N-1)*np.sum((np.log(data[1:]/data[:-1])-v_hat_d)**2)
    sigma2 = sigma2_hat_d*days
    sigma = np.sqrt(sigma2)
    q_exact = 1/2 +1/2*1/np.sqrt(sigma2/(v_hat**2*dt)+1)
    upper_exact = np.exp(np.sqrt(sigma2*dt+(v_hat*dt)*2))
    lower_exact = 1/upper_exact
    q_approx = 1/2+1/2*v_hat/sigma*np.sqrt(dt)
    upper_approx = np.exp(np.sqrt(dt*sigma2))
    lower_approx = 1/upper_approx
    u_diff = np.abs(upper_exact-upper_approx)
    d_diff = np.abs(lower_approx-lower_exact)
    print(q_exact)
    print(upper_exact,u_diff, u_diff/upper_exact)
    print(lower_exact,d_diff,d_diff/lower_exact)
    return upper_exact, lower_exact, upper_approx, lower_approx


def option_price_with_data(file_name, file_key, strike, risk_free_interest, time_to_maturity):
    data = read_data(file_name, file_key)
    data = np.array(data)
    Ue, De, Ua, Da = get_UD(data)
    sum_approx = 0
    sum_exact = 0
    m = len(data)
    dt = time_to_maturity/m
    diskont = np.exp(-risk_free_interest*time_to_maturity)
    pa = (np.exp(risk_free_interest*dt)-Da)/(Ua-Da)
    pe = (np.exp(risk_free_interest*dt)-De)/(Ue-De)
    for k in range(m):
        sum_approx += payoff(strike - data[0]*Ua**k *Da**(m-k))
        sum_exact += payoff(strike - data[0]*Ue**k *De**(m-k))

    put_approx = sum_approx*diskont
    put_exact = sum_exact*diskont
    return put_approx, put_exact

def main():
    K = 100
    file = "IBM.csv"
    key = "Adj Close"
    rf = .03
    T = 1
    put_approx, put_exact = option_price_with_data(file, key, K, rf,T)
    print(put_approx)
    print(put_exact)

main()

# Ue =1.0146 
# Daily = .007
# Annual .185