import numpy as np
from scipy.optimize import minimize

def loss_function_mishra(param):
    """
    Computes Mishra's Bird function with a constraint.
    """
    x, y = param
    
    # Define the function
    f = lambda x, y: np.sin(y) * np.exp((1 - np.cos(x))**2) + np.cos(x) * np.exp((1 - np.sin(y))**2) + (x - y)**2
    
    # Constraint check
    in_bounds = (x + 5)**2 + (y + 5)**2 < 25
    
    return f(x, y) if in_bounds else 1e10

def loss_function_beale(param):
    x, y = param
    f = lambda x,y: (1.5-x+x*y)**2 + (2.25 -x +x*y**2)**2 + (2.625-x + x*y**3)**2
    if x >= -4.5 and y <= 4.5:
        inbounds = True
    else:
        inbounds = False
    return f(x,y) if inbounds else 1e10

# Optimization settings
opt_settings = {
    'maxiter': 200,
    'disp': True,
    'fatol': 1e-10,
    'xatol': 1e-10
}

# Function to be minimized
fun = loss_function_beale

# Initial values 
param0 = [4, -4]

# Optimization
res = minimize(fun, param0, method='Nelder-Mead', options=opt_settings)

# Results
print(f"Starting values: {param0}")
print(f"Optimized values: {res.x}")
