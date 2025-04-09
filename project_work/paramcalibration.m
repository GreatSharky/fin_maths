load("empVolatilitySurfaceData.mat");
settings = calibrationSettings;
kappa = settings.parameters0(1);
theta = settings.parameters0(2);
eta = settings.parameters0(3);
rho = settings.parameters0(4);
V0 = settings.parameters0(5);
param = [V0, kappa, theta, eta, rho];


options = settings.calibrOptions;
fun = @(x) lossfunction(x, data, settings);
for i = 1:100
    param = create_ini_params(param, settings);
    [param_final, fFinal, exitFlag] = fminsearch(fun, param, options);
    loss_matrix(i) = fFinal;
    param_matrix(i,:) = [fFinal, param_final];
    break
end
% Best_param = 0.020884593943414   5.347473692942446   0.000682258718764   0.286788948704852  -0.602970312709901
std = calculate_std(param_final, data, settings.model)

function std = calculate_std(param, data, model)
    T = data.T;
    S0 = data.S0;
    K = data.K;
    r = data.r;
    params = {param(1), param(2), param(3), param(4), param(5)};
    f =@(x) CallPricingFFT(model, 13, S0, K(4), T(4), r, 0, x{:});
    model_iv = ones(length(data.T), length(data.K));
    for i = 1 : length(data.T)
        call = abs(S0.*CallPricingFFT(model, 13, S0, K, T(i), r, 0, params{:}));
        model_iv(i,:) = blsimpv(S0, K, data.r, T(i), call);
    end
    error = (data.IVolSurf(:) - model_iv(:)).^2;
    J = jacobianest(f,params);
    sigma2 = error/(length(error)-length(params));
    std = sqrt(sigma2*diag(J/J));
end


function hestonIV = ivmse(model, S0, K, T, r, parameters)
    hestonIV = ones(length(T), length(K));
    for i = 1:length(T)
        call = abs(S0.*CallPricingFFT(model, 13, S0, K, T(i), r, 0, parameters{:}));
        implied_volatility = blsimpv(S0,K,r,T(i),call);
        hestonIV(i,:) = implied_volatility;
    end
end

function loss = lossfunction(param, data, settings)
    inside = inbounds(param, settings);
    if inside
        r = data.r;
        S0 = data.S0;
        K = data.K;
        T = data.T;
        parameters = {param(1), param(2), param(3), param(4), param(5)};
        heston = ivmse(settings.model, S0, K, T, r, parameters);
        loss = mse(data.IVolSurf, heston);
    else
        loss = 1E6;
    end
end

function inside_bounds = inbounds(param, settings)
    if param(1) < settings.minV0 || param(1) > settings.maxV0
        inside_bounds = false;
    elseif param(2) < settings.minKappa || param(2) > settings.maxKappa
        inside_bounds = false;
    elseif param(3) < settings.minTheta || param(3) > settings.maxKappa
        inside_bounds = false;
    elseif param(4) < settings.minEta || param(4) > settings.maxEta
        inside_bounds = false;
    elseif param(5) < settings.minRho || param(5) > settings.maxRho
        inside_bounds = false;
    else
        inside_bounds = true;
    end
end

function params = create_ini_params(params, settings)
    random_params = rand(5);
    v0 = settings.minV0 + abs(settings.maxV0-settings.minV0)*random_params(1);
    kappa = settings.minKappa + abs(settings.maxKappa-settings.minKappa)*random_params(2);
    theta = settings.minTheta + abs(settings.maxTheta-settings.minTheta)*random_params(3);
    eta = settings.minEta + abs(settings.maxEta-settings.minEta)*random_params(4);
    rho = settings.minRho + abs(settings.maxRho-settings.minRho)*random_params(5);
    params = [v0, kappa, theta, eta, rho];
end
