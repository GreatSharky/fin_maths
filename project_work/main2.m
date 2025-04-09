
S0 = 1;
K = [0.7:0.1:1.3];
r = 0;
T = 1;

model = 'Heston';

V0 = 0.04;
theta = 0.05;                   % long term variance
kappa = 6;                      % mean reversion speed of variance
lambda = 0;
kappaQ = kappa + lambda;
thetaQ = kappa*theta/(kappa + lambda);
eta =  0.5;                     % volatility of variance
rho = -0.8;                     % correlation between returns and variance

parameters = {V0, thetaQ kappaQ, eta, rho};

n = 13;
callPrices = ...
    S0.*CallPricingFFT(model, n, ...
    1, K./S0, T, r, 0, ...
    parameters{:});

putPrices = callPrices - S0 + K.*exp(-r*T);

disp(['Heston call: ', num2str(callPrices)]);
disp(['Heston put: ', num2str(putPrices)]);

disp(['BS call: ', num2str(blsprice(S0, K, r, T, sqrt(V0)))]);
disp(['BS put: ', num2str(blsprice(S0, K, r, T, sqrt(V0)) - S0 + K.*exp(-r*T))]);
