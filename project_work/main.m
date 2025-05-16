% First calibrate model
initial_params = [0.2206    0.0138    1.0521    1.8432   -0.5638]; % Best results
iterations = 1; % Add iterations interested in finding better parameters
[params, std, loss, p, initial_params] = model_calibration(initial_params, "empVolatilitySurfaceData.mat",iterations, true);
disp(['Heston model params: ', num2str(params)]);
disp(['STD: ', num2str(std)])
disp(['Model loss: ', num2str(loss)])

% Price Down and In Asian call option using calibrated model
S0 = 1;
barrier = .85;
dt = 1/252;
T = 1;
iterations = 1E6;
r = .0466;

[option_price, avg_stock, stocks] = option_pricing(S0, dt, barrier, T, iterations, r, params);

disp(['Down and In Asian call option value: ', num2str(option_price)]);

