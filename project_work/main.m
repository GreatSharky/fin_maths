% First calibrate model
initial_params = [0.2206    0.0138    1.0521    1.8432   -0.5638]; % Best results
iterations = 1; % Add iterations interested in finding better parameters
[params, std, loss] = model_calibration(initial_params, "empVolatilitySurfaceData.mat",iterations, true);
% params = [0.0094    0.0240    6.2693    1.1382   -0.5700]
% std = [0.0017    0.0003    0.3438    0.0423    0.0079]
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

