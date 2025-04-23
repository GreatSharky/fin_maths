% First calibrate model

[params, std, loss, p] = model_calibration(1, "empVolatilitySurfaceData.mat", true);
disp(['Heston model params: ', num2str(params)]);
disp(['STD: ', num2str(std)])
disp(['Model loss: ', num2str(loss)])

% Price Down and In Asian call option
S0 = 1;
barrier = .85;
dt = 1/252;
T = 1;
iterations = 1E6;
r = .0466;

[option_price, avg_stock, stocks] = option_pricing(S0, dt, barrier, T, iterations, r, params);

disp(['Heston call: ', num2str(option_price)]);

