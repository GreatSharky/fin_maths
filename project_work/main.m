

% Optimization settings
OptSettings.MaxFunEvals = 500;
OptSettings.MaxIter = 200;
OptSettings.TolFun = 1e-10;
OptSettings.TolX = 1e-10;
OptSettings.Display = 'iter';
OptSettings.FunValCheck = 'on';

% Function to be minimized
fun = @(x)lossFunction(x);

% Initial values 
param0 = [-4, -2];

% Optimization
[param_final, fFinal, exitFlag] = fminsearch(fun, param0, OptSettings);
%[param_final, fFinal, exitFlag] = fminunc(fun, param0, OptSettings);


disp(['Starting values: ', num2str(param0)]);
disp(['Optimized values: ', num2str(param_final)]);
