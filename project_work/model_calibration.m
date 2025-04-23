function [model_parameters, std, loss, param_matrix, best_initial_params] = model_calibration(initial_params, data_file, iterations, draw_results)
    d = load(data_file);
    data = d.data;
    settings = calibrationSettings;
    param = initial_params;
    
    options = settings.calibrOptions;
    %options.PlotFcns = @optimplotfval;
    marketsurf = data.IVolSurf;
    T = data.T;
    K = data.K;
    [K_grid, T_grid] = meshgrid(K,T);
    Kq = linspace(K(1), K(end),length(K)*2);
    Tq = linspace(T(1), T(end),length(T)*2);
    T_train = Tq(1:3:end);
    K_train = Kq(1:3:end);
    [Ktm, Ttm] = meshgrid(K_train,T_train);
    trainsurf = interp2(K_grid,T_grid,marketsurf, Ktm, Ttm, "spline");
    
    
    fun1 = @(x) lossfunction(x, trainsurf, T_train, K_train, data.S0, data.r, settings);
    ini_matrix = ones(iterations,5);
    param_matrix = ones(iterations,6);
    for j = 1:iterations
        [param_final, fFinal, exitFlag] = fminsearch(fun1, param, options);
        ini_matrix(j,:) = param;
        param_matrix(j,:) = [fFinal, param_final];
        param = create_ini_params(settings);
    end
    [best, index] = min(param_matrix(:,1));
    best_initial_params = ini_matrix(index,:);
    best_param = param_matrix(index, 2:end);
    % best_param = [0.009415411834850   0.024027732885822   6.270622086572633   1.138496430434848  -0.570134584083674]
    std_calc = @(x) calculate_std(x, data, settings.model);
    [std, loss] = std_calc(best_param);
    std = std';
    model = @(x) model_ivs(data.T, data.S0, data.K, data.r, settings.model, x);
    iv = model(best_param);
    if draw_results
        figure
        surf(K_grid,T_grid, marketsurf)
        hold on
        surf(K_grid,T_grid, iv, "EdgeColor", "red")
        hold off
    end
    model_parameters = best_param;

end

    function [std, loss] = calculate_std(param, data, model)
        fun = @(x) model_ivs(data.T, data.S0, data.K, data.r, model, x);
        iv = fun(param);
        mask_iv = ~isnan(iv(:));
        stderror = sum((data.IVolSurf(mask_iv)-iv(mask_iv)).^2);
        J = jacobianest(fun,param);
        J = J(mask_iv,:);
        shape = size(J);
        sigma2 = stderror/(shape(1)-shape(2));
        Jnew = J.'*J;
        diagJ = diag(inv(Jnew));
        loss = mse(data.IVolSurf, iv);
        std = sqrt(sigma2*diagJ);
    end

    function ivs = model_ivs(T,S0, K, r, model, param)
        params = {param(1), param(2), param(3), param(4), param(5)};
        ivs = ones(length(T), length(K));
        for i = 1: length(T)
            call_price = S0.*CallPricingFFT(model,14,S0, K, T(i), r, 0, params{:});
            call_price = max(call_price, 1E-8);
            iv = blsimpv(S0, K, r, T(i), call_price);
            ivs(i,:) = iv;
        end
    end

    function loss = lossfunction(param, surface, T,K,S0, r, settings)
        inside = inbounds(param, settings);
        if inside
            heston = model_ivs(T, S0, K, r, settings.model, param);
            if any(isnan(heston), "all")
                loss = 1E3;
            else
                loss = mse(surface, heston);
            end
        else
            loss = 1E3;
        end
    end

    function inside_bounds = inbounds(param, settings)
        if param(1) < settings.minV0 || param(1) > settings.maxV0
            inside_bounds = false;
        elseif param(3) < settings.minKappa || param(3) > settings.maxKappa
            inside_bounds = false;
        elseif param(2) < settings.minTheta || param(2) > settings.maxKappa
            inside_bounds = false;
        elseif param(4) < settings.minEta || param(4) > settings.maxEta
            inside_bounds = false;
        elseif param(5) < settings.minRho || param(5) > settings.maxRho
            inside_bounds = false;
        else
            inside_bounds = true;
        end
    end

    function params = create_ini_params(settings)
        random_params = rand(5);
        v0 = settings.minV0 + abs(settings.maxV0-settings.minV0)*random_params(1);
        kappa = settings.minKappa + abs(settings.maxKappa-settings.minKappa)*random_params(3);
        theta = settings.minTheta + abs(settings.maxTheta-settings.minTheta)*random_params(2);
        eta = settings.minEta + abs(settings.maxEta-settings.minEta)*random_params(4);
        rho = settings.minRho + abs(settings.maxRho-settings.minRho)*random_params(5);
        params = [v0, theta, kappa, eta, rho];
    end

