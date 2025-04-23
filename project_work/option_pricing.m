
function [option_price, avg_stock, stokcs] = option_pricing(S0, dt, H, T, M, r, model_parameters)
    N = T/dt;
    V0 = model_parameters(1);
    theta = model_parameters(2);
    kappa = model_parameters(3);
    eta = model_parameters(4);
    rho = model_parameters(5);
    
    stokcs = ones(2*M, N);
    options = ones(2*M,1);
    
    for i = 1:M
        stock_price = stock_sim(S0, r, V0, theta, kappa, eta, rho, dt, N);
        payoff = exp(-r*T)*option_value(stock_price(:,1),H);
        payoff_a = exp(-r*T)*option_value(stock_price(:,2),H);
        stokcs(i,:) = stock_price(:,1)';
        stokcs(i+M, :) = stock_price(:,2)';
        options(i,:) = payoff;
        options(i+M,:) = payoff_a;
    end
    avg_stock = mean(stokcs(:,end));
    option_price = mean(options(:,end));

    function payoff = option_value(S,H)
        if min(S) > H
            payoff = 0;
        else
            K = mean(S);
            payoff = max(S(end)-K,0);
        end
    end
    
    function stock_prices = stock_sim(S0, r, V0, theta, kappa, eta, rho, dt, days)
        S = ones(days,2);
        epsilon = randn(days,2);
        V = ones(days,2);
        V(1,1) = heston_volatility(V0, theta, kappa, eta, rho, epsilon(1,1), dt, epsilon(1,2));
        V(1,2) = heston_volatility(V0, theta, kappa, eta, rho, -epsilon(1,1), dt, -epsilon(1,2));
        % V(:,:) =.3;
        S(1,1) = S0*exp((r-1/2*V(1,1))*dt + epsilon(1,1)*V(1,1)*sqrt(dt));
        S(1,2) = S0*exp((r-1/2*V(1,2))*dt - epsilon(1,1)*V(1,2)*sqrt(dt));
        for d = 2:days
            V(d,1) = heston_volatility(V(d-1,1), theta, kappa, eta, rho, epsilon(d,1), dt, epsilon(d,2));
            V(d,2) = heston_volatility(V(d-1,2), theta, kappa, eta, rho, -epsilon(d,1), dt, -epsilon(d,2));
            % V(:,:) = .3;
            S(d,1) = S(d-1,1)*exp((r-1/2*V(d,1))*dt + epsilon(d,1)*V(d,1)*sqrt(dt));
            S(d,2) = S(d-1,2)*exp((r-1/2*V(d,2))*dt - epsilon(d,1)*V(d,2)*sqrt(dt));
        end
        stock_prices = S;
    end
    
    function volatility = heston_volatility(V, theta, kappa, eta, rho, epsilon, dt, e)
        epsilon2 = rho*epsilon + sqrt(1-rho^2)*e;
        volatility = V + kappa*(theta-V)*dt + eta*sqrt(V*dt)*epsilon2 + 1/4*eta^2*dt*(epsilon2^2-1);
        volatility = max(volatility,0);
    end

end