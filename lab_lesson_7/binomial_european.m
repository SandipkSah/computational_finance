function [price] = binomial_european(S0, K, r, T, sigma, M)
    % Binomial tree method for European option pricing
    % S0: initial stock price
    % K: strike price
    % r: risk-free interest rate
    % T: time to maturity (in years)
    % sigma: volatility
    % N: number of time steps

    dt = T / M;  % time step
    beta = 0.5 * (exp(-r * dt) + exp((r + sigma^2) * dt));  % β calculation
    u = beta + sqrt(beta^2 - 1);  % up factor
    d = beta - sqrt(beta^2 - 1);  % down factor
    p = (exp(r * dt) - d) / (u - d);  % risk-neutral probability

    % Initialize asset prices at maturity
    ST = zeros(M+1, 1);
    for i = 0:M
        ST(i+1) = S0 * (u^i) * (d^(M-i));
    end

    % Initialize option values at maturity
    C = max(0, ST - K);  % Call option payoff

    % Backward induction to calculate option price at t=0
    for j = M:-1:1
        for i = 0:j-1
            C(i+1) = exp(-r * dt) * (p * C(i+2) + (1-p) * C(i+1));
        end
    end

    price = C(1);  % Option price at t=0
end