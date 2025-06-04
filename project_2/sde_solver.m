function [t, X_em, X_milstein] = sde_solver(a, b, dbdx, x0, T, N)
    % Inputs:
    % a     - function handle a(t,x)
    % b     - function handle b(t,x)
    % dbdx  - function handle for ∂b/∂x(t,x) (required for Milstein)
    % x0    - initial condition
    % T     - final time
    % N     - number of steps

    % Time grid
    h = T / N;
    t = linspace(0, T, N+1);

    % Pre-allocate solution arrays
    X_em = zeros(1, N+1);     % Euler-Maruyama
    X_milstein = zeros(1, N+1); % Milstein

    % Initial conditions
    X_em(1) = x0;
    X_milstein(1) = x0;

    % Generate Brownian increments
    dW = sqrt(h) * randn(1, N);

    % Time stepping
    for n = 1:N
        tn = t(n);
        Xn_em = X_em(n);
        Xn_mil = X_milstein(n);
        dWn = dW(n);

        % Euler-Maruyama
        X_em(n+1) = Xn_em + a(tn, Xn_em) * h + b(tn, Xn_em) * dWn;

        % Milstein
        X_milstein(n+1) = Xn_mil + a(tn, Xn_mil) * h + b(tn, Xn_mil) * dWn ...
                          + 0.5 * b(tn, Xn_mil) * dbdx(tn, Xn_mil) * (dWn^2 - h);
    end
end
