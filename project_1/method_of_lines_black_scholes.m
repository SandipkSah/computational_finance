function method_of_lines_black_scholes(NS, Nt)
    % Parameters
    r = 0.06;      % Risk-free interest rate
    sigma = 0.3;   % Volatility
    T = 1;         % Time to maturity
    K = 10;        % Strike price
    S_max = 15;    % Maximum stock price
    
    % Discretization
    dS = S_max / NS;            % Spatial step size
    dt = T / Nt;                % Time step size
    S = linspace(0, S_max, NS+1); % Stock price grid
    t = linspace(0, T, Nt+1);    % Time grid
    
    % Initial condition (European put option payoff)
    % U(:, 1) = max(K - S, 0); % For European Put Option
    U(:, 1) = max(S - K, 0); % For European Call Option (uncomment this line for call)
    
    
    % Solve using RK4 method
    for j = 1:Nt
        U(:, j+1) = RK4(@(t, U) function_f(t, U, S, r, sigma, dS, dt), dt, U(:, j), (Nt-j)*dt);
    end
    
    % Plot results
    figure;
    surf(t, S, U);
    xlabel('Time to Maturity t');
    ylabel('Stock Price S');
    zlabel('Option Price U');
    title('Method of Lines with RK4 for Option');
    shading interp;
end

% Function for Black-Scholes ODE system
function dUdt = function_f(t, U, S, r, sigma, dS, dt)
    NS = length(S);
    dUdt = zeros(size(U));
    
    % Finite difference approximations for derivatives
    dU_dS = zeros(size(U));
    d2U_dS2 = zeros(size(U));
    
    for i = 2:NS-1
        dU_dS(i) = (U(i+1) - U(i-1)) / (2 * dS);
        d2U_dS2(i) = (U(i+1) - 2*U(i) + U(i-1)) / (dS^2);
    end
    
    % Boundary conditions
    dU_dS(1) = -1;  % At S = 0, dU/dS = -1 (for put option)
    dU_dS(NS) = 0;  % At S = S_max, dU/dS = 0
    
    % Black-Scholes ODE system
    for i = 2:NS-1
        dUdt(i) = - (0.5 * sigma^2 * S(i)^2 * d2U_dS2(i) + r * S(i) * dU_dS(i) - r * U(i));
    end
end

% Fourth-Order Runge-Kutta Method
function U_new = RK4(f, dt, U, t)
    % this is probably wrong
    k1 = f(t, U);
    k2 = f(t + dt/2, U + dt/2 * k1);
    k3 = f(t + dt/2, U + dt/2 * k2);
    k4 = f(t + dt, U + dt * k3);
    U_new = U + dt * (k1 + 2*k2 + 2*k3 + k4) / 6;
end
