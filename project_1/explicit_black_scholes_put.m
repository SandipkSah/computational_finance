function explicit_black_scholes_put(NS, Nt)
    % EXPLICIT_BLACK_SCHOLES_PUT Solves the Black-Scholes equation for a European put option
    % using the explicit finite difference method.
    %
    %   explicit_black_scholes_put(NS, Nt) solves the Black-Scholes PDE for a European put option
    %   using NS spatial steps and Nt time steps.
    %
    %   Parameters:
    %       r  = 0.06 (Risk-free interest rate)
    %       sigma = 0.3 (Volatility)
    %       T = 1 (Time to maturity)
    %       K = 10 (Strike price)
    %       S_max = 15 (Upper limit for stock price domain)
    %
    %   The solution is visualized as a surface plot.

    % Define parameters
    r = 0.06;
    sigma = 0.3;
    T = 1;
    K = 10;
    S_max = 15;

    % Discretization
    dS = S_max / NS;
    dt = T / Nt;
    S = linspace(0, S_max, NS+1); % Stock price grid
    t = linspace(0, T, Nt+1); % Time grid
    
    % Stability condition check
    if dt > (dS^2 / (sigma^2 * S_max^2))
        error('Stability condition violated: Choose smaller dt or larger dS.');
    end
    
    % Initialize solution matrix
    U = zeros(NS+1, Nt+1);
    
    % Set initial condition: Payoff at maturity for put option
    U(:, 1) = max(K - S, 0);
    
    % Boundary conditions
    U(1, :) = K * exp(-r * (T - t)); % U(0,t) = K * exp(-r(T-t))
    U(end, :) = 0; % U(S_max, t) = 0
    
    % Time-stepping loop (Explicit Finite Difference Method)
    for j = 1:Nt 
        for i = 2:NS
            Si = S(i);
            i_param = Si/dS;
            a = 0.5 * dt * (sigma^2 * i_param^2 - r * i_param);
            b = 1 - sigma^2 * i_param^2 * dt - r * dt;
            c = 0.5 * dt * (sigma^2 * i_param^2 + r * i_param);
            U(i, j+1) = a * U(i-1, j) + b * U(i, j) + c * U(i+1, j);
        end
    end

    % Plot the computed solution
    figure;
    surf(t, S, U);
    xlabel('Time to Maturity t');
    ylabel('Stock Price S');
    zlabel('Option Price U');
    title('Explicit Finite Difference Method for European Put Option');
    shading interp;
end
