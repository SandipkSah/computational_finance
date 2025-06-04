function task4a_single_simulation()
    rng(123);  % Fixed seed

    %% Create output directory if it doesn't exist
    if ~exist('output', 'dir')
        mkdir('output');
    end

    % Parameters
    mu = 0.5;
    sigma = 0.3;
    S0 = 1;
    T = 1;
    h = 0.001;
    N = T / h;
    t = linspace(0, T, N+1);

    % Generate Brownian path
    dW = sqrt(h) * randn(1, N);
    W = [0, cumsum(dW)];

    % Exact solution
    S_exact = S0 * exp((mu - 0.5 * sigma^2) * t + sigma * W);

    % Define drift and diffusion
    a = @(t, x) mu * x;
    b = @(t, x) sigma * x;
    dbdx = @(t, x) sigma;

    % Simulate numerically
    [~, S_em, S_milstein] = sde_solver_given_path(a, b, dbdx, S0, T, N, dW);

    % Plot
    figure;
    plot(t, S_exact, 'k', 'DisplayName', 'Exact Solution'); hold on;
    plot(t, S_em, 'b--', 'DisplayName', 'Euler-Maruyama');
    plot(t, S_milstein, 'r-.', 'DisplayName', 'Milstein');
    xlabel('t'); ylabel('S(t)');
    title('Geometric Brownian Motion: Numerical vs Exact');
    legend show; grid on;

    % Save plot
    saveas(gcf, fullfile('output', 'task4a_gbm_single_path.png'));
end
