function task4b_convergence_analysis()
    rng(42);
    mu = 0.5; sigma = 0.3; S0 = 1; T = 1;
    dbdx = @(t, x) sigma;
    a = @(t, x) mu * x;
    b = @(t, x) sigma * x;

    hs = 0.005 * (0.5).^(0:3);   % h = 0.005, 0.0025, 0.00125, 0.000625
    M = 500000;                 % Number of simulations

    strong_em = zeros(1, length(hs));
    weak_em = zeros(1, length(hs));
    strong_mil = zeros(1, length(hs));
    weak_mil = zeros(1, length(hs));

    for i = 1:length(hs)
        h = hs(i);
        N = T / h;

        % Preallocate arrays
        err_em = zeros(1, M);
        err_mil = zeros(1, M);
        exact_vals = zeros(1, M);
        em_vals = zeros(1, M);
        mil_vals = zeros(1, M);

        fprintf("Running simulations for h = %.5f...\n", h);

        parfor j = 1:M
            dW = sqrt(h) * randn(1, N);
            W = sum(dW);
            S_exact = S0 * exp((mu - 0.5 * sigma^2) * T + sigma * W);

            [~, S_em, S_milstein] = sde_solver_given_path(a, b, dbdx, S0, T, N, dW);

            S_em_final = S_em(length(S_em));             % avoid using `end`
            S_mil_final = S_milstein(length(S_milstein));

            em_vals(j) = S_em_final;
            mil_vals(j) = S_mil_final;
            exact_vals(j) = S_exact;

            err_em(j) = abs(S_em_final - S_exact);
            err_mil(j) = abs(S_mil_final - S_exact);
        end

        % Compute strong and weak errors
        strong_em(i) = mean(err_em);
        strong_mil(i) = mean(err_mil);
        weak_em(i) = abs(mean(em_vals) - mean(exact_vals));
        weak_mil(i) = abs(mean(mil_vals) - mean(exact_vals));

        fprintf("h = %.5f | Strong EM: %.5e | Weak EM: %.5e | Strong Milstein: %.5e | Weak Milstein: %.5e\n", ...
            h, strong_em(i), weak_em(i), strong_mil(i), weak_mil(i));
    end

    % Plot convergence rates
    figure;
    loglog(hs, strong_em, 'o-', 'DisplayName', 'Strong EM'); hold on;
    loglog(hs, strong_mil, 's--', 'DisplayName', 'Strong Milstein');
    loglog(hs, weak_em, 'd-', 'DisplayName', 'Weak EM');
    loglog(hs, weak_mil, '^--', 'DisplayName', 'Weak Milstein');
    xlabel('Step size h'); ylabel('Error');
    legend show; grid on;
    title('Strong and Weak Convergence of SDE Solvers');

    % Save figure to output directory
    if ~exist('output', 'dir')
        mkdir('output');
    end
    saveas(gcf, fullfile('output', 'task4b_convergence_plot.png'));
end
