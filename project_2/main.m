% MAIN SCRIPT TO RUN PROJECT TASKS 1–3
clc; clear; close all;

if ~exist('output', 'dir')
    mkdir('output');
end

% %% Task 1: Random number generators
% fprintf('--- Running Random Generators ---\n');
% random_generators();

% %% Task 2: Estimate Mandelbrot area using MC and QMC
% fprintf('\n--- Estimating Mandelbrot Area ---\n');
% estimate_mandelbrot_area();

%% Task 3: Solve SDE using Euler-Maruyama and Milstein methods
fprintf('\n--- Solving SDE using Euler-Maruyama and Milstein ---\n');

% Example SDE: dX = 0.5*X dt + 0.3*X dB
a = @(t, x) 0.5 * x;
b = @(t, x) 0.3 * x;
dbdx = @(t, x) 0.3;  % ∂b/∂x

x0 = 1;       % Initial condition
T = 1;        % Final time
N = 1000;     % Number of time steps

[t, X_em, X_milstein] = sde_solver(a, b, dbdx, x0, T, N);

% Plot the solutions
fig = figure;
hold on;
plot(t, X_em, 'b-', 'DisplayName', 'Euler-Maruyama');
plot(t, X_milstein, 'r--', 'DisplayName', 'Milstein');
xlabel('t'); ylabel('X(t)');
title('SDE Solution: Euler-Maruyama vs Milstein');
legend show;
grid on;
hold off;

% Save the plot
saveas(fig, fullfile('output', 'sde_solution.png'));

% %% Task 4: Convergence analysis of SDE solvers
% fprintf('\n--- Task 4(a): Plot Single Path ---\n');
% task4a_single_simulation();

% fprintf('\n--- Task 4(b): Run Convergence Study (500k paths) ---\n');
% task4b_convergence_analysis();