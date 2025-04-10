% main.m
% Computational Project 1 — Computational Methods in Finance
% Author: Sandip Sah, Ugnius Braun, Joao Maria Relvas 
% Date: [April 2025]
%
% This script solves:
% 1. A European put option using the Method of Lines with 4th order Runge-Kutta (Part 1)
% 2. An American put option using the Crank-Nicolson scheme with PSOR (Part 2)
% All numerical results and figures will be generated here and saved to ./graph

clc; clear;

output_folder = 'graph'; % Folder to save plots

if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end


%% ------------------------------------------------------------------------
%  PART 1 — European Put Option using Method of Lines + RK4
% -------------------------------------------------------------------------
disp('=== Part 1: Solving European Put Option using Method of Lines + RK4 ===');

% Parameters for Part 1
r = 0.06;           % Risk-free interest rate
sigma = 0.3;        % Volatility
K = 10;             % Strike price
T = 1;              % Maturity
Smax = 15;          % Maximum asset price
Ns = 100;           % Spatial steps
Nt = 1000;          % Time steps

% Solve using function
[S1, tau1, U_all] = solve_european_put_MOL(r, sigma, K, T, Smax, Ns, Nt);

% Plot result
fig1 = figure;
mesh(S1, tau1, U_all')
title('Part 1: European Put Option U(S, t) — Method of Lines + RK4')
xlabel('Asset Price S'), ylabel('Time t'), zlabel('Option Value U')
colormap turbo
view(135, 30)
saveas(fig1, fullfile(output_folder, 'european_put_surface.png'));

%% ------------------------------------------------------------------------
%  PART 2 — American Put Option using Crank-Nicolson + PSOR
% -------------------------------------------------------------------------
disp('=== Part 2: Solving American Put Option using Crank-Nicolson + PSOR ===');

% Parameters for Part 2 (can reuse some from Part 1)
Ns = 200;
Nt = 1000;
omega = 1.2;         % Relaxation factor for PSOR
tol = 1e-6;          % Convergence tolerance
max_iter = 10000;    % Max iterations for PSOR

% Solve using function
[S2, tau2, U_eu, U_am] = solve_american_put_CN_PSOR(r, sigma, K, T, Smax, Ns, Nt, omega, tol, max_iter);

% Plot European put result
fig2 = figure;
mesh(S2, tau2, U_eu')
title('Part 2: European Put Option U(S, t) — Crank-Nicolson')
xlabel('Asset Price S'), ylabel('Time t'), zlabel('Option Value U')
colormap turbo
view(135, 30)
saveas(fig2, fullfile(output_folder, 'american_european_surface.png'));

% Plot American put result
fig3 = figure;
mesh(S2, tau2, U_am')
title('Part 2: American Put Option U(S, t) — Crank-Nicolson + PSOR')
xlabel('Asset Price S'), ylabel('Time t'), zlabel('Option Value U')
colormap turbo
view(135, 30)
saveas(fig3, fullfile(output_folder, 'american_put_surface.png'));

% Plot early exercise premium
fig4 = figure;
mesh(S2, tau2, U_am' - U_eu')
title('Part 2: Early Exercise Premium (American - European)')
xlabel('Asset Price S'), ylabel('Time t'), zlabel('Premium Value')
colormap turbo
view(135, 30)
saveas(fig4, fullfile(output_folder, 'early_exercise_premium.png'));

disp('=== All computations completed. Figures saved in "graph" directory. ===');
