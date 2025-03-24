clc; clear;

% Choose values of N
% N_values = [10, 20, 40, 80, 160];
N_values = [10, 20];
errors = zeros(size(N_values));

for k = 1:length(N_values)
    N = N_values(k);
    x = linspace(0, 1, N+1); % Grid points
    u_exact = (1/9) * (13 + 4*x.*(-1 + cos(3)) - 4*cos(3*x)); % Exact solution
    u_num = solve_bvp(N); % Numerical solution
    
    % Compute max norm error
    errors(k) = max(abs(u_num - u_exact'));
end

% Estimate order of convergence using polyfit (log-log regression)
p = polyfit(log(N_values), log(errors), 1);
order_of_convergence = -p(1);

% Display results
fprintf('Estimated order of convergence: %.2f\n', order_of_convergence);

% Plot error vs N in log-log scale
figure;
loglog(N_values, errors, 'r-o', 'LineWidth', 1.5);
xlabel('N');
ylabel('Error');
title('Convergence Study');
grid on;


