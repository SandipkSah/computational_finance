clc; clear;

% Define the exact derivative of f(x) = sin(x) at x = 0
exact_derivative = cos(0); % cos(0) = 1

% Define the range of i values and compute corresponding h values
i_values = 10:10:100;       % Range of i values
h_values = 1 ./ i_values;   % Compute h_i = 1/i

% Initialize arrays to store errors for each finite difference method
forward_errors = zeros(size(h_values));   % Errors for forward difference
backward_errors = zeros(size(h_values));  % Errors for backward difference
centered_errors = zeros(size(h_values));  % Errors for centered difference

% Loop through each h value to compute approximations and errors
for k = 1:length(h_values)
    h = h_values(k); % Current step size
    
    % Compute finite difference approximations
    [forward, backward, centered] = finite_differences(h);
    
    % Compute absolute errors for each method
    forward_errors(k) = abs(forward - exact_derivative);
    backward_errors(k) = abs(backward - exact_derivative);
    centered_errors(k) = abs(centered - exact_derivative);
end

% Display the results in a formatted table
disp(' i      h      Forward Error   Backward Error   Centered Error');
disp('------------------------------------------------------------------');
for k = 1:length(h_values)
    fprintf('%3d  %8.5f  %12.8f  %12.8f  %12.8f\n', ...
            i_values(k), h_values(k), ...
            forward_errors(k), backward_errors(k), centered_errors(k));
end

% Estimate the order of convergence using linear regression on log-log scale
log_h = log(h_values); % Log of step sizes

% Fit linear models to log(error) vs log(h)
p_forward = polyfit(log_h, log(forward_errors), 1);   % Forward difference slope
p_backward = polyfit(log_h, log(backward_errors), 1); % Backward difference slope
p_centered = polyfit(log_h, log(centered_errors), 1); % Centered difference slope

% Display the estimated orders of convergence
fprintf('\nEstimated order of convergence:\n');
fprintf('Forward Difference: %.2f\n', p_forward(1));
fprintf('Backward Difference: %.2f\n', p_backward(1));
fprintf('Centered Difference: %.2f\n', p_centered(1));

% Plot the errors in log-log scale for visualization
figure;
loglog(h_values, forward_errors, 'r-o', 'DisplayName', 'Forward');
hold on;
loglog(h_values, backward_errors, 'b-s', 'DisplayName', 'Backward');
loglog(h_values, centered_errors, 'g-d', 'DisplayName', 'Centered');
xlabel('Step Size (h)');
ylabel('Absolute Error');
title('Finite Difference Errors vs Step Size (h)');
legend;
grid on;

%% Helper Function: Compute Finite Differences
function [forward, backward, centered] = finite_differences(h)
    % Compute finite difference approximations for f(x) = sin(x) at x=0.
    %
    % Input:
    %   h - Step size.
    %
    % Output:
    %   forward - Forward difference approximation.
    %   backward - Backward difference approximation.
    %   centered - Centered difference approximation.
    
    forward = (sin(h) - sin(0)) / h;          % Forward difference formula
    backward = (sin(0) - sin(-h)) / h;        % Backward difference formula
    centered = (sin(h) - sin(-h)) / (2 * h);  % Centered difference formula
end
