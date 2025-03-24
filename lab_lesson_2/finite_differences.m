function [forward, backward, centered] = finite_differences(h)
    % Function to compute finite difference approximations of f'(0)
    
    % Define f(x) = sin(x)
    f = @(x) sin(x);

    % Compute forward, backward, and centered difference approximations
    forward = (f(h) - f(0)) / h;
    backward = (f(0) - f(-h)) / h;
    centered = (f(h) - f(-h)) / (2 * h);
end
