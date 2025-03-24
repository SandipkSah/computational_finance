function u = solve_bvp(N)
    % Solve u''(x) = f(x) in (0,1) with given boundary conditions
    % Input: N - Number of sub-intervals
    % Output: u - Numerical solution at grid points

    % Step size
    h = 1 / (N);
    x = linspace(0, 1, N+1); % Grid points including boundaries

    % Exact boundary conditions
    u0 = (1/9) * (13 - 4 * cos(3));  % u(0)
    u1 = (1/9) * (13 + 4 * (-1 + cos(3)) - 4 * cos(3)); % u(1)

    % Compute the source term f(x) correctly based on the given exact solution
    f = 4 * cos(3 * x(2:end-1));  % Evaluate f at interior points

    % Construct the tridiagonal matrix A
    A = (1/h^2) * (diag(-2 * ones(N-1, 1)) + diag(ones(N-2, 1), 1) + diag(ones(N-2, 1), -1));
    disp("check if it is tridiagnol matrix")
    disp(A)

    % Construct the right-hand side vector b
    b = f(:); % Ensure column vector
    b(1) = b(1) - u0 / h^2; % Incorporate boundary conditions
    b(end) = b(end) - u1 / h^2;

    % Solve the system A * u_interior = b
    u_interior = A \ b;

    % Construct full solution vector including boundary values
    u = [u0; u_interior; u1];

    % Plot numerical vs exact solution
    figure;
    plot(x, u, 'ro-', 'DisplayName', 'Numerical Solution');
    hold on;
    plot(x, (1/9) * (13 + 4*x.*(-1 + cos(3)) - 4*cos(3*x)), 'b-', 'DisplayName', 'Exact Solution');
    xlabel('x');
    ylabel('u(x)');
    title(sprintf('Numerical vs Exact Solution for N = %d', N));
    legend;
    grid on;
end
