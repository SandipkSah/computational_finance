function [u, x] = solveObstaclePSOR(N, omega, tol, maxIter)
    % Parameters
    h = 2/N;               % Grid spacing
    x = linspace(-1, 1, N+1)';  % Grid points
    
    % Obstacle function g(x)
    g = @(x) 15/16 + (3/8)*x - (25/16)*x.^2;
    G = g(x);              % Obstacle values at grid points
    
    % Initialize solution (u ≥ g)
    u = max(0, G);         % Initial guess (u = max(0,g))
    W = u;                 % Alias for solution vector
    
    % Finite difference matrix A (tridiagonal)
    A = gallery('tridiag', N+1, -1, 2, -1) / h^2;
    A = A(2:N, 2:N);       % Remove boundary points (u(-1)=u(1)=0)
    
    % PSOR iteration
    for iter = 1:maxIter
        u_old = u;
        
        for i = 2:N
            % Gauss-Seidel update with projection
            sigma = A(i,1:i-1)*u(1:i-1) + A(i,i+1:N)*u_old(i+1:N);
            u(i) = max(G(i), u_old(i) + omega*( -A(i,i)*u_old(i) - sigma ) / A(i,i));
        end
        
        % Check convergence
        if norm(u - u_old, 'inf') < tol
            fprintf('Converged in %d iterations\n', iter);
            break;
        end
    end
    
    % Set boundary conditions
    u(1) = 0;   % u(-1) = 0
    u(end) = 0; % u(1) = 0
    
    % Plot results
    figure;
    plot(x, u, 'b-', 'LineWidth', 2); hold on;
    plot(x, G, 'r--', 'LineWidth', 2);
    xlabel('x'); ylabel('u(x)');
    legend('Solution u(x)', 'Obstacle g(x)');
    title('PSOR Solution of Obstacle Problem');
    grid on;
end


% Run with:
N = 100;       % Number of grid points
omega = 1.5;   % Relaxation parameter (1 < omega < 2)
tol = 1e-6;    % Convergence tolerance
maxIter = 1000;% Maximum iterations

[u, x] = solveObstaclePSOR(N, omega, tol, maxIter);