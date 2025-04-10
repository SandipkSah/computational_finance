function [S, tau, U_all] = solve_european_put_MOL(r, sigma, K, T, Smax, Ns, Nt)
    % Solve European put using Method of Lines and RK4
    % Returns stock grid S, time grid tau, and full solution matrix U_all
    
    dS = Smax / Ns;
    dt = T / Nt;
    S = linspace(0, Smax, Ns+1)';
    tau = linspace(0, T, Nt+1);
    
    % Initial condition
    U0 = max(K - S, 0);
    U = U0;
    
    % Boundary conditions
    bc_low = @(t) K * exp(-r * (T - t));  % S = 0
    bc_high = @(t) 0;                     % S = Smax
    
    % Derivative matrices
    e = ones(Ns+1,1);
    M = spdiags([-0.5*e 0*e 0.5*e], [-1 0 1], Ns+1, Ns+1);
    L = spdiags([e -2*e e], [-1 0 1], Ns+1, Ns+1);
    M(1,:) = 0; M(end,:) = 0;
    L(1,:) = 0; L(end,:) = 0;
    
    F = @(U, t) 0.5 * sigma^2 * (S.^2) .* (L * U) / dS^2 ...
              + r * S .* (M * U) / (2*dS) ...
              - r * U;
    
    U_all = zeros(Ns+1, Nt+1);
    U_all(:,1) = U0;
    
    for n = 1:Nt
        t = tau(n);
        U(1) = bc_low(t);
        U(end) = bc_high(t);
        
        k1 = dt*F(U, t);
        k2 = dt*F(U + 0.5*k1, t + 0.5*dt);
        k3 = dt*F(U + 0.5*k2, t + 0.5*dt);
        k4 = dt*F(U + dt*k3, t + dt);
        
        U = U + 1/6 * (k1 + 2*k2 + 2*k3 + k4);
        U(1) = bc_low(t + dt);
        U(end) = bc_high(t + dt);
        U_all(:,n+1) = U;
    end
    end
    