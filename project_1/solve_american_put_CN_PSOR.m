function [S, tau, U_eu, U_am] = solve_american_put_CN_PSOR(r, sigma, K, T, Smax, Ns, Nt, omega, tol, max_iter)

    dS = Smax / Ns;
    dt = T / Nt;
    S = linspace(0, Smax, Ns+1)';
    tau = linspace(0, T, Nt+1);
    
    % Terminal condition
    U = max(K - S, 0);
    
    % Coefficients
    i = (0:Ns)';
    alpha = 0.25 * dt * (sigma^2 * i.^2 - r * i);
    beta  = -0.5 * dt * (sigma^2 * i.^2 + r);
    gamma = 0.25 * dt * (sigma^2 * i.^2 + r * i);
    
    A = sparse(Ns+1, Ns+1);
    B = sparse(Ns+1, Ns+1);
    for j = 2:Ns
        A(j,j-1) = -alpha(j);
        A(j,j)   = 1 - beta(j);
        A(j,j+1) = -gamma(j);
        
        B(j,j-1) = alpha(j);
        B(j,j)   = 1 + beta(j);
        B(j,j+1) = gamma(j);
    end
    A(1,1) = 1; A(end,end) = 1;
    B(1,1) = 1; B(end,end) = 1;
    
    % European Put
    U_eu = zeros(Ns+1, Nt+1);
    U_eu(:,1) = U;
    
    for n = 1:Nt
        b = B * U;
        b(1) = K * exp(-r * (T - tau(n+1)));
        b(end) = 0;
        U = A \ b;
        U_eu(:,n+1) = U;
    end
    
    % American Put using PSOR
    U = max(K - S, 0);
    U_am = zeros(Ns+1, Nt+1);
    U_am(:,1) = U;
    
    for n = 1:Nt
        b = B * U;
        b(1) = K * exp(-r * (T - tau(n+1)));
        b(end) = 0;
    
        U_old = U;
        for iter = 1:max_iter
            U_new = U_old;
            for j = 2:Ns
                g = max(K - S(j), 0);
                rhs = b(j) - A(j,1:j-1)*U_new(1:j-1) - A(j,j+1:end)*U_old(j+1:end);
                U_temp = (1 - omega) * U_old(j) + omega / A(j,j) * rhs;
                U_new(j) = max(g, U_temp);
            end
            if norm(U_new - U_old, inf) < tol
                break
            end
            U_old = U_new;
        end
        U = U_new;
        U_am(:,n+1) = U;
    end
    end
    