function implicit_heat(Nx, Nt)
    L = 1;
    T = 1;

    dx = L / Nx;
    dt = T / Nt;
    r = dt / dx^2;

    x = linspace(0, L, Nx+1);
    t = linspace(0, T, Nt+1);
    U = zeros(Nx+1, Nt+1);

    % Initial condition
    U(:,1) = sin(pi*x);

    % Boundary conditions
    U(1,:) = 0;
    U(end,:) = 0;

    % Matrix setup
    A = (1 + 2*r) * eye(Nx-1) - r * diag(ones(Nx-2,1), 1) - r * diag(ones(Nx-2,1), -1);

    for n = 1:Nt
        U(2:Nx, n+1) = A \ U(2:Nx, n);
    end

    % Plot results
    figure;
    surf(t, x, U);
    xlabel('Time t');
    ylabel('Position x');
    zlabel('Temperature U');
    title('Implicit Method Solution');
end
