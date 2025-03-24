function explicit_heat(Nx, Nt)
    % EXPLICIT_HEAT Solves the 1D heat equation using the explicit finite difference method.
    %
    %   explicit_heat(Nx, Nt) solves the heat equation on a spatial domain [0,1] and
    %   a time domain [0,1] using Nx spatial steps and Nt time steps.
    %
    %   The function enforces stability by ensuring that r = dt / dx^2 <= 0.5.
    %   The initial condition is U(x,0) = sin(pi*x), and Dirichlet boundary
    %   conditions U(0,t) = 0 and U(1,t) = 0 are applied.
    %
    %   The computed solution is visualized as a surface plot.

    % Define domain limits
    L = 1;  % Length of the spatial domain
    T = 1;  % Total simulation time

    % Compute step sizes
    dx = L / Nx;       % Spatial step size
    dt = T / Nt;       % Time step size
    r = dt / dx^2;     % Stability parameter

    % Check stability condition for explicit method
    if r > 0.5
        error('Stability condition violated: Choose smaller dt or larger dx.');
    end

    % Discretize space and time
    x = linspace(0, L, Nx+1);  % Spatial grid
    t = linspace(0, T, Nt+1);  % Time grid
    U = zeros(Nx+1, Nt+1);     % Solution matrix initialization

    % Set initial condition: U(x,0) = sin(pi*x)
    U(:,1) = sin(pi*x);

    % Apply boundary conditions: U(0,t) = 0 and U(1,t) = 0
    U(1,:) = 0;
    U(end,:) = 0;

    % Time-stepping loop using explicit finite difference method
    for n = 1:Nt
        for j = 2:Nx
            U(j, n+1) = U(j, n) + r * (U(j+1, n) - 2*U(j, n) + U(j-1, n));
        end
    end

    % Plot the computed solution as a surface plot
    figure;
    surf(t, x, U);
    xlabel('Time t');
    ylabel('Position x');
    zlabel('Temperature U');
    title('Explicit Finite Difference Method for 1D Heat Equation');
    shading interp; % Smooth shading for better visualization
end