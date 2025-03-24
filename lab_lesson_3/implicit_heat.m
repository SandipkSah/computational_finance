function implicit_heat(Nx, Nt)
    % implicit_heat - Solves the 1D heat equation using the implicit method.
    %
    % Syntax: implicit_heat(Nx, Nt)
    %
    % Inputs:
    %   Nx - Number of spatial grid points.
    %   Nt - Number of time steps.
    %
    % Outputs:
    %   A 3D surface plot of the temperature distribution over space and time.

    % Define the length of the rod and the total time
    L = 1;  % Length of the rod
    T = 1;  % Total time

    % Calculate spatial and temporal step sizes
    dx = L / Nx;  % Spatial step size
    dt = T / Nt;  % Temporal step size
    r = dt / dx^2;  % Ratio for stability

    % Create spatial and temporal grids
    x = linspace(0, L, Nx+1);  % Spatial grid points
    t = linspace(0, T, Nt+1);  % Temporal grid points

    % Initialize the solution matrix U(x,t)
    U = zeros(Nx+1, Nt+1);  % Matrix to store temperature values

    %% Initial Condition
    U(:,1) = sin(pi*x);  % Initial temperature distribution (sinusoidal)

    %% Boundary Conditions
    U(1,:) = 0;     % Temperature at x=0 is always 0
    U(end,:) = 0;   % Temperature at x=L is always 0

    %% Matrix Setup for Implicit Method
    A = (1 + 2*r) * eye(Nx-1) - r * diag(ones(Nx-2,1), 1) - r * diag(ones(Nx-2,1), -1);
    
        %% Time-stepping Loop
        for n = 1:Nt
            U(2:Nx, n+1) = A \ U(2:Nx, n);  
        end
    
        %% Plot Results
        figure;
        surf(t, x, U);  
        xlabel('Time t');  
        ylabel('Position x');  
        zlabel('Temperature U');  
        title('Implicit Method Solution to the Heat Equation');
end
    
