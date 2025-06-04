function random_generators()
    clc;
    rng(12345);  % Set seed for reproducibility
    N = 10000;

    %% Create output folder if not present
    if ~exist('output', 'dir')
        mkdir('output');
    end

    %% Uniform U([a,b])
    a = 2; b = 5;
    U_ab = (b - a) * rand(1, N) + a;
    figure;
    histogram(U_ab, 50, 'Normalization', 'pdf');
    title('U([a,b]) Distribution'); xlabel('x'); ylabel('Density');
    saveas(gcf, fullfile('output', 'uniform_ab_distribution.png'));

    %% Exponential Exp(theta)
    theta = 2;
    U = rand(1, N);
    Exp_theta = -log(1 - U) / theta;
    figure;
    histogram(Exp_theta, 50, 'Normalization', 'pdf');
    title('Exponential(\theta) Distribution'); xlabel('x'); ylabel('Density');
    saveas(gcf, fullfile('output', 'exponential_theta_distribution.png'));

    %% Normal N(0,1) using Box-Muller
    U1 = rand(1, N/2);
    U2 = rand(1, N/2);
    R = sqrt(-2 * log(U1));
    Theta = 2 * pi * U2;
    Z1 = R .* cos(Theta);
    Z2 = R .* sin(Theta);
    Normal_samples = [Z1, Z2];
    figure;
    histogram(Normal_samples, 50, 'Normalization', 'pdf');
    title('Standard Normal N(0,1) Distribution'); xlabel('x'); ylabel('Density');
    saveas(gcf, fullfile('output', 'normal_box_muller_distribution.png'));

    %% Halton nodes in [0,1]x[0,1]
    halton_pts = halton2d(1000, [2, 3]);
    figure;
    scatter(halton_pts(:,1), halton_pts(:,2), '.');
    title('Halton Nodes in [0,1]x[0,1]'); xlabel('x'); ylabel('y');
    axis equal;
    saveas(gcf, fullfile('output', 'halton_nodes_scatter.png'));
end

function H = halton2d(N, bases)
    H = zeros(N, 2);
    for i = 1:N
        H(i, 1) = radical_inverse(i, bases(1));
        H(i, 2) = radical_inverse(i, bases(2));
    end
end

function r = radical_inverse(n, base)
    r = 0; f = 1 / base;
    while n > 0
        digit = mod(n, base);
        r = r + digit * f;
        n = floor(n / base);
        f = f / base;
    end
end
