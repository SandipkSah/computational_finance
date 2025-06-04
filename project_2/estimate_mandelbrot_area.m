function estimate_mandelbrot_area()
    rng(12345);  % Reproducibility

    %% Create output directory if it doesn't exist
    if ~exist('output', 'dir')
        mkdir('output');
    end

    %% Load Mandelbrot characteristic mask
    load('mandelbrot.mat');  % Expects 'mandelbrot'
    M = mandelbrot;
    [ny, nx] = size(M);

    %% Monte Carlo Method
    N = 100000;
    Ux = rand(N, 1);
    Uy = rand(N, 1);
    monte_idx = mask_indices(Ux, Uy, M);
    area_mc = sum(monte_idx) / N;

    %% Quasi-Monte Carlo Method (Halton)
    halton_pts = halton2d(N, [2, 3]);
    Hx = halton_pts(:,1);
    Hy = halton_pts(:,2);
    halton_idx = mask_indices(Hx, Hy, M);
    area_qmc = sum(halton_idx) / N;

    %% Output
    fprintf('Estimated Area (Monte Carlo): %.5f\n', area_mc);
    fprintf('Estimated Area (Quasi-Monte Carlo): %.5f\n', area_qmc);

    %% Plot and save
    figure;
    subplot(1,2,1);
    scatter(Ux(monte_idx), Uy(monte_idx), 1, 'b'); title('Monte Carlo Points Inside');
    axis equal; xlim([0 1]); ylim([0 1]);

    subplot(1,2,2);
    scatter(Hx(halton_idx), Hy(halton_idx), 1, 'r'); title('Halton Points Inside');
    axis equal; xlim([0 1]); ylim([0 1]);

    % Save the figure
    saveas(gcf, fullfile('output', 'mandelbrot_area_estimation.png'));
end


function idx = mask_indices(x, y, M)
    [ny, nx] = size(M);
    xi = max(min(floor(x * nx) + 1, nx), 1);
    yi = max(min(floor(y * ny) + 1, ny), 1);
    idx = M(sub2ind([ny, nx], yi, xi)) > 0;
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
