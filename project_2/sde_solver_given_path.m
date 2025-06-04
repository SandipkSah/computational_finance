function [t, X_em, X_milstein] = sde_solver_given_path(a, b, dbdx, x0, T, N, dW)
    h = T / N;
    t = linspace(0, T, N+1);
    X_em = zeros(1, N+1);
    X_milstein = zeros(1, N+1);
    X_em(1) = x0;
    X_milstein(1) = x0;

    for n = 1:N
        tn = t(n);
        Xn_em = X_em(n);
        Xn_mil = X_milstein(n);
        dWn = dW(n);

        X_em(n+1) = Xn_em + a(tn, Xn_em)*h + b(tn, Xn_em)*dWn;
        X_milstein(n+1) = Xn_mil + a(tn, Xn_mil)*h + b(tn, Xn_mil)*dWn + ...
                          0.5 * b(tn, Xn_mil) * dbdx(tn, Xn_mil) * (dWn^2 - h);
    end
end
