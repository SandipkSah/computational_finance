%% Helper Functions for Finite Differences
function W_T = func_F(T, Y)
    % W = 1 + T;
    
end

function y_j_1=next_step(y_j, h, t_j)
    f1 = func_F(t_j, y_j)
    f2 = func_F(t_j + h, y_j + f1)
    y_j_1 = y_j + h/2 * (f1 + f2)
end 


function y_t=RK_Heun()

    y_t="return value"
end 



% Parameters and exact solution
r = 0.06;
sigma = 0.3;
T = 1;
K = 10;
S_* = 15;
