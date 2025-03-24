function [der] = forward(f,x0,h)
der=(f(x0+h)-f(x0))./h;
end