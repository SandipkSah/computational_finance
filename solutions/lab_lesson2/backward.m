function [der] = backward(f,x0,h)
der=(f(x0)-f(x0-h))./h;
end