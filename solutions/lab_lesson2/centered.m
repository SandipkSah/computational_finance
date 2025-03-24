function [der] = centered(f,x0,h)
der=(f(x0+h)-f(x0-h))./(2*h);
end