function [c,b,s] = equation(x,t,u,DuDx)
%a PDE in time and one space dimension.
c = pi^2;
b = DuDx;
s = exp(-t)+exp(-2*t)*cos(3*pi*x/4);
end