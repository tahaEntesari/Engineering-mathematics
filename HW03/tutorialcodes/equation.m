function [c,b,s] = equation(x,t,u,DuDx)
%a PDE in time and one space dimension.
c = 1;
b = DuDx;
s = 0;
end