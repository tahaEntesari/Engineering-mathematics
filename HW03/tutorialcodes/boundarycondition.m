function [pl,ql,pr,qr] = boundarycondition(xl,ul,xr,ur,t)
%for a PDE in time and one space dimension.
pl = ul ;
ql = 0;
pr = ur - exp(-t) ;
qr = 0;
end