 syms x k n l

evalin(symengine,'assume(k,Type::Integer)');

a = @(f,x,k,l) int(f.*cos(k*pi*x/l)/l,x,-l,l);

b = @(f,x,k,l) int(f.*sin(k*pi*x/l)/l,x,-l,l);

fs = @(f,x,n,l) a(f,x,0,l)/2 + symsum(a(f,x,k,l)*cos(k*pi*x/l) + b(f,x,k,l)*sin(k*pi*x/l),k,1,n);

f =@(x) sign(x);

pretty(fs(f,x,3,1));

ezplot(f,-1,1)
hold on
ezplot(fs(f,x,5,1),-1,1)
pause(1.5)
ezplot(fs(f,x,15,1),-1,1)
pause(1.5)
ezplot(fs(f,x,25,1),-1,1)