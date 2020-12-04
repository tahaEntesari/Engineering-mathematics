
%% write Greek letters in plot 
theta = linspace(0,2*pi,100);
alpha = theta.^2 + theta ;
plot(theta,alpha,'linewidth',2)
xlabel('\theta')
ylabel('\alpha (\theta)') 

%% Define a symbolic function 
syms x t f g

f = @(x) (x.^2 + 2*x ) ;
pretty(f(x));

g = @(x,t) x.*exp(t)./t;
pretty(g(x,t));

%% calculate Integral 
Intf = int(f,x);
pretty(Intf);

Intg1 = int(g,x);
pretty(Intg1);

Intg2 = int(g,t);
pretty(Intg2);

%% calculate derivative
DIFFf = diff(f,x);
pretty(DIFFf);

DIFFg1 = diff(g,x);
pretty(DIFFg1);

DIFFg2 = diff(g,t);
pretty(DIFFg2);

%% calculate Definite Integral 
int (f,x,-1,1) 

int (g,x,-1,1) 

int (g,t,1,2)   

%% plot symbolic function 
fplot(f,[-1 1])
title('f(x)')
figure
subplot(2,2,1)
fplot(g(x,3),[-1 1])
title('g(x,3)')
subplot(2,2,2)
fplot(g(1,t),[-1 1])
title('g(1,t)')
subplot(2,2,3)
fsurf(g)
title('fsurf(g)')
subplot(2,2,4)
fcontour(g)
title('fcontour(g)')


