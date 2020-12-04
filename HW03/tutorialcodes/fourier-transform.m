% using matlab function 'fourier'
syms x 
g = @(x) sinc(x);
Fg = fourier(g(x));

figure
fplot(Fg)
xlabel('\omega')
ylabel('G(j\omega)')
axis([-2*pi 2*pi -0.5 1.5])
%%
% using matlab function 'ifourier'
syms z 
f = @(z) sinc(z);
Ff = ifourier(f(z));

figure
fplot(Ff)
xlabel('t')
ylabel('f(t)')
axis([-2*pi 2*pi -0.5 0.5])
%%
% calculate fourier transform
syms xx w
h = @(xx) exp(-xx.^2);
fh = @(xx,w) h(xx).*exp(-2*j.*w.*xx);
Fh  = int(fh,xx,-inf,inf);
figure
fplot(Fh)
xlabel('\omega')
ylabel('H(j\omega)')
axis([-2 2 -0.5 2])
%%
% calculate inverse fourier transform
syms ww t
l = @(ww) exp(-ww.^2);
fl = @(ww,t) l(ww).*exp(2*j*t.*ww);
Fl  = int(fl,ww,-inf,inf);
figure
fplot(Fl)
xlabel('t')
ylabel('l(t)')