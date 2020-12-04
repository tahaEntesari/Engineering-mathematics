%% Question 1
syms x y;
f=@(x,y)(sin(x)*y^2);
u=@(x,y)(y*log(x)+x*sin(y));
g=@(x,y)(cos(y)/(1+x^2));
v=@(x,y)(besselj(0,x*y));

ft3=taylor(f,[x,y],[1,0],'order',3);
ft6=taylor(f,[x,y],[1,0],'order',6);
ft8=taylor(f,[x,y],[1,0],'order',8);
fprintf('f(x,y) taylor series expantion at [1,0]:\n');
fprintf('third order taylor:\n');
disp(ft3);
fprintf('sixth order taylor:\n');
disp(ft6);
fprintf('8''th order taylor:\n');
disp(ft8);

ut3=taylor(u,[x,y],[1,0],'order',3);
ut6=taylor(u,[x,y],[1,0],'order',6);
ut8=taylor(u,[x,y],[1,0],'order',8);
fprintf('u(x,y) taylor series expantion at [1,0]:\n');
fprintf('third order taylor:\n');
disp(ut3);
fprintf('sixth order taylor:\n');
disp(ut6);
fprintf('8''th order taylor:\n');
disp(ut8);

gt3=taylor(g,[x,y],[0,1],'order',3);
gt6=taylor(g,[x,y],[0,1],'order',8);
gt8=taylor(g,[x,y],[0,1],'order',8);
fprintf('g(x,y) taylor series expantion at [1,0]:\n');
fprintf('third order taylor:\n');
disp(gt3);
fprintf('sixth order taylor:\n');
disp(gt6);
fprintf('8''th order taylor:\n');
disp(gt8);

vt3=taylor(v,[x,y],[0,0],'order',3);
vt6=taylor(v,[x,y],[0,0],'order',6);
vt8=taylor(v,[x,y],[0,0],'order',8);
fprintf('v(x,y) taylor series expantion at [1,0]:\n');
fprintf('third order taylor:\n');
disp(vt3);
fprintf('sixth order taylor:\n');
disp(vt6);
fprintf('8''th order taylor:\n');
disp(vt8);
%% part b
figure(1);
fsurf(f,'r');
axis ([-3 3 -3 3 -25 25])
hold on
title('f(x,y) plot');
xlabel('x');
ylabel('y');
zlabel('f(x,y)');
fsurf(ft3,'b');
fsurf(ft6,'g');
fsurf(ft8,'y');
legend('main function ','3''rd order taylor','6''th order taylor','8''th order taylor')

figure(2);
fsurf(u,'r');
axis ([-3 3 -3 3 -25 25])
hold on
title('u(x,y) plot');
xlabel('x');
ylabel('y');
zlabel('u(x,y)');
fsurf(ut3,'b');
fsurf(ut6,'g');
fsurf(ut8,'y');
legend('main function ','3''rd order taylor','6''th order taylor','8''th order taylor')

figure(3);
fsurf(g,'r');
axis ([-3 3 -3 3 -25 25])
hold on
title('g(x,y) plot');
xlabel('x');
ylabel('y');
zlabel('g(x,y)');
fsurf(gt3,'b');
fsurf(gt6,'g');
fsurf(gt8,'y');
legend('main function','3''rd order taylor','6''th order taylor','8''th order taylor')

figure(4);
fsurf(v,'r');
axis ([-3 3 -3 3 -25 25])
hold on
title('v(x,y) plot');
xlabel('x');
ylabel('y');
zlabel('v(x,y)');
fsurf(vt3,'b');
fsurf(vt6,'g');
fsurf(vt8,'y');
legend('main function ','3''rd order taylor','6''th order taylor','8''th order taylor')

%% Question 2

clear
format long
syms x y t;
assume (x,'real')
assume(y,'real')
assume(t,'real')
z=x+1i*y;
f=@(z)(3*z.^4+z)./((z.^4).*(z+1));
fprintf('F(z) is:\n');
pretty(f(z));
fzt1=@(t) 2+2*exp(1i*t);
fzt2=@(t) exp(1i*t);
fzt1p=@(t) 2i*exp(1i*t);
fzt2p=@(t) 1i*exp(1i*t);
foverc=integral(@(t)f(fzt1(t)).*fzt1p(t),0,2*pi);
fovercp=integral(@(t)f(fzt2(t)).*fzt2p(t),0,2*pi);
fprintf('The integral of f over C is:\n%.16f\t+%.16fi\n Over C'' is:\n%.16f\t+%.16fi\n',real(foverc),imag(foverc),real(fovercp),imag(fovercp))

%%
g=@(z)(z.^3.*exp(-(z.^2+.64)))./(z.^6+z.^2);
fprintf('G(z) is:\n');
pretty(g(z))
% first let's calculate over C
gzt1=@(t).5*(1+1i+exp(1i*t));
gzt2=@(t)(t+1i);
gzt3=@(t).5*(-1+1i+exp(1i*t));
gzt4=@(t)t;
gzt1p=@(t)(1i/2*exp(1i*t));

goverc=integral(@(t)g(gzt1(t)).*gzt1p(t),-pi/2,pi/2)+integral(@(t)g(gzt3(t)).*gzt1p(t),pi/2,3*pi/2)+integral(@(t)g(gzt2(t)),.5,-.5)+integral(@(t)g(gzt4(t)),-.5,.5);
% now c'
gcpzt1=@(t).5*(1-1i+exp(1i*t));
gcpzt2=@(t)(t);
gcpzt3=@(t).5*(-1-1i+exp(1i*t));
gcpzt4=@(t)t-1i;
gcpzt1p=@(t)(1i/2*exp(1i*t));

govercp=integral(@(t)g(gcpzt1(t)).*gcpzt1p(t),-pi/2,pi/2)+integral(@(t)g(gcpzt3(t)).*gcpzt1p(t),pi/2,3*pi/2)+integral(@(t)g(gcpzt2(t)),.5,-.5)+integral(@(t)g(gcpzt4(t)),-.5,.5);

fprintf('The integral of g over C is:\n%.16f\t+%.16fi\n Over C'' is:\n%.16f\t+%.16fi\n',real(goverc),imag(goverc),real(govercp),imag(govercp))

assume (x,'clear')
assume (y,'clear')
assume (t,'clear')
%% Question 3

clear
format long

syms x y t;
assume(x,'real');
assume(y,'real');
assume(t,'real');
z=x+1i*y;
f=@(z)(z.^2+5).*(z+exp(-z))./((exp(z+.5)).*(z-.5i));
fprintf('f(z) is :\n');
pretty(f(z));
g=@(z)(z.^2+5).*(z+exp(-z))./((exp(z+.5)));
fz=@(t).5i+exp(1i*t);
fzp=@(t)1i*exp(1i*t);   
I=integral(@(t)f(fz(t)).*fzp(t),0,2*pi);
fprintf('the integral by parameterizing or cauchy''s integral formula is:\n%.16f\t+%.16fi\n',real(I),imag(I));
fprintf('we take g(z) to be f(z)*(z-0.5j),that is g(z)=:\n')
pretty(g(z));
fprintf('by cauchy''s integral formula we know that g(0.5j)*2pi*i\n is the integral that we want to calculate.\nthus by dividing the calculated integral by (g(0.5j)*2i) we will get pi\n');
dt=0.00001;
w=0:dt:2*pi;
ant=f(fz(w)).*fzp(w);
res=sum(ant)*dt;
caledpi=res/(g(.5j)*2i);
fprintf('the calculated integral by taking delta to be %f is :\n%.16f\t+%.16fi\n',dt,real(res),imag(res));
fprintf('the calculated pi is:%.16f\n',caledpi);
fprintf('the error is :%.16f\n',abs(pi-caledpi));
fprintf('the delta is %f\n so for intervals smaller that this we will get a better result\n',dt);


assume (x,'clear')
assume (y,'clear')
assume (t,'clear')
%% Question 4
clear
format long
syms x y t;
assume(x,'real');
assume(y,'real');
assume(t,'real');
z=x+1i*y;
f=@(z)exp(z+1./z);
fz1=@(t)exp(1i*t);
fz1p=@(t)1i*exp(1i*t);
I=integral(@(t)f(fz1(t)).*fz1p(t),0,2*pi);

% by testing I found that n=50 is an appropriate form of infinity.
n=1:50;
an=1./(factorial(n).*factorial(n-1));
sn=cumsum(an);
% In order to show that a series is convergant an essential condition is
% that the an's approach zero as n approaches infinity
fprintf('The an''s approach zero and thus the essential condition that an is convergent holds.\nan for n=50 is as follows:\n');
disp(an(50));
format short
fprintf('To show that the series is convergent we must show that a(n+1)/a(n) is smaller than unity:\n');
L=an(2:50)./an(1:49);
disp(L)
fprintf('by the calculation done above it is clear that the series is convergent\n');
fprintf('limit of the series is as follows:\n');
disp(sn(50));   
fprintf('I is %.16f \t 2*pi*j*l is %.16f \n the error of the two mentioned is %.16f\n',imag(I),imag(2j*pi*sn(50)),abs(imag(I-2j*pi*sn(50))));

assume (x,'clear')
assume (y,'clear')
assume (t,'clear')
