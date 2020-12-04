%% question 1 part a

clear
syms T f n x
a=@(T ,f ,n, x)(2/T*int(f.*cos(2*pi*n*x/T),x,-T/2,T/2));
b=@(T ,f ,n, x)(2/T*int(f.*sin(2*pi*n*x/T),x,-T/2,T/2));
fouriers=@(f , x ,k ,T)(a(T,f,0,x)/2+symsum(a(T,f,n,x)*cos(2*pi*n*x/T)+b(T,f,n,x)*sin(2*pi*n*x/T),n,1,k));


f=@(x)(x*(sign(x)+1)/2);
ffs=fouriers(f,x,7,2);
disp(f(x))
pretty(f(x));
pretty(ffs)


g=@(x)(sin(pi*x/2));
gfs=eval(['@(x)',char(fouriers(g,x,7,2))]);
pretty(g(x));
pretty(gfs(x))



h=@(x)(x^2);
pretty(h(x));   
hfs=fouriers(h,x,7,2);
pretty(hfs);


l=@(x)(sin(x)*(sign(x)+1)/2);
lfs=fouriers(l,x,7,2*pi);
pretty(l(x));
pretty(lfs)

%% question 1 part b
figure(1)
cla;
fplot(f(x),[-1,1]);
hold on
fplot(fouriers(f,x,15,2),[-1,1]);
pause(1)
fplot(fouriers(f,x,50,2),[-1,1]);
pause(1)
fplot(fouriers(f,x,150,2),[-1,1]);
title('f(x) with its fourier series')
xlabel('x');
ylabel('f(x)');
legend('f(x)','15 term fourier','50 term fourier','150 term fourier');

figure(2)
cla;
fplot(g(x),[-1,1]);
hold on
fplot(fouriers(g,x,15,2),[-1,1]);
pause(1)
fplot(fouriers(g,x,50,2),[-1,1]);
pause(1)
fplot(fouriers(g,x,150,2),[-1,1]);
title('g(x) with its fourier series')
xlabel('x');
ylabel('g(x)');
legend('g(x)','15 term fourier','50 term fourier','150 term fourier');

figure(3)
cla;
fplot(h(x),[-1,1]);
hold on
fplot(fouriers(h,x,15,2),[-1,1]);
pause(1)
fplot(fouriers(h,x,50,2),[-1,1]);
pause(1)
fplot(fouriers(h,x,150,2),[-1,1]);
title('h(x) with its fourier series')
xlabel('x');
ylabel('h(x)');
legend('h(x)','15 term fourier','50 term fourier','150 term fourier');


figure(4)
cla;
fplot(l(x),[-1,1]);
hold on
fplot(fouriers(l,x,15,2*pi),[-pi,pi]);
pause(1)
fplot(fouriers(l,x,50,2*pi),[-pi,pi]);
pause(1)
fplot(fouriers(l,x,150,2*pi),[-pi,pi]);
title('l(x) with its fourier series')
xlabel('x');
ylabel('l(x)');
legend('l(x)','15 term fourier','50 term fourier','150 term fourier');

%% question 1 part c
figure(5)
cla;
title('absolute Error of f(x) with its fourier series')
xlabel('x');
ylabel('|fs(x)-f(x)|');
hold on
fplot(abs(fouriers(f,x,15,2)-f(x)),[-1,1]);
pause(1)
fplot(abs(fouriers(f,x,50,2)-f(x)),[-1,1]);
pause(1)
fplot(abs(fouriers(f,x,150,2)-f(x)),[-1,1]);
legend('15 term fourier Error','50 term fourier Error','150 term fourier Error');

figure(6)
cla;
title('absolute Error of g(x) with its fourier series')
xlabel('x');
ylabel('|gs(x)-g(x)|');
hold on
fplot(abs(fouriers(g,x,15,2)-g(x)),[-1,1]);
pause(1)
fplot(abs(fouriers(g,x,50,2)-g(x)),[-1,1]);
pause(1)
fplot(abs(fouriers(g,x,150,2)-g(x)),[-1,1]);
legend('15 term fourier Error','50 term fourier Error','150 term fourier Error');

figure(7)
cla;
title('absolute Error of h(x) with its fourier series')
xlabel('x');
ylabel('|hs(x)-h(x)|');
hold on
fplot(abs(fouriers(h,x,15,2)-h(x)),[-1,1]);
pause(1)
fplot(abs(fouriers(h,x,50,2)-h(x)),[-1,1]);
pause(1)
fplot(abs(fouriers(h,x,150,2)-h(x)),[-1,1]);
legend('15 term fourier Error','50 term fourier Error','150 term fourier Error');


figure(8)
cla;
title('absolute Error of l(x) with its fourier series')
xlabel('x');
ylabel('|ls(x)-l(x)|');
hold on
fplot(abs(fouriers(l,x,15,2*pi)-l(x)),[-pi,pi]);
pause(1)
fplot(abs(fouriers(l,x,50,2*pi)-l(x)),[-pi,pi]);
pause(1)
fplot(abs(fouriers(l,x,150,2*pi)-l(x)),[-pi,pi]);
legend('15 term fourier Error','50 term fourier Error','150 term fourier Error');

%% Question 2 part a

clear
syms T f n x
a=@(T ,f ,n, x)(2/T*int(f.*cos(2*pi*n*x/T),x,-T/2,T/2));
b=@(T ,f ,n, x)(2/T*int(f.*sin(2*pi*n*x/T),x,-T/2,T/2));
fouriers=@(f , x ,k ,T)(a(T,f,0,x)/2+symsum(a(T,f,n,x)*cos(2*pi*n*x/T)+b(T,f,n,x)*sin(2*pi*n*x/T),n,1,k));
E=@(T,x,f)(1/T*int(f(x)^2,x,-T/2,T/2));

f=@(x)(abs(x));
Ef=E(2,x,f);
pretty(f(x));
pretty(Ef);

g=@(x)(x);
Eg=E(1,x,g);
pretty(g(x));
pretty(Eg);

u=@(x)(sinh(x));
Eu=E(4,x,u);
pretty(u(x));
pretty(Eu);

%% Question 2 part b
tic
N=1:100;
Efs=ones(1,100)*(a(2,f,0,x)/2)^2;
for i=1:100
        for j=i:100
        Efs(j)=Efs(j)+1/2*((a(2,f,i,x))^2+(b(2,f,i,x))^2);
    end
end
figure(9);

plot(N,Efs);
title('f(x) Energy summation using fourier series')
xlabel('x');
ylabel('Energy');



Egs=ones(1,100)*(a(1,g,0,x)/2)^2;
for i=1:100
        for j=i:100
        Egs(j)=Egs(j)+1/2*((a(1,g,i,x))^2+(b(1,g,i,x))^2);
    end
end
figure(10);
plot(N,Egs);
title('g(x) Energy summation using fourier series')
xlabel('x');
ylabel('Energy');

Eus=ones(1,100)*(a(4,u,0,x)/2)^2;
for i=1:100
        for j=i:100
        Eus(j)=Eus(j)+1/2*((a(4,u,i,x))^2+(b(4,u,i,x))^2);
    end
end
figure(11);
plot(N,Eus);
title('u(x) Energy summation using fourier series')
xlabel('x');
ylabel('Energy');
toc

%% Question 2 part c


ferr=eval(abs(Ef-Efs));
ffirstN=find((ferr./Ef)<.1);
figure(12);
cla;
plot(N,ferr);
title('f(x) energy difference');
xlabel('x');
ylabel('|E(f)-E(fourier(f))|');

fprintf('the least N for which energy using fourier series has 0.1percent err for f(x): %d  \n',N(ffirstN(1)));


gerr=eval(abs(Eg-Egs));
gfirstN=find(gerr./Eg<0.1);
figure(13);
cla;
plot(N,gerr);
title('g(x) energy difference');
xlabel('x');
ylabel('|E(g)-E(fourier(g))|');

fprintf('the least N for which energy using fourier series has 0.1percent err for g(x): %d  \n',N(gfirstN(1)));



uerr=eval(abs(Eu-Eus));
ufirstN=find(uerr/Eu<.1);
figure(14);
cla;
plot(N,uerr);
title('u(x) energy difference');
xlabel('x');
ylabel('|E(u)-E(fourier(u))|');

fprintf('the least N for which energy using fourier series has 0.1percent err for u(x): %d\n',N(ufirstN(1)));


%% Question 3 part a

clear 
syms x T n
a=@(T ,f ,n, x)(2/T*int(f.*cos(2*pi*n*x/T),x,-T/2,T/2));
b=@(T ,f ,n, x)(2/T*int(f.*sin(2*pi*n*x/T),x,-T/2,T/2));
fouriers=@(f , x ,k ,T)(a(T,f,0,x)/2+symsum(a(T,f,n,x)*cos(2*pi*n*x/T)+b(T,f,n,x)*sin(2*pi*n*x/T),n,1,k));

p=@(x)(x*(sign(x)+1)/2);
pfs=fouriers(p,x,300,2);
figure(15);
cla;
fplot(p(x),[-1,1]);
hold on
fplot(pfs,[-1,1]);
fplot(pfs-p(x),[-1,1]);
legend('p(x)','fourier of p','F(p(x))-p(x)');
title('overshoot calculation for p(x)');
xlabel('x');
x=(-1-eps):.00001:(1+eps);
pfsevaluatednumerically=(eval(pfs));
maxima=(max(pfsevaluatednumerically));
minima=(min(pfsevaluatednumerically));
fprintf('overshoot of p(x) fourier calculation at the end of the boundry: %f\n',maxima-1);
fprintf('overshoot of p(x) fourier calculation at the beginning of the boundry: %f\n',minima);



clear x
syms x

q=@(x)(sign(-x+.5)-sign(-x-.5));
qfs=fouriers(q,x,300,2);
figure(16);
cla;
fplot(q(x),[-1,1]);
hold on
fplot(qfs,[-1,1]);
fplot(qfs-q(x),[-1,1]);
legend('q(x)','fourier of q','F(q(x))-q(x)');
title('overshoot calculation for q(x)');
xlabel('x');
x=(-1-eps):.00001:(1+eps);
qfsevaluatednumerically=(eval(qfs));
maxima=(max(qfsevaluatednumerically));
fprintf('overshoot of q(x) fourier calculation: %f\n',maxima-2);

%% Question 3 part b
clear x
syms x


loncozfouriers=@(f , x ,k ,T)(a(T,f,0,x)/2+symsum((a(T,f,n,x)*cos(2*pi*n*x/T)+b(T,f,n,x)*sin(2*pi*n*x/T))*sinc(n*pi/2/k),n,1,k));
figure(17);
cla;
hold on
fplot(loncozfouriers(p,x,20,2),[-1,1]);
fplot(fouriers(p,x,20,2),[-1,1]);
pause(1);
fplot(loncozfouriers(p,x,50,2),[-1,1]);
fplot(fouriers(p,x,50,2),[-1,1]);
pause(1);
fplot(loncozfouriers(p,x,400,2),[-1,1]);
fplot(fouriers(p,x,400,2),[-1,1]);
title('lonczos fourier series coparrison with normal fourier series for p(x)')
xlabel('x');
legend('20 term lonczos','20 term fourier','50 term lonczos','50 term fourier','400 term lonczos','400 term fourier','location','northwest');


figure(18);
cla;
hold on
fplot(loncozfouriers(q,x,20,2),[-1,1]);
fplot(fouriers(q,x,20,2),[-1,1]);
pause(1);
fplot(loncozfouriers(q,x,50,2),[-1,1]);
fplot(fouriers(q,x,50,2),[-1,1]);
pause(1);
fplot(loncozfouriers(q,x,400,2),[-1,1]);
fplot(fouriers(q,x,400,2),[-1,1]);
title('lonczos fourier series coparrison with normal fourier series for q(x)')
xlabel('x');
legend('20 term lonczos','20 term fourier','50 term lonczos','50 term fourier','400 term lonczos','400 term fourier','location','northwest');




%% Question 4 part a

clear
load x;
load y;
a=zeros(1,5);
b=zeros(1,5);
for i=1:5
a(i)=2/2001*sum(y.*cos(pi*i*x));
b(i)=2/2001*sum(y.*sin(pi*i*x));
end
fprintf('The coefficients are:\n \t\ta \t\t b\n');
disp([a' b'])
% Is is clear that the most significant coefficients are a2,b1 and b3.
%thus we estimate the main function using these coefficients
f=a(2)*cos(2*pi*x)+b(1)*sin(pi*x)+b(3)*sin(3*pi*x);


%% part b
F=fit(x',y','fourier5')
%% part c
figure(19)
cla;
plot(F,x,y);
hold on
plot(x,f,'r.')
legend('The initial data','using fit function','using only major fourier coefficient','location','northwest');
