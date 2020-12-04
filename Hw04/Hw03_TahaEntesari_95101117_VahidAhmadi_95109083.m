%% Qustion 1
clear
syms z;

f=@(z)(sinh(z)-sqrt(2*i));
fans=solve(f(z)==0)
figure(1);
plot(real(fans),imag(fans),'*')
title('f(z) roots');
xlabel('real f(z)=0');
ylabel('imag f(z)=0');

g=@(z)(sin(z^6)-2*i);
gans=solve(g(z)==0)
figure(2);
plot(real(gans),imag(gans),'*')
title('g(z) roots');
xlabel('real g(z)=0');
ylabel('imag g(z)=0');

h=@(z)(z^5+2*z^3+1);
hans=solve(h(z)==0)
figure(3);
plot(real(hans),imag(hans),'*')
title('h(z) roots');
xlabel('real h(z)=0');
ylabel('imag h(z)=0');

l=@(z)((sin(z^2)-i)^3);
lans=solve(l(z)==0)
figure(4);
plot(real(lans),imag(lans),'*')
title('l(z) roots');
xlabel('real l(z)=0');
ylabel('imag l(z)=0');

w=@(z)(z^6-log(j));
wans=solve(w(z)==0)
figure(5);
plot(real(wans),imag(wans),'*')
title('w(z) roots');
xlabel('real w(z)=0');
ylabel('imag w(z)=0');

v=@(z)(cos(z^(1/4))-j^j+2j);
vans=solve(v(z)==0)
figure(6);
plot(real(vans),imag(vans),'*')
title('v(z) roots');
xlabel('real v(z)=0');
ylabel('imag v(z)=0');
%% Question 2
clear
syms x y
assume(x,'real');
assume(y,'real');
z=x+1i*y;

f=@(z)(sinc(z));
realf=real(f(z));
imagf=imag(f(z));
pretty(realf);
pretty(imagf);


g=@(z)(z^3+3*z^2);
realg=real(g(z));
imagg=imag(g(z));
pretty(realg);
pretty(imagg);

h=@(z)(cos(1/z))^2;
realh=real(h(z));
imagh=imag(h(z));
pretty(realh);
pretty(imagh);


l=@(z)(sinh(z^2));
reall=real(l(z));
imagl=imag(l(z));
pretty(reall);
pretty(imagl);


w=@(z)(imag(z^3));
realw=real(w(z));
imagw=imag(w(z));
pretty(realw);
pretty(imagw);


v=@(z)(1/z);
realv=real(v(z));
imagv=imag(v(z));
pretty(realv);
pretty(imagv);


%% part b
fulaplacian=eval(['@(x,y)' char(diff(realf,x,2) +diff(realf,y,2))]);
%
%fulaplacian=eval(['@(x)' char(diff(realf,x,2))])+ eval(['@(y)' char(diff(realf,y,2))]);
%fvlaplacian=eval(['@(x)' char(diff(imagf,x,2))])+eval(['@(y)' char(diff(imagf,y,2))]);
fvlaplacian=eval(['@(x,y)' char(diff(imagf,x,2) + diff(imagf,y,2))]);
ftest=isequal(fulaplacian,fvlaplacian);
if ftest ,fprintf('f(z) is analytic and it''s real and imaginary parts are harmonic conjugates\n');
else fprintf('f(z) is not analytic\n');
end
%%
gulaplacian=diff(realg,x,2)+diff(realg,y,2);
gvlaplacian=diff(imagg,x,2)+diff(imagg,y,2);
gtest=isequal(gulaplacian,gvlaplacian);
if gtest ,fprintf('g(z) is analytic and it''s real and imaginary parts are harmonic conjugates\n');
else fprintf('g(z) is not analytic\n');
end



hulaplacian=diff(realh,x,2)+diff(realh,y,2);
hvlaplacian=diff(imagh,x,2)+diff(imagh,y,2);
htest=isequal(hulaplacian,hvlaplacian);
if htest ,fprintf('h(z) is analytic and it''s real and imaginary parts are harmonic conjugates\n');
else fprintf('h(z) is not analytic\n');
end


lulaplacian=diff(reall,x,2)+diff(reall,y,2);
lvlaplacian=diff(imagl,x,2)+diff(imagl,y,2);
ltest=isequal(lulaplacian,lvlaplacian);
if ltest ,fprintf('l(z) is analytic and it''s real and imaginary parts are harmonic conjugates\n');
else fprintf('l(z) is not analytic\n');
end


wulaplacian=diff(realw,x,2)+diff(realw,y,2);
wvlaplacian=diff(imagw,x,2)+diff(imagw,y,2);
wtest=isequal(wulaplacian,wvlaplacian);
if wtest ,fprintf('w(z) is analytic and it''s real and imaginary parts are harmonic conjugates\n');
else fprintf('w(z) is not analytic\n');
end


vulaplacian=diff(realv,x,2)+diff(realv,y,2);
vvlaplacian=diff(imagv,x,2)+diff(imagv,y,2);
vtest=isequal(vulaplacian,vvlaplacian);
if vtest ,fprintf('v(z) is analytic and it''s real and imaginary parts are harmonic conjugates\n');
else fprintf('v(z) is not analytic\n');
end
    
%% Question 3
clear
syms x y
assume(x,'real');
assume(y,'real');
z=x+1i*y;

f=@(z)(z^4+z^2);

g=@(z)(imag(z)+2*z);

h=@(z)(1/(1+z^2));

l=@(z)(cosh(z^3+3*z));

w=@(z)(exp(z^2+z));

v=@(z)(conj(z)+2*log(z));

fu=real(f(z));
fv=imag(f(z));
f1stCR=isequal(diff(fu,x),diff(fv,y));
f2ndCR=isequal(diff(fu,y),-diff(fv,x));
if(f1stCR&&f2ndCR)
    fprintf('f(z) is analytic\nthe derivative is:\n');
    disp(diff(f(z),x));
else fprintf('f(z) is not analytic\n');
end


gu=real(g(z));
gv=imag(g(z));
g1stCR=isequal(diff(gu,x),diff(gv,y));
g2ndCR=isequal(diff(gu,y),-diff(gv,x));
if(g1stCR&&g2ndCR)
    fprintf('g(z) is analytic\nthe derivative is:\n');
    disp(diff(g(z),x));
else fprintf('g(z) is not analytic\n');
end

hu=real(h(z));
hv=imag(h(z));
h1stCR=isequal(diff(hu,x),diff(hv,y));
h2ndCR=isequal(diff(hu,y),-diff(hv,x));
if(h1stCR&&h2ndCR)
    fprintf('h(z) is analytic\nthe derivative is:\n');
    disp(diff(h(z),x));
else fprintf('h(z) is not analytic\n');
end

lu=real(l(z));
lv=imag(l(z));
l1stCR=isequal(diff(lu,x),diff(lv,y));
l2ndCR=isequal(diff(lu,y),-diff(lv,x));
if(l1stCR&&l2ndCR)
    fprintf('l(z) is analytic\nthe derivative is:\n');
    disp(diff(l(z),x));
else fprintf('l(z) is not analytic\n');
end

wu=real(w(z));
wv=imag(w(z));
w1stCR=isequal(diff(wu,x),diff(wv,y));
w2ndCR=isequal(diff(wu,y),-diff(wv,x));
if(w1stCR&&w2ndCR)
    fprintf('w(z) is analytic\nthe derivative is:\n');
    disp(diff(w(z),x));
else fprintf('w(z) is not analytic\n');
end

vu=real(v(z));
vv=imag(v(z));
v1stCR=isequal(diff(vu,x),diff(vv,y));
v2ndCR=isequal(diff(vu,y),-diff(vv,x));
if(v1stCR&&v2ndCR)
    fprintf('v(z) is analytic\nthe derivative is:\n');
    disp(diff(v(z),x));
else fprintf('v(z) is not analytic\n');
end


%% Question 4
clear
[x y] = meshgrid(-1:0.002:1,-1:0.002:1);
 
au=(x+exp(y).*cos(x));
bu=(exp(-2*x.*y).*sin(x.^2-y.^2));
cu=(x.*exp(x).*cos(y)-y.*exp(x).*sin(y));
du=(exp(-x).*(x.*sin(y)-y.*cos(y)));

av=y-exp(y).*sin(x);
bv=-exp(-2*x.*y).*cos(x.^2-y.^2);
cv=exp(x).*(y.*cos(y)+x.*sin(y));
dv=exp(-x).*(x.*cos(y)+y.*sin(y));


figure(7)
contour(x,y,au,'g')
hold on
contour(x,y,av,'r')
legend('real((potential))','imaginary((field))');
title('A');
xlabel('x');
ylabel('y');

figure(8)
contour(x,y,bu,'g')
hold on
contour(x,y,bv,'r')
legend('real((potential))','imaginary((field))');
title('B');
xlabel('x');
ylabel('y');

figure(9)
contour(x,y,cu,'g')
hold on
contour(x,y,cv,'r')
legend('real((potential))','imaginary((field))');
title('C');
xlabel('x');
ylabel('y');

figure(10)
contour(x,y,du,'g')
hold on
contour(x,y,dv,'r')
legend('real((potential))','imaginary((field))');
title('D');
xlabel('x');
ylabel('y');


