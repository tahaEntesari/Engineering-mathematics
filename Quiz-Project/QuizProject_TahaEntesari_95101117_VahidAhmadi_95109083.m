%% 1 
clear
x=zeros(1,1000);
x(1)=1;
H1=System01(x);
aaa=max(H1);
nnn=find(H1==max(H1));
fprintf('the maximum output is %f and it occutrs for n=%d\n',aaa,nnn);
H2=System02(x);
n=1:1000;
figure(1);
plot(n,H1);
grid;
hold on
plot(n,H2);
title('answer of both systems to Delta');
xlabel('n');
ylabel('H');
legend('System01','System02');
fprintf('the plot shows no difference between the two systems.\nAnd by calculating H(system01)/H(system02) \nwe will see that all terms are 1.\n');
figure(2);
plot(n,H1./H2);
title('H1/H2');
%% 2
syms n z;
w=pi/50;
nmax=find(H1==max(H1));
%nmax=22;
%max(H1)=.7450;
%%a=1/nmax*log(sin(nmax*pi/50)/.7450);
a=0.012568;
h=@(n)sin(w*n).*exp(-a*n);
Hbymatlab=ztrans(h(n));
Hbyhand=@(z)1/2i*(z/(z-exp(w*1i-a))-z/(z-exp(-w*1j-a)));
b=[0,exp(w*1j-a)-exp(-w*1j-a),0];
a=[1,-(exp(w*1j-a)+exp(-w*1j-a)),exp(-2*a)];
a=a*2i;
Z=[1,1/z,1/z^2];
fprintf('the delta answer of both systems is:\n');
disp(h(n));
fprintf('the Z-transform is calculated:\n');
disp(Hbyhand);
fprintf('the z-transform calculated by matlab is:\n');
disp(Hbymatlab);
%clear Hbyhand;
%Hbyhand=@(z)sum(b.*Z)/sym(a.*Z);
[r, p k]=residuez(b,a);
r=roots(b);
zplane(r,p);

%% 4 & 5
format long;
n=0:999;
np=800:1000;
clear x w;
x=ones(20,1000);
for l=1:20
    w=l^4;
for k=0:999
    x(l,k+1)=sin(w*k);
end
end
y1=zeros(20,1000);
y2=zeros(20,1000);
for l=1:20
    figure(l+2);
    y1(l,:)=System01(x(l,:));
    y2(l,:)=System02(x(l,:));
    plot(n,y1(l,:));
    title(['output plot at frequency ' num2str(l^4)]);
    xlabel('n');
    ylabel('y');
    hold on;
    plot(n,y2(l,:));
    legend('System01','System02');
end

% difference plots
clear x y1 y2 n;
n=1:500;
x=ones(1,500);
y1to1=System01(x);
y2to1=System02(x);
x=-x;
y1tom1=System01(x);
y2tom1=System02(x);
figure(300);
plot(n,System01(0)-(y1to1+y1tom1));
axis([1,500,-.2,.2]);
title('linearity test for x1=1 and x2=-1 for System01');
figure(301);
plot(n,System02(0)-(y2to1+y2tom1));
title('linearity test for x1=1 and x2=-1 for System02');
fprintf('just by this test it is obvious that system02 is non linear.\nso further tests will only be run on System01');


for k=0:499
    x(k+1)=sin(w*k);
end
y1to1=System01(x);

y2to1=System02(x);
for k=0:499
    x(k+1)=log(1+k);
end
y1tom1=System01(x);
y2tom1=System02(x);

for k=0:499
    x(k+1)=log(1+k)+sin(100*k);
end
y1t=System01(x);
y2t=System02(x);
figure(302);
plot(n,y1to1+y1tom1-y1t);
title('System01 linearity test with inputs log(1+n) and sin(100n)');
figure(303);
plot(n,y2to1+y2tom1-y2t);
title('System02 linearity test with inputs log(1+n) and sin(100n)');
fprintf('for System01 though this output is not ideal but we can still count the system Linear\n');




%%  6 
%let's try u(n)
clear x y1 y2 n ;
n=1:1000;
x=ones(1,1000);
y1=System01(x);
y2=System02(x);
figure(23);
plot(n,y1);
hold on;
plot(n,y2);
title('answer of both systems to u(n)');
legend('System01','System02');
fprintf(' The outputs are the same\n');

% 2*u(n)
x=x*2;
y1=System01(x);
y2=System02(x);
figure(30);
plot(n,y1);
hold on;
plot(n,y2);
title('answer of both systems to 2*u(n)');
legend('System01','System02');
fprintf(' The outputs are the same\n');

fprintf('the systems are linear towards the above 2 inputs.but we can not judge just by a pair of inputs');



% now (-1)^n
for i=1:1000
    x(i)=(-1)^(i-1);
end
y1=System01(x);
y2=System02(x);
figure(24);
plot(n,y1);
hold on;
plot(n,y2);
title('answer of both systems to (-1)^nu(n)');
legend('System01','System02');
fprintf(' The outputs are different\n');
% now- n
x=-1*ones(1,1000);
y1=System01(x);
y2=System02(x);
figure(31);
plot(n,y1);
hold on;
plot(n,y2);
title('answer of both systems to -1u(n)');
legend('System01','System02');
fprintf(' The outputs are different\n');


fprintf('By comparing these 4 tests we can surely state that System02 is not linear.\nbecause as it can be seen in the figure,it''s answer to n and -n are the same.\nbut we can not yet generalize System01.');

% now n
x=0:999;
y1=System01(x);
y2=System02(x);
figure(25);
plot(n,y1);
hold on;
plot(n,y2);
title('answer of both systems to n');
legend('System01','System02');
fprintf(' The outputs are the same\n');


% now sigme(j) from j to n
x=cumsum(x);
y1=System01(x);
y2=System02(x);
figure(26);
plot(n,y1);
hold on;
plot(n,y2);
title('answer of both systems to sigma(j) from 0 to n');
legend('System01','System02');
fprintf(' The outputs are the same\n');

% now 1.001^n
for i=1:1000
    x(i)=1.001^i;
end
y1=System01(x);
y2=System02(x);
figure(27);
plot(n,y1);
hold on;
plot(n,y2);
title('answer of both systems to 1.001^n');
legend('System01','System02');
fprintf(' The outputs are the same\n');


% now nchoosek(n,2)
x(1)=0;
for i=2:1000
    x(i)=nchoosek(i,2);
end

y1=System01(x);
y2=System02(x);
figure(28);
plot(n,y1);
hold on;
plot(n,y2);
title('answer of both systems to (n,2)');
legend('System01','System02');
fprintf(' The outputs are the same\n');


fprintf('one distinc difference between the two systems is that System01(-x)=-System01(x) but System02(-x)=System02(x)\n');




