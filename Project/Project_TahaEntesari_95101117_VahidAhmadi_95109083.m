%% Question 1 part a
clear
syms a n;
% Since matlab itself calculates the one-sided Z-transform, there is no need
% for the heaviside function to be written.
x=@(a,n) a^n;
X=ztrans(x(a,n));
pretty(X);
fprintf('the Region of convergence is |z|>a\n');


%% Question 1 part b

syms k c;
y=@(a,n)n*x(a,n);
w=@(a,n,c)c^n*x(a,n);
t=@(a,n)symsum(x(a,k),k,0,n);
r=@(a,n)x(a,n+1);


Y=ztrans(y(a,n));
W=ztrans(w(a,n,c));
T=ztrans(t(a,n));
R=ztrans(r(a,n));
fprintf('y=nx(n) and its ztransform are:\n');
pretty(y(a,n));
pretty(Y);
fprintf('this transform is (-z) time X(z)''\n\n\n');

fprintf('w=c^nx(n) and its ztransform are:\n');
pretty(w(a,n,c));
pretty(W);
fprintf('this transform is X(z/c)\n\n\n');


fprintf('t=(sum(x(k))from zero to n) and its ztransform are:\n');
pretty(t(a,n));
pretty(T);
fprintf('the transform of the given serires is equal to:\nX(z)/(z-1)+z/(z-a).\n had the summation been up to (n-1) terms,the second term\n i.e. z/(z-a) would have been omitted.\n\n');


fprintf('r=x(n+1) and its ztransform are:\n');
pretty(r(a,n));
pretty(R);
fprintf('Had this transform been taken completely and two-sided,it would have been\n equal z times X(z).but matlab only calculates onesided Z transform\n thus in this case it results in z times (X(z)-a) \n\n');


%% Question 2 part 1

% matlab's Heaviside function has a slight difference with what we already
% know.the difference is the point n=0 where instead of 1,it returns 0.5
% by simple tests I found that instead of heaviside(n),it should be written
% heaviside(n+1),and it is obvious why.
clear
syms n
x1=@(n)(1/3^n+3^n)*heaviside(n+1);
X1=ztrans(x1(n));
pretty(X1);
fprintf('with ROC :|z|>3\n');
%% part 2
% the requested transform is calculated and is as follows:
%     X2(z)=z/(z-1/3)+z/(z-3)
syms z
Xx2=@(z)z/(z-1/3)+z/(z-3);
fprintf('the transform is:\n');
pretty(Xx2(z));
fprintf('the transform is the same as X1(z) but with a different ROC\n');
fprintf('with ROC 1/3<|z|<3\n');
%% part 3
fprintf('we can make up the following series that would seem to give the \n same Z transform but the ROC is obviusly diffenrent.\nwhere as in the first series mentioned in the question\n the ROC is a disk, the latter has no region of convergence\n');
syms n
x3=@(n)3^n*heaviside(n)-1/3^n*heaviside(-n-1);
pretty(x3(n));
fprintf('or we can make this one:\n');
x3p=@(n)-3^n*heaviside(-n-1)-1/3^n*heaviside(-n-1);
pretty(x3p(n));
%% part 4
b=[2 -10/3 0];
a=[1 -10/3 1];
Z=[1 1/z 1/z^2];
X2=@(z)sum(b.*Z)/sum(a.*Z);
y=iztrans(X2(z));
fprintf('the inverse z transform of X2(z) is:\n');
disp(y);
fprintf('the above series is equal to x1(n)\n');
%% part 5 & 6
[r p k]=residuez(b,a);
figure();
clear r;
r=roots(b);
zplane(r,p);
title('zero-pole plot of X3(z)');
fprintf('the ROC of X1(z) is |z|>3\n');
fprintf('the ROC of X2(z) is 1/3<|z|<3\n');
fprintf('the transform of X3(z) does not converge\n\n\n\n');

fprintf('in two sided Z-transform given only the z transform without the ROC\n we cannot determine the initial series.\nIn the one sided Z-transform we can determine a unique series\nwhich its transform would be the given transform.\n');

%% Question 3 part 1 & 2
clear 
syms z n t;
X=@(z)1./(z-2)-1./(z-1);
pretty(X(z));
fprintf('Since it is assumed that all Z transform are one-sided\nthus the ROC is |z|>2\n');
fprintf('the inverse z transform by seperation of fractions:\n');
g=@(t) 5*exp(1i*t);
gp=@(t) 5i*exp(1i*t);
integrand=@(n,z)z.^(n-1)./((z-1).*(z-2));
I=@(n)integral(@(t)integrand(n,g(t)).*gp(t),0,2*pi)/(pi*2i);
%
%
x=iztrans(X(z));
pretty(x);
[r,p,k]=residuez([1 0 0],[1,-3,2]);
%zplane([],[1;2])
figure();
zplane([],p)
title('zero-pole plot of X(z)');

%% part 4 & 5
fprintf('since the evaluating contour must be inside the ROC\nthus any circle of radius higher than 2 is sufficient\n');
xns=zeros(1,15);
for k=0:14
fprintf('x(%d)=',k);
    xns(k+1)=I(k);
    disp(xns(k+1));
end
n=0:14;
xns=real(xns);
plot(n,xns,'*');




