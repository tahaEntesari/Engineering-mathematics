%% Question 1 Fourier Transform
clear
syms x w;
f=@(x)(1/(1+x^2));
g=@(x)(exp(-abs(x)));
h=@(x)(x*exp(-x)*heaviside(x));
l=@(x)(x*exp(x)*heaviside(-x));
fft=eval(['@(w)' char(fourier(f(x)))]);
gft=eval(['@(w)' char(fourier(g(x)))]);
hft=eval(['@(w)' char(fourier(h(x)))]);
lft=eval(['@(w)' char(fourier(l(x)))]);
%fft=@(w)(fourier(f,x,w));
%gft=@(w)(fourier(g,x,w));
%hft=@(w)(fourier(h,x,w));
%lft=@(w)(fourier(l,x,w));

fprintf('f(x) fourier transform:\n');
pretty(fft(w));
fprintf('g(x) fourier transform:\n');
pretty(gft(w));
fprintf('h(x) fourier transform:\n');
pretty(hft(w));
fprintf('l(x) fourier transform:\n');
pretty(lft(w));

figure(1);
fplot(abs(fft(w)));
title('f fourier transform');
xlabel('w');
ylabel('fourier transform');

figure(2);
fplot(abs(gft(w)));
title('g fourier transform');
xlabel('w');
ylabel('fourier transform');

figure(3);
fplot(abs(hft(w)));
title('h fourier transform');
xlabel('w');
ylabel('fourier transform');

figure(4);
fplot(abs(lft(w)));
title('l fourier transform');
xlabel('w');
ylabel('fourier transform');
%% Part b
FcrossG=eval(['@(w)' char(fft(w)*gft(w))]);
HcrossL=eval(['@(w)' char(hft(w)*lft(w))]);
fprintf('F(w)*G(w) :\n');
pretty(FcrossG(w));
fprintf('H(w)*L(w) :\n');
pretty(HcrossL(w));
FcrossGinvf=eval(['@(x)' char((ifourier(FcrossG(w))))]);
HcrossLinvf=eval(['@(x)' char((ifourier(HcrossL(w))))]);

fprintf('F cross G inverse fourier transform:\n');
pretty(FcrossGinvf(x));
fprintf('According to matlab''s documentation it is said that "If ifourier cannot transform the input, then it returns an unevaluated call to fourier." which is the case that is happening with F cross G');
fprintf('H cross L inverse fourier transform:\n');
pretty(HcrossLinvf(x));
figure(5)
fplot(FcrossGinvf(x));
title('F cross G inverse fourier transform');
xlabel('x');
figure(6)
fplot(HcrossLinvf(x));
title('H cross L inverse fourier transform');
xlabel('x');

%% Part c
syms y
fconvg=eval(['@(x)' char(int(f(y)*g(x-y),y,-inf,inf))]);
fprintf('f convolve g is :\n');
pretty(fconvg(x))
hconvl=eval(['@(x)' char(int(h(y)*l(x-y),y,-inf,inf))]);
fprintf('h convolve l is :\n');
pretty(hconvl(x));
figure(8)
fplot(fconvg(x));
title('convolution of f and  g');
xlabel('x');
ylabel('f*g');
figure(9)
fplot(hconvl(x));
title('convolution of h and  l');
xlabel('x');
ylabel('h*l');

%% Question 2 part a
%The function defintions in Question 2 are listed at the end of the file

clear

x=linspace(0,1,100);
t=0:.5:30;
%t=linspace(0,30,100)
%t=[5,15,25];
u=pdepe(0,@equation,@initialcondition,@boundarycondition,x,t);
figure(10);
surf(x,t,u) 
title('Numerical solution computed with 100 mesh points.')
xlabel('Distance x')
ylabel('Time t')

figure(11)
plot(x,u(find(t==5),:))
title('t=5');
xlabel('x');
ylabel('tempreture');
axis([0 1 -0.1 1])

figure(12)
plot(x,u(find(t==15),:))
title('t=15');
xlabel('x');
ylabel('tempreture');
axis([0 1 -0.1 1])

figure(13)
plot(x,u(find(t==25),:))
title('t=25');
xlabel('x');
ylabel('tempreture');
axis([0 1 -0.1 1])




%% Part b

clear

x=linspace(0,1,11);
t=0:.5:15;
%t=linspace(0,30,100)
%t=[5,15,25];
u=pdepe(0,@equation2,@initialcondition2,@boundarycondition2,x,t);
figure(14);
surf(x,t,u) 
title('Numerical solution computed with 100 mesh points.')
xlabel('Distance x')
ylabel('Time t')

figure(15)

plot(x,u(find(t==2),:))
title('t=2');
xlabel('x');
ylabel('tempreture');
axis([0 1 -0.1 1])

figure(16)
plot(x,u(find(t==5),:))
title('t=5');
xlabel('x');
ylabel('tempreture');
axis([0 1 -0.1 1])

figure(17)
plot(x,u(find(t==15),:))
title('t=15');
xlabel('x');
ylabel('tempreture');
axis([0 1 -0.1 1])

%% Question 3 part a
dz = 0.1;                           % each depth step is 0.1 meter
Nz = ceil(1 / dz) ;               % Choose the number of depth steps 
Nt = 300;                          % Choose the number of time steps
dt = 15/Nt;                          % Length of each time step in seconds 
K = 1/10;                           
z = linspace(0,1,Nz+1); 

T = ones(Nz+1,Nt+1);        % Create temperature matrix with Nz+1 rows, and Nt+1 columns
                       
time = [0:dt:dt*Nt];
T(1,:) = 1;  
T(end,:) = 0;  
T(:,1) = 1 - z.^2  ;

for i=2:Nt+1
        depth_2D = (T(1:end-2,i-1)-2*T(2:end-1,i-1)+T(3:end,i-1))/dz^2;
        time_1D = K*depth_2D ;
        T(2:end-1,i) = time_1D*dt + T(2:end-1,i-1);
end
figure(18)
plot(z,T(:,find(time==2)))
title('t=2');
xlabel('x');
ylabel('tempreture');


figure(19)
plot(z,T(:,find(time==5)))
title('t=5');
xlabel('x');
ylabel('tempreture');


figure(20)
plot(z,T(:,find(time==15)))
title('t=15');
xlabel('x');
ylabel('tempreture');
%% part b

figure(21)
plot(z,T(:,find(time==2))-(u(find(t==2),:))')
title('err t=2');
figure(22)
plot(z,T(:,find(time==5))-u(find(t==5),:)')
title('err t=5');
figure(23)
plot(z,T(:,find(time==15))-u(find(t==15),:)')
title('err t=15');
%% Part c

dz = 0.1;                           % each depth step is 0.1 meter
Nz = ceil(1 / dz);                % Choose the number of depth steps 
Nt = 150    ;                          % Choose the number of time steps
dt = 15/Nt;                          % Length of each time step in seconds 
K = 1/10;                           
z = linspace(0,1,Nz+1); 

d=K*dt/dz^2

T = ones(Nz+1,Nt+1);        % Create temperature matrix with Nz+1 rows, and Nt+1 columns
                       
time = [0:dt:dt*Nt];
T(1,:) = 1;  
T(end,:) = 0;  
T(:,1) = 1 - z.^2  ;

for i=2:Nt+1
        depth_2D = (T(1:end-2,i-1)-2*T(2:end-1,i-1)+T(3:end,i-1))/dz^2;
        time_1D = K*depth_2D ;
        T(2:end-1,i) = time_1D*dt + T(2:end-1,i-1);
end

figure(20)
plot(z,T(:,find(time==15)))
title('t=15');
xlabel('x');
ylabel('tempreture');

%% Question 4
clear
%syms locations
%defining geometry
gd = [1;0;0;1];
ns = [83;81;49];
sf = 'SQ1';
dl = decsg(gd,sf,ns);
model = createpde(1);
pg = geometryFromEdges(model,dl);
pdegplot(model,'EdgeLabels','on');
axis([-1.2 1.2 -1.2 1.2])

% Boundary condition
applyBoundaryCondition(model,'dirichlet','edge',1,'u',0);
applyBoundaryCondition(model,'dirichlet','edge',2,'u',0);
applyBoundaryCondition(model,'dirichlet','edge',3,'u',0);
applyBoundaryCondition(model,'dirichlet','edge',4,'u',0);

% initial condition
%ut0 = @(locations) 0;
u0  = @(locations) 1-(locations.x).^2-(locations.y).^2;

% solving PDE
%generateMesh(model,'Hmax',0.05);
generateMesh(model);
setInitialConditions(model,u0,0);
specifyCoefficients(model,'m',1,'d',0,'c',1,'a',0,'f',2);
tlist = 0:0.5:15;
results = solvepde(model,tlist);
u = results.NodalSolution;

% plot 
for k=1:length(tlist)
    figure(24)
    if(.5*(k-1)==2)
        figure(25)
    end
    if(.5*(k-1)==8)
        figure(26)
    end
    if(.5*(k-1)==15)
        figure(27)
    end
    
a = pdeplot(model,'XYData',u(:,k),'ZData',u(:,k));
title(.5*(k - 1))
axis([-1 1 -1 1])
pause(.2)

end

%%
function [c,b,s] = equation(x,t,u,DuDx)
%a PDE in time and one space dimension.
c = 100;
b = DuDx;
s = exp(-t)+exp(-2*t)*cos(3*pi*x/4);
end

function value = initialcondition(x)
value = 0.5-abs(x-.5) ;
end

function [pl,ql,pr,qr] = boundarycondition(xl,ul,xr,ur,t)
%for a PDE in time and one space dimension.
pl = ul ;
ql = 0;
pr = -exp(-t) ;
qr = 1;
end


function [c,b,s] = equation2(x,t,u,DuDx)
%a PDE in time and one space dimension.
c = 40;
b = DuDx;
s = -4*u^2;
end

function value = initialcondition2(x)
value = 1-x^2 ;
end

function [pl,ql,pr,qr] = boundarycondition2(xl,ul,xr,ur,t)
%for a PDE in time and one space dimension.
pl = ul-1 ;
ql = 0;
pr = ur ;
qr = 0;
end


