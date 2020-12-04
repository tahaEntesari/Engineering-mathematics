m = 0;
%Define the solution mesh
x = linspace(0,1,11);
t = linspace(0,5,101);
%Solve the PDE
u = pdepe(m,@equation,@initialcondition,@boundarycondition,x,t);

%Plot solution
figure
for k=1:length(t)
plot(x,u(k,:))
axis([0 1 -0.1 1])
pause(0.1)
end


