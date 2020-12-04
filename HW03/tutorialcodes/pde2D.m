%defining geometry
gd = [1;0;0;1];
ns = [67:49];
sf = 'C1';
dl = decsg(gd,sf,ns);
model = createpde(1);
pg = geometryFromEdges(model,dl);
pdegplot(model,'EdgeLabels','on')
axis([-1.2 1.2 -1.2 1.2])

%% Boundary condition
applyBoundaryCondition(model,'dirichlet','edge',1,'u',0);
%applyBoundaryCondition(model,'dirichlet','edge',2,'u',0);
%applyBoundaryCondition(model,'dirichlet','edge',3,'u',0);
%applyBoundaryCondition(model,'dirichlet','edge',4,'u',0);

%% initial condition
ut0 = @(locations) 0;
u0  = @(locations) 1-locations^2;

%% solving PDE
generateMesh(model,'Hmax',0.05);
setInitialConditions(model,u0,ut0);
specifyCoefficients(model,'m',1,'d',0,'c',1,'a',0,'f',0);
tlist = 0:0.25:5;
results = solvepde(model,tlist);
u = results.NodalSolution;

%% plot 
for k=1:length(tlist)
a = pdeplot(model,'XYData',u(:,k),'ZData',u(:,k))
axis([0 1 0 1 -1 1])
pause(0.5)
end