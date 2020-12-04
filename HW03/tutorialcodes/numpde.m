dz = 0.1;                           % each depth step is 0.1 meter
Nz = ceil(1 / dz);                % Choose the number of depth steps 
Nt = 200;                          % Choose the number of time steps
dt = 5/Nt;                          % Length of each time step in seconds 
K = 1/8;                           
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
        plot(z,T(:,i))
        pause(0.1)
end