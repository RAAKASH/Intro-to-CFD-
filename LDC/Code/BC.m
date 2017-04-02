function [ psi,w,u,err,err1,ERR1,ERR2,iter] = BC( psi,psi_1,u0,j1,j2,j3,j4,dx)
%% Imposes Boundary Conditions
%
[Ny,Nx] = size(psi);
psi = zeros(Ny,Nx);
w   = zeros(Ny,Nx);
u   = zeros(Ny,Nx);
psi(1,:) = psi_1;  % Bottom
psi(1:(j1-1),1) = psi_1; % Left
psi(1:(j3-1),end) = psi_1; %Right
psi(j1:j2,1) = psi_1-psi_1*(1:(j2-j1+1))/(j2-j1+1); %Project specialized
psi(j3:j4,end) = psi_1-psi_1*(1:(j2-j1+1))/(j2-j1+1); %Project specialized
w(:,end) = -2*u0/dx;  %Vorticity
u(end,:) = u0;
fprintf('Boundary condition imposed \n \n');

err = 10;  
err1 = err;
ERR1 = err1; %Error in streamfunction convergence
ERR2 = err1; %Error in navierstrokes statisfying
iter = 1;

end

