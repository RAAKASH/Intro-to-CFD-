%% Project 
% Stream vorticity approach
% Equation 1 : psi(i,j) = (w*(dx)^2 + psi(i-1,j)+ psi(i+1,j) + psi(i,j-1) + psi(i,j+1))/4
% Equation 2 : u(i,j) = (psi(i,j+1) - psi(i,j-1))/(2*dy)
% Equation 3 : w(i,j) = (w(i+1,j)*[1 - 2*dx*u(i,j)/(2*gamma)] + w(i-1,j)*[1 + 2*dx*u(i,j)/(2*gamma)] 
%                        ... w(i,j+1)*[1 - 2*dx*v(i,j)/(2*gamma)] + w(i,j-1)*[1 + 2*dx*v(i,j)/(2*gamma)])/4

close all;
clear; 
clc;
%% Variable initialization
x = 1; % X length
y = 1; % Y length
dx = 0.0100; % Grid Size for X direction
dy = 0.0100; % Grid Size for Y direction

Nx = x/dx + 1; % No of Grid Points in the X direction
Ny = y/dy + 1; % No of Grid Points in the Y direction

j1 = 20;  % Specialized for Project
j2 = 25;

j3 = 5;
j4 = 10;

psi_1 = 0;  % Stream function Value specialized for project 
u0 = 1;     % Velocity of the Lid (m/s)

psi = zeros(Ny,Nx);
w   = zeros(Ny,Nx);
u   = zeros(Ny,Nx);
v   = zeros(Ny,Nx);

Re = 800;
gamma = u0/Re;
alpha = 1; % Relaxation parameter for 


% dt = 0.0001; % Time step
dt = 0.2/gamma/(1/dx^2 + 1/dy^2); % Minimum time step for least computational expense
%% Error variables initialization
err = 10;  
err1 = err;
ERR1 = err1;
ERR2 = err1;
iter = 1;

fprintf('Variables Initialized \n');
%% Boundary conditions
psi(1,:) = psi_1;
psi(1:(j1-1),1) = psi_1;
psi(1:(j3-1),end) = psi_1;
%psi(j1:j2,1) = psi_1-psi_1*(1:(j2-j1+1))/(j2-j1+1); %Project specialized
w(:,end) = -2*u0/dx;
u(end,:) = u0;
fprintf('Boundary condition imposed \n');
%% 
%  while((ERR2>10^-3))
  while(dt*iter<100)
     W = w;
    PSI = psi;
    psi = streamfunc( w ,psi ,x,j1,j2,j3,j4,alpha);
   
    [u,v] = velocity( u,v,psi,x,j1,j2,j3,j4);
     %u
   
    [ w,iter1] = omega( u,v,psi,u0,w,x,gamma,j1,j2,j3,j4,Re,dt);
     %w
   
%Different Error calculations 
err = rms(rms((W - w)));
err1 = rms(rms((PSI - psi)));
iter = iter+1;
[ERR1,ERR2] = rmse_psi( psi,w,dx,gamma,v,u );

% figure(1); clf
%     subplot(2,1,1),  pcolor(0:dx:x,0:dy:y,psi);   caxis([-0.1 0.01]);  hold on; colorbar; shading interp;  axis square; 
%     subplot(2,1,2),  pcolor(0:dx:x,0:dy:y,w); caxis([-50  50]);  hold on; colorbar; shading interp;  axis square; drawnow

fprintf('Iter in main');
iter

    %Non Convergent Checks
    if((err>100))
        fprintf('Omega Non convergent');
        err
    break;
    end

 end

%% Plotting
close all;
Nx = 1/dx+1;
Ny = 1/dx+1;
contourf(0:dx:x,0:dy:y,psi,[-linspace(-0.0005,0.0013,20),-logspace(-6.64,-1,30)],'ShowText','off');
caxis([-0.05 0.005])
colorbar;
title('Stream Function Contours');
axis equal
pause;
figure(1); clf
    subplot(2,1,1),  pcolor(0:dx:x,0:dy:y,psi);   caxis([-0.0005 0.0014]);  hold on; colorbar; shading interp;  axis square; 
%     subplot(2,1,2),  pcolor(0:dx:x,0:dy:y,w); caxis([-50  50]);  hold on; colorbar; shading interp;  axis square; drawnow
 
pause;
%% 
close all;
vlevels= linspace(-0.0005,0.0013,20);
contourf(0:dx:x,0:dy:.2,psi(1:length(0:dy:.2),length(0:dx:0):Nx),vlevels,'ShowText','off');
title('Stream Function Contours');
colorbar;
axis equal
pause
figure(1); clf
    subplot(2,1,1),  pcolor(.6:dx:x,0:dy:.3,psi(1:length(0:dy:.3),length(0:dx:0.6):Nx));   caxis([-0.0005 0.0014]);  hold on; colorbar; shading interp;  axis square; 
%     subplot(2,1,2),  pcolor(0:dx:x,0:dy:y,w); caxis([-50  50]);  hold on; colorbar; shading interp;  axis square; drawnow
 
pause;