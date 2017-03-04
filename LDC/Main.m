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
x = 1;
y = 1;
dx = 0.025;
dy = 0.025;

Nx = x/dx + 1;
Ny = y/dy + 1;

j1 = 20;
j2 = 25;

j3 = 5;
j4 = 10;

psi_1 = -0.1;
u0 = 1;

psi = zeros(Ny,Nx);
w   = zeros(Ny,Nx);
u   = zeros(Ny,Nx);
v   = zeros(Ny,Nx);

Re = 500;
gamma = u0/Re;
alpha = 1.5;
err = 10;
err1 = err;
iter = 1;


dt = 0.02;
%% Boundary conditions
psi(1,:) = psi_1;
psi(1:(j1-1),1) = psi_1;
psi(1:(j3-1),end) = psi_1;
%psi(j1:j2,1) = psi_1-psi_1*(1:(j2-j1+1))/(j2-j1+1);

u(end,:) = u0;

%% 
 while((err>10^-8))
% while(15/dt/iter>1) % Uncomment for time wise solution
     W = w;
   PSI = psi;
    psi = streamfunc( w ,psi ,x,j1,j2,j3,j4,alpha);
   
    [u,v] = velocity( u,v,psi,x,j1,j2,j3,j4);
     %u
   
    [ w,iter1] = omega( u,v,psi,u0,w,x,gamma,j1,j2,j3,j4,Re,dt);
     %w
   
 
err = rms(rms((W - w)));
err1 = rms(rms((PSI - psi)));
iter = iter+1;

% figure(1); clf
%     subplot(2,1,1),  pcolor(0:dx:x,0:dy:y,psi);   caxis([-0.1 0.01]);  hold on; colorbar; shading interp;  axis square; 
%     subplot(2,1,2),  pcolor(0:dx:x,0:dy:y,w); caxis([-50  50]);  hold on; colorbar; shading interp;  axis square; drawnow

%  if((err>10^4))
%      
%  fprintf('Non convergent');
%  break;
%  end
%  end
end
%% check
close all;
contour(0:dx:x,0:dy:y,psi,'ShowText','on');
title('Stream Function Contours');

figure(1); clf
    subplot(2,1,1),  pcolor(0:dx:x,0:dy:y,psi);   caxis([-0.115 0.01]);  hold on; colorbar; shading interp;  axis square; 
    subplot(2,1,2),  pcolor(0:dx:x,0:dy:y,w); caxis([-50  50]);  hold on; colorbar; shading interp;  axis square; drawnow
 
% contour(0:dx:x,0:dy:y,u,'ShowText','on');
% title('Stream Function Contours');
% contour(0:dx:x,0:dy:y,w,'ShowText','on');
% title('Stream Function Contours');
pause;
%% Testing solver
clc 
clear all
close all;
i = 1;
f=1;
 for u0 = 1:100
for alpha = 0.01:0.01:2
Input
[w,u,v,psi,f] = Solver( w,psi,u,v,j1,j2,j3,j4,x,u0,Re,dt,alpha ,psi_1);
if(f==0)
    alpha
    u0
    f=0;
    break;
end
if(f==0)
break
end
i=i+1
end
 end