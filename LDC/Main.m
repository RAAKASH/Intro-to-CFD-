%% Project 
% Stream vorticity approach
% Equation 1 : psi(i,j) = (w*(dx)^2 + psi(i-1,j)+ psi(i+1,j) + psi(i,j-1) + psi(i,j+1))/4
% Equation 2 : u(i,j) = (psi(i,j+1) - psi(i,j-1))/(2*dy)
% Equation 3 : w(i,j) = (w(i+1,j)*[1 - 2*dx*u(i,j)/(2*gamma)] + w(i-1,j)*[1 + 2*dx*u(i,j)/(2*gamma)] 
%                        ... w(i,j+1)*[1 - 2*dx*v(i,j)/(2*gamma)] + w(i,j-1)*[1 + 2*dx*v(i,j)/(2*gamma)])/4

close all;
clear; 
clc;
% Note to Get results for Lid Driven Cavity Comment project specialized statements in functions
%% Variable initialization
x = 1; % X length
y = 1; % Y length
dx = 0.0100; % Grid Size for X direction
dy = 0.0100; % Grid Size for Y direction
dx = 1/128;  %For Ghia Ghia and Shin reference
dy = 1/128;
Nx = x/dx + 1; % No of Grid Points in the X direction
Ny = y/dy + 1; % No of Grid Points in the Y direction

j1 = 20;  % Specialized for Project
j2 = 25;  %Specialized for project

j3 = 5;   %Specialized for project
j4 = 10;  %Specialized for project

psi_1 = 0;  % Stream function Value specialized for project 
u0 = 1;     % Velocity of the Lid (m/s)

%General Variable Declaration
psi = zeros(Ny,Nx);
w   = zeros(Ny,Nx);
u   = zeros(Ny,Nx);
v   = zeros(Ny,Nx);

Re = 400;
gamma = u0/Re;
alpha = 1.5; % Relaxation parameter for stream function (to Increase speed of convergence)
alpha1 = 0.3; % Relaxation parameter for navierstokes function convergence

% dt = 0.0001; % Time step
dt = 0.4/gamma/(1/dx^2 + 1/dy^2); % Minimum time step for least computational expense
t = 1000; %Total time for computation
%% Error variables initialization
err = 10;  
err1 = err;
ERR1 = err1; %Error in streamfunction convergence
ERR2 = err1; %Error in navierstrokes statisfying
iter = 1;

fprintf('Variables Initialized \n');
%% Boundary conditions
psi(1,:) = psi_1;  
psi(1:(j1-1),1) = psi_1;
psi(1:(j3-1),end) = psi_1;
psi(j1:j2,1) = psi_1-psi_1*(1:(j2-j1+1))/(j2-j1+1); %Project specialized
psi(j3:j4,end) = psi_1-psi_1*(1:(j2-j1+1))/(j2-j1+1); %Project specialized
w(:,end) = -2*u0/dx;  
u(end,:) = u0;
fprintf('Boundary condition imposed \n');
%% 
%  while((ERR2>10^-3))
  while(dt*iter<t)
     W = w;
    PSI = psi;
    psi = streamfunc( w ,psi ,x,j1,j2,j3,j4,alpha);
   
    [u,v] = velocity( u,v,psi,x,j1,j2,j3,j4);
     %u
   
    [ w,iter1] = omega( u,v,psi,u0,w,x,gamma,j1,j2,j3,j4,Re,dt,alpha1);
     %w
   
%Different Error calculations 
err = rms(rms((W - w)))/rms(rms(W));
err1 = rms(rms((PSI - psi)))/rms(rms(PSI));
iter = iter+1;
 [ERR1,ERR2] = rmse_psi( psi,w,dx,gamma,v,u );

    if(mod(iter,200)==0)
      fprintf('Iter in main');
      iter
    end
    %Non Convergent Checks
    if((err>100))
       fprintf('Omega Non convergent');
       err
       alpha1 = alpha1 - 0.05;
       fprintf('New Value of relaxation parameter : %d',alpha1); %Need to reinitialize values
       err = 10;
       err1 = 10;
       if(alpha1 ==0)
        break;
       end
    end

  end
  
%% Check -Re 400,1000,3200 (Lid Driven Cavity)
if(Nx == 129)
if(Re ==400)
close all;
u1 = [1,0.75837,0.68439,0.61756,0.55892,0.29093,0.16256,0.02135,-0.11477,-0.17119,-0.32726 ,-0.24299 ,-0.14612,-0.10338,-0.09266 ,-0.08186,0.00000];
y1 = [129,126,125,124,123,110,95,80,65,59,37,23,14,10,9,8,1];
plot(u1,y1,'r',u(y1,65),y1,'b');
legend('Ghia Ghia and Shin','My Code');
xlabel('U - Velocity at Mid plane');
ylabel('Ny');
title('Results Comparison');
pause;
v1 = [0.00000,-0.12146,-0.15663,-0.19254,-0.22847,-0.23827,-0.44993,-0.38598,0.05188,0.30174,0.30203,0.28124,0.22965,0.20920,0.19713,0.18360,0.00000];
Nx = [129.0000  125.0000  124.0000  123.0000  122.0000  117.0000  111.0000  104.0000   65.0000   31.0000   30.0000   21.0000   13.0000 11.0000   10.0000    9.0000    1.0000
    1.0000    0.9688    0.9609    0.9531    0.9453    0.9063    0.8594    0.8047    0.5000    0.2344    0.2266    0.1563  0.0781    0.0703    0.0625   0  0.0938  ];
  close all
plot(v1,Nx(1,:),'r',v(65,Nx(1,:)),Nx(1,:),'b');
legend('Ghia Ghia and Shin','My Code');
xlabel('V - Velocity at Mid plane');
ylabel('Nx');
title('Results Comparison');
pause;
end

if(Re ==1000)
close all;
u1 = [1.0000,0.65928,0.57492,0.51117,0.46604,0.33304,0.18719,0.05702,-0.06080,-0.10648,-0.27805,-0.38289,-0.29730,-0.22220,-0.20196,-0.18109,0.00000];
y1 = [129,126,125,124,123,110,95,80,65,59,37,23,14,10,9,8,1];
plot(u1,y1,'r',u(y1,65),y1,'b');
legend('Ghia Ghia and Shin','My Code');
xlabel('U - Velocity at Mid plane');
ylabel('Ny');
title('Results Comparison');
pause;
v1 = [0.0,-0.21388,-0.27669,-0.33714,-0.39188,-0.51550,-0.42665,-0.31966,0.02526,0.32235,0.33075,0.37095,0.32627,0.30353,0.29012,0.27485,0.00000];
Nx = [129.0000  125.0000  124.0000  123.0000  122.0000  117.0000  111.0000  104.0000   65.0000   31.0000   30.0000   21.0000   13.0000 11.0000   10.0000    9.0000    1.0000
    1.0000    0.9688    0.9609    0.9531    0.9453    0.9063    0.8594    0.8047    0.5000    0.2344    0.2266    0.1563  0.0781    0.0703    0.0625   0  0.0938  ];
  close all
plot(v1,Nx(1,:),'r',v(65,Nx(1,:)),Nx(1,:),'b');
legend('Ghia Ghia and Shin','My Code');
xlabel('V - Velocity at Mid plane');
ylabel('Nx');
title('Results Comparison');
pause;
end
if(Re ==3200)
    close all;
v1 = [0.00000,-0.39017,-0.47425,-0.52357,-0.54053,-0.44307,-0.37401,-0.31184,0.00999,0.28188,0.29030,0.37119,0.42768,0.41906,0.40917,0.39560,0.00000]; 
x1 = [129,126,125,124,123,110,95,80,65,59,37,23,14,10,9,8,1];
plot(v1,x1,'r',v(65,y1),y1,'b');
legend('Ghia Ghia and Shin','My Code');
xlabel('V - Velocity at Mid plane');
ylabel('Nx');
title('Results Comparison');
pause;
close all;
u1 =[1,0.53236,0.48296,0.46547,0.46101,0.34682,0.19791,0.07156,-0.04272,-0.86636,-0.24427,-0.34323,-0.41933,-0.37827,-0.35344,-0.32407,0];
y1 = [129,126,125,124,123,110,95,80,65,59,37,23,14,10,9,8,1];
plot(u1,y1,'r',u(y1,65),y1,'b');
legend('Ghia Ghia and Shin','My Code');
xlabel('U - Velocity at Mid plane');
ylabel('Ny');   
title('Results Comparison');
end
end
%% Plotting -General stream Function Plot
close all;
Nx = 1/dx+1;
Ny = 1/dx+1;
contourf(0:dx:x,0:dy:y,psi,[-linspace(-0.0005,0.0013,10),-logspace(-6.64,-0.5,30)],'ShowText','off');
caxis([-0.15 0.02])
colorbar;
title('Stream Function Contours');
axis equal
pause;
figure(1); clf
    subplot(1,1,1),  pcolor(0:dx:x,0:dy:y,psi);   caxis([-0.0005 0.0014]);  hold on; colorbar; shading interp;  axis square; 
%     subplot(2,1,2),  pcolor(0:dx:x,0:dy:y,w); caxis([-50  50]);  hold on; colorbar; shading interp;  axis square; drawnow
 
pause;
 subplot(1,1,1),  pcolor(0:dx:x,0:dy:y,psi);   caxis([-0.12 0.00]);  hold on; colorbar; shading interp;  axis square; 

%%  Secondary Vortices Plot
close all;
vlevels= linspace(0.05,0.1,10);
contourf(0:dx:x,0:dy:.5,psi(1:length(0:dy:.5),length(0:dx:0):Nx),vlevels,'ShowText','off');
title('Stream Function Contours');
colorbar;
axis equal
pause
close all;
vlevels= linspace(-0.0005,0.0013,20);
contourf(0:dx:0.5,0.5:dy:1,psi((length(0:dy:0.5)):end,1:length(0:dx:0.5)),vlevels,'ShowText','off');
title('Stream Function Contours');
colorbar;
axis equal
pause
close all;
figure(1); clf
    subplot(1,1,1),  pcolor(.6:dx:x,0:dy:.3,psi(1:length(0:dy:.3),(length(0:dx:0.6)+1):Nx));   caxis([-0.0005 0.0014]);  hold on; colorbar; shading interp;  axis square; 
pause;
