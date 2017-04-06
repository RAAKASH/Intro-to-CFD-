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
tic; %CPU time 
x = 0.3; % X length
y = 0.3; % Y length
dx = 0.00200; % Grid Size for X direction
dy = 0.00200; % Grid Size for Y direction
% dx = x/128;  %For Ghia Ghia and Shin reference
% dy = x/128;
% dx = x/160;  %For O. BOTELLA reference
% dy = x/160;
Nx = x/dx + 1; % No of Grid Points in the X direction
Ny = y/dy + 1; % No of Grid Points in the Y direction

% j1 =  20;  % Specialized for Project
% j2 = 25;  %Specialized for project
%j3 = 5;   %Specialized for project
%j4 = 10;  %Specialized for project
j1 = round((20-1)*0.01/0.3/(dy/y)+1);% Specialized for Project
j2 = round((25-1)*0.01/0.3/(dy/y)+1);% Specialized for Project
j3 = round((5-1)*0.01/0.3/(dy/y)+1);% Specialized for Project
j4 = round((10-1)*0.01/0.3/(dy/y)+1);% Specialized for Project

psi_1 = -0.1;  % Stream function Value specialized for project ,Non Zero Positive Number for Project
u0 = 1;     % Velocity of the Lid (m/s)

%General Variable Declaration
psi = zeros(Ny,Nx);
w   = zeros(Ny,Nx);
u   = zeros(Ny,Nx);
v   = zeros(Ny,Nx);

Re = 100;
gamma = u0/Re;
alpha = 1.5; % Relaxation parameter for stream function (to Increase speed of convergence)
alpha1 = 1.9; % Relaxation parameter for navierstokes function convergence

% dt = 0.00001; % Time step
dt = 0.2/gamma/(1/dx^2 + 1/dy^2); % Minimum time step for least computational expense
t = 10000; %Total time for computation

f = 0 ;% flag
%% Error variables initialization - parameters
err = 10;    % Percent error of omega
err1 = err;  % Peercent error of psi
ERR1 = err1; % Error in streamfunction convergence
ERR2 = err1; % Error in navierstrokes statisfying
iter = 1;

fprintf('Variables Initialized \n');
%% Boundary conditions
psi(1,:) = psi_1;  % Bottom
psi(1:(j1-1),1) = psi_1; % Left
psi(1:(j3-1),end) = psi_1; %Right
psi(j1:j2,1) = psi_1-psi_1*(1:(j2-j1+1))/(j2-j1+1); %Project specialized
psi(j3:j4,end) = psi_1-psi_1*(1:(j2-j1+1))/(j2-j1+1); %Project specialized
w(end,:) = -2*u0/dx;  %Vorticity
u(end,:) = u0;
fprintf('Boundary condition imposed \n');
%% Solving Using Stream Vorticity Approach
  while((ERR2>9*10^-8))
%  while(dt*iter<t)
    W = w; % Old Values of vorticity
    PSI = psi; % Old value of psi
    [psi,f] = streamfunc( w ,psi ,x,j1,j2,j3,j4,alpha,u,v);
   % Checking Break at StreamFunction  
    if(f==1)
     fprintf('Error in Streamfunc \n');
     alpha = alpha -0.01;
     [ psi,w,u,err,err1,ERR1,ERR2,iter] = BC( psi,psi_1,u0,j1,j2,j3,j4,dx);
     continue
    end
    
    [u,v] = velocity( u,v,psi,x,y,j1,j2,j3,j4);
   
    [w,iter1] = omega( u,v,psi,u0,w,x,y,gamma,j1,j2,j3,j4,Re,dt,alpha1);
   
   %Different Error calculations 
    [ERR1,ERR2] = rmse_psi( psi,w,dx,gamma,v,u,x );
    err = rms(rms((W - w)))/rms(rms(W));
    err1 = rms(rms((PSI - psi)))/rms(rms(PSI));


    if(mod(iter,200)==0)
      fprintf('***********Iteration in main**************** \nIteration Number : %d \n Error : %d \n',iter,ERR2);
      Error(iter/200) = ERR2;
      
    end
    %Non Convergent Checks ,optimizing Alpha
    if((err>100)||(isnan(err)))
       fprintf('Omega Non convergent,Error: %d - \n *********Wait Testing New Parameter************* \n',err);
       alpha1 = alpha1 - 0.01;
       fprintf('New Value of relaxation parameter : %d  \n',alpha1); %Need to reinitialize values
       [ psi,w,u,err,err1,ERR1,ERR2,iter] = BC( psi,psi_1,u0,j1,j2,j3,j4,dx);
       
       if(alpha1 <= 0)
        clc;
        fprintf('**********No Solution******************');
        break;
       end
       continue    
    end
    iter = iter+1;
    
    
    
  end
  TotalCPUtime = toc;
%% Check -Re 100,400,1000,3200 (Lid Driven Cavity) -With Ghia
if(psi_1==0)
if(Re ==100)
close all;
u1 = [1,0.84123,0.78871,0.73722,0.68717,0.23151,0.00332,-0.13641,-0.20581,-0.21090,-0.15662,-0.10150,-0.06434 ,-0.04775,-0.04192 ,-0.03717,0];
y1 = [129,126,125,124,123,110,95,80,65,59,37,23,14,10,9,8,1];
plot(u1,(y1-1)/128,'r +',u(:,(Ny+1)/2),(0:(Ny-1))*dy,'b');
legend('Ghia Ghia and Shin','My Code');
xlabel('U - Velocity at Mid plane');
ylabel('Ny');
title('Results Comparison');
pause;
v1 = [0.00000,-0.05906,-0.07391,-0.08864,-0.10313,-0.16914,-0.22445,-0.24533,0.05454,0.17527,0.17507,0.16077,0.12317,0.10890,0.10091,0.09233,0.00000];
Nx1 = [129.0000  125.0000  124.0000  123.0000  122.0000  117.0000  111.0000  104.0000   65.0000   31.0000   30.0000   21.0000   13.0000 11.0000   10.0000    9.0000    1.0000
    1.0000    0.9688    0.9609    0.9531    0.9453    0.9063    0.8594    0.8047    0.5000    0.2344    0.2266    0.1563  0.0781    0.0703    0.0625   0  0.0938  ];
  close all
plot(v1,(Nx1(1,:)-1)/128,'r +',v((Nx+1)/2,:),(0:(Nx-1)*dx),'b');
legend('Ghia Ghia and Shin','My Code');
xlabel('V - Velocity at Mid plane');
ylabel('Nx');
title('Results Comparison');
pause;
end

if(Re ==400)
close all;
u1 = [1,0.75837,0.68439,0.61756,0.55892,0.29093,0.16256,0.02135,-0.11477,-0.17119,-0.32726 ,-0.24299 ,-0.14612,-0.10338,-0.09266 ,-0.08186,0.00000];
y1 = [129,126,125,124,123,110,95,80,65,59,37,23,14,10,9,8,1];
plot(u1,(y1-1)/128,'r +',u(:,(Ny+1)/2),(0:(Ny-1))*dy,'b');
legend('Ghia Ghia and Shin','My Code');
xlabel('U - Velocity at Mid plane');
ylabel('Ny');
title('Results Comparison');
pause;
v1 = [0.00000,-0.12146,-0.15663,-0.19254,-0.22847,-0.23827,-0.44993,-0.38598,0.05188,0.30174,0.30203,0.28124,0.22965,0.20920,0.19713,0.18360,0.00000];
Nx1 = [129.0000  125.0000  124.0000  123.0000  122.0000  117.0000  111.0000  104.0000   65.0000   31.0000   30.0000   21.0000   13.0000 11.0000   10.0000    9.0000    1.0000
    1.0000    0.9688    0.9609    0.9531    0.9453    0.9063    0.8594    0.8047    0.5000    0.2344    0.2266    0.1563  0.0781    0.0703    0.0625   0  0.0938  ];
  close all
plot(v1,(Nx1(1,:)-1)/128,'r +',v((Nx+1)/2,:),(0:(Ny-1))*dy,'b');
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
u2 = -[-1;-0.664422700000000;-0.580835900000000;-0.516927700000000;-0.472332900000000;-0.337221200000000;-0.188674700000000;-0.0570178000000000;0.0620561000000000;0.108199900000000;0.280369600000000;0.388569100000000;0.300456100000000;0.222895500000000;0.202330000000000;0.181288100000000;0];
plot(u1,(y1-1)/128,'r +',u2,(y1-1)/128,'g x',u(:,(Ny+1)/2),(0:(Ny-1))*dy,'b');
legend('Ghia Ghia and Shin','O. BOTELLA','My Code');
xlabel('U - Velocity at Mid plane');
ylabel('Ny');
title('Results Comparison');
pause;
v1 = [0.0,-0.21388,-0.27669,-0.33714,-0.39188,-0.51550,-0.42665,-0.31966,0.02526,0.32235,0.33075,0.37095,0.32627,0.30353,0.29012,0.27485,0.00000];
v2 = [0;-0.227922500000000;-0.293686900000000;-0.355321300000000;-0.410375400000000;-0.526439200000000;-0.426454500000000;-0.320213700000000;0.0257995000000000;0.325359200000000;0.333992400000000;0.376918900000000;0.333044200000000;0.309909700000000;0.296270300000000;0.280705600000000;0];
Nx1 = [129.0000  125.0000  124.0000  123.0000  122.0000  117.0000  111.0000  104.0000   65.0000   31.0000   30.0000   21.0000   13.0000 11.0000   10.0000    9.0000    1.0000
    1.0000    0.9688    0.9609    0.9531    0.9453    0.9063    0.8594    0.8047    0.5000    0.2344    0.2266    0.1563  0.0781    0.0703    0.0625   0  0.0938  ];
  close all
plot(v1,(Nx1(1,:)-1)/128,'r +',v2,(Nx1(1,:)-1)/128,'g x',v((Nx+1)/2,:),(0:(Ny-1))*dx,'b');
legend('Ghia Ghia and Shin','O. BOTELLA','My Code');
xlabel('V - Velocity at Mid plane');
ylabel('Nx');
title('Results Comparison');
pause;
end
if(Re ==3200)
    close all;
v1 = [0.00000,-0.39017,-0.47425,-0.52357,-0.54053,-0.44307,-0.37401,-0.31184,0.00999,0.28188,0.29030,0.37119,0.42768,0.41906,0.40917,0.39560,0.00000]; 
x1 = [129,126,125,124,123,110,95,80,65,59,37,23,14,10,9,8,1];
plot(v1,(x1-1)/128,'r +',v((Nx+1)/2,:),(0:(Nx-1))*dx,'b');
legend('Ghia Ghia and Shin','My Code');
xlabel('V - Velocity at Mid plane');
ylabel('Nx');
title('Results Comparison');
pause;
close all;
u1 =[1,0.53236,0.48296,0.46547,0.46101,0.34682,0.19791,0.07156,-0.04272,-0.86636,-0.24427,-0.34323,-0.41933,-0.37827,-0.35344,-0.32407,0];
y1 = [129,126,125,124,123,110,95,80,65,59,37,23,14,10,9,8,1];
plot(u1,(y1-1)/128,'r +',u(:,(Nx+1)/2),(0:(Ny-1))*dy,'b');
legend('Ghia Ghia and Shin','My Code');
xlabel('U - Velocity at Mid plane');
ylabel('Ny');   
title('Results Comparison');
end
end
%% Plotting -General stream Function Plot
close all;
Nx = x/dx+1;
Ny = y/dx+1;
contourf(0:dx:x,0:dy:y,psi,linspace(-1,1,500),'ShowText','off');
caxis([-0.15 0.02])
colorbar;
title('Stream Function Contours');
axis equal
pause;

close all;
contourf(0:dx:x,0:dy:y,w,'ShowText','on');
colorbar;
title('Vorticity Contours');
axis equal
pause;
figure(1); clf
    subplot(2,1,1),  pcolor(0:dx:x,0:dy:y,psi);   caxis([-0.0005 0.0014]);  hold on; colorbar; shading interp;  axis square; 
    subplot(2,1,2),  pcolor(0:dx:x,0:dy:y,w); caxis([-2  1]);  hold on; colorbar; shading interp;  axis square; drawnow
 
pause;
 subplot(1,1,1),  pcolor(0:dx:x,0:dy:y,psi);   caxis([-0.12 0.00]);  hold on; colorbar; shading interp;  axis square; 

%% U ,V velocity plots
 close all;
 plot(u(:,(Nx+1)/2),0:dx:x);
 xlabel('U velocity')
 ylabel('y - Displacement');
 title('U vs y');
 pause;
 close all;
  plot(v((Ny+1)/2,:),0:dy:y);
 xlabel('V velocity')
 ylabel('x -Displacement');
 title('V vs x');
pause
%%  Secondary Vortices Plot
close all;
vlevels= linspace(-0.0005,0.0013,10);
contourf(0:dx:x,0:dy:(y/2),psi(1:length(0:dy:(y/2)),length(0:dx:0):Nx),vlevels,'ShowText','off');
title('Stream Function Contours');
colorbar;
axis equal
pause
close all;
vlevels= linspace(-0.0005,0.0013,20);
contourf(0:dx:x/2,y/2:dy:1,psi((length(0:dy:y/2)-1):end,1:length(0:dx:x/2)),vlevels,'ShowText','off');
title('Stream Function Contours');
colorbar;
axis equal
pause
close all;
figure(1); clf
    subplot(1,1,1),  pcolor(.6:dx:x,0:dy:.3,psi(1:length(0:dy:.3),(length(0:dx:0.6)):Nx));   caxis([-0.0005 0.0014]);  hold on; colorbar; shading interp;  axis square; 
pause;
%% Report
x1 = [129,126,125,124,123,110,95,80,65,59,37,23,14,10,9,8,1];
x11 = (x1-1)/128;
N = round(x11/dx + 1);
a = u(N,(Nx+1)/2);
b = v((Ny+1)/2,N)';