%% Input Variables
 close ;
 clear ;
 clc;


 L = 10 ; % Length of rod in meters
 t = 25;  % Max time of observation in seconds
 alpha = 0.5; % SI units
 

 Nx = 11:10:101;  % No of grid points in space
 dt = [0.001,0.01,0.1];% Time differential in seconds.
 %% Running both schemes , comparing time
 for i=1:1:5
      for j=1:3
    e(i,j)= Assignment1( L,t,alpha,Nx(i),dt(j)); % Explicit
    %f(i,j)= Assignment2( L,t,alpha,Nx(i),dt(j)); % Implicit
    f(i,j)= BTCS( L,t,alpha,Nx(i),dt(j)); % Implicit
     end
 end
 
