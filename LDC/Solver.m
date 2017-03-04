function [ w,u,v,psi,f] = Solver( w,psi,u,v,j1,j2,j3,j4,x,u0,Re,dt,alpha,psi_1 )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% Solver
[Ny,Nx] = size(w);
dx = x/(Nx-1);
dy=dx;
y = x;


psi(1,:) = psi_1;
psi(1:(j1-1),1) = psi_1;
psi(1:(j3-1),end) = psi_1;
%psi(j1:j2,1) = psi_1-psi_1*(1:(j2-j1+1))/(j2-j1+1);

u(end,:) = u0;


iter = 1;
err = 10;

f=0;
gamma = u0/Re;
while((err>10^-4))
    W = w;
    [psi,g] = streamfunc( w ,psi ,x,j1,j2,j3,j4,alpha);
   
    [u,v] = velocity( u,v,psi,x,j1,j2,j3,j4);
    
   
    [w] = omega( u,v,psi,u0,w,x,gamma,j1,j2,j3,j4,Re,dt);
    
   
 
err = rms(rms((W - w)));

iter = iter+1;

   if((err>10^3)||(sum(sum((w>10^3)))>1)||(g==1))
     f = 1 ;
     fprintf('Non convergent,Broke at solver \n');
     break;
   end
end
 
%% Plotting
if(f~=1)
figure(1); clf
    subplot(2,1,1),  pcolor(0:dx:x,0:dy:y,psi);   caxis([-0.1 0.01]);  hold on; colorbar; shading interp;  axis square; 
    subplot(2,1,2),  pcolor(0:dx:x,0:dy:y,w); caxis([-50  50]);  hold on; colorbar; shading interp;  axis square; drawnow
 
end
end

