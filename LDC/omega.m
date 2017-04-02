function [ w ,iter] = omega( u,v,psi,u0,w,x,gamma,j1,j2,j3,j4,Re,dt,alpha1)
%% Solving Navier Stokes Equation 
%% Central difference Scheme
[Ny,Nx] = size(psi);
dx = x/(Nx-1);
dy = dx;
err = 10;
iter = 1;
alpha1 ; %To check if PSOR has any stability difference compared to Gauss seidel

%% Boundary Conditions -Edges
for j =1:(j1-1) %left
% w(j,1) = 2*(psi(j,1)-psi(j,2))/dx^2 +2*(psi(j,2)-psi(j,1))/dx; 
% w(j,end) = 2*(psi(j,end)-psi(j,end-1))/dx^2 -2*(psi(j,end)-psi(j,end-1))/dx; 
w(j,1) = 2*(psi(j,1)-psi(j,2))/dx^2; 
end

for j =(j1):j2 %left
%w(j,end) = -2*(psi(j+1,1)-2*psi(j,1)+psi(j-1,1))/dx^2 ;
w(j,1) = 2*(psi(j,1)-psi(j,2))/dx^2; 
end

for j =(j2+1):Ny %left
% w(j,1) = 2*(psi(j,1)-psi(j,2))/dx^2 +2*(psi(j,2)-psi(j,1))/dx; 
% w(j,end) = 2*(psi(j,end)-psi(j,end-1))/dx^2 -2*(psi(j,end)-psi(j,end-1))/dx; 
w(j,1) = 2*(psi(j,1)-psi(j,2))/dx^2; 
end


for j =1:(j3-1) %right
w(j,end) = 2*(psi(j,end)-psi(j,end-1))/dx^2 ; 
end

for j =(j3):j4%right
%w(j,end) = -2*(psi(j,end)-psi(j,end-1))/dx^2 ;
w(j,end) = 2*(psi(j,end)-psi(j,end-1))/dx^2; 

end


for j =(j4+1):Ny %right
w(j,end) = 2*(psi(j,end)-psi(j,end-1))/dx^2 ;
end

for i =2:(Nx-1) %bottom
% w(1,i) = 2*(psi(1,i)-psi(2,i))/dy^2 + 2*(psi(2,i)-psi(1,i))/dy; 
% w(end,i) = 2*(psi(end,i)-psi(end-1,i))/dy^2 - 2*(psi(end,i)-psi(end-1,i))/dy; 
 w(1,i) = 2*(psi(1,i)-psi(2,i))/dy^2;  
end
for i =2:(Nx-1) %Top
w(end,i) = 2*(psi(end,i)-psi(end-1,i))/dy^2 - 2*u0/dy; 
end

%% Boundary Conditions Corners
w(1,1) = 0;
w(end,1)=0;
w(end,end)= w(end,end-1);
w(end,1) = w(end,2);
%while(err>0.0001)
W = w;
temp = w;
%% Solving for w 
for i = 2:(Nx-1)
    for j = 2:(Ny-1)
        
    %--------------------- Done initially[Gauss seidel method](M1)(Converged only for Re<60)-------------------------------
%       w(j,i) = (w(j,i+1)*(1 - dx*u(j,i)/(2*gamma)) + w(j,i-1)*(1 + dx*u(j,i)/(2*gamma)) ...
%                + w(j+1,i)*(1 - dx*v(j,i)/(2*gamma)) + w(j-1,i)*(1 + dx*v(j,i)/(2*gamma)))/4 ;

    %--------------------- Done at the 3rd attempt[Jacobi method(same as the time step (almost))](Converged for Re<250)(M2)---------------
%      temp(j,i) =    (w(j,i+1)*(1 - dx*u(j,i)/(2*gamma)) + w(j,i-1)*(1 + dx*u(j,i)/(2*gamma)) ...
%                     + w(j+1,i)*(1 - dx*v(j,i)/(2*gamma)) + w(j-1,i)*(1 + dx*v(j,i)/(2*gamma)))/4 ;
     
    % ---------------------- Time step (reference -IISC)( Converges based on the FTCS condition )(M3) ------------------            
%     w(j,i) =  w(j,i) -dt*(Upwind(w,u,v,i,j,dx,dy,0.5*((i<(Nx-1))&&(j<(Ny-1))&&(j>2)&&(i>2)))+( u(j,i)*(w(j,i+1)-w(j,i-1)))/(2*dx) ...
%                           +(v(j,i)*(w(j+1,i)-w(j-1,i)))/(2*dx) ...
%                           - gamma*( w(j,i+1)+w(j,i-1)+ w(j+1,i)+ w(j-1,i) - 4*w(j,i))/(dx^2)   );

 %--------------------- Done 4th[PSOR Gauss seidel method](M4)(Converges for Re<1500 depending on alpha1)-------------------------------
    w(j,i) = (1-alpha1)*w(j,i)+ alpha1*(w(j,i+1)*(1 - dx*u(j,i)/(2*gamma)) + w(j,i-1)*(1 + dx*u(j,i)/(2*gamma)) ...
               + w(j+1,i)*(1 - dx*v(j,i)/(2*gamma)) + w(j-1,i)*(1 + dx*v(j,i)/(2*gamma)))/4 ;



    end
        
end
% Uncomment below for (M3)
% w(2:(Ny-1),2:(Nx-1)) = w(2:(Ny-1),2:(Nx-1)) +temp(2:(Ny-1),2:(Nx-1));

% Uncomment below for (M2)
%  w(2:(Ny-1),2:(Nx-1)) = temp(2:(Ny-1),2:(Nx-1));

iter = iter+1;
err = rms(rms(W-w));

%end
end
