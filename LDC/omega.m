function [ w ,iter] = omega( u,v,psi,u0,w,x,gamma,j1,j2,j3,j4,Re,dt)
%% Solving for u,v using w
%% Central difference
[Ny,Nx] = size(psi);
dx = x/(Nx+1);
dy = dx;
err = 10;
iter = 1;
alpha1 = 0.45; %To check if PSOR has any stability difference compared to Gauss seidel
for j =1:(j1-1)
% w(j,1) = 2*(psi(j,1)-psi(j,2))/dx^2 +2*(psi(j,2)-psi(j,1))/dx; 
% w(j,end) = 2*(psi(j,end)-psi(j,end-1))/dx^2 -2*(psi(j,end)-psi(j,end-1))/dx; 
w(j,1) = 2*(psi(j,1)-psi(j,2))/dx^2; 
end

for j =(j1):j2
%w(j,end) = -2*(psi(j+1,1)-2*psi(j,1)+psi(j-1,1))/dx^2 ;
w(j,1) = 2*(psi(j,1)-psi(j,2))/dx^2; 
end

for j =(j2+1):Ny
% w(j,1) = 2*(psi(j,1)-psi(j,2))/dx^2 +2*(psi(j,2)-psi(j,1))/dx; 
% w(j,end) = 2*(psi(j,end)-psi(j,end-1))/dx^2 -2*(psi(j,end)-psi(j,end-1))/dx; 
w(j,1) = 2*(psi(j,1)-psi(j,2))/dx^2; 
end


for j =1:(j3-1)
w(j,end) = 2*(psi(j,end)-psi(j,end-1))/dx^2 ; 
end

for j =(j3):j4
%w(j,end) = -2*(psi(j,end)-psi(j,end-1))/dx^2 ;
w(j,end) = 2*(psi(j,end)-psi(j,end-1))/dx^2; 

end


for j =(j4+1):Ny
w(j,end) = 2*(psi(j,end)-psi(j,end-1))/dx^2 ;
end

for i =2:(Nx-1)
% w(1,i) = 2*(psi(1,i)-psi(2,i))/dy^2 + 2*(psi(2,i)-psi(1,i))/dy; 
% w(end,i) = 2*(psi(end,i)-psi(end-1,i))/dy^2 - 2*(psi(end,i)-psi(end-1,i))/dy; 
 w(1,i) = 2*(psi(1,i)-psi(2,i))/dy^2; 
 w(end,i) = 2*(psi(end,i)-psi(end-1,i))/dy^2 - 2*u0/dy; 
end

%while(err>0.0001)
W = w;
temp = w;
for i = 2:(Nx-1)
    for j = 2:(Ny-1)
        
    %--------------------- Done initially[Gauss seidel method](M1)(Converged only for Re<60)-------------------------------
%       w(j,i) = (w(j,i+1)*(1 - dx*u(j,i)/(2*gamma)) + w(j,i-1)*(1 + dx*u(j,i)/(2*gamma)) ...
%                + w(j+1,i)*(1 - dx*v(j,i)/(2*gamma)) + w(j-1,i)*(1 + dx*v(j,i)/(2*gamma)))/4 ;

    %--------------------- Done at the 3rd attempt[Jacobi method](Converged for Re<250)(M2)---------------
%      temp(j,i) =    (w(j,i+1)*(1 - dx*u(j,i)/(2*gamma)) + w(j,i-1)*(1 + dx*u(j,i)/(2*gamma)) ...
%                     + w(j+1,i)*(1 - dx*v(j,i)/(2*gamma)) + w(j-1,i)*(1 + dx*v(j,i)/(2*gamma)))/4 ;
     
    % ---------------------- Done 2nd (reference)( Converges well at  Re=500 )(M3) ------------------            
%     temp(j,i) =   -dt*(( u(j,i)*(w(j,i+1)-w(j,i-1)))/(2*dx) ...
%                           +(v(j,i)*(w(j+1,i)-w(j-1,i)))/(2*dx) ...
%                           - ( w(j,i+1)+w(j,i-1)+ w(j+1,i)+ w(j-1,i) - 4*w(j,i))/(Re*dx^2)   );

 %--------------------- Done 4th[PSOR Gauss seidel method](M4)(Converges for Re<1500 depending on alpha1)-------------------------------
    w(j,i) = (1-alpha1)*w(j,i)+ alpha1*(w(j,i+1)*(1 - dx*u(j,i)/(2*gamma)) + w(j,i-1)*(1 + dx*u(j,i)/(2*gamma)) ...
               + w(j+1,i)*(1 - dx*v(j,i)/(2*gamma)) + w(j-1,i)*(1 + dx*v(j,i)/(2*gamma)))/4 ;



    end
        
end
% Uncomment below for (M3)
% w(2:(Ny-1),2:(Nx-1)) = w(2:(Ny-1),2:(Nx-1)) +temp(2:(Ny-1),2:(Nx-1));

% Uncomment below for (M2)
% w(2:(Ny-1),2:(Nx-1)) = temp(2:(Ny-1),2:(Nx-1));

iter = iter+1;
err = rms(rms(W-w));

%end
end
