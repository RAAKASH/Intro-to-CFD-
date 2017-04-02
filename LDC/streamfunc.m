function [ psi,f] = streamfunc( w,psi,x,j1,j2,j3,j4,alpha)
%% Solving laplace equation
%% Point Gauss seidel - With relaxation parameter alpha
err = 10;
iter = 0;
[Ny,Nx] = size(w);
dx = x/(Nx-1);
f=0;
 while((err>0.01))
    iter = iter +1;
    PSI =psi;
    
    
    %psi(j3:j4,end) = 2*psi(j3:j4,end-1) - psi(j3:j4,end-2);
    
    for i= 2:(Nx-1)
    for j= 2:(Ny-1)
         psi(j,i) = (1-alpha)*psi(j,i)+alpha*(psi(j-1,i)+psi(j+1,i)+psi(j,i-1)+psi(j,i+1)+w(j,i)*dx^2 )/4; 
    end
    end
    %  Some Project Specialized Boundary Conditions
    %   psi(j1:j2,1) = psi(j1:j2,2);       %Project specialized
       psi(j3:j4,end) = psi(j3:j4,end-1); %Project Specialized

    err = rms(rms((PSI - psi)));
    
    if(iter>100)
      fprintf('Broken in stream function at \n Iteration no =  %d, Error : %d \n',iter,err);
      f=1; % Flag for Breaking
      break
    end
 end

