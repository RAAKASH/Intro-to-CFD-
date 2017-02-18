function [ TotalTime ] = Assignment1( L,t,alpha,Nx,dt)
 %% CFD Assignment1 -Intro
 % One dimentional unsteady heat conduction equation
 close all;
 %% Variable initialization -1
  dx = (L/(Nx-1)); % Distance differential in m
 
 m = round(t/dt+1); % No of grid points in time
 T = zeros(m,Nx); % Grid generation ,Initial condition
  
 T(:,1) = 1;  %Boundary condition
 T(:,Nx) = 0; %Boundary condition
 
 %% CSFT scheme 
 % dT/dt = alpha d2T/d2x
 % Tn+1(i) = Tn(i) + alpha*dt(Tn(i+1)+2*Tn(i)+Tn(i-1))/dx2
 t = cputime;
 for n = 2:m
     for i = 2:(Nx-1)
      T(n,i) =  T(n-1,i) + alpha*dt*(T(n-1,i+1)-2*T(n-1,i)+T(n-1,i-1))/dx/dx;
     
     end
 end
 TotalTime = cputime - t;
 %% Plotting data for t = 0.1,0.5,1,5,10,15,20 s
    
    plot(0:dx:L , T(0.1/dt+1,:),0:dx:L , T(0.5/dt+1,:),0:dx:L , T(1/dt+1,:),0:dx:L , T(5/dt+1,:),0:dx:L , T(10/dt+1,:),0:dx:L , T(15/dt+1,:),0:dx:L , T(20/dt+1,:));
     xlabel('Length along rod')
     ylabel('Temperatures')
     legend('At 0.1s','At 0.5s','At 1s','At 5s','At 10s','At 15s','At 20s');
     s1 = num2str(dt);
     s2 = 'For dt =' ;
     s4= num2str(Nx);
     s3 = strcat(s2,s1,'s','Nx=',s4,'- Explicit scheme');
     title(s3);
% pause


print(strcat(s3,'.jpg'),'-dpng')

end

