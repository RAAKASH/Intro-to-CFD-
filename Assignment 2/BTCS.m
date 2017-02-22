function [ TotalTime,T ] = BTCS(L,t,alpha,Nx,dt )
 
%% CFD Assignment 2 -Intro (BTCS scheme) (Time wise optimized )
 % One dimentional unsteady heat conduction equation
close all;

 %% Variable initialization -1
 dx = (L/(Nx-1)); % Distance differential in m
 
 m = round(t/dt); % No of grid points in time
 T = zeros(m,Nx); % Grid generation ,Initial condition
  
 T(:,1) = 1;  %Boundary condition
 T(:,Nx) = 0; %Boundary condition
 gamma = alpha*dt/(dx^2);

 %% BTCS scheme 
 % dT/dt = alpha d2T/d2x
 % Tn+1(i) = Tn(i) + alpha*dt(Tn+1(i+1)-2*Tn+1(i)+Tn+1(i-1))/dx2
 % Tn+1(i)-gamma*Tn+1(i+1)+2*gamma*Tn+1(i)-gamma*Tn+1(i-1)= Tn(i) 
 % -gamma*Tn+1(i-1)+Tn+1(i)(2*gamma+1) + -gamma*Tn+1(i+1) =  Tn(i)
 
    t = cputime; % Calculating Time
    % Matrix Construction
     M = zeros(Nx);
     M(1,1) =1;
     M(Nx,Nx) =1;
     
      for i = 2:(Nx-1)
      M(i,(i-1):(i+1)) =  [-gamma,1+2*gamma,-gamma];
      end
      N = M;
    % Upper triangular Matrix Conversion
       Factors = zeros(Nx,1);
      
    for i = 2:(Nx-1)
        Factors(i) = (N(i,i-1)/N(i-1,i-1));
        N(i,:) = N(i,:) - N(i-1,:)*(N(i,i-1)/N(i-1,i-1));
      
    end
      
    
   %% Computing Grid Values   

for n = 2:m
    
      X=T(n-1,:);
      for i = 2:(Nx-1)
      X(i) = (X(i) - X(i-1)*Factors(i));  
      end
    
      % Solving
     T(n,1)=1;
     T(n,Nx)=0;
      for i = (Nx-1):-1:(2)
      T(n,i) = ( X(i) - N(i,i+1)*T(n,i+1))/N(i,i);
      end
      
 end
  
 TotalTime = cputime - t; % Computational time
 %% Plotting data for t = 0.1,0.5,1,5,10,15,20 s
     plot(0:dx:L , T(0.1/dt+1,:),0:dx:L , T(0.5/dt+1,:),0:dx:L , T(1/dt+1,:),0:dx:L , T(5/dt+1,:),0:dx:L , T(10/dt+1,:),0:dx:L , T(15/dt+1,:),0:dx:L , T(20/dt+1,:));
     xlabel('Length along rod')
     ylabel('Temperatures')
     legend('At 0.1s','At 0.5s','At 1s','At 5s','At 10s','At 15s','At 20s');
     s1 = num2str(dt);
     s2 = 'For dt =' ;
     s4= num2str(Nx);
     s3 = strcat(s2,s1,'s','Nx=',s4,'- Implicit scheme');
     title(s3);
      %pause;
print(strcat(s3,'.jpg'),'-dpng')

end

