function [ TotalTime,T ] = Assignment2(L,t,alpha,Nx,dt )
 
%% CFD Assignment 2 -Intro (BTCS scheme)
 % One dimentional unsteady heat conduction equation
close all;

 %% Variable initialization -1
 dx = (L/(Nx-1)); % Distance differential in m
 
 m = round(t/dt); % No of grid points in time
 T = zeros(m,Nx); % Grid generation ,Initial condition
  
 T(:,1) = 1;  %Boundary condition
 T(:,Nx) = 0; %Boundary condition

 %% BTCS scheme 
 % dT/dt = alpha d2T/d2x
 % Tn+1(i) = Tn(i) + alpha*dt(Tn+1(i+1)-2*Tn+1(i)+Tn+1(i-1))/dx2
 % Tn+1(i)-gamma*Tn+1(i+1)+2*gamma*Tn+1(i)-gamma*Tn+1(i-1)= Tn(i) 
 % -gamma*Tn+1(i-1)+Tn+1(i)(2*gamma+1) + -gamma*Tn+1(i+1) =  Tn(i)
 
  t = cputime; % Calculating Time
 for n = 2:m
     M = zeros(Nx);
     M(1,1) =1/dt;
     M(Nx,Nx) =1/dt;
     
      % Matrix Construction
      for i = 2:(Nx-1)
      M(i,(i-1):(i+1)) =  [-alpha/(dx^2),1/dt+2*alpha/(dx^2),-alpha/(dx^2)];
      end
      % A = M\(T(n-1,:)'/dt); %CHECK 1
      
      % Upper triangular Matrix Conversion
      M = [M,(T(n-1,:)'/dt)];
      for i = 2:(Nx-1)
      M(i,:) = M(i,:) - M(i-1,:)*(M(i,i-1)/M(i-1,i-1));    
      end
     
      %B = M(:,1:(end-1))\M(:,end);%CHECK 2
    
   
      % Solving
      X = M(:,1:(end-1));
      T(n,1)=1;
      T(n,Nx)=0;
      for i = (Nx-1):-1:(2)
      T(n,i) = ( M(i,end) - X(i,i+1)*T(n,i+1))/M(i,i);
      end
           
      %D(:,:,n)=M; %CHECK 3
 end
 % T(n,:)- B' %CHECK 4
  
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

