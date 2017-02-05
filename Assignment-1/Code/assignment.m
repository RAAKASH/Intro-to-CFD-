 %% CFD Assignment -Intro
 % One dimentional unsteady heat conduction equation
 close ;
 clear ;
 clc;

 %% Variable initialization -1
 Input  % Running input file
 
 m = round(t/dt); % No of grid points in time
 T = zeros(m,Nx); % Grid generation ,Initial condition
  
 T(:,1) = 1;  %Boundary condition
 T(:,Nx) = 0; %Boundary condition
 
 %% CSFT scheme 
 % dT/dt = alpha d2T/d2x
 % Tn+1(i) = Tn(i) + alpha*dt(Tn(i+1)+2*Tn(i)+Tn(i-1))/dx2
 
 for n = 2:m
     for i = 2:(Nx-1)
      T(n,i) =  T(n-1,i) + alpha*dt*(T(n-1,i+1)-2*T(n-1,i)+T(n-1,i-1))/dx/dx;
      %pause;
     end
 end
 
 %% Plotting data for t = 0.1,0.5,1,5,10,15,20 s
     total = 20;
     step = 5;
     j = (0:step:total)/dt;
     plot(0:dx:L , T(0.1/dt+1,:),0:dx:L , T(0.5/dt+1,:),0:dx:L , T(1/dt+1,:),0:dx:L , T(5/dt+1,:),0:dx:L , T(10/dt+1,:),0:dx:L , T(15/dt+1,:),0:dx:L , T(20/dt+1,:));
     xlabel('Length along rod')
     ylabel('Temperatures')
     legend('At 0.1s','At 0.5s','At 1s','At 5s','At 10s','At 15s','At 20s');
     s1 = num2str(dt);
     s2 = 'For dt =' ;
     s3 = strcat(s2,s1,'s');
     title(s3);


