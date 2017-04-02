function [ERR1,ERR2 ] = rmse_psi( psi,w,dx,gamma,v,u,x)
% Finds error matrix
[m,n]=size(psi);
dy=dx;
Nx=x/dx+1;
Ny=Nx;
err1 = zeros(m,n);
err2 =err1;
for i=2:n-1
    for j=2:m-1
     err1(j,i) = (psi(j-1,i)+psi(j+1,i)-4*psi(j,i)+psi(j,i-1)+psi(j,i+1)+w(j,i)*dx^2 )/4;
     err2(j,i)= ( u(j,i)*(w(j,i+1)-w(j,i-1)))/(2*dx) ...
                          +(v(j,i)*(w(j+1,i)-w(j-1,i)))/(2*dx) ...
                          - gamma*( w(j,i+1)+w(j,i-1)+ w(j+1,i)+ w(j-1,i) - 4*w(j,i))/(dx^2)  ;
    end
end
ERR1 = rms(rms(err1));
ERR2 = rms(rms(err2));
end