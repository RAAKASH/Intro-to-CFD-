function [ u,v ] = velocity( u,v,psi,x,y,j1,j2,j3,j4)
%% Solving for u,v using w
%% Central difference
[Ny,Nx] = size(psi);
dx = x/(Nx-1);
dy = y/(Ny-1);
for j= j1:j2
 u(j,1) = (psi(j+1,1) - psi(j,1))/(dy);
 v(j,1) = -(-3*psi(j,1)+4*psi(j,2)-psi(j,3))/(2*dx); 
% Since velocity is zero at inlet.

end

for j= j3:j4
 u(j,end) = (psi(j+1,end) - psi(j,end))/(dy);
 v(j,end) = -(-3*psi(j,end)+4*psi(j,end-1)-psi(j,end-2))/(2*dx);
end


for i= 2:(Nx-1)
    for j = 2:(Ny-1)
    u(j,i) =  (psi(j+1,i) - psi(j-1,i))/(2*dy);
    v(j,i) = -(psi(j,i+1) - psi(j,i-1))/(2*dx);
    end
end


end

