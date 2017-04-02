function [ u,v ] = velocity( u,v,psi,x,j1,j2,j3,j4)
%% Solving for u,v using w
%% Central difference
[Ny,Nx] = size(psi);
dx = x/(Nx-1);
dy = dx;
for j= j1:j2
 u(j,1) = (psi(j+1,1) - psi(j,1))/(dy);
end

for j= j3:j4
 u(j,end) = (psi(j+1,end) - psi(j,end))/(dy);
end


for i= 2:(Nx-1)
    for j = 2:(Ny-1)
    u(j,i) =  (psi(j+1,i) - psi(j-1,i))/(2*dy);
    v(j,i) = -(psi(j,i+1) - psi(j,i-1))/(2*dx);
    end
end


end

