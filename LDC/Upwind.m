function [ fin ] = Upwind(w,u,v,i,j,dx,dy,q)
% Upwind scheme
% q = 0.5;
if(q~=0)
w_x1 = (w(j,i-2)-3*w(j,i-1)+3*w(j,i)-w(j,i+1))/(3*dx); 
w_y1 = (w(j-2,i)-3*w(j-1,i)+3*w(j,i)-w(j+1,i))/(3*dy);
w_x2 = (w(j,i-1)-3*w(j,i)+3*w(j,i+1)-w(j,i+2))/(3*dx);
w_y2 = (w(j-1,i)-3*w(j,i)+3*w(j+1,i)-w(j+2,i))/(3*dy);
u1 = min(u(j,i),0);
u2 = max(u(j,i),0);
v1 = min(v(j,i),0);
v2 = max(v(j,i),0);
fin  = q*(u2*w_x1+u1*w_x2 +v2*w_y1 +v1*w_y2);
else
    fin=0;
end

end

