dx = 0.01;

[xx,yy] = meshgrid([0:0.01:1]);

u = xx*0 + 1;
dt = dx*dx/4;
for t = 1:50000
    
    u(2:end-1,2:end-1) = u(2:end-1,2:end-1) + (dt/(dx*dx))*(u(3:end,2:end-1) + u(1:end-2,2:end-1) + u(2:end-1,3:end) + u(2:end-1,1:end-2) - 4*u(2:end-1,2:end-1)) + dt*u(2:end-1,2:end-1);
%    surf(u);
%     pause(0.01);
    u(:,1) = 0;
    u(:,end) = 0;
    u(1,:) = 1;
    u(end,:) = 0;
end
    surf(u)
    
    %gradient operators:
    
    
    ux = (u(3:end,2:end-1) - u(1:end-2,2:end-1))/(2*dx);
    uy = (u(2:end-1,3:end) - u(2:end-1,1:end-2))/(2*dx);
    uxx = (u(3:end,2:end-1) + u(1:end-2,2:end-1) - 2*u(2:end-1,2:end-1))/dx^2;
    uyy = (u(2:end-1,1:end-2) + u(2:end-1,3:end) - 2*u(2:end-1,2:end-1))/dx^2;