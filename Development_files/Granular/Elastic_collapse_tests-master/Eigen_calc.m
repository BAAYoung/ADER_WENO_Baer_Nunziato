clear all

dx = 0.01;
Lx = 4;
Ly = 0.25;


lambda = 1;
zeta = pi*50/(180);
%zeta = 0;
delta = 23*pi/180;
%delta = 0;
neta = 0;

[xx,yy] = meshgrid([0:dx:Lx],[0:dx:Ly]);
Nx = size(xx,2);
Ny = size(xx,1);

%h = xx*0 + 0.2;
%h = 0.21 -0.19*sign(xx-2.5);
h = xx*0 + 0.2;
%h = 0.5 + 0.1*real(sqrt(0.2-(xx-0.5).^2 - (yy-0.5).^2));
%h = xx*0 + 0.2;
u = xx*0 ;
v = xx*0;
T = xx*0 + 1;
omega = 200;
eta = 0;
U(:,:,1) = h;
U(:,:,2) = u.*h;
U(:,:,3) = v.*h;
U(:,:,4) = T.*h;

phi = sqrt(abs(0.501.^2 - (xx-3).^2 - (yy-2).^2)).* sign(0.501.^2 - (xx-3).^2 - (yy-2).^2); 
phi = xx.*0 - 1;
phix = phi*0;
phiy = phi*0;

phix(2:end-1,2:end-1) = (phi(2:end-1,3:end) - phi(2:end-1,1:end-2))/(2*dx);
phiy(2:end-1,2:end-1) = (phi(3:end,2:end-1) - phi(1:end-2,2:end-1))/(2*dx);

nx = phix./sqrt(phix.^2 + phiy.^2);
ny = phiy./sqrt(phix.^2 + phiy.^2);
nx(isnan(nx)) = 0;
ny(isnan(ny)) = 0;


X_circ2D = ceil((sign ( 0.501.^2 - (xx-3).^2 - (yy-2).^2) +1)/2);

for i = 2:3
    X_circ(:,:,i) = X_circ2D;
end
X_circ(:,:,1) = X_circ2D*0.5;

t_accu = 0;
for t = 1:10000
    
    h = U(:,:,1);
    u = U(:,:,2)./h;
    v = U(:,:,3)./h;
    c_max = max(max(sqrt( h*cos(zeta)  ) +  sqrt(u.^2 + v.^2)));
    dt = 0.45*dx/c_max;
    
    UL = U(:,1:end-1,:);
    UR = U(:,2:end,:);
    
    ULW = 0.5*(UL+UR) - 0.5*(dt/dx)*(Fx(UR,zeta,lambda)-Fx(UL,zeta,lambda));
    
    FLW = Fx(ULW,zeta,lambda);
    Fforce = 0.25*(Fx(UL,zeta,lambda) + 2*FLW + Fx(UR,zeta,lambda) - (dx/dt)*(UR - UL));
    U(:,2:end-1,:) = U(:,2:end-1,:) -(dt/dx)*(Fforce(:,2:end,:)- Fforce(:,1:end-1,:));
    
    UL = U(1:end-1,:,:);
    UR = U(2:end,:,:);
    
    ULW = 0.5*(UL+UR) - 0.5*(dt/dx)*(Fz(UR,zeta,lambda)-Fz(UL,zeta,lambda));
    
    FLW = Fz(ULW,zeta,lambda);
    Fforce = 0.25*(Fz(UL,zeta,lambda) + 2*FLW + Fz(UR,zeta,lambda) - (dx/dt)*(UR - UL));
    U(2:end-1,:,:) = U(2:end-1,:,:) -(dt/dx)*(Fforce(2:end,:,:)- Fforce(1:end-1,:,:));
    
    
    
    
    
    
    u_n = u./sqrt(u.^2 + v.^2);
    v_n = v./sqrt(u.^2 + v.^2);
    u_n(isnan(u_n))=1;
    v_n(isnan(v_n))=1;
    
   U(:,:,2) = U(:,:,2) + dt*h.*(sin(zeta) - cos(zeta)*tan(delta).*u_n) ;
   U(2:end-1,2:end-1,2) = U(2:end-1,2:end-1,2) + eta*( U(3:end,2:end-1,2) + U(1:end-2,2:end-1,2) - 2*U(2:end-1,2:end-1,2) )*dt/(dx^2);
   U(:,:,3) = U(:,:,3) - dt*h.*cos(zeta).*tan(delta).*v_n;
   U(2:end-1,2:end-1,4) = U(2:end-1,2:end-1,4)  + U(2:end-1,2:end-1,2).*eta.*( U(3:end,2:end-1,2) + U(1:end-2,2:end-1,2) - 2*U(2:end-1,2:end-1,2) )*dt/(dx^2) + 0.5*(dt/dx)*(U(3:end,2:end-1,2)-U(1:end-2,2:end-1,2)).^2;
   U(:,:,4) = U(:,:,4) + dt*U(:,:,4)*neta;
    
    U(1,:,:) = U(2,:,:);
    U(end,:,:) = U(end-1,:,:);
    
    U(:,1,:) = U(:,2,:);
    U(:,end,:) = U(:,end-1,:);
    
    
    U(1,:,3) = 0;
    U(end,:,3) = 0;
    U(1,:,2) = 0;
    U(end,:,2) = 0;
    
    
    
%     U(1,:,4) = 1*U(1,:,1);
%     U(end,:,4) = 1*U(end,:,1);
    
    
    %U(X_circ==1) = 0;
    %U(X_circ==0.5) = 0.01;
    
    dtau = dt/10;
    for tau = 1:100
        for i = 2:Nx-1
            for j = 2:Ny-1
                if phi(j,i)>=0
                    if nx(j,i)>0
                        dhdx = (U(j,i,1) - U(j,i-1,1))/dx;
                        dTdx = (U(j,i,4) - U(j,i-1,4))/dx;
                    elseif nx(j,i)<0
                        dhdx = ( U(j,i+1,1) - U(j,i,1))/dx;
                        dTdx = ( U(j,i+1,4) - U(j,i,4))/dx;
                    else 
                        dhdx = 0;
                        dTdx = 0;
                    end
                    
                    if ny(j,i)>0
                        dhdy = (U(j,i,1) - U(j-1,i,1))/dx;
                        dTdy = (U(j,i,4) - U(j-1,i,4))/dx;
                    elseif ny(j,i)<0
                        dhdy = ( U(j+1,i,1) - U(j,i,1))/dx;
                        dTdy = ( U(j+1,i,4) - U(j,i,4))/dx;
                    else 
                        dhdy = 0;
                        dTdy = 0;
                    end
                    
                    U(j,i,1) = U(j,i,1) - dtau*(dhdy.*ny(j,i) + dhdx.*nx(j,i));
                    U(j,i,4) = U(j,i,4) - dtau*(dTdy.*ny(j,i) + dTdx.*nx(j,i));
                end
            end
        end
    end
                        
    
    
    
%     U(36:45,46:55,2) = 0;
%     U(36:45,46:55,3) = 0;
%     U(36:45,46:55,1) = 0.01;
    T = U(:,:,4)./U(:,:,1);
        
    surf(h);
    %axis equal;
    %axis([0 150 0 100 0.4 0.6]);
    %axis equal
    %axis([0 21 0 101 0 10])
%     if isnan(sum(sum(sum(U))))
%         disp('stop')
%         pause(100);
%     end
    pause(0.01);
    t_accu = t_accu + dt;
end