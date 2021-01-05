clear all

dx = 0.01;
Lx = 3;
Ly = 0.5;
gamm = 1.4;
mu = 1e-5;
[xx,yy] = meshgrid([0:dx:Lx],[0:dx:Ly]);
Nx = size(xx,2);
Ny = size(xx,1);

rho = xx*0 + 1.1;
%p = xx*0 + 1 + exp(-100*( (xx-0.5).^2 + (yy-0.5).^2) );
p = xx*0 + 1e5;
u = xx*0 + 4 ;
v = xx*0 ;

U(:,:,1) = rho;
U(:,:,2) = rho.*u;
U(:,:,3) = rho.*v;
U(:,:,4) = 0.5*rho.*(u.^2+v.^2) + p/(gamm-1);

t = 0;

while t<1.0
    rho = U(:,:,1);
    u = U(:,:,2)./rho;
    v = U(:,:,3)./rho;
    p = (U(:,:,4) - 0.5*rho.*(u.^2 + v.^2))*(gamm-1);
    a = sqrt(gamm*p./rho);
    c_max = max(max(a + sqrt(u.^2 + v.^2)));
    dt = 0.1*dx/c_max;
    
    UL = U(:,1:end-1,:);
    UR = U(:,2:end,:);
    
    ULW = 0.5*(UL+UR) - 0.5*(dt/dx)*(Fx(UR,gamm,mu,dx)-Fx(UL,gamm,mu,dx));
    
    FLW = Fx(ULW,gamm,mu,dx);
    Fforce = 0.25*(Fx(UL,gamm,mu,dx) + 2*FLW + Fx(UR,gamm,mu,dx) - (dx/dt)*(UR - UL));
    U(:,2:end-1,:) = U(:,2:end-1,:) -(dt/dx)*(Fforce(:,2:end,:)- Fforce(:,1:end-1,:));
    
    UL = U(1:end-1,:,:);
    UR = U(2:end,:,:);
    
    ULW = 0.5*(UL+UR) - 0.5*(dt/dx)*(Fz(UR,gamm,mu,dx)-Fz(UL,gamm,mu,dx));
    
    FLW = Fz(ULW,gamm,mu,dx);
    Fforce = 0.25*(Fz(UL,gamm,mu,dx) + 2*FLW + Fz(UR,gamm,mu,dx) - (dx/dt)*(UR - UL));
    U(2:end-1,:,:) = U(2:end-1,:,:) -(dt/dx)*(Fforce(2:end,:,:)- Fforce(1:end-1,:,:));
    
    
    U(1,:,:) = U(2,:,:);
    U(end,:,:) = U(end-1,:,:);
    U(:,1,:) = U(:,2,:);
    U(:,end,:) = U(:,end-1,:);
    
    U(:,1,2) = 4;
    U(24:26,24:28,2) = 0;
    U(24:26,24:28,3) = 0;
    
%     surf(p);
%     pause(0.01);
    t = t+dt

end
