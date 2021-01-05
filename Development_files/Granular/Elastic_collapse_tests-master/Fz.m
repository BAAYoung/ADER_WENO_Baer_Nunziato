function [ F ] = Fz( U,gamm,mu,dx )
%UNTITLED14 Summary of this function goes here
%   Detailed explanation goes here

    lambda = -2*mu/3;
    rho = U(:,:,1);
    u = U(:,:,2)./rho;
    v = U(:,:,3)./rho;
    p = (U(:,:,4) - 0.5*rho.*(u.^2 + v.^2))*(gamm-1);

    
    ux = u*0;
    uy = u*0;
    vx = u*0;
    vy = u*0;
    %Stress tensors:
    
    ux(:,2:end-1) = (u(:,3:end) - u(:,1:end-2))/(2*dx);
    uy(2:end-1,:) = ( u(3:end,:) - u(1:end-2,:))/(2*dx);
    vx(:,2:end-1) = (v(:,3:end) - v(:,1:end-2))/(2*dx);
    vy(2:end-1,:) = ( v(3:end,:) - v(1:end-2,:))/(2*dx);
    
    
    ux(:,1) = (u(:,2)-u(:,1))/dx;
    ux(:,end) = (u(:,end)-u(:,end-1))/dx;
    
    vx(:,1) = (v(:,2)-v(:,1))/dx;
    vx(:,end) = (v(:,end)-v(:,end-1))/dx;
    
    uy(1,:) = ( u(2,:) - u(1,:) )/dx;
    uy(end,:) = ( u(end,:) - u(end-1,:) )/dx;
    
    vy(1,:) = ( v(2,:) - v(1,:) )/dx;
    vy(end,:) = ( v(end,:) - v(end-1,:) )/dx;
    
    
    
    
    
    
    tau_yy = 2*mu*vy + lambda*(ux + vy);
    tau_xy = mu*(uy + vx);
    
    
    
    F(:,:,1) = rho.*v;
    F(:,:,3) = rho.*v.*v + p - tau_yy;
    F(:,:,2) = rho.*u.*v - tau_xy;
    F(:,:,4) = v.*(U(:,:,4)+p) - u.*tau_xy - v.*tau_yy;

end

