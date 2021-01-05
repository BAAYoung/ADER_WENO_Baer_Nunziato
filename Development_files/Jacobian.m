syms S1 S0 S4 S8 alpha1 rho1 u1 v1 p1 rho2 u2 v2 p2 e1 e2 Fd Ft q Q gamma1 gamma2 pi1 mu



e1 = (p1 + gamma1*pi1)/((gamma1-1)) ;
e2 = p2/((gamma2-1));

alpha2 = 1-alpha1;

U = [alpha1;
    alpha1*rho1;
    alpha1*rho1*u1;
    alpha1*rho1*v1;
    alpha1*rho1*(0.5*u1*u1 + 0.5*v1*v1 + e1);
    alpha2*rho2;
    alpha2*rho2*u2;
    alpha2*rho2*v2;
    alpha2*rho2*(0.5*u2*u2 + 0.5*v2*v2 + e2)];

JAC = jacobian(U,[alpha1 rho1 u1 v1 p1 rho2 u2 v2 p2]); 

SU = [S1;
    0;
    Ft;
    0;
    u1*Ft + Q;
    0;
    -Ft;
    0;
    -u1*Ft - Q];

SS = [mu(p2-p1);
    0;
    0;
    0;
    S4;
    0;
    0;
    0;
    S8];
   