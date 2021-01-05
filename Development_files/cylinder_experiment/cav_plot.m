[xx,yy] = meshgrid([0:0.002:0.2],[0:0.002:0.14]);
[xp,yp] = meshgrid([0:0.0002:0.2],[0:0.0002:0.14]);
R = (1- sign(0.038*0.038 - (xp-0.1).^2 - (yp-0.07).^2))/2;
Vq = interp2(xx,yy,alpha1830,xp,yp);

surf(xp,yp,sign((Vq-0.601).*R));
shading interp
box on

colormap bone