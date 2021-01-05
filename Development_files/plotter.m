
x_handle = xlabel('Cell','Fontsize',24);
y_handle = ylabel('Velocity, u','Fontsize',24);

set(x_handle,'Fontname','Lucida bright');
set(y_handle,'Fontname','Lucida bright');
axis square

[xx,yy] = meshgrid([0:0.2:20]);
R = (sign( (xx-10).^2 + (yy-10).^2 - 4) +1)/2;