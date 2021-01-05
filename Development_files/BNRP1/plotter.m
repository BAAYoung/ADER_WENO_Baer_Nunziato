%clear axis
%figure

xi = [1:0.2:201];
xj = [1:1:201];
x_handle = xlabel('Cell','Fontsize',24);
y_handle = ylabel('Velocity, u_1','Fontsize',24);

set(x_handle,'Fontname','Lucida bright');
set(y_handle,'Fontname','Lucida bright');
axis square

[xx,yy] = meshgrid([0:0.1:20]);
R = (sign( (xx-10).^2 + (yy-10).^2 - 16) +1)/2;
box on
%imagesc(alpha110.*R);
%hold on
% plot(rho28(3,11:end),'k')
% plot(rho216(3,11:end),'k')
% plot(rho224(3,11:end),'k')
% plot(rho232(3,11:end),'k')
% plot(rho240(3,11:end),'k')
% plot(rho248(3,11:end),'k')
%caxis([0.5 1]);

colormap bone