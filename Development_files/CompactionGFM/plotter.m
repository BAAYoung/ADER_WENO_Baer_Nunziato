%clear axis
figure
x_handle = xlabel('Cell','Fontsize',24);
y_handle = ylabel('Gas Density, \rho_2','Fontsize',24);

set(x_handle,'Fontname','Lucida bright');
set(y_handle,'Fontname','Lucida bright');
axis square

[xx,yy] = meshgrid([0:0.2:20]);
R = (sign( (xx-10).^2 + (yy-10).^2 - 4) +1)/2;
hold on
% plot(rho28(3,11:end),'k')
% plot(rho216(3,11:end),'k')
% plot(rho224(3,11:end),'k')
% plot(rho232(3,11:end),'k')
% plot(rho240(3,11:end),'k')
% plot(rho248(3,11:end),'k')
