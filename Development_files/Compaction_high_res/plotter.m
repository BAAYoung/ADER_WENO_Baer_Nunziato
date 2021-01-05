%clear axis
%figure

[xx,yy] = meshgrid([0:0.16:20]);
R = (sign( (xx-10).^2 + (yy-10).^2 - 16) +1)/2;
imagesc(m130.*R + (1-R).*mean(m130));
x_handle = xlabel('Cell','Fontsize',24);
y_handle = ylabel('Cell','Fontsize',24);
hcb=colorbar;

title(hcb,'\rho','Fontsize',16,'Fontname','Lucida bright');
set(x_handle,'Fontname','Lucida bright');
set(y_handle,'Fontname','Lucida bright');
axis square

%hold on
% plot(rho28(3,11:end),'k')
% plot(rho216(3,11:end),'k')
% plot(rho224(3,11:end),'k')
% plot(rho232(3,11:end),'k')
% plot(rho240(3,11:end),'k')
% plot(rho248(3,11:end),'k')
caxis([2 5]);

colormap bone