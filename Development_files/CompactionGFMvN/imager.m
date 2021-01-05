[xx,yy] = meshgrid([0:0.16:20]);
u_piston = 0.013318;
R = (sign( (xx-10).^2 + (yy-10).^2 - 16) +1)/2;
imagesc(umag20.*R + mean(mean(umag20)).*(1-R));
x_handle = xlabel('Cell','Fontsize',24);
y_handle = ylabel('Cell','Fontsize',24);
hcb=colorbar;

title(hcb,'|u|','Fontsize',16,'Fontname','Lucida bright');
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
caxis([0 14e-3]);

colormap bone