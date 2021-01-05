%clear axis
%figure

[xx,yy] = meshgrid([0:0.16:20]);
R = (sign( (xx-10).^2 + (yy-10).^2 - 16) +1)/2;
%imagesc(alpha120);
x_handle = xlabel('Cell','Fontsize',24);
y_handle = ylabel('Normalized Density','Fontsize',24);
%hcb=colorbar;

%title(hcb,'\alpha_1','Fontsize',16,'Fontname','Lucida bright');
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
%caxis([0.4 1]);

colormap bone