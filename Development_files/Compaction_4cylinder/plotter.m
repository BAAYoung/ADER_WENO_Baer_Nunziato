%clear axis
%figure

[xx,yy] = meshgrid([0:0.2:20]);
R = (sign( (xx-10).^2 + (yy-10).^2 - 16) +1)/2;
imagesc(u1_mag.*R + (1-R).*mean(u1_mag));
x_handle = xlabel('Cell','Fontsize',24);
y_handle = ylabel('Cell','Fontsize',24);
hcb=colorbar;

title(hcb,'|u_1|','Fontsize',16,'Fontname','Lucida bright');
set(x_handle,'Fontname','Lucida bright');
set(y_handle,'Fontname','Lucida bright');
axis square
grad_rho = exp(-20*(abs(rho18(2:end-1,3:end) - rho18(2:end-1,1:end-2))  -20*abs(rho18(3:end,2:end-1) - rho18(1:end-2,2:end-1)))./(0.2*2*1000*rho18(2:end-1,2:end-1))); 

%hold on
% plot(rho28(3,11:end),'k')
% plot(rho216(3,11:end),'k')
% plot(rho224(3,11:end),'k')
% plot(rho232(3,11:end),'k')
% plot(rho240(3,11:end),'k')
% plot(rho248(3,11:end),'k')
%caxis([0 16e-4]);

colormap bone