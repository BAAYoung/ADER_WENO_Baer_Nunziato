%clear axis
x_handle = xlabel('Cell','Fontsize',24);
y_handle = ylabel('Solid Volume Fraction, \alpha_1','Fontsize',24);

set(x_handle,'Fontname','Lucida bright');
set(y_handle,'Fontname','Lucida bright');
axis square

[xx,yy] = meshgrid([0:0.2:20]);
R = (sign( (xx-10).^2 + (yy-10).^2 - 4) +1)/2;
hold on
plot(alpha18(3,:),'k')
plot(alpha116(3,:),'k')
plot(alpha124(3,:),'k')
plot(alpha132(3,:),'k')
plot(alpha140(3,:),'k')
plot(alpha148(3,:),'k')
