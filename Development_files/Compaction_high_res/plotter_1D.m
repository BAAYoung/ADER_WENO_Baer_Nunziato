[xx,yy] = meshgrid([0:0.16:20]);
R = (sign( (xx-10).^2 + (yy-10).^2 - 16) +1)/2;



hold on
h1 = plot(alpha120(63,:),'color',[0.25 0.25 0.25]);
h2 = plot(alpha110(63,:),'color',[0.5 0.5 0.5]);
plot(alpha130(63,:),'k')
%axis([1 126 0.3 1])


xi = [1:0.2:201];
xj = [1:1:201];
x_handle = xlabel('Cell','Fontsize',24);
y_handle = ylabel('Volume fraction, \alpha_1','Fontsize',24);

set(x_handle,'Fontname','Lucida bright');
set(y_handle,'Fontname','Lucida bright');
%axis square

[xx,yy] = meshgrid([0:0.16:20]);
R = (sign( (xx-10).^2 + (yy-10).^2 - 16) +1)/2;

box on