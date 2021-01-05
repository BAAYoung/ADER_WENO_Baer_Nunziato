[xx,yy] = meshgrid([0:0.002:0.2],[0:0.002:0.14]);

R = (1- sign(0.038*0.038 - (xx-0.1).^2 - (yy-0.07).^2))/2;
hold on
%for i = [29 40 51 60 72]
i = 51;
    plot((yy(:,1)-0.07)/0.038, R(:,i).*umag(:,i)/mean(mean(umag)),'k -- ');

%end

x_handle = xlabel('y/D','Fontsize',24);
y_handle = ylabel('|u|/U_p','Fontsize',24);

set(x_handle,'Fontname','Lucida bright');
set(y_handle,'Fontname','Lucida bright');
axis square
axis([-2 2 -0.2 1.4]);