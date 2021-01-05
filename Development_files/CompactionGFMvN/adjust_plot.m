%close all

slic = 63;

u_piston = 0.013318;

save_t_inc = 5;
dx = 0.16;

Lx = 40;
Ly = 20;
xx = [0:0.16:40];
yy = [0:0.16:20];
t = 30;
t_end = 30;
y = yy(slic);
x_pos =(t_end- t)*save_t_inc*u_piston;

r = sqrt(16-(y-10)^2);
x_1 = 30+r-x_pos;
x_2 = 30-r -x_pos;
xx_gfm = xx-x_pos;

plot(xx,u210(slic,:)+u_piston,'color',[0.6 0.6 0.6]);
hold on
%lot(xx_gfm,g_alpha130(slic,:),'k .');
%plot([x_1 x_1+1e-5],[0.4 1],'k','LineWidth',4);
%plot([x_2 x_2+1e-5],[0.4 1],'k','LineWidth',4');
axis([15 40 0 u_piston])

x_handle = xlabel('Distance','Fontsize',24);
y_handle = ylabel('Horizontal velocity, u','Fontsize',24);

set(x_handle,'Fontname','Lucida bright');
set(y_handle,'Fontname','Lucida bright');

