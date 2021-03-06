
%importing data:

N_p =52;
for i = 1:N_p
    n = (i-1)*5;
var_name = ['p2_' num2str(n) '.dat'];

p2(:,:,i) = importdata(var_name);
end


for i = 1:N_p
    px50(i) = p2(3,50,i);
end

for i = 1:N_p
    px180(i) = p2(3,180,i);
end

tt = [-0.4e-3 :0.00005:-0.4e-3 + 0.00005*(N_p-1)];
ttdown = tt+1.6e-3;
plot(tt,px50,'k');
hold on
plot(ttdown,px180,'color',[0.5 0.5 0.5]);

x_handle = xlabel('Time','Fontsize',24);
y_handle = ylabel('Gas Pressure, p_2','Fontsize',24);

set(x_handle,'Fontname','Lucida bright');
set(y_handle,'Fontname','Lucida bright');