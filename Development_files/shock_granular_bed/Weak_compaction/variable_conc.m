
%importing data:

N_p =40;
for i = 1:N_p
    n = (i-1)*5;
var_name = ['p2_' num2str(n) '.dat'];

p2(:,:,i) = importdata(var_name);
end

for i = 1:N_p
    px50(i) = p2(3,50,i);
end

for i = 1:N_p
    px199(i) = p2(3,199,i);
end

tt = [-0.4e-3 :0.00005:-0.4e-3 + 0.00005*(N_p-1)];
%ttdown = [
plot(tt,px50);
hold on
plot(tt,px199);