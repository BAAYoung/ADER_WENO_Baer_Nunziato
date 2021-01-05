clear all


[xx,yy] = meshgrid([0:0.1:1]);
Vert_num = xx*0;

X = xx(:);
Y = yy(:);

Y = Y+2*X.*(X-1).*(Y-1);

for i = 1:11
    for j = 1:11
        Vert_num(i,j) = (i-1)*11 + j;
    end
end
Count = 1;
for i = 1:10
    for j = 1:10
        Vert_sq((i-1)*10 + j,:) = [Vert_num(i,j) Vert_num(i,j+1) Vert_num(i+1,j+1) Vert_num(i+1,j)];
    end
end


for i = 1:100;
    Vert_tri((i-1)*2 + 1:(i-1)*2+2,:) = [Vert_sq(i,1:3);Vert_sq(i,3:4) Vert_sq(i,1)];
end


%Calculating normals:

for i = 1:200
    Nx(i,:) = [Y(Vert_tri(i,2))-Y(Vert_tri(i,1)) Y(Vert_tri(i,3))-Y(Vert_tri(i,2)) Y(Vert_tri(i,1))-Y(Vert_tri(i,3))];
    Ny(i,:) = [X(Vert_tri(i,1))-X(Vert_tri(i,2)) X(Vert_tri(i,2))-X(Vert_tri(i,3)) X(Vert_tri(i,3))-X(Vert_tri(i,1))];
    CMx(i,:) = [X(Vert_tri(i,1))+X(Vert_tri(i,2)) X(Vert_tri(i,2))+X(Vert_tri(i,3)) X(Vert_tri(i,3))+X(Vert_tri(i,1))]/2;
    CMy(i,:) = [Y(Vert_tri(i,1))+Y(Vert_tri(i,2)) Y(Vert_tri(i,2))+Y(Vert_tri(i,3)) Y(Vert_tri(i,3))+Y(Vert_tri(i,1))]/2;
end
N_mag = sqrt(Nx.^2 + Ny.^2);
Nx = -Nx./N_mag;
Ny = -Ny./N_mag;


hold on
for i = 1:200
    fill(X(Vert_tri(i,:)),Y(Vert_tri(i,:)),'r');
    
end
for i = 1:200
    quiver(CMx(i,:),CMy(i,:),Nx(i,:)*0.1,Ny(i,:)*0.1,'k');
end
axis equal
