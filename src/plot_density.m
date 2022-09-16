function plot_density(h,psi_angle,elevation,azimuth,zfactor,view_flag)
% plot seismicity density on the GUI axes using the viewing elevation and
% azimuth from the slider scales

%  rotation angles:
%     x axis - psi_angle
%     y axis - elevation
%     x axis - azimuth

%  zfactor is the zoom factor

%  view_flag =
%              0, no rotation
%              1, rotate about x axis
%              2, rotate about y axis
%              3, rotate about z axis

global data1 xc yc zc vec_plane xb_old yb_old zb_old xs ys zs N Nc
global PLOT_Latest x_poly y_poly x_orig y_orig n_orig DI D
global data2 xd yd zd
global xclus yclus zclus Cluster_number N_cluster L W Strike Dip xv yv zv lambda3

psi=deg2rad(psi_angle);
phi=deg2rad(azimuth);
theta=deg2rad(elevation);

%  calculate the rotation matrixes
ct=cos(theta);
st=sin(theta);

cp=cos(phi);
sp=sin(phi);

cs=cos(psi);
ss=sin(psi);

%  determine what angle change is occurring
if view_flag == 0
    %  no angle change
    DC=[1.0  0.0  0.0; 0.0  1.0  0.0; 0.0 0.0 1.0];
end
if view_flag == 1
    %  rotation about x axis
    DC=[1.0  0.0  0.0; 0.0  cs  -ss; 0.0 ss cs];
end
if view_flag == 2
    %  rotation about y axis
    DC=[ct  0.0  -st; 0.0  1.0  0.0; st 0.0 ct];
end
if view_flag == 3
    %  rotation about z axis
    DC=[cp  sp  0.0; -sp  cp  0.0 ; 0.0 0.0 1.0];
end

D=DC*D;

DI=inv(D);

data=data1';
data2(1:N,1)=data(1:N,1) - x_orig;
data2(1:N,2)=data(1:N,2) - y_orig;
data2(1:N,3)=data(1:N,3);

data2=D*data2';
data2=data2';

%  hypocentral locations, N=number of hypocenters
xd(1:N)=data2(1:N,1);
yd(1:N)=data2(1:N,2);
zd(1:N)=data2(1:N,3);

hold on;

plot(h,xd(1:N),yd(1:N),'o','MarkerEdgeColor','k','MarkerFaceColor',[0 0.75 0.75]);


if PLOT_Latest == 1; plot(h,x_poly,y_poly,'-k'); end

hold off;

zoom(h,zfactor);
axis(h,'equal');
grid(h,'on');
if PLOT_Latest == 1
    cdum=strcat('rotated hypocenters ','Cluster ',num2str(Cluster_number));
else
    cdum='rotated hypocenters';
end
title(h,cdum);
xlabel(h,'X km');
ylabel(h,'Y km');
zlabel(h,'Z km');

end
