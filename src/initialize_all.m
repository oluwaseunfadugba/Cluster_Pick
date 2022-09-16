function initialize_all
%  Initialize the global variables
global data1 xc yc zc vec_plane xb_old yb_old zb_old xs ys zs N Nc
global PLOT_Latest x_poly y_poly x_orig y_orig n_orig DI D
global data2 xd yd zd
global xclus yclus zclus Cluster_number N_cluster L W Strike Dip xv yv zv lambda3
global elevation_old azimuth_old

Nm=5000; % maximum number of events in a catalog
Nn=30;   % maximum number of clusters
Np=10;  %maximum number of polygon vertices

N=0;
Nc=0;

data1=zeros(3,Nm);
xs=zeros(Nm,1);
ys=zeros(Nm,1);
zs=zeros(Nm,1);
xc=zeros(Nm,1);
yc=zeros(Nm,1);
zc=zeros(Nm,1);
vec_plane=zeros(Nn,3);
xb_old=zeros(Nn,1);
yb_old=zeros(Nn,1);
zb_old=zeros(Nn,1);

PLOT_Latest=0;
x_poly=zeros(Np,1);
y_poly=zeros(Np,1);
x_orig=0.0;
y_orig=0.0;
n_orig=0;
DI=zeros(3,3);
D=[1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0];

data2=zeros(Nm,3);
xd=zeros(Nm,1);
yd=zeros(Nm,1);
zd=zeros(Nm,1);

xclus=zeros(Nn,Nm);
yclus=zeros(Nn,Nm);
zclus=zeros(Nn,Nm);
Cluster_number=0;
N_cluster=zeros(Nn,1);
L=zeros(Nn,1);
W=zeros(Nn,1);
Strike=zeros(Nn,1);
Dip=zeros(Nn,1);
xv=zeros(Nn,4);
yv=zeros(Nn,4);
zv=zeros(Nn,4);
lambda3=zeros(Nn,1);

elevation_old=0.0;
azimuth_old=0.0;
end
