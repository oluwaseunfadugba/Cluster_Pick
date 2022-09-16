% OADC_3D will split the thickest fault N_loop times to find the
% configuration with the best fit.           
% load up arrays with good cluster parameters   

N_loop = 1; dip_threshold = 10;

xb_tmp=xb; yb_tmp=yb; zb_tmp=zb; % Barycenters
xv_tmp=xv; yv_tmp=yv; zv_tmp=zv; % fault plane vertices
xt_tmp=xt; yt_tmp=yt; zt_tmp=zt; % hypocenter location in a cluster
vec_plane_tmp=vec_plane; % eigenvector that describes each plane
Nt_tmp=Nt; % number of events in each trial cluster
lambda3_tmp=lambda3; % minimum eigenvalue
L_tmp=L; W_tmp=W; Strike_tmp=Strike; Dip_tmp=Dip; % fault plane parameters

% Initialize loop temporary arrays
[m,n] = size(xb); xb_tmp_i = zeros(m,n,N_loop); 
yb_tmp_i = zeros(m,n,N_loop); zb_tmp_i = zeros(m,n,N_loop);
[m,n] = size(xv); xv_tmp_i = zeros(m,n,N_loop); 
    yv_tmp_i = zeros(m,n,N_loop); zv_tmp_i = zeros(m,n,N_loop);
[~,n] = size(xt); xt_tmp_i = zeros(Kfaults+1,n,N_loop); yt_tmp_i = ...
    zeros(Kfaults+1,n,N_loop); zt_tmp_i = zeros(Kfaults+1,n,N_loop);
[m,n] = size(vec_plane); vec_plane_tmp_i = zeros(m,n,N_loop);
[m,n] = size(Nt); Nt_tmp_i = zeros(m,n,N_loop);
[m,n] = size(lambda3); lambda3_tmp_i = zeros(m,n,N_loop);
[m,n] = size(L); L_tmp_i = zeros(m,n,N_loop);
[m,n] = size(W); W_tmp_i = zeros(m,n,N_loop);
[m,n] = size(Strike); Strike_tmp_i = zeros(m,n,N_loop);
[m,n] = size(Dip); Dip_tmp_i = zeros(m,n,N_loop);

%textprogressbar('Determining the best random fault configurations: ');

splitfault_Nloop

%textprogressbar('done');

JFINAL;

% Determine the iteration with the minimum value of JFINAL, and
% store its corresponding arrays
[~,index] = sort(JFINAL);

%  print the good news
fprintf('Lambda3 of all configs: ');
fprintf('[');fprintf('%8.4f', JFINAL);fprintf(']')
fprintf('\n')
fprintf('Best configuration found!\n');

% load up arrays with good cluster parameters    
xb = xb_tmp_i(:,:,index(1)); yb = yb_tmp_i(:,:,index(1)); 
zb=zb_tmp_i(:,:,index(1)); % Barycenters
xv = xv_tmp_i(:,:,index(1)); yv=yv_tmp_i(:,:,index(1)); 
zv=zv_tmp_i(:,:,index(1)); % fault plane vertices
xt=xt_tmp_i(:,:,index(1)); yt=yt_tmp_i(:,:,index(1)); 
zt=zt_tmp_i(:,:,index(1)); % hypocenter location in a cluster
vec_plane=vec_plane_tmp_i(:,:,index(1)); % eigenvector that describes each plane
Nt=Nt_tmp_i(:,:,index(1)); % number of events in each trial cluster
lambda3=lambda3_tmp_i(:,:,index(1)); % minimum eigenvalue
L=L_tmp_i(:,:,index(1)); W=W_tmp_i(:,:,index(1)); 
Strike=Strike_tmp_i(:,:,index(1)); Dip=Dip_tmp_i(:,:,index(1)); % fault plane parameters
