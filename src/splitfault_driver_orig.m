function splitfault_driver_orig(Kfaults)

% OADC_3D will split the thickest fault N_loop times to find the
% configuration with the best fit.           
% load up arrays with good cluster parameters   

global N_loop dip_threshold
%N_loop = 1; dip_threshold = 10;

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

i = 0; ii = 0; skip_model = 0; N_max = 30; clus_mineqs = 4;

while (ii  <= N_max) && (i< N_loop) 
    ii = ii + 1;
    
    % Put back the temporary arrays into the original arrays
    % load up arrays with good cluster parameters
    xb = xb_tmp; yb = yb_tmp; zb = zb_tmp; % Barycenters
    xv=xv_tmp; yv=yv_tmp; zv=zv_tmp; % fault plane vertices
    xt=xt_tmp; yt=yt_tmp; zt=zt_tmp; % hypocenter location in a cluster
    vec_plane=vec_plane_tmp; % eigenvector that describes each plane
    Nt=Nt_tmp; % number of events in each trial cluster
    lambda3=lambda3_tmp; % minimum eigenvalue
    L=L_tmp; W=W_tmp; Strike=Strike_tmp; Dip=Dip_tmp; % fault plane parameters

    splitfault(Kfaults)

    % Store each parameter to know which to advance to the next
    % iteration. load up arrays with good cluster parameters
    xb_tmp_ii=xb; yb_tmp_ii=yb; zb_tmp_ii=zb; % Barycenters
    xv_tmp_ii=xv; yv_tmp_ii=yv; zv_tmp_ii=zv; % fault plane vertices
    vec_plane_tmp_ii=vec_plane; % eigenvector that describes each plane
    Nt_tmp_ii=Nt; % number of events in each trial cluster
    lambda3_tmp_ii=lambda3; % minimum   
    L_tmp_ii=L; W_tmp_ii=W; Strike_tmp_ii=Strike;
    Dip_tmp_ii=Dip; % fault plane parameters

    % increase the fault number
    Kfaults=Kfaults+1;

    JFINALL=faultcluster(con_tol,Kfaults);

    % Reduce the number of fault because faultcluster.m have increased it by 1.
    Kfaults=Kfaults-1;          

    minDip = min(Dip(1:(Kfaults+1)));
    minNt = min(Nt(1:(Kfaults+1)));

    minDip; minNt;

    if (ii  == N_max)
        fprintf('Terminating model: Subhorizontal Planes can not be avoided!\n');
        Kfaults = kmax+1;
        skip_model = 1;
        break;
    end

    if (minDip >= dip_threshold) && (i < N_loop) && (minNt >= clus_mineqs) 
        % minDip >= 0 is normal as if no minDip constraint and no while loop
        i = i+1;

        JFINAL(i) = JFINALL;

        % Store each parameter to know which to advance to the next
        % iteration. load up arrays with good cluster parameters
        xb_tmp_i(:,:,i)=xb_tmp_ii; yb_tmp_i(:,:,i)=yb_tmp_ii; 
        zb_tmp_i(:,:,i)=zb_tmp_ii; % Barycenters
        xv_tmp_i(:,:,i)=xv_tmp_ii; yv_tmp_i(:,:,i)=yv_tmp_ii; 
        zv_tmp_i(:,:,i)=zv_tmp_ii; % fault plane vertices
        vec_plane_tmp_i(:,:,i)=vec_plane_tmp_ii; % eigenvector that describes each plane
        Nt_tmp_i(:,:,i)=Nt_tmp_ii; % number of events in each trial cluster
        lambda3_tmp_i(:,:,i)=lambda3_tmp_ii; % minimum eigenvalue
        L_tmp_i(:,:,i)=L_tmp_ii; W_tmp_i(:,:,i)=W_tmp_ii; Strike_tmp_i(:,:,i)=Strike_tmp_ii;
        Dip_tmp_i(:,:,i)=Dip_tmp_ii; % fault plane parameters

    end                    
end

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

end