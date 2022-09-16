function splitfault_driver(Kfaults)

% OADC_3D will split the thickest fault N_loop times to find the
% configuration with the best fit.           
% load up arrays with good cluster parameters   

global N_loop dip_threshold skip_model

global xb yb zb; % Barycenters
global xv yv zv; % fault plane vertices
global xt yt zt; % hypocenter location in a cluster
global vec_plane; % eigenvector that describes each plane
global Nt; % number of events in each trial cluster
global lambda3; % minimum eigenvalue
global L W Strike Dip con_tol clus_mineqs

global xb_tmp_i yb_tmp_i zb_tmp_i
global xv_tmp_i yv_tmp_i zv_tmp_i
global xt_tmp_i yt_tmp_i zt_tmp_i
global vec_plane_tmp_i Nt_tmp_i lambda3_tmp_i
global L_tmp_i W_tmp_i Strike_tmp_i Dip_tmp_i

xb_tmp=xb; yb_tmp=yb; zb_tmp=zb; % Barycenters
xv_tmp=xv; yv_tmp=yv; zv_tmp=zv; % fault plane vertices
xt_tmp=xt; yt_tmp=yt; zt_tmp=zt; % hypocenter location in a cluster
vec_plane_tmp=vec_plane; % eigenvector that describes each plane
Nt_tmp=Nt; % number of events in each trial cluster
lambda3_tmp=lambda3; % minimum eigenvalue
L_tmp=L; W_tmp=W; Strike_tmp=Strike; Dip_tmp=Dip; % fault plane parameters

i = 0; ii = 0; skip_model = 0; N_max = 30;
JFINAL_Nloop = zeros(1,N_loop);

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

    % split the thickest fault into two new random fault planes
    splitfault(Kfaults)

    % increase the fault number
    Kfaults=Kfaults+1;

    JFINALL=faultcluster(con_tol,Kfaults);
    
    % Reduce the number of fault bcos faultcluster.m have increased it by 1.
    Kfaults=Kfaults-1;          

    minDip = min(Dip(1:(Kfaults+1)));
    minNt = min(Nt(1:(Kfaults+1)));

    if (ii  == N_max)
        fprintf('Terminating model: Subhorizontal Planes can not be avoided!\n');
        %Kfaults = kmax+1;
        skip_model = skip_model+1;
        break;
    end

    if (minDip >= dip_threshold) && (i < N_loop) && (minNt >= clus_mineqs) 
        % minDip >= 0 is normal as if no minDip constraint and no while loop
        i = i+1

        JFINAL_Nloop(i) = JFINALL;
        
        xb_tmp_i
        xb
        
        % Store each parameter to know which to advance to the next
        % iteration. load up arrays with good cluster parameters
        xb_tmp_i(:,:,i)=xb; yb_tmp_i(:,:,i)=yb; zb_tmp_i(:,:,i)=zb; % Barycenters
        xv_tmp_i(:,:,i)=xv; yv_tmp_i(:,:,i)=yv; zv_tmp_i(:,:,i)=zv; % fault plane vertices
        vec_plane_tmp_i(:,:,i)=vec_plane; % eigenvector that describes each plane
        Nt_tmp_i(:,:,i)=Nt; % number of events in each trial cluster
        lambda3_tmp_i(:,:,i)=lambda3; % minimum eigenvalue
        L_tmp_i(:,:,i)=L; W_tmp_i(:,:,i)=W; 
        Strike_tmp_i(:,:,i)=Strike; Dip_tmp_i(:,:,i)=Dip; % fault plane parameters

    end                    
end

% Determine the iteration with the minimum value of JFINAL, and
% store its corresponding arrays
[~,index] = sort(JFINAL_Nloop);

%  print the good news
fprintf('Global variance of the fit of all configs: ');
fprintf('[');fprintf('%8.4f', JFINAL_Nloop);fprintf(']')
fprintf('\n')
fprintf('Best configuration found! Config %3s\n',num2str(index(1)));

% load up arrays with good cluster parameters    
xb = xb_tmp_i(:,:,index(1)); yb = yb_tmp_i(:,:,index(1)); zb=zb_tmp_i(:,:,index(1)); % Barycenters
xv = xv_tmp_i(:,:,index(1)); yv=yv_tmp_i(:,:,index(1)); zv=zv_tmp_i(:,:,index(1)); % fault plane vertices
xt=xt_tmp_i(:,:,index(1)); yt=yt_tmp_i(:,:,index(1)); zt=zt_tmp_i(:,:,index(1)); % hypocenter location in a cluster
vec_plane=vec_plane_tmp_i(:,:,index(1)); % eigenvector that describes each plane
Nt=Nt_tmp_i(:,:,index(1)); % number of events in each trial cluster
lambda3=lambda3_tmp_i(:,:,index(1)); % minimum eigenvalue
L=L_tmp_i(:,:,index(1)); W=W_tmp_i(:,:,index(1)); 
Strike=Strike_tmp_i(:,:,index(1)); Dip=Dip_tmp_i(:,:,index(1)); % fault plane parameters

end