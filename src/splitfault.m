function splitfault(Kfaults)
%  splitfault - split the fault with the largest lambda3 eigenvalue into 2.

% n0 = number of fault clusters to determine
% J  = global variance is output

% Nt = number of event hypocenters in each cluster
% N = total number of hypocenters over all clusters

% xs,ys,zs = location of each hypocenter in original data
% xb,yb,zb = location of cluster barycenter
% xt,yt,zt = location of hypocenter in a cluster

global xc yc zc vec_plane xb_old yb_old zb_old xs ys zs N Nc
global xt yt zt Nt xb yb zb lambda3
global L W xv yv zv L_old W_old xv_old yv_old zv_old fscale
global xt_old yt_old zt_old vec_plane_old lambda3_old
global Strike Dip Strike_old Dip_old Nt_old

%  find the index of the cluster with the largest lambda3
[lambda3max,kthick]=max(lambda3);

%  Save the fault models for the better fitting faults
kg=0;
for k=1:Kfaults
    
    if k ~= kthick
        kg=kg+1;
        % load up arrays with good cluster parameters
        
        % Barycenters
        xb_old(kg)=xb(k);
        yb_old(kg)=yb(k);
        zb_old(kg)=zb(k);
        
        % fault plane vertices
        xv_old(kg,:)=xv(k,:);
        yv_old(kg,:)=yv(k,:);
        zv_old(kg,:)=zv(k,:);
        
        % hypocenter location in a cluster
        xt_old(kg,:)=xt(k,:);
        yt_old(kg,:)=yt(k,:);
        zt_old(kg,:)=zt(k,:);
        
        % eigenvector that describes each plane
        vec_plane_old(kg,1:3)=vec_plane(k,:);

        % number of events in each trial cluster
        Nt_old(kg)=Nt(k);

        % minimum eigenvalue
        lambda3_old(kg)=lambda3(k);

        % fault plane parameters
        L_old(kg)=L(k);
        W_old(kg)=W(k);
        Strike_old(kg)=Strike(k);
        Dip_old(kg)=Dip(k);
    end
end

% The cluster with the greatest fault thickness based on the minimum
% eigenvalue will now be split into two parts using focal-mechanism-seeded 
% planes if FM_file exist. If FM_file does not exist, the thick fault is 
% split into two random faults. This is the "kthick" cluster

% saving temp variables
xt_tmp=xt; yt_tmp=yt; zt_tmp=zt;
Nt_tmp=Nt; lambda3_tmp=lambda3;
L_tmp=L; W_tmp=W;
        
FM_seeded_planes(2,kthick);

while length(xb)<=1
    fprintf('One of the planes is not stable!\n');
    fprintf('Trying different positions\n');
    
    xt=xt_tmp; yt=yt_tmp; zt=zt_tmp;
    Nt=Nt_tmp; lambda3=lambda3_tmp;
    L=L_tmp; W=W_tmp;

    FM_seeded_planes(2,kthick); 
end 

%  Now add these new faults to the other clusters in the "old" storage
for k=1:2
    
        kg=kg+1;
        % load up old arrays with new cluster parameters
            
        % Barycenters
        xb_old(kg)=xb(k);
        yb_old(kg)=yb(k);
        zb_old(kg)=zb(k);
        
        % fault plane vertices
        xv_old(kg,:)=xv(k,:);
        yv_old(kg,:)=yv(k,:);
        zv_old(kg,:)=zv(k,:);
        
        % eigenvector that describes each plane
        vec_plane_old(kg,1:3)=vec_plane(k,:);

        % fault plane parameters
        L_old(kg)=L(k);
        W_old(kg)=W(k);
        
end

%  Load up working arrays with all data
xb=xb_old;
yb=yb_old;
zb=zb_old;

xv=xv_old;
yv=yv_old;
zv=zv_old;

vec_plane=vec_plane_old;
L=L_old;
W=W_old;

%  All done.  Start the analysis again

end
