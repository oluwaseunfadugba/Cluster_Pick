function combine_coplananr_faults(savevar_filename)
%clear all; close all; clc;

load(savevar_filename)

Kfaults_ini = Kfaults

% ------------------------------------------------------------------------
mincopl_lambda3 = err_av;
counter = 0;
mindist_clus = err_av;%10;%err_av;%1;

k=0; index = 1; 

for each_clus = 1:Kfaults-1 % for each cluster less 1     
    [~,sorted_Nt]= sort(Nt,'descend');

    i = sorted_Nt(index);

    [bool,coplanar] = is_consistent_cluster...
        (i,k,index,Kfaults,xt,yt,zt,Nt,sorted_Nt,mincopl_lambda3,mindist_clus,xv,yv,zv);

    if bool == 1
        counter = counter + 1;
        [Kfaults,xb,yb,zb,xv,yv,zv,xt,yt,zt,vec_plane,Nt,lambda3,...
            Strike,Dip, L,W] = combine_coplanar_clusters(Kfaults,coplanar,...
            xb,yb,zb,xv,yv,zv,xt,yt,zt,vec_plane,Nt,lambda3,Strike,...
            Dip, L,W,counter,xs,ys,zs,PLOT_FLAG1);

    else
        index = index+1;
    end

end

if PLOT_FLAG0 == 1
    %  plot final planes after combining coplanar faults
    picname= ['Final Planes (Combining Coplanar Faults: ' num2str(Kfaults_ini) ' to ' num2str(Kfaults) ')'];
    datplot(xs,ys,zs,Kfaults,xv,yv,zv,picname);
    datplot_with_colors(xs,ys,zs,Kfaults,xv,yv,zv,picname);

end

Kfaults

% Save all vaiables
!rm -r *savevar_filename_coplanar.mat
savevar_filename = [simul_tag '.saved_variables_coplanar.mat'];
save(savevar_filename)

end

%% Helper functions
function [bool,coplanar] = is_consistent_cluster...
    (i,k,index,Kfaults,xt,yt,zt,Nt,sorted_Nt,mincopl_lambda3,mindist_clus,xv,yv,zv)

    cluserss = [];
    
    % Compare a cluster to another cluster that are within a certain
    % distance from it. This will avoid combining clusters that are very
    % far from each other thereby crossing the entire seismic zone.
    
    j = sorted_Nt(index+1:Kfaults);    
    qclusts = find_close_clusters(i,j,xv,yv,zv,mindist_clus);

    if ~isempty(qclusts)
        for j=qclusts %sorted_Nt(index+1:Kfaults)
            k = k+1;

            cluster1 = [xt(i,1:Nt(i))' yt(i,1:Nt(i))' zt(i,1:Nt(i))'];
            cluster2 = [xt(j,1:Nt(j))' yt(j,1:Nt(j))' zt(j,1:Nt(j))'];

            comb_cluster = [cluster1; cluster2];

            % compute the covariance matrix for this cluster
            Cxy=cov(comb_cluster,0);

            % compute the eigenvalues and eigenvectors for this cluster
            [~,D]=eig(Cxy);

            lambda3n=sqrt(12.*D(1,1));
            cluserss(k,1:3) = [i j lambda3n];
        end
        
        if min(cluserss(:,3)) < mincopl_lambda3
            bool = 1;
            coplanar = cluserss(cluserss(:,3) == min(cluserss(:,3)),:);

        else
            bool = 0;
            coplanar = [];
        end
    
    else
        bool = 0;
        coplanar = [];
    end
    
end

function qclusts = find_close_clusters(i,j,xv,yv,zv,mindist_clus)

    qclusts = [];

    for jj = j
        
        min_dist = dmin_btw2clus(xv(i,:),yv(i,:),zv(i,:),xv(jj,:),yv(jj,:),zv(jj,:),0.03);

        if min_dist <= mindist_clus

            qclusts = [qclusts jj];
        end
    end
end

function qclusts = find_close_clusters_old(i,j,xt,yt,zt,Nt,mindist_clus)

    qclusts = [];

    for jj = j
        cluster1 = [xt(i,1:Nt(i))' yt(i,1:Nt(i))' zt(i,1:Nt(i))'];
        cluster2 = [xt(jj,1:Nt(jj))' yt(jj,1:Nt(jj))' zt(jj,1:Nt(jj))'];

        dist =[];

        for ii = 1: Nt(i)
            ref_pt = cluster1(ii,:);
            diff = cluster2 - ref_pt;

            C =  sqrt(diff(:,1).^2 + diff(:,2).^2 + diff(:,3).^2);
            dist = [dist; C];
        end

        min(dist);
        
        if min(dist) <= mindist_clus

            qclusts = [qclusts jj];
        end

    end
    
end

function [Kfaults,xb,yb,zb,xv,yv,zv,xt,yt,zt,vec_plane,Nt,lambda3,Strike,Dip, L,W]=...
    combine_coplanar_clusters(Kfaults,coplanar,xb,yb,zb,xv,yv,zv,xt,yt,zt,...
    vec_plane,Nt,lambda3,Strike,Dip,L,W,counter,xs,ys,zs,PLOT_FLAG1)

%  Save the fault models for the better fitting faults
kg=0;
for k=1:Kfaults

    if k ~= coplanar(1:2)
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

% The coplanar clusters will be added and the fault parameters recalculated.
i = coplanar(1);
j = coplanar(2);

cluster1 = [xt(i,1:Nt(i))' yt(i,1:Nt(i))' zt(i,1:Nt(i))'];
cluster2 = [xt(j,1:Nt(j))' yt(j,1:Nt(j))' zt(j,1:Nt(j))'];

comb_cluster = [cluster1; cluster2];

% compute the covariance matrix for this cluster
Cxy=cov(comb_cluster,0);

% compute the eigenvalues and eigenvectors for this cluster
[V,D]=eig(Cxy);

% calculate fault plane parameters from the eigen results
% and calculate the vertices of the fault plane
k=1;
xb(k)=mean(comb_cluster(:,1));
yb(k)=mean(comb_cluster(:,2));
zb(k)=mean(comb_cluster(:,3));

[L(k),W(k),Strike(k),Dip(k),xv(k,:),yv(k,:),zv(k,:)] = fltplane(V,D,xb(k),yb(k),zb(k));

% save the plane unit normal vector and eigenvalue
vec_plane(k,1:3)=V(1:3,1);
lambda3(k)=sqrt(12.*D(1,1));

%  Now add the new fault to the other clusters in the "old" storage
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
Strike_old(kg)=Strike(k);
Dip_old(kg)=Dip(k);

% number of events in each trial cluster
Nt_old(kg)=length(comb_cluster(:,1));

% hypocenter location in a cluster
xt_old(kg,1:Nt_old(kg))=comb_cluster(:,1)';
yt_old(kg,1:Nt_old(kg))=comb_cluster(:,2)';
zt_old(kg,1:Nt_old(kg))=comb_cluster(:,3)';

% minimum eigenvalue
lambda3_old(kg)=lambda3(k);

% Update number of faults
Kfaults = Kfaults-1;

%  Load up working arrays with all data
% Barycenters
xb=xb_old(1:Kfaults);
yb=yb_old(1:Kfaults);
zb=zb_old(1:Kfaults);

% fault plane vertices
xv=xv_old(1:Kfaults,:);
yv=yv_old(1:Kfaults,:);
zv=zv_old(1:Kfaults,:);

% eigenvector that describes each plane
vec_plane=vec_plane_old(1:Kfaults,:);

% fault plane parameters
L=L_old(1:Kfaults);
W=W_old(1:Kfaults);
Strike=Strike_old(1:Kfaults);
Dip=Dip_old(1:Kfaults);

% hypocenter location in a cluster
xt = xt_old(1:Kfaults,:);
yt = yt_old(1:Kfaults,:);
zt = zt_old(1:Kfaults,:);

% number of events in each trial cluster
Nt = Nt_old(1:Kfaults);

% minimum eigenvalue
lambda3 = lambda3_old(1:Kfaults);

if PLOT_FLAG1 == 1 
    
    %  plot final planes
    picname=['Combining Coplanar Faults ' num2str(counter)];
    datplot(xs,ys,zs,Kfaults,xv,yv,zv,picname);

    hold on;
    plot3(cluster1(:,1),cluster1(:,2),cluster1(:,3),'o','MarkerEdgeColor','b','MarkerFaceColor','b');
    hold on;
    plot3(cluster2(:,1),cluster2(:,2),cluster2(:,3),'o','MarkerEdgeColor','r','MarkerFaceColor','r');

end

end
