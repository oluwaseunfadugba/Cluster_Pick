function [Nt,xt,yt,zt]=new_pcluster(xs,ys,zs,xv,yv,zv,L,W,vec_plane,Kfaults)
%  pcluster - form n0 clusters of seismicity around the input fault planes
%  using the modified k-means method.
%  Also compute the global variance of the fit.

% n0 = number of fault clusters to determine
% J  = global variance is output

% Nt = number of event hypocenters in each cluster
% N = total number of hypocenters over all clusters

% xs,ys,zs = location of each hypocenter in original data
% xb,yb,zb = location of cluster barycenter
% xt,yt,zt = location of hypocenter in a cluster

% Now go through the entire catalog of hypocenters to determine new
% clusters.
    
% J=0.0;
Nt(1:Kfaults)=0;
N = length(xs);

xt = zeros(Kfaults,N);
yt = zeros(Kfaults,N);
zt = zeros(Kfaults,N);

for k=1:N   %per hypocenter
    for m=1:Kfaults  %per fault plane
        xsn = xs(k); ysn = ys(k); zsn = zs(k);
        xvn = xv(m,1:4); yvn = yv(m,1:4); zvn = zv(m,1:4); 
        
        dst(m)=New_rectdist(m,xsn,ysn,zsn,xvn,yvn,zvn,L,W,vec_plane);         
    end
 
    %  find the closest fault plane
    [~,index]=min(dst);
    
    Nt(index)=Nt(index) + 1;
    xt(index,Nt(index))=xs(k);
    yt(index,Nt(index))=ys(k);
    zt(index,Nt(index))=zs(k);
        
end

% figure
% plot3(xt(1,:),yt(1,:),zt(1,:),'o','MarkerEdgeColor','r','MarkerFaceColor','r'); hold on;
% plot3(xt(2,:),yt(2,:),zt(2,:),'o','MarkerEdgeColor','b','MarkerFaceColor','b'); hold on;
% fill3(xv(1,1:4),yv(1,1:4),zv(1,1:4),'w','FaceAlpha',0.7,'FaceColor','r'); hold on;
% fill3(xv(2,1:4),yv(2,1:4),zv(2,1:4),'w','FaceAlpha',0.7,'FaceColor','b'); hold on;

return;
 