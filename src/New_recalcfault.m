function [Ln,Wn,Striken,Dipn,xvn,yvn,zvn,vec_planen,sum_lambda3qn] = New_recalcfault(Nt,xt,yt,zt,Kfaults)
%  recalcfault - Compute the covariance matrix Cxy of each cluster, perform
%  principal components analysis, create best fault planes around each
%  cluster barycenter.

% n0 = number of fault planes
% Nt = number of event hypocenters in each cluster
% N = total number of hypocenters over all clusters

% xs,ys,zs = location of each hypocenter in original data
% xb,yb,zb = location of cluster barycenter
% xt,yt,zt = location of hypocenter in a cluster

%  Analyze each cluster
for k=1:Kfaults
    
    % compute the covariance matrix for this cluster
    Cxy=cov( [xt(k,1:Nt(k))' yt(k,1:Nt(k))' zt(k,1:Nt(k))'],0);
    
    % Seun checks if Cxy contains NaN. This happens if no hypocenter is close to
    % one of the splitted faults.
    NrNaN = sum(isnan(Cxy(:)));
  
    if (NrNaN > 0) || (length(Cxy)==1); continue; end
    
    % compute the eigenvalues and eigenvectors for this cluster
    [V,D]=eig(Cxy);
            
    % calculate fault plane parameters from the eigen results
    % and calculate the vertices of the fault plane
    X = [xt(k,1:Nt(k))' yt(k,1:Nt(k))' zt(k,1:Nt(k))'];
    
    xb = mean(X(:,1)); yb = mean(X(:,2)); zb = mean(X(:,3));
    
    % [Ln(k),Wn(k),Striken(k),Dipn(k),xvn(k,:),yvn(k,:),zvn(k,:)] = fltplane(X,V,D,xb,yb,zb);
    [Ln(k),Wn(k),Striken(k),Dipn(k),xvn(k,:),yvn(k,:),zvn(k,:)] = fltplane(V,D,xb,yb,zb);
            
    % save the plane unit normal vector and eigenvalue
    vec_planen(k,1:3)=V(1:3,1);
    lambda3qn(k)=sqrt(12.*D(1,1));
            
end

sum_lambda3qn = sum(lambda3qn);

end
