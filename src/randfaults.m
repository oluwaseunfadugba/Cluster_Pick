function randfaults(n0,FAULT_FLAG)

%  Construct n0 randomn faults
%  just need the vertice locations

%  FAULT_FLAG = 0,  Use all hypocenters.  Initialization of the Process
%  FAULT_FLAG = kthick (from splitfault.m),  Split the thickest fault into
%               two random planes of length L/2

% Nt = number of event hypocenters in each cluster
% N = total number of hypocenters over all clusters

% xs,ys,zs = location of each hypocenter in original data
% xb,yb,zb = location of cluster barycenter
% xt,yt,zt = location of hypocenter in a cluster

global xc yc zc vec_plane xb_old yb_old zb_old xs ys zs N Nc
global L W xv yv zv L_old W_old xv_old yv_old zv_old fscale
global xt yt zt Nt xb yb zb lambda3

%fprintf('From randfaults\n');

for k=1:n0
    
    % generate a random fault using a random square 3x3 matrix.  The
    % smallest eigenvalue will be the fault normal.
    A=rand(3);
    B=A'*A;
    [V,~]=eig(B);
    
    % choose 2 random numbers between 0 and 1
    L1=rand(1,2);
    % choose a permutation of 1-10
    a=randperm(10);
    
    if FAULT_FLAG == 0
        L2=L1(1).*fscale;
        W2=L1(2).*fscale;
    else
        L2=L(FAULT_FLAG)./2.0;
        W2=W(FAULT_FLAG);
    end
    
    %  make L2 >= W2
    if L2 <= W2
        L(k)=W2;
        W(k)=L2;
    else
        L(k)=L2;
        W(k)=W2;
    end
    
    %  find the plane unit vector
    vec_plane(k,1:3)=V(1:3,1);
    
    % now randomly pick a hypocenter location to be the center of the fault
    if FAULT_FLAG == 0
        nb=randperm(N);
        xb(k)=xs(nb(1));
        yb(k)=ys(nb(1));
        zb(k)=zs(nb(1));
    else
        nb=randperm(Nt(FAULT_FLAG));
        xb(k)=xt(FAULT_FLAG,nb(1));
        yb(k)=yt(FAULT_FLAG,nb(1));
        zb(k)=zt(FAULT_FLAG,nb(1));
    end

    
    L2=L(k)./2.0;
    W2=W(k)./2.0;

    % compute vertice locations
    xv(k,1)=W2.*V(1,2)+L2.*V(1,3) + xb(k);
    yv(k,1)=W2.*V(2,2)+L2.*V(2,3) + yb(k);
    zv(k,1)=W2.*V(3,2)+L2.*V(3,3) + zb(k);

    xv(k,2)=W2.*V(1,2)-L2.*V(1,3) + xb(k);
    yv(k,2)=W2.*V(2,2)-L2.*V(2,3) + yb(k);
    zv(k,2)=W2.*V(3,2)-L2.*V(3,3) + zb(k);

    xv(k,3)=-W2.*V(1,2)-L2.*V(1,3) + xb(k);
    yv(k,3)=-W2.*V(2,2)-L2.*V(2,3) + yb(k);
    zv(k,3)=-W2.*V(3,2)-L2.*V(3,3) + zb(k);

    xv(k,4)=-W2.*V(1,2)+L2.*V(1,3) + xb(k);
    yv(k,4)=-W2.*V(2,2)+L2.*V(2,3) + yb(k);
    zv(k,4)=-W2.*V(3,2)+L2.*V(3,3) + zb(k);
    
end

% Force the random faults to be within the thick cluster
% We will locally partition the earthquakes within the "thick" cluster and
% determine the associated faults. These associated faults will be the
% desired "random" faults within the "thick" cluster.

if FAULT_FLAG ~= 0
    
    xs_thick = xt(FAULT_FLAG,1:Nt(FAULT_FLAG));
    ys_thick = yt(FAULT_FLAG,1:Nt(FAULT_FLAG));
    zs_thick = zt(FAULT_FLAG,1:Nt(FAULT_FLAG));

    [Nt_n,xt_n,yt_n,zt_n] = new_pcluster...
        (xs_thick,ys_thick,zs_thick,xv,yv,zv,L,W,vec_plane,n0);    
    [L,W,~,~,xv_n,yv_n,zv_n,vec_plane_n,~] = New_recalcfault...
        (Nt_n,xt_n,yt_n,zt_n,n0);
    
    xv=xv_n; yv=yv_n; zv=zv_n;
    
    vec_plane = vec_plane_n; 
    xb = mean(xv_n,2)'; yb = mean(yv_n,2)'; zb = mean(zv_n,2)';

end

return;

%
%     for k=1:n0
%         L2=L(k)./2.0;
%         W2=W(k)./2.0;
% 
%         V(1:3,1) = vec_plane(k,1:3);
% 
%         % compute vertice locations
%         xv(k,1)=W2.*V(1,2)+L2.*V(1,3) + xb(k);
%         yv(k,1)=W2.*V(2,2)+L2.*V(2,3) + yb(k);
%         zv(k,1)=W2.*V(3,2)+L2.*V(3,3) + zb(k);
% 
%         xv(k,2)=W2.*V(1,2)-L2.*V(1,3) + xb(k);
%         yv(k,2)=W2.*V(2,2)-L2.*V(2,3) + yb(k);
%         zv(k,2)=W2.*V(3,2)-L2.*V(3,3) + zb(k);
% 
%         xv(k,3)=-W2.*V(1,2)-L2.*V(1,3) + xb(k);
%         yv(k,3)=-W2.*V(2,2)-L2.*V(2,3) + yb(k);
%         zv(k,3)=-W2.*V(3,2)-L2.*V(3,3) + zb(k);
% 
%         xv(k,4)=-W2.*V(1,2)+L2.*V(1,3) + xb(k);
%         yv(k,4)=-W2.*V(2,2)+L2.*V(2,3) + yb(k);
%         zv(k,4)=-W2.*V(3,2)+L2.*V(3,3) + zb(k);
% 
%     end