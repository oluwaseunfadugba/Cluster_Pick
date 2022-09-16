function generate_plane_using_FM(n0,FAULT_FLAG,strike_pr,dip_pr,kk,xs_fm,ys_fm,zs_fm)

%  Construct n0 faults using nearby focal mechanisms when we are splitting
%  the thickest fault in OADC_3D
%  just need the vertice locations

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

L22=L(FAULT_FLAG)./2.0;
W22=W(FAULT_FLAG);
    
for k=1:n0
    if L22 <= W22
        L(k)=W22; 
        W(k)=L22;
    else
        L(k)=L22;
        W(k)=W22;
    end
   
    L2=L(k)./2.0;
    W2=W(k)./2.0;
end

%
strikes_fm_now = strike_pr;
dips_fm_now = dip_pr;

% Get unit vector from strike and dip.
%  find the plane unit vector
V(1) = -cosd(strikes_fm_now)*sind(dips_fm_now);
V(2) =  sind(strikes_fm_now)*sind(dips_fm_now);
V(3) = -cosd(dips_fm_now);

% compute vertice locations
% First, for horizontal fault
xvt = [-W2 -W2 W2 W2];
yvt = [L2 -L2 -L2 L2];
zvt = [0 0 0 0];

% Fault corner coordinates have been created - Rotate the data into the
% geographical coordinate system in order of rake, dip, strike
% rotate into rake direction
R=[xvt ; yvt ;zvt];
[~,n] = size(R); rake =0;

Drake=[cosd(rake) -sind(rake) 0 ;sind(rake) cosd(rake) 0 ;0 0 1];
Rrake=Drake*R;

xvt(1:n)=Rrake(1,1:n);
yvt(1:n)=Rrake(2,1:n);
zvt(1:n)=Rrake(3,1:n);

% rotate into dip direction
Ddip=[cosd(dips_fm_now) 0 sind(dips_fm_now); 0 1  0 ;...
     -sind(dips_fm_now) 0 cosd(dips_fm_now)];
Rdip=Ddip*Rrake;

xvt(1:n)=Rdip(1,1:n);
yvt(1:n)=Rdip(2,1:n);
zvt(1:n)=Rdip(3,1:n);

% rotate into strike direction
Dstrike=[cosd(strikes_fm_now) sind(strikes_fm_now) 0 ; ...
        -sind(strikes_fm_now) cosd(strikes_fm_now) 0 ; 0 0 1];
Rstrike=Dstrike*Rdip;

xvt(1:n)=Rstrike(1,1:n);
yvt(1:n)=Rstrike(2,1:n);
zvt(1:n)=Rstrike(3,1:n);

%
for k = 1:n0
    
    % now randomly pick a hypocenter location to be the center of the fault
    %FAULT_FLAG
    Nt(FAULT_FLAG);
    nb=randperm(Nt(FAULT_FLAG)); 
    xb(k)=xt(FAULT_FLAG,nb(1));
    yb(k)=yt(FAULT_FLAG,nb(1));
    zb(k)=zt(FAULT_FLAG,nb(1));
   
    vec_plane(k,1:3)= V'; 

    xv(k,1:n)=xvt(1:n) + xb(k);
    yv(k,1:n)=yvt(1:n) + yb(k);
    zv(k,1:n)=zvt(1:n) + zb(k);
    
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