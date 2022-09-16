function dmin=New_rectdist(m,xsn,ysn,zsn,xvn,yvn,zvn,L,W,vec_plane)

% global xc yc zc vec_plane xb_old yb_old zb_old xs ys zs N Nc
% global xt yt zt Nt xb yb zb lambda3
% global L W xv yv zv Lold Wold xv_old yv_old zv_old fscale


%global L W vec_plane

%fprintf('From rectdist %g %g\n',[k m]);

% k = index of hypocenter wanted
% m = index of fault plane wanted

% Choose hypocenter vector
u=[xsn ysn zsn;xsn ysn zsn;xsn ysn zsn;xsn ysn zsn];
u1=[xsn ysn zsn];

% Construct vertices vector
v=[xvn' yvn' zvn'];

% barycenter
xb = mean(xvn)'; yb = mean(yvn)'; zb = mean(zvn)';

% Construct distance vector from each vertex
d=u-v;
d1=[d(1,1) d(1,2) d(1,3)];
d2=[d(2,1) d(2,2) d(2,3)];
d3=[d(3,1) d(3,2) d(3,3)];
d4=[d(4,1) d(4,2) d(4,3)];

d1n=norm(d1);
d2n=norm(d2);
d3n=norm(d3);
d4n=norm(d4);

d12=dot(d1,d2);
d23=dot(d2,d3);
d34=dot(d3,d4);
d41=dot(d4,d1);

L1=L(m);
W1=W(m);

% Find all angles and perpendicular distances to edge
% 12 Face
%gam12=acos(d12./(d1n.*d2n))
[alph12,beta12,gam12]=trisol(d2n,d1n,L1,'r');
a12=d1n.*d2n.*sin(gam12)./L1;

% 23 Face
%gam23=acos(d23./(d2n.*d3n))
[alph23,beta23,gam23]=trisol(d3n,d2n,W1,'r');
a23=d2n.*d3n.*sin(gam23)./W1;

% 34 Face
%gam34=acos(d34./(d3n.*d4n))
[alph34,beta34,gam34]=trisol(d4n,d3n,L1,'r');
a34=d3n.*d4n.*sin(gam34)./L1;

% 41 Face
%gam41=acos(d41./(d4n.*d1n))
[alph41,beta41,gam41]=trisol(d1n,d4n,W1,'r');
a41=d4n.*d1n.*sin(gam41)./W1;

% Find the perpendicular distance to the infinite plane
vec(1:3)=vec_plane(m,1:3);

u2=u1-[xb yb zb];
dperp=abs(dot(u2,vec));

% Find the position field for the hypocenter relative to the plane and find
% the minimum distance.
p2=pi/2.0;
dmin=dperp;
if (beta41 > p2) && (alph23 > p2); dmin=a12; end
if (beta12 > p2) && (alph23 > p2); dmin=d2n; end
if (beta12 > p2) && (alph34 > p2); dmin=a23; end
if (beta23 > p2) && (alph34 > p2); dmin=d3n; end
if (beta23 > p2) && (alph41 > p2); dmin=a34; end
if (beta34 > p2) && (alph41 > p2); dmin=d4n; end
if (beta34 > p2) && (alph12 > p2); dmin=a41; end
if (beta41 > p2) && (alph12 > p2); dmin=d1n; end

%dmin;
return;
