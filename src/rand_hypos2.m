function [rxp, ryp, rzp] = rand_hypos2(L,W,zerr_av,nhypos,strike,dip,rake,rt)

%  construct set of random hypocenters on the x-y plane with an additional
%  random error in z.  The x-y plane is specified by -L/2 < x < L/2 and
%  -W/2 < y W/2.  nhypos is the number of hypocenters desired.

%  zerr_av is the standard deviation of z error
%  strike, dip, and rake are in the geographical coordinate system
%  (degrees)

%  rt is a constant translation vector of the origin.

%  outfile is the name of the output file of hypocenter locations
% L=10; W=5; zerr_av=0.5; nhypos=200;strike=45;dip=90;rake=0;rt=[0 0 -5];
% outfile = 'rand_hypo_seun.txt';

con=pi/180.;
strike=strike.*con;
dip=dip.*con;
rake=rake.*con;

n=nhypos/4;
rx=[ (randperm(n) + randn(1,n))  (randperm(n) + randn(1,n)) -(randperm(n) + randn(1,n)) -(randperm(n) + randn(1,n))] ;
ry=[ (randperm(n) + randn(1,n)) -(randperm(n) + randn(1,n)) -(randperm(n) + randn(1,n))  (randperm(n) + randn(1,n))] ;

% scale to L and W
rx=rx.*L./(2.0.*n);
ry=ry.*W./(2.0.*n);

xlen=length(rx);
ylen=length(ry);

% now construct an error in z
rz=randn(1,nhypos).*zerr_av;

% Random hypocenters have been created - Rotate the data into the
% geographical coordinate system in order of rake, dip, strike

% rotate into rake direction

R=[rx ; ry ;rz];
Drake=[ cos(rake)  -sin(rake) 0 ; sin(rake) cos(rake) 0 ; 0 0 1];

Rrake=Drake*R;

rxp(1:nhypos)=Rrake(1,1:nhypos);
ryp(1:nhypos)=Rrake(2,1:nhypos);
rzp(1:nhypos)=Rrake(3,1:nhypos);

% rotate into dip direction

Ddip=[ 1 0 0; 0 cos(dip)  -sin(dip) ; 0 sin(dip) cos(dip)];

Rdip=Ddip*Rrake;

rxp(1:nhypos)=Rdip(1,1:nhypos);
ryp(1:nhypos)=Rdip(2,1:nhypos);
rzp(1:nhypos)=Rdip(3,1:nhypos);

% rotate into strike direction

Dstrike=[ sin(strike)  -cos(strike) 0 ; cos(strike) sin(strike) 0 ; 0 0 1];

Rstrike=Dstrike*Rdip;

rxp(1:nhypos)=Rstrike(1,1:nhypos);
ryp(1:nhypos)=Rstrike(2,1:nhypos);
rzp(1:nhypos)=Rstrike(3,1:nhypos);

%  translate the final rotated collection of points
rxp=rxp+rt(1);
ryp=ryp+rt(2);
rzp=rzp+rt(3);
