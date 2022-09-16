function [m]=dismom(strike,dip,rake)
%
%   function dismom(strike,dip,rake)
%
%   compute the moment tensor for a point dislocation
%   
%
%  Input:   strike = strike of fault plane (clockwise from north, dip to the
%                    right)
%           dip    = fault dip from horizontal
%           rake   = fault rake (slip vector on fault showing slip of upper
%                    plane
%           all angles in degrees
%
%  Output:  M11=m(1)
%           M22=m(2)
%           M33=m(3)
%           M12=m(4)
%           M13=m(5)
%           M23=m(6)
%*************************************************************************
con=pi/180.;
s=strike*con;
d=dip*con;
r=rake*con;
%
m(1)=sin(s)*sin(s)*sin(r)*sin(2*d) + sin(2*s)*cos(r)*sin(d);
m(2)=cos(s)*cos(s)*sin(r)*sin(2*d) - sin(2*s)*cos(r)*sin(d);
m(3)=-sin(r)*sin(2*d);
m(4)=-cos(2*s)*cos(r)*sin(d) - 0.5*sin(2*s)*sin(r)*sin(2*d);
m(5)=cos(s)*cos(r)*cos(d) + sin(s)*sin(r)*cos(2*d);
m(6)=sin(s)*cos(r)*cos(d) - cos(s)*sin(r)*cos(2*d);
%
return;
