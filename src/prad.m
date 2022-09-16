function [pamp]=prad(az,inc,vp,vs,dens,eps,m, z_displ)
%
%   function [pamp]=prad(az,inc,vp,vs,dens,eps,m)
%
%   calculate the radial pwave radiation pattern
%   given an azimuth and incidence angle
%
%   Input:   az = array of azimuths (degrees)
%            inc = array of incidence angles (degrees)
%            eps = array of ray direction values
%                = +1, downgoing ray
%                = -1, upgoing ray
%            vp = source p wave velocity
%            vs = source s wave velocity
%            dens = source density
%            m = moment tensor elements for source (e.g., computed using
%                                                   dismom);
%   Output:   pamp = array of amplitudes
%             pamp(k,j) where k=length(inc); j=length(az);
%
%**************************************************************************
%
%  zero the p_sph array
%
%vp=1000*vp; % m/s
%dens=1000*dens; %kg/m3
%
con=pi/180.;
p=sin(inc*con)/vp;              %ray parameter array
a=az*con;                       %azimuth array
eba=real(sqrt(1/(vp*vp) - p.^2));     %vertical slowness array
%
ileng=length(inc);
aleng=length(az);
%
%  compute spherical P wave source Green's functions
%
%  The resulting amplitude is in cm at 1 km distance
%  from a source with a seismic moment of 10^25 dyn cm.
%
con1=(10^5)/(4.0*pi*dens*vp);
hr0=con1*1.0/(vp*vp);                       %isotropic source
hr1=con1*p.*p;                              %vertical strike-slip
hr2=-2*con1.*p.*eba;  %.*eps                %vertical dip-slip
hr3=-con1*(p.*p -2.*eba.*eba);              %clvd
%
%   compute displacements
%
for k=1:ileng
    for j=1:aleng
        p_sph(k,j)=((m(1)+m(2)+m(3))/3)*hr0 ...
        + (0.5*(m(2)-m(1))*cos(2*a(j)) - m(4)*sin(2*a(j)))*hr1(k) ...
        + (m(5)*cos(a(j)) + m(6)*sin(a(j)))*hr2(k) ...
        + ((m(1)+m(2)-2*m(3))/6)*hr3(k);
    end
end
pamp=p_sph;

if z_displ == 1
    pamp = eba'.*p_sph; % positive upward  .*eps
end

return;