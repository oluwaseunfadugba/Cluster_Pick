function [xn,yn,zn,vol]=tetvol(x,y,z,N)
% tetvol - calculate the volume of 4 closest neighbor tetrahedra for every
% hypocenter in the catalog.  Then sort catalog in ascending order by
% tetrahedron volume

v(1:N)=0.0;
for k=1:N;
    
    dist(1:N)=0.0;
    
    x0=x(k);
    y0=y(k);
    z0=z(k);
    % calculate distance from reference source location
    dist=sqrt((x-x0).^2 + (y-y0).^2 + (z-z0).^2);
    
    % sort by ascending distance
    [~,I]=sort(dist);
    
    % check to make sure minimum distance is for the reference source
    %ichck=I(1);
    %if ichck == k;
        %  use first 4 locations to compute volume of the tetrahedron
        %  associated with this source point
       
        x0=x(I(1));
        x1=x(I(2));
        x2=x(I(3));
        x3=x(I(4));
        
        y0=y(I(1));
        y1=y(I(2));
        y2=y(I(3));
        y3=y(I(4));
        
        z0=z(I(1));
        z1=z(I(2));
        z2=z(I(3));
        z3=z(I(4));
        
        a1=(x0-x1).*(y0-y2).*(z0-z3);
        a2=(y0-y1).*(z0-z2).*(x0-x3);
        a3=(z0-z1).*(x0-x2).*(y0-y3);
        a4=-(x0-x3).*(y0-y2).*(z0-z1);
        a5=-(y0-y3).*(z0-z2).*(x0-x1);
        a6=-(z0-z3).*(x0-x2).*(y0-y1);
        
        v(k)=abs((a1+a2+a3+a4+a5+a6)./6.0);
        
end

% sort catalog by ascending tetrahedron volume
[~,II]=sort(v);

xn(1:N)=0.0;
yn(1:N)=0.0;
zn(1:N)=0.0;
vol(1:N)=0.0;

for kk=1:N
    xn(kk)=x(II(kk));
    yn(kk)=y(II(kk));
    zn(kk)=z(II(kk));
    vol(kk)=v(II(kk));
end
    
end
