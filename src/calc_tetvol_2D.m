function [xn, yn, Area]=calc_tetvol_2D(x,y)
% tetvol - calculate the area of 3 closest neighbor tetrahedra for every
% hypocenter in the catalog.  Then sort catalog in ascending order by
% tetrahedron volume

%x = xs; y= ys; z= zs;

N = length(x);

Ao(1:N)=0.0; 
for k=1:N
    
    x0=x(k);
    y0=y(k);
    
    % calculate distance from reference source location
    dist=sqrt((x-x0).^2 + (y-y0).^2);
    
    % sort by ascending distance
    [~,I]=sort(dist);
    
    %  use first 3 locations to compute area of the tetrahedron
    %  associated with this source point
    x0=x(I(1)); x1=x(I(2)); x2=x(I(3)); 
    y0=y(I(1)); y1=y(I(2)); y2=y(I(3)); 
        
    matr = [x0 y0 1; x1 y1 1; x2 y2 1];
    Ao(k) = abs(det(matr))./2;
    
end

% sort catalog by ascending tetrahedron area
[~,II]=sort(Ao);

xn(1:N)=0.0;
yn(1:N)=0.0;
Area(1:N)=0.0;

for kk=1:N
    xn(kk)=x(II(kk));
    yn(kk)=y(II(kk));
    Area(kk)=Ao(II(kk));
end
    
end