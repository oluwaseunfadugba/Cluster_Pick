function [xr,yr]=randomcat_2D(xs,ys,xp,yp)
% randomcat - produce a randomized catalog of epicenters based on an input
% catalog.  xs and ys will be randomized. The data
% will be randomized within the bounded polygon specified by (xp,yp).

%  find x and y bounds of the seismicity
N=length(xs);
N_poly=length(xp);

xmin=min(min(xs),min(xp));
xmax=max(max(xs),max(xp));
dx=xmax-xmin;

ymin=min(min(ys),min(yp));
ymax=max(max(ys),max(yp));
dy=ymax-ymin;

for k=1:N
    
    IN_POLY=0;  %  assume the hypocenter is out of the polygon
    while IN_POLY == 0
        xr_test=xmin + rand(1).*dx;
        yr_test=ymin + rand(1).*dy;

        % determine if a hypocenter is inside or outside the polygon
        IN=inpolygon(xr_test,yr_test,xp,yp);
        if IN == 1
            IN_POLY=1;
            xr(k)=xr_test;
            yr(k)=yr_test;            
        end
    end   % end of while loop 
end
end

