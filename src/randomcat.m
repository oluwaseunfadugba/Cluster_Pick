function [xr,yr,zr]=randomcat(xs,ys,zs,xp,yp,rnd_in_z)
% randomcat - produce a randomized catalog of epicenters based on an input
% catalog.  xs and ys will be randomized, zs will remain the same. The data
% will be randomized within the bounded polygon specified by (xp,yp).

%rnd_in_z = 1;

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
            zr(k)=zs(k);
        end
    end   % end of while loop
    
end

if rnd_in_z == 1
    % Randomized the hypocenters in depth as well

    zmin=min(zs);
    zmax=max(zs);
    dz=zmax-zmin;

    for k=1:N

        IN_DEP_RANGE=0;  %  assume the hypocenter is out of the polygon
        while IN_DEP_RANGE == 0
            zr_test= zmin + rand(1).*dz;
            % determine if a hypocenter is within the depth range (zmin, zmax)
            IN= (zr_test >= zmin & zr_test <= zmax);
            if IN == 1
                IN_DEP_RANGE=1;
                zr(k)=zr_test;
            end
        end   % end of while loop

    end

end



% plot the random catalog and bounding box
%figure;

% hold on
% plot3(xr,yr,zr,'o','MarkerEdgeColor','k','MarkerFaceColor',[0 0.75 0.75]);
% plot3(xp,yp,zp,'-k');
% hold off
% 
% axis equal;
% title('randomized hypocenters');
% xlabel('X km');
% ylabel('Y km');
% zlabel('Z km');

end

