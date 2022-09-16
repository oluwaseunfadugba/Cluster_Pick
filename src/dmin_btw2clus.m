function dmin = dmin_btw2clus(xv1,yv1,zv1,xv2,yv2,zv2,resolution)

    % The function determines the minimum distance between two finite planes 
    % defined by the corners. e.g., xv1 has 1x4 for the x-cord of plane 1

    N = round(1/resolution); N2 = N*N;

    pts = get_pts_on_plane(xv1,yv1,zv1,N);
    pts2 = get_pts_on_plane(xv2,yv2,zv2,N);

    C = zeros(1,N2);
    for i=1:N2
        
        dd = pts2 - pts(i,:);
        C(i) = min(sqrt(dd(:,1).^2 + dd(:,2).^2 + dd(:,3).^2));
        
    end
    
    dmin = min(C);    
 
end



function pts = get_pts_on_plane(x,y,z,N)

    x_now1 = linspace(x(1), x(2), N);
    y_now1 = linspace(y(1), y(2), N);
    z_now1 = linspace(z(1), z(2), N);

    x_now2 = linspace(x(4), x(3), N);
    y_now2 = linspace(y(4), y(3), N);
    z_now2 = linspace(z(4), z(3), N);

    xx_now = []; yy_now = []; zz_now = []; 
    
    for i=1:N

        xx_now = [xx_now linspace(x_now1(i), x_now2(i), N)];
        yy_now = [yy_now linspace(y_now1(i), y_now2(i), N)];
        zz_now = [zz_now linspace(z_now1(i), z_now2(i), N)];

    end
    
    pts = [xx_now' yy_now' zz_now'];
    
end
