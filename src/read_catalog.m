function read_catalog(infile)
%  read_catalog - read the hypocenters file to analyze

global xc yc zc vec_plane xb_old yb_old zb_old xs ys zs N Nc PLOT_FLAG0

fid=fopen(infile,'r');

[data,count]=fscanf(fid,'%g %g %g',[3,inf]);

fclose(fid);

data=data';
[N,m]=size(data);

%  hypocentral locations, N=number of hypocenters
xs(1:N)=data(1:N,1);
ys(1:N)=data(1:N,2);
zs(1:N)=data(1:N,3);

% plot the input data

if PLOT_FLAG0 == 1
    figure;
    plot3(xs,ys,zs,'o');
    axis equal;
    grid on;
    title('Input hypocenters');
    xlabel('X km');
    ylabel('Y km');
    zlabel('Z km');
end

end

