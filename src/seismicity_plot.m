function seismicity_plot( infile, outfile, polygonfile )

%  Read in and plot ETSZ seismicity

%   infile = file containing seismicity parameters
%   outfile = file containing chosen block of seismicity
%   polygonfile = file containing the bounding polygon vertices for
%                 seismicity

close all;
%****************  Input hypocentral locations ****************************
fid=fopen(infile,'r');

[data,count]=fscanf(fid,'%g %g %g %g %g %g %g %g %g %g %g',[11,inf]);
fclose(fid);

data=data';
[N,m]=size(data);

%  hypocentral locations, N=number of hypocenters
xs(1:N)=data(1:N,5);
ys(1:N)=data(1:N,6);
zs(1:N)=-data(1:N,7);

% plot the input data
h1=figure;
plot3(xs,ys,zs,'o');
axis equal;
grid on;
title('Input hypocenters');
xlabel('X km');
ylabel('Y km');
zlabel('Z km');

% plot the input data on the plane
h2=figure;
plot(xs,ys,'o');
axis equal;
grid on;
title('Input hypocenters');
xlabel('X km');
ylabel('Y km');
zlabel('Z km');

%rect=getrect(h2)
[x,y]=getline(h2,'closed');

% now plot the input data and chosen polygon on the plane
h3=figure;
plot(xs,ys,'o',x,y,'-k');
axis equal;
grid on;
title('Input hypocenters');
xlabel('X km');
ylabel('Y km');

% determine if a hypocenter is inside or outside the polygon
IN=inpolygon(xs,ys,x,y);

% subset the seismicity in this polygon
k=0;
for kk=1:N;
    
    if IN(kk) == 1;
        k=k+1;
        xn(k)=xs(kk);
        yn(k)=ys(kk);
        zn(k)=zs(kk);
    end;
    
end;

nblock=k;

% plot the chosen data and polygon alone
h4=figure;
plot(xn,yn,'o',x,y,'-k');
axis equal;
grid on;
title('Chosen Hypocenters');
xlabel('X km');
ylabel('Y km');
    
% write the block of seismicity to a file
fid=fopen(outfile,'w');
for kk=1:nblock;
    fprintf(fid,'%12.5f %12.5f %12.5f\n',[xn(kk) yn(kk) zn(kk)]);
end
fclose(fid);

% write the vertices for the bounding polygon to a file
fid=fopen(polygonfile,'w');
npoly=length(x);
for kk=1:npoly;
    fprintf(fid,'%12.5f %12.5f\n',[x(kk) y(kk)]);
end
fclose(fid);
end

