function seis_convert

%  Convert lat ong depth yr mo day hr min msec seismicity file from Chris
%  into x (E), y (N), z (up) relative to some reference lat and lon

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Parameters
%  Use NVT Array location - station 10

latc=36.3137;      % reference latitude
longc=-89.5328;    % reference longitude

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;

%  Read in seismicity file

fid=fopen('src10_latlon');
[data,count]=fscanf(fid,'%g  %g  %g %g  %g  %g %g  %g  %g',[9,inf]);
fclose(fid);
data=data';
[m,n]=size(data)

long(1:m)=data(1:m,2);
lat(1:m)=data(1:m,1);
z(1:m)=-data(1:m,3);

% plot the data
figure;
plot(long,lat,'o');
axis equal;
grid on;
title('Hypocenters');
xlabel('Longitude, degrees');
ylabel('Latitude, degrees');


for k=1:m; fprintf('%g %g %g\n',lat(k),long(k),z(k));end;

%  produce an equal distance azimuth projection of the data from the
%  reference pole

deg2rad=pi/180;
r_earth=6371.0;

for j=1:m;
        
    [dist(j),az(j)]=distance(latc,longc,lat(j),long(j));
    x=r_earth.*deg2rad.*dist.*sin(az.*deg2rad);
    y=r_earth.*deg2rad.*dist.*cos(az*deg2rad);
    
end;

for k=1:m; fprintf('%g %g %g\n',dist(k),az(k),z(k));end;

% plot the converted data on a 3D plot
figure;
plot3(x,y,z,'o');
axis equal;
grid on;
title('Hypocenters');
xlabel('X km');
ylabel('Y km');
zlabel('Z km');

% plot the converted data on a 2D plot
figure;
plot(x,y,'o');
axis equal;
grid on;
title('Hypocenters Relative to NVT Array');
xlabel('X km');
ylabel('Y km');
    
% write the block of seismicity to a file

fid=fopen('NM_powell_src10.txt','w');

for kk=1:m;
    fprintf(fid,'%12.5f %12.5f %12.5f\n',[x(kk) y(kk) z(kk)]);
end

fclose(fid);
