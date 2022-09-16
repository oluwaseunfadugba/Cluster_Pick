function cluster_sort( catfilein, polyfilein, catfileout_clustered, catfileout_diffuse )
% cluster_sort: separate diffuse from clustered seismicity in a catalog
% based on an algorithm suggested by Ouillon and Sornette (2001),
% Segmentation of fault networks determined from spatial clustering of
% earthquakes, JGR, 116, doi:10.1029/2010JB007752.

%  Input:
%        catfilein - catalog file subset from seismicity_plot.m
%        polyfilein - bounding polygon file from seismicity_plot.m

%  Ouput:
%        catfileout_clustered - clustered seismicity
%        catfileout_diffuse - diffuse seismicity

close all;
%**************** Assorted Parameters *************************************
PROB=0.05;

%****************  Input hypocentral locations ****************************
fid=fopen(catfilein,'r');

[data,count]=fscanf(fid,'%g %g %g',[3,inf]);

fclose(fid);

data=data';
[N,m]=size(data);

%  hypocentral locations, N=number of hypocenters
xs(1:N)=data(1:N,1);
ys(1:N)=data(1:N,2);
zs(1:N)=data(1:N,3);

%****************  Input polygon vertices *********************************
fid2=fopen(polyfilein,'r');

[data1,count]=fscanf(fid2,'%g %g',[2,inf]);

fclose(fid2);

data1=data1';
[N_poly,m]=size(data1);

%  hypocentral locations, N=number of hypocenters
x_poly(1:N_poly)=data1(1:N_poly,1);
y_poly(1:N_poly)=data1(1:N_poly,2);
z_poly(1:N_poly)=0.0;

%**************************************************************************
% plot the input data with bounding polygon
figure;
hold on
plot3(xs,ys,zs,'o','MarkerEdgeColor','k','MarkerFaceColor',[0 0.75 0.75]);
plot3(x_poly,y_poly,z_poly,'-k');
hold off
axis equal;
grid on
title('Input hypocenters');
xlabel('X km');
ylabel('Y km');
zlabel('Z km');

%***************  Produce a Randomized Catalog ****************************
[xr,yr,zr]=randomcat(xs,ys,zs,x_poly,y_poly,z_poly);

%**************************************************************************
%  Calculate tetrahedra volumes for each event location in random catalog
%  and sort in ascending order of volume
[xvr,yvr,zvr,volr]=tetvol(xr,yr,zr,N);

%plot volume of tetrahedra vs index
figure;
plot(volr);

%**************************************************************************
%  Calculate tetrahedra volumes for each event location in original catalog
%  and sort in ascending order of volume
[xvs,yvs,zvs,vols]=tetvol(xs,ys,zs,N);
%zvs

%plot volume of tetrahedra vs index
figure;
plot(vols);

%**************************************************************************
%  Compute the normalized cumulative probability density function of the
%  random catalog tetrahedra volumes
[cdf_rnd,vr]=ecdf(volr);

nkh=length(vr);
% for kh=1:nkh;
%     fprintf('cdf_rnd %g  vr %g\n',cdf_rnd(kh),vr(kh));
% end

%**************************************************************************
%  Compute the normalized cumulative probability density function of the
%  original catalog tetrahedra volumes0
[cdf_cat,vs]=ecdf(vols);
 
%**************************************************************************
%  Calculate the volume of tetrahedra at the 5% probability level (0.05) 
%  from the random catalog
for kh=1:nkh-1
    if cdf_rnd(kh) < PROB & cdf_rnd(kh+1) >= PROB
        ncalc=kh;
    end
end
ncalc;
V05=vr(ncalc);
fprintf('ncalc = %i, Volume at 5 percent probability (V05, km cubed) = %g\n',ncalc,V05);

%V05=3.0;
%**************************************************************************
%  Calculate N(V(0.05)) from original catalog
nks=length(vs);
for ks=1:nks-1
    if vs(ks) < V05 & vs(ks+1) >= V05
        ncat=ks;
    end
end
ncat;
NV05_cat=cdf_cat(ncat);
fprintf('ncat = %i, Probability at target volume (NV05) = %g\n',ncalc,NV05_cat);

% plot the cdfs with determined values
figure;
hold on;
plot(log(vr),cdf_rnd,'-k',log(vs),cdf_cat);
plot([log(V05) log(V05)],[0 1],'-r');
hold off;

title('Cumulative Distribution Functions');
xlabel('ln(volume) km^3');
ylabel('Probability');
legend('Random','Catalog');

%**************************************************************************
%  Separate diffuse and clustered events in original catalog
%zvs
[x_dif,y_dif,z_dif,x_cls,y_cls,z_cls]=sep_cat(xvs,yvs,zvs,ncat,N);
%z_cls

N_clustered=length(x_cls);
N_diffuse=length(x_dif);
fprintf('N_clustered = %i  N_diffuse = %i\n',N_clustered,N_diffuse);

% plot the separated catalogs
figure;
plot3(x_cls,y_cls,z_cls,'o','MarkerEdgeColor','k','MarkerFaceColor',[0 0.75 0.75]);
axis equal;
grid on
title('clustered hypocenters');
xlabel('X km');
ylabel('Y km');
zlabel('Z km');

figure;
plot3(x_dif,y_dif,z_dif,'o','MarkerEdgeColor','k','MarkerFaceColor',[0 0.75 0.75]);
axis equal;
grid on
title('diffuse hypocenters');
xlabel('X km');
ylabel('Y km');
zlabel('Z km');

%**************************************************************************
%  Save new catalogs to file
% write the block of clustered seismicity to a file

fid=fopen(catfileout_clustered,'w');
nblock=length(x_cls);
for kk=1:nblock
    fprintf(fid,'%12.5f %12.5f %12.5f\n',[x_cls(kk) y_cls(kk) z_cls(kk)]);
end

fclose(fid);

% write the block of diffuse seismicity to a file

fid=fopen(catfileout_diffuse,'w');
nblock=length(x_dif);
for kk=1:nblock
    fprintf(fid,'%12.5f %12.5f %12.5f\n',[x_dif(kk) y_dif(kk) z_dif(kk)]);
end

fclose(fid);

end

