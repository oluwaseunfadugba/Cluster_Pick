function [x_dif,y_dif,z_dif,x_cls,y_cls,z_cls]=sep_cat(xvs,yvs,zvs,ncat,N)
% sep_cat - separate the catalog into diffuse and clustered seismicity


ndif=N-ncat;

x_dif(1:ndif)=0.0;
y_dif(1:ndif)=0.0;
z_dif(1:ndif)=0.0;

x_cls(1:ncat)=0.0;
y_cls(1:ncat)=0.0;
z_cls(1:ncat)=0.0;

x_cls(1:ncat)=xvs(1:ncat);
y_cls(1:ncat)=yvs(1:ncat);
z_cls(1:ncat)=zvs(1:ncat);

x_dif(1:ndif)=xvs(1+ncat:N);
y_dif(1:ndif)=yvs(1+ncat:N);
z_dif(1:ndif)=zvs(1+ncat:N);

end

