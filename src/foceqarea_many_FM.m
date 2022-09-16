function foceqarea_many_FM(c_all,icontour_all,numstart_all,numfin_all,p_pol, data_flag,line_color)
%
% numstart_all
% numfin_all
%   produce a figure of the zero amplitude contour
%   on an equal area, lower hemisphere, projection
%
%   Input:  c = contour array of zero amplitude contour
%           icontour = number of contours in plot
%           num(icontour)=array of data numbers per contour
%           numstart(icontour)=array of starting positions
%           numfin(icontour)=array of finishing positions for
%                            each contour
%           p_pol = array of P wave polarity data
%                  arranged in polarity, azimuth, incidence
%                  angle 3-tuples
%                  polarity varies as -1,0,+1
%           data_flag = 1 for plot P-wave polarity
%                       0 otherwise
%*********************************************************
%
%   create the outside circle of the projection
%
con=pi/180.;
radius=1.0;
a=linspace(0,360,1000)*con;
yc=radius*cos(a);
xc=radius*sin(a);

% convert Pwave amplitudes to polarities (i,e. +1 or -1)
p_pol(:,1) = sign(p_pol(:,1));

%figure;
%figure('color','w','name','P Wave Nodal Surfaces');
%hold on;
daspect([1 1 1]);
%axis off;
plot(xc,yc,'linewidth',1.5,'color','r')%[0.5,0.5,0.5]);

[no_of_FM,~] = size(icontour_all);

for j = 1:no_of_FM
    j;
    icontour = icontour_all(j,:);
    a = 2*j - 1; b = 2*j;
        
    for m=1:icontour
        
        c = c_all(a:b,:); c( :, ~any(c,1) ) = []; % Removing columns containing zeros
        numstart = numstart_all(j,:); numstart( :, ~any(numstart,1) ) = [];
        numfin = numfin_all(j,:); numfin( :, ~any(numfin,1) ) = [];

        azn=c(1,numstart(m):numfin(m)-1); 
        incn=c(2,numstart(m):numfin(m)-1);
        r=sqrt(2)*radius*sin(incn*con/2);
        ycnod=r.*cos(azn*con);
        xcnod=r.*sin(azn*con);
        plot(xcnod,ycnod,'linewidth',1.5,'color', line_color);        
    end
end
%
%   plot N,E,S,W and center tick marks
%
for it=1:4
    az=90*(it-1)*con;
    r1=radius*0.97;
    r2=radius*1.03;
    xt=radius*1.1*sin(az);
    yt=radius*1.1*cos(az);
    plot([r1*sin(az) r2*sin(az)],[r1*cos(az) r2*cos(az)]);
    if it == 1
        text(xt,yt,'N','horizontalalignment','center', 'FontSize',14);
    end
    if it == 2
        text(xt,yt,'E','horizontalalignment','center', 'FontSize',14);
    end
    if it == 3
        text(xt,yt,'S','horizontalalignment','center', 'FontSize',14);
    end
    if it == 4
        text(xt,yt,'W','horizontalalignment','center', 'FontSize',14);
    end
end
plot([-0.03 0.03],[0 0]);
plot([0 0],[-0.03 0.03]);
%
%   plot wave polarities
%
if data_flag == 1

    [pix,~]=size(p_pol);
    ik=1;
    while ik <= pix
        pol=p_pol(ik,1);
        az=p_pol(ik,2);
        inc=p_pol(ik,3);
        r=sqrt(2)*radius*sin(inc*con/2);
        xp=r*sin(az*con);
        yp=r*cos(az*con);
        if pol == -1.0
            mark='o';
            fcol='w';
        end
        if pol == 0.0
            mark='x';
            fcol='k';
        end
        if pol == 1.0
            mark='o';
            fcol='k';
        end
        plot(xp,yp,'marker',mark,'markeredgecolor','k','markerfacecolor',fcol);
        ik=ik+1;
    end
end
%
hold off;
return;