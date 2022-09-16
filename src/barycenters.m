function barycenters

close all;

% make a plot of the Ouillion et al (2008) clustering figure

% Area ranges from 0 to 10 in x and y

% produce a x,y vectors of random integers between 1 and 10

rx=randperm(9) + randn(1,9);
ry=randperm(9) + randn(1,9);

% produce 2 seed locations

sxr=randperm(9) + randn(1,9);
syr=randperm(9) + randn(1,9);

sx(1:2)=sxr(1:2);
sy(1:2)=sxr(1:2);

figure;

subplot(3,2,1);
plot(rx,ry,'ok',sx,sy,'sk','Markerfacecolor','k');
axis square;


%  Iteration (use 3)

nplot=1;

for kit=1:6


% calculate distances and sort to each seed location

nc1=0;  %cluster 1 counter
nc2=0;  %cluster 2 counter

J=0;

for m=1:9
    
    dist1=sqrt((rx(m)-sx(1))^2 + (ry(m)-sy(1))^2);
    dist2=sqrt((rx(m)-sx(2))^2 + (ry(m)-sy(2))^2);
    
    % sort
    if dist1 <= dist2
        nc1=nc1+1;
        xc1(nc1)=rx(m);
        yc1(nc1)=ry(m);
        dc1(nc1)=dist1;
    else
        nc2=nc2+1;
        xc2(nc2)=rx(m);
        yc2(nc2)=ry(m);
        dc2(nc2)=dist2;
    end
    
end

% calculate norm

jnorm=0;
if nc1 > 0
    for k=1:nc1; jnorm=jnorm + dc1(k).*dc1(k); end
end
if nc2 > 0
    for k=1:nc2; jnorm=jnorm + dc2(k).*dc2(k); end
end

kit
jnorm=jnorm./9.0


if kit <= 2 || kit == 6
    
% plot line from each seed location to cluster events
nplot=nplot+1;
subplot(3,2,nplot);
hold on;
if nc1 > 0
    plot(sx(1),sy(1),'sk',xc1,yc1,'ok','Markerfacecolor','k');
    for k=1:nc1
        plot([sx(1) xc1(k)],[sy(1) yc1(k)],'-k');
    end
end
if nc2 > 0
    plot(sx(2),sy(2),'sk',xc2,yc2,'ok','Markerfacecolor','k');
    for k=1:nc2
        plot([sx(2) xc2(k)],[sy(2) yc2(k)],'-k');
    end
end
hold off;
axis square;
norm_nam=strcat('J=',num2str(jnorm));
title(norm_nam);

end

% find each barycenter
bxc1=0.;
byc1=0.;
bxc2=0.;
byc2=0.;

if nc1 > 0
    for k=1:nc1; bxc1=bxc1+xc1(k); byc1=byc1+yc1(k); end
    bxc1=bxc1./nc1;
    byc1=byc1./nc1;
end
if nc2 > 0
    for k=1:nc2; bxc2=bxc2+xc2(k); byc2=byc2+yc2(k); end
    bxc2=bxc2./nc2;
    byc2=byc2./nc2;
end

if kit <= 2
    
% plot this new information
nplot=nplot+1;
subplot(3,2,nplot);
hold on;
if nc1 > 0
    plot(bxc1,byc1,'sk',xc1,yc1,'ok','Markerfacecolor','k');
    for k=1:nc1
        plot([bxc1 xc1(k)],[byc1 yc1(k)],'-k');
    end
end
if nc2 > 0
    plot(bxc2,byc2,'sk',xc2,yc2,'ok','Markerfacecolor','k');
    for k=1:nc2
        plot([bxc2 xc2(k)],[byc2 yc2(k)],'-k');
    end
end
hold off;
axis square;

end
        
% replace seed locations for next iteration
sx(1)=bxc1;
sy(1)=byc1;
sx(2)=bxc2;
sy(2)=byc2;

end
