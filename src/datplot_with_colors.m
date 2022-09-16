function datplot_with_colors(xs,ys,zs,n0,xv,yv,zv,picname)

global xt yt zt Nt xb yb zb lambda3

hold on;

plot3(xt(1,1:Nt(1)),yt(1,1:Nt(1)),zt(1,1:Nt(1)),'o','MarkerEdgeColor',...
    'k','MarkerFaceColor','g'); hold on;

colors=[0. 1 1; 
    0.5 0.5 0.8;
    0.7   0   0;
    1   0.8 0.8;
    0.5 0.8 0.5;
    0.5 0.8 0.5;
    0.7 0.9 0.5;
    0.7   0   0;
    1   0.8 0.8;
    0.5 0.8 0.5;
    0.5 0.8 0.5;
    0.7 0.9 0.5];

for k=1:n0
    fill3(xv(k,1:4),yv(k,1:4),zv(k,1:4),'w','FaceAlpha',0.6,'FaceColor',...
        colors(mod(k,11)+2,:));
end

hold on;

for k=2:n0
    plot3(xt(k,1:Nt(k)),yt(k,1:Nt(k)),zt(k,1:Nt(k)),'o','MarkerEdgeColor',...
        'k','MarkerFaceColor',colors(mod(k,11)+1,:)); hold on;%'k'
end

hold off;
axis equal;
grid on;
%title('Input hypocenters');
xlabel('X km');
ylabel('Y km');
zlabel('Z km');
grid on;
view(3);

title(picname);
set(gca, 'fontsize', 18);

return;