function datplot2(xs,ys,zs,n0,xv,yv,zv,picname)

% plot a rendition of the data and the planes that fit the data
figure('Name',picname);
hold on;
plot3(xs,ys,zs,'o','MarkerEdgeColor','k','MarkerFaceColor','k');
%for k=1:n0;
    fill3(xv(n0,1:4),yv(n0,1:4),zv(n0,1:4),'w','FaceAlpha',0.7,'FaceColor',[0.5 0.5 0.5]);
%end;
hold off;
axis equal;
grid on;
%title('Input hypocenters');
xlabel('X km');
ylabel('Y km');
zlabel('Z km');
grid on;
view(3);
return;