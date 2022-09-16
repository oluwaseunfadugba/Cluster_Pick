function plot_faults(h,xs,ys,zs,xv,yv,zv,flt_no)
% plot a rendition of the data and the planes that fit the data
%figure('Name','picname');
%hold on;

% global tag 

for i = 1:length(flt_no)
   
    plot3(h,xs,ys,zs,...
        'o','MarkerEdgeColor','k','MarkerFaceColor',[0 0.75 0.75]);hold on;

end
hold on;

for k=flt_no
    fill3(h,xv(k,1:4),yv(k,1:4),zv(k,1:4),'w','FaceAlpha',0.7,'FaceColor',[0.5 0.5 0.5]);
end
hold off;

axis equal;
grid on;
%title('Input hypocenters');
xlabel('X km');
ylabel('Y km');
zlabel('Z km');
set(h, 'fontsize', 18);
grid on;

if length(xs) >1
xlim([min(xs)-std(xs) max(xs)+std(xs)])
ylim([min(ys)-std(ys) max(ys)+std(ys)])
zlim([min(zs)-std(zs) max(zs)+std(zs)])
end
% AZ 
% EL

% view(AZ,EL);