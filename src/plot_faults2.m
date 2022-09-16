function plot_faults2(h)

global tag_3D xs_3D ys_3D zs_3D flt_no_3D is_3D
global tag_2D xs_2D ys_2D  flt_no_2D is_2D

if is_3D == 1
    for i = 1:length(flt_no_3D)
        plot3(h,xs_3D(tag_3D==flt_no_3D(i)),ys_3D(tag_3D==flt_no_3D(i)),...
            zs_3D(tag_3D==flt_no_3D(i)),'o','MarkerEdgeColor','k',...
            'MarkerFaceColor',[0 0.75 0.75]);hold on;
    end
    hold off;

    axis equal;
    xlabel('X km');
    ylabel('Y km');
    zlabel('Z km');
    set(h, 'fontsize', 18);
    grid MINOR;

    if length(xs_3D) >1
        xlim([min(xs_3D)-std(xs_3D) max(xs_3D)+std(xs_3D)])
        ylim([min(ys_3D)-std(ys_3D) max(ys_3D)+std(ys_3D)])
        zlim([min(zs_3D)-std(zs_3D) max(zs_3D)+std(zs_3D)])
    end
    
elseif is_2D == 1
     for i = 1:length(flt_no_2D)
        plot(h,xs_2D(tag_2D==flt_no_2D(i)),ys_2D(tag_2D==flt_no_2D(i)),...
            'o','MarkerEdgeColor','k','MarkerFaceColor',[0 0.75 0.75]);
        hold on;
    end
    hold off;

    axis equal;
    xlabel('X km');
    ylabel('Y km');
    set(h, 'fontsize', 18);
    grid MINOR;

    if length(xs_2D) >1
        xlim([min(xs_2D)-std(xs_2D) max(xs_2D)+std(xs_2D)])
        ylim([min(ys_2D)-std(ys_2D) max(ys_2D)+std(ys_2D)])
    end
  
end