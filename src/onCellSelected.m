function onCellSelected(hObject, g)

        global picks h_table h1 Kfaults
        global xs ys zs xv yv zv

        % fprintf('Editting cell (%d,%d)\n', g.Indices(1), g.Indices(2));        
    
        picks(g.Indices(1), 1) = num2cell(g.NewData);
            
        viewpt = h1.View;
        
        cla(h1,'reset');
        set(h1,'NextPlot','add');
 
        flt_no = 1:Kfaults;
        flt_no = flt_no(cell2mat(picks(:,1)));
        plot_faults(h1,xs,ys,zs,xv,yv,zv,flt_no)  
        
        xlim([min(xs)-1 max(xs)+1])
        ylim([min(ys)-1 max(ys)+1])
        zlim([min(zs)-1 max(zs)+1])

        view(h1,viewpt(1),viewpt(2));
        rotate3d on
end