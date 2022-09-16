function onCellSelected_gen_synthetic(hObject, g)

        global h_table h1 
        global picks_3D Kfaults_3D xs_3D ys_3D zs_3D tag_3D flt_no_3D is_3D
        global picks_2D Kfaults_2D xs_2D ys_2D  tag_2D flt_no_2D is_2D
        
        if is_3D == 1

            if g.Indices(2) == 1
                picks_3D(g.Indices(1), 1) = num2cell(g.NewData);

                viewpt = h1.View;

                cla(h1,'reset');
                set(h1,'NextPlot','add');

                % remove datasets and plot the remaining datasets            
                flt_no_3D = 1:Kfaults_3D;
                flt_no_3D = flt_no_3D(cell2mat(picks_3D(:,1)));

                plot_faults2(h1)  

                view(h1,viewpt(1),viewpt(2));
            else

                picks_3D(g.Indices(1), g.Indices(2)) = num2cell(g.NewData);

            end
            
        elseif is_2D == 1
            
            if g.Indices(2) == 1
                picks_2D(g.Indices(1), 1) = num2cell(g.NewData);

                viewpt = h1.View;
                
                cla(h1,'reset');
                set(h1,'NextPlot','add');

                % remove datasets and plot the remaining datasets            
                flt_no_2D = 1:Kfaults_2D;
                flt_no_2D = flt_no_2D(cell2mat(picks_2D(:,1)));

                plot_faults2(h1)  

                view(h1,viewpt(1),viewpt(2));
            else

                picks_2D(g.Indices(1), g.Indices(2)) = num2cell(g.NewData);

            end             
        end
end