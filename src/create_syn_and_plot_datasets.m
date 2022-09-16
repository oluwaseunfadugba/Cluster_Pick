function create_syn_and_plot_datasets()

global  h_table h1 
global picks_3D Kfaults_3D xs_3D ys_3D zs_3D tag_3D flt_no_3D is_3D 
global picks_2D Kfaults_2D xs_2D ys_2D  tag_2D flt_no_2D is_2D 

if is_3D == 1
    xs_3D = []; ys_3D = []; zs_3D = []; tag_3D = [];

    Kfaults_3D = length(picks_3D(:,1));
    flt_no_3D = 1:Kfaults_3D;
    flt_no_3D = flt_no_3D(cell2mat(picks_3D(:,1)));
    picks_mat_3D = cell2mat(picks_3D(:,2:end));

    for i=flt_no_3D

            strike_3D=picks_mat_3D(i,1); dip_3D=picks_mat_3D(i,2); rake_3D = 0;
            L_3D=picks_mat_3D(i,3); W_3D=picks_mat_3D(i,4);  
            zerr_av_3D=picks_mat_3D(i,5); nhypos_3D=picks_mat_3D(i,6); 
            rt_3D=[picks_mat_3D(i,7) picks_mat_3D(i,8) picks_mat_3D(i,9)];
            [xs_i_3D,ys_i_3D,zs_i_3D] = rand_hypos2...
                (L_3D,W_3D,zerr_av_3D,nhypos_3D,strike_3D,dip_3D,rake_3D,rt_3D);
            tag_i_3D = i*ones(1,nhypos_3D);

            xs_3D = [xs_3D xs_i_3D];
            ys_3D = [ys_3D ys_i_3D];
            zs_3D = [zs_3D zs_i_3D];
            tag_3D = [tag_3D tag_i_3D];

    end

    plot_faults2(h1) 
    rotate3d on
    
elseif is_2D == 1
    xs_2D = []; ys_2D = []; tag_2D = [];

    Kfaults_2D = length(picks_2D(:,1));
    flt_no_2D = 1:Kfaults_2D;
    flt_no_2D = flt_no_2D(cell2mat(picks_2D(:,1)));
    picks_mat_2D = cell2mat(picks_2D(:,2:end));

    for i=flt_no_2D

            strike_2D=picks_mat_2D(i,1); 
            L_2D=picks_mat_2D(i,2); 
            zerr_av_2D=picks_mat_2D(i,3);
            nhypos_2D=picks_mat_2D(i,4);
            rt_2D=[picks_mat_2D(i,5) picks_mat_2D(i,6)];
          
            [xs_i_2D,ys_i_2D] = rand_hypos2D(L_2D,zerr_av_2D,nhypos_2D,strike_2D,rt_2D);
            tag_i_2D = i*ones(1,nhypos_2D);

            xs_2D = [xs_2D xs_i_2D];
            ys_2D = [ys_2D ys_i_2D];
            tag_2D = [tag_2D tag_i_2D];

    end

    plot_faults2(h1) 
    
end
