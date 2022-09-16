function cluster_pick_Driver_HPC(i)
tic
% Add the path to the source codes
addpath('src')

arg_kmin = 1; arg_kmax=100; arg_err_av = 2.5; 
arg_hypo_infile = 'testdata.txt';
arg_clus_mineqs = 4; arg_N_loop = 5; 
arg_dip_threshold = 10; arg_FM_file = 'All_write_mechanisms.txt'; %'';%
arg_PLOT_FLAG0 = 1;   % =0, no plots at all
arg_PLOT_FLAG1 = 0;   % =0, no intermediate loop plots of data and planes
arg_comb_coplanar = 1; % =0, Do not check or combine coplananr faults
arg_plot_avg_FM = 0; % =1, Plot intermediate average FM focalsphere

arg_simul_tag = ['Simul3_Faults.i.' num2str(i)];

OADC_3D_now_12_4_20(i,arg_kmin,arg_kmax,arg_err_av,...
            arg_hypo_infile,arg_clus_mineqs,arg_N_loop,arg_simul_tag,...
            arg_dip_threshold,arg_FM_file,arg_PLOT_FLAG0,arg_PLOT_FLAG1,...
            arg_comb_coplanar, arg_plot_avg_FM);
        
close all;

toc