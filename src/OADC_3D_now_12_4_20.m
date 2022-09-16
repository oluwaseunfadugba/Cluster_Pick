function OADC_3D_now_12_4_20(i,arg_kmin,arg_kmax,arg_err_av,...
    arg_hypo_infile,arg_clus_mineqs,arg_N_loop,arg_simul_tag,...
    arg_dip_threshold,arg_FM_file,arg_PLOT_FLAG0,arg_PLOT_FLAG1,...
    arg_comb_coplanar, arg_plot_avg_FM)

%  Implementation of 3-D Optimal Anisotropic Dynamic Clustering from

%  Ouillon, Ducorbier, and Sornette (2008).  Automatic reconstruction of
%  fault networks from seismicity catalogs: Three-dimensional optimal
%  anisotropic dynamic clustering, JGR, 113, B01306,
%  doi:10.1029/2007JB00503.
%
%  specify:
%
%       i = index number for random number generator, e.g., rng(i)
%       kmin = starting number of fault planes
%       kmax = maximum number of fault planes analyzed
%       err_av = average hypocentral error in km
%       hypo_infile = file containing (x,y,z) positions of hypocenters
%       clus_mineqs = minimum number of eqs per cluster
%       Nloop_per_iter = Number of loop per each iteration
%       dip_thresh = dip angle threshold
%       PLOT_FLAG1 = 1;   % =0, no intermediate loop plots of data and planes
%       comb_coplanar = 1; % =0, Do not check or combine coplananr faults 
%       FM_file = file containing focal mechanism containing (xsf,ysf,zsf,
%                   strike, dip, rake)
%       arg_plot_avg_FM = 0; % =1, Plot intermediate average FM focalsphere
%       simul_tag = simulation tag for folder name.
%
% <command execution_time="8810">OADC_3D(1,7,0.01,'testdata.txt')</command>
% <command execution_time="4459">OADC_3D(1,3,0.5,'testdata.txt')</command>
% <command execution_time="5385">OADC_3D(1,4,0.5,'testdata.txt')</command>

clear GLOBAL; 
%close all; 
clc; tic

global xc yc zc vec_plane xb_old yb_old zb_old xs ys zs N Nc
global xt yt zt Nt xb yb zb lambda3 Strike Dip
global L W xv yv zv L_old W_old xv_old yv_old zv_old fscale
global con_tol skip_model 

global kmin kmax err_av hypo_infile clus_mineqs 
global N_loop simul_tag dip_threshold FM_file 
global PLOT_FLAG0 PLOT_FLAG1 comb_coplanar plot_avg_FM

 
% % Test Input parameters
% i = 1; arg_kmin = 1; arg_kmax=10; arg_err_av = 0.35; 
% arg_hypo_infile = 'testdata.txt';%'Simul_hypos.txt';
% arg_clus_mineqs = 4; arg_N_loop = 1; arg_simul_tag = ['OADC.Model.i' num2str(i)];
% arg_dip_threshold = 10; arg_FM_file = '';%'All_write_mechanisms.txt'; 
% arg_PLOT_FLAG1 = 0;   % =0, no intermediate loop plots of data and planes
% arg_comb_coplanar = 1; % =0, Do not check or combine coplananr faults
% arg_plot_avg_FM = 0; % =1, Plot intermediate average FM focalsphere
% 

kmin = arg_kmin; kmax = arg_kmax; err_av = arg_err_av; hypo_infile = arg_hypo_infile;
clus_mineqs = arg_clus_mineqs; N_loop = arg_N_loop; simul_tag = arg_simul_tag; 
dip_threshold = arg_dip_threshold; FM_file = arg_FM_file;
PLOT_FLAG0 = arg_PLOT_FLAG0;PLOT_FLAG1 = arg_PLOT_FLAG1; 
comb_coplanar = arg_comb_coplanar; plot_avg_FM = arg_plot_avg_FM;

% Checking for consistency 
if PLOT_FLAG0 == 0
    PLOT_FLAG1 = 0;
end

fprintf('******************** i = %s *******************************\n',num2str(i));

% remove previous simulations with the same simul_tag
eval(sprintf('%s%s%s %s','! rm -rf ',simul_tag, '*', '*~'))

% Save diary
diaryname = [simul_tag '.Diary.txt'];
diary(diaryname)

%********************** Set Parameters ************************************
%   Fault length scale for random faults. Will be between 0 and fscale in km
fscale=50.0;
%   Convergence tolerance value for the clustering algorithm in 'faultcluster'.  
%   Represents the smallest change in global variance with hypocenter 
%   clustering iteration.  The clustering process will stop once the change
%   in global variance with iteration drops to this value or smaller.
con_tol=0.1;  %0.00001;%  units usually in km
%PLOT_FLAG1=1;   % =0, no intermediate loop plots of data and planes

rng(i); %rng('default'); % Setting random number generator

%***************** Read Catalog of Hypocenters ****************************
read_catalog(hypo_infile);

%********************** Initialize Space **********************************
init_space(kmax);

%******************* Initialize random faults *****************************
FAULT_FLAG=0;   % Initialization, use all hypocenters
randfaults(kmin,FAULT_FLAG);

if PLOT_FLAG1 == 1
    %  plot initial planes
    picname='Initial Model';
    datplot(xs,ys,zs,kmin,xv,yv,zv,picname);
end

SOL_FLAG=0; Kfaults=kmin;

%******************** Big Loop over Kfaults *******************************
while Kfaults <= kmax
    Kfaults

    %  form clusters of seismicity using present number of random faults.
    %  Much of the work is done here
    JFINAL=faultcluster(con_tol,Kfaults)
    
    %  plot initial planes
    if PLOT_FLAG1 == 1
        picname=strcat('Iteration',num2str(Kfaults),'Model');
        datplot(xs,ys,zs,Kfaults,xv,yv,zv,picname);
    end

    %  test to see if fit is within the error. Look at largest lambda3 eigenvalue
    lambda3max=max(lambda3)
    if lambda3max <= err_av
        %  print the good news
        fprintf('Fault model converged to within error!\n');
        SOL_FLAG=1;
        Kfaults_good=Kfaults;
        Kfaults=kmax+1;
         
    else
        % split the thickest fault into two new random fault planes
        if Kfaults < kmax
            fprintf('Splitting thickest fault, Kfaults= %i +1\n',Kfaults);
            
            % split the thickest fault
            splitfault_driver(Kfaults);
            
            % increase the fault number
            Kfaults=Kfaults+1;
         
            %  plot new planes with data
            if skip_model == 0
                if PLOT_FLAG1 == 1 
                    picname='Model from fault split';
                    datplot(xs,ys,zs,Kfaults,xv,yv,zv,picname);
                end
            end
            
        else
             Kfaults=Kfaults+1;

        end        
    end
end

if skip_model == 0
    %  Output the final fault model (knowing it is not optimal)
    if SOL_FLAG == 0
        fprintf('Analysis Complete, fault model is not optimal\n');
        Kfaults=Kfaults-1;
    else
        Kfaults=Kfaults_good;
    end
    
    if comb_coplanar == 0 % Checking and combining coplananr faults 
        if PLOT_FLAG0 == 1 
            %  plot final planes
            picname='Final Model';
            datplot(xs,ys,zs,Kfaults,xv,yv,zv,picname);
            datplot_with_colors(xs,ys,zs,Kfaults,xv,yv,zv,picname);   
        end
    end

    % Postprocessing and Saving files
    % saving all variables to file
    !rm -r *savevar_filename.mat
    savevar_filename = [simul_tag '.saved_variables.mat'];
    save(savevar_filename)

    if comb_coplanar == 1
        fprintf('\n');
        fprintf('Checking and combining coplananr faults!\n');
        
        % Checking and combining coplananr faults
        combine_coplananr_faults(savevar_filename)
    end
    
    if PLOT_FLAG0 ==1
        % saving all figures
        FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
        savefig(FigList(end:-1:1),[simul_tag '.figures.fig'])
    end
    
    eval(sprintf('%s%s%s','! mkdir ',simul_tag,'_results'))
    eval(sprintf('%s%s%s %s%s','! mv ',simul_tag, '*',simul_tag,'_results'))
end

toc
%return