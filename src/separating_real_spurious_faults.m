clear all; close all; clc; tic

% Main Results
% 0. All Models
% All 44 modesls
%
% 1. To be included in the paper
% 9 14 18 25 28 29 58 61 72 96 
%
% 2. Favorite Models
% 9 14 18 24 25 28 29 33 50 55 58 61 72 96 100   
%
% 3. Second Models 
% 2 3 26 33 60 65 70 71 74 83   
%
% 4. Thrid Models
% 45 64
%
% 5. Other Models
%
Main_Results=1;

fid = fopen('OADC_results_database.txt','w'); % open the original vtk file 
fprintf(fid, 'OADC_results_database: \n');
fprintf(fid, 'Model_no Real_faults_no Total_fault_no Fault_no Strike Dip Length Width Lambda_3 xv1 xv2 xv3 xv4 yv1 yv2 yv3 yv4 zv1 zv2 zv3 zv4\n');

if Main_Results==0; Main_Results=1; end
if Main_Results==1
%To_be_included
% 9 14 18 25 28 29 58 61 72 96 
model_no = 9; real_fault_index = [1 2 3 11 12 13 14 15 16 17 18];            just_do_it(fid,model_no,real_fault_index); %8
model_no = 14; real_fault_index = [1 2 4 5 8 9 11 12 13 14 16 17];             just_do_it(fid,model_no,real_fault_index);
model_no = 18; real_fault_index = [2 3 4 5 8 9 11 12 14 15 16];                just_do_it(fid,model_no,real_fault_index);
model_no = 25; real_fault_index = [1 2 3 6 8 11 13 15 16 17 18];               just_do_it(fid,model_no,real_fault_index);
model_no = 28; real_fault_index = [1 2 3 5 6 8  9 11 14 15 16 17];             just_do_it(fid,model_no,real_fault_index);
model_no = 29; real_fault_index = [1 2 3 4 8 10 11 13 14 16 17];               just_do_it(fid,model_no,real_fault_index);
model_no = 58; real_fault_index = [1 3 4 5 6 10 12 13 15 16 17 18];          just_do_it(fid,model_no,real_fault_index); %7
model_no = 61; real_fault_index = [1 2 3 4 6 7 8 9 10 11 12];               just_do_it(fid,model_no,real_fault_index);%16
model_no = 72; real_fault_index = [1 3 4 5 6 9 12 13 15 16 17];                just_do_it(fid,model_no,real_fault_index);%5
model_no = 96; real_fault_index = [1 3 4 5  9 12 13 15 16 17 14];           just_do_it(fid,model_no,real_fault_index); % 8
end

if Main_Results==0; Main_Results=2; end
if Main_Results==2
% Favorite Models
% 9 14 18 24 25 28 29 33 50 55 58 61 72 96 100 
model_no = 9; real_fault_index = [1 2 3 11 12 13 14 15 16 17 18];            just_do_it(fid,model_no,real_fault_index); %8
model_no = 14; real_fault_index = [1 2 4 5 8 9 11 12 13 14 16 17];             just_do_it(fid,model_no,real_fault_index);
model_no = 18; real_fault_index = [2 3 4 5 8 9 11 12 14 15 16];                just_do_it(fid,model_no,real_fault_index);
model_no = 25; real_fault_index = [1 2 3 6 8 11 13 15 16 17 18];               just_do_it(fid,model_no,real_fault_index);
model_no = 28; real_fault_index = [1 2 3 5 6 8  9 11 14 15 16 17];             just_do_it(fid,model_no,real_fault_index);
model_no = 29; real_fault_index = [1 2 3 4 8 10 11 13 14 16 17];               just_do_it(fid,model_no,real_fault_index);
model_no = 58; real_fault_index = [1 3 4 5 6 10 12 13 15 16 17 18];          just_do_it(fid,model_no,real_fault_index); %7
model_no = 61; real_fault_index = [1 2 3 4 6 7 8 9 10 11 12];               just_do_it(fid,model_no,real_fault_index);% 16
model_no = 72; real_fault_index = [1 3 4 5 6 9 12 13 15 16 17];                just_do_it(fid,model_no,real_fault_index);%5
model_no = 96; real_fault_index = [1 3 4 5  9 12 13 15 16 17 14];           just_do_it(fid,model_no,real_fault_index); %  8
model_no = 24; real_fault_index = [2 4 5 6 7 8 9 10 11 13 14 15 16];           just_do_it(fid,model_no,real_fault_index);
model_no = 50; real_fault_index = [3 5 6 7 8 10 12 14 15 16 17 18];            just_do_it(fid,model_no,real_fault_index);
model_no = 55; real_fault_index = [2 4 5 6 10 11 12 13 15 16 18];              just_do_it(fid,model_no,real_fault_index);
model_no = 100;real_fault_index = [1 2 5 8 9 11 13 14 16];                     just_do_it(fid,model_no,real_fault_index);
end
 
if Main_Results==0; Main_Results=3; end
if Main_Results==3
% Second Models 
% 2 3 26 33 60 65 70 71 74 83 
model_no = 2; real_fault_index = [1 3 4 6 7 12 15 16];                         just_do_it(fid,model_no,real_fault_index);
model_no = 3; real_fault_index = [1 2 3 5 6 7 8 11 12 13 14 15 16 17 18 19];   just_do_it(fid,model_no,real_fault_index);
model_no = 26; real_fault_index = [4 5 6 7 8 9 10 11 12 15 16 17 19];          just_do_it(fid,model_no,real_fault_index);
model_no = 33; real_fault_index = [1 2 3 6 10 11 12 13 14];                    just_do_it(fid,model_no,real_fault_index);
model_no = 60; real_fault_index = [1 2 4 5 8 10 11 12 13 14 16 17 18 19];      just_do_it(fid,model_no,real_fault_index);
model_no = 65; real_fault_index = [1 5 7 8 11 12 15 16];                       just_do_it(fid,model_no,real_fault_index);
model_no = 70; real_fault_index = [1 3 5 6 9 10 12 14 16 17 18];               just_do_it(fid,model_no,real_fault_index); %11 12
model_no = 71; real_fault_index = [1 2 3 4 5 9 10 11 13 14 15 19];             just_do_it(fid,model_no,real_fault_index);%9
model_no = 74; real_fault_index = [2 5 6 8 9 10 11 13 15 17   12];             just_do_it(fid,model_no,real_fault_index);
model_no = 83; real_fault_index = [1 2 3 4 5 6 7 9 10 12 13 14 17 19];         just_do_it(fid,model_no,real_fault_index);
end

if Main_Results==0; Main_Results=4; end
if Main_Results==4
% Third models
% 45 64
model_no = 45; real_fault_index = [1 2 3 6 10 11 12 13];                       just_do_it(fid,model_no,real_fault_index);
model_no = 64; real_fault_index = [1 3 4 6 7 8 9 10 12 13 16 17];              just_do_it(fid,model_no,real_fault_index); %6
end

if Main_Results==0; Main_Results=5; end
if Main_Results==5
% Other models
model_no = 1; real_fault_index = [1 6 7 8 10 11 12 13 14 16 18 19 21];         just_do_it(fid,model_no,real_fault_index);
model_no = 46; real_fault_index = [2 3 4 5 6 8 9 10 11 12 13 14 15];           just_do_it(fid,model_no,real_fault_index);
model_no = 73; real_fault_index = [1 3 4 6 7 8 10 11 12 ];                     just_do_it(fid,model_no,real_fault_index);
model_no = 84; real_fault_index = [2 3 6 9 10 11 12 13 15 16 19 20 21];        just_do_it(fid,model_no,real_fault_index);
model_no = 85; real_fault_index = [1 2 4 5 8 9 11 14 12];                      just_do_it(fid,model_no,real_fault_index);
model_no = 93; real_fault_index = [1 2 5 6 7 8 9 12 14 15 13];                 just_do_it(fid,model_no,real_fault_index);
model_no = 98; real_fault_index = [1 2 3 9 10 11 13 14 15 16];                 just_do_it(fid,model_no,real_fault_index);
model_no = 99; real_fault_index = [1 4 5 6 7 8 11 12 13 15 17 18 19];          just_do_it(fid,model_no,real_fault_index);
end


fclose(fid);

!rm -r OADC_Results_now
!mkdir OADC_Results_now
!mv fig_OADC_Model_No_* OADC_Results_now
!mv OADC_results_database.txt OADC_Results_now

toc


% ----------------- DO NOT EDIT BELOW THIS LINE ---------------------------
function just_do_it(fid,model_no,real_fault_index)
data=load(['/Users/oluwaseunfadugba/Documents/OADC_project/' ...
    'OADC_orig_start_afresh_with_new_modifications/'...
    'OADC.Err_3.Model.' num2str(model_no) '_results/'...
    'OADC.Err_3.Model.' num2str(model_no) '.saved_variables_coplanar.mat']);

xmin = 0;   xmax = 70;
ymin = 0;   ymax = 70;
zmin = -30; zmax = 0;

xs=data.xs; ys=data.ys; zs=data.zs;
xv=data.xv; yv=data.yv; zv=data.zv;
xt=data.xt; yt=data.yt; zt=data.zt;
Kfaults=data.Kfaults; Nt=data.Nt;
Strike=data.Strike; Dip=data.Dip;
L=data.L; W=data.W; lambda3=data.lambda3;

% Writing fault information to file
write_results_to_file(fid,model_no,Kfaults,real_fault_index,Strike,Dip,L,W,lambda3,xv,yv,zv)

picname = ['OADC Model No = ' num2str(model_no)];
Fig1 = figure('Name',picname,'Position', get(0, 'Screensize'));
n_flts = length(real_fault_index);

ax1 = subplot(1,2,1);

Ntr = Nt(real_fault_index);
xvr=xv(real_fault_index,:); yvr=yv(real_fault_index,:); zvr=zv(real_fault_index,:);
xtr=xt(real_fault_index,:); ytr=yt(real_fault_index,:); ztr=zt(real_fault_index,:);

datplot_with_colors_sep_spurious_faults...
    (xs,ys,zs,n_flts,xvr,yvr,zvr,xtr, ytr, ztr, Ntr,picname)

hold on; plot_structures(xmin, xmax, ymin,ymax,zmin,zmax);


ax2 = subplot(1,2,2);

Nt(real_fault_index)=[];
xv(real_fault_index,:)=[]; yv(real_fault_index,:)=[]; zv(real_fault_index,:)=[];
xt(real_fault_index,:)=[]; yt(real_fault_index,:)=[]; zt(real_fault_index,:)=[];

datplot_with_colors_sep_spurious_faults...
    (xs,ys,zs,Kfaults-n_flts,xv,yv,zv,xt, yt, zt, Nt,'Spurious Faults')

hold on; plot_structures(xmin, xmax, ymin,ymax,zmin,zmax);

%camproj('perspective')
% campos([0,0,2.5])
% camtarget([50,0,2.5])

Link = linkprop([ax1, ax2],{'CameraUpVector', 'CameraPosition', 'CameraTarget'});
setappdata(gcf, 'StoreTheLink', Link);
view(-75,50)
rotate3d on

%CaptureFigVid(vw, 'CSZ3',OptionZ)


savefig(Fig1,['fig_OADC_Model_No_' num2str(model_no) '.fig'])


end

function write_results_to_file(fid,model_no,Kfaults,real_fault_index,Strike,Dip,L,W,lambda3,xv,yv,zv)
    Real_faults_no = length(real_fault_index);
    for i=1:Real_faults_no
        fault_no = real_fault_index(i);

        if Dip(fault_no)>90
            Dip(fault_no)=180-Dip(fault_no);

            %Strike(fault_no)=wrapTo360(Strike(fault_no)+180);
        end

        fprintf(fid,'%21f\t',[model_no Real_faults_no Kfaults fault_no Strike(fault_no)...
            Dip(fault_no)   L(fault_no)     W(fault_no)     lambda3(fault_no)...
            xv(fault_no,1)  xv(fault_no,2)  xv(fault_no,3)  xv(fault_no,4)  ...
            yv(fault_no,1)  yv(fault_no,2)  yv(fault_no,3)  yv(fault_no,4)  ...
            zv(fault_no,1)  zv(fault_no,2)  zv(fault_no,3)  zv(fault_no,4)]);

        fprintf(fid,'\n');

    end
    fprintf(fid,'\n');

end

function plot_structures(xmin, xmax, ymin,ymax,zmin,zmax)
plot_box(xmin, xmax, ymin,ymax,zmin,zmax) ;
linewidth = 1.5;


plot3([(xmax-10) (xmax-10) (xmax-15) (xmax-10) (xmax-5)],...
    [5 20 15 20 15]+ymin,[0 0 0 0 0],'LineWidth',linewidth,'color','r'); hold on;

xt = [0 0 10 30 0];
yt = [30 5 5 5 70];
zt = [0 0 0 0 0];
str = {'GNF','SLF','CHF','SSF','Impact Structure'};
text(xt,yt,zt,str,'color','b','FontWeight','bold','FontSize',16)

text((xmax-10),25,0,{'N'},'color','b','FontWeight','bold','FontSize',20)


[x,y,z] = extract_crater_data('Crater_dataset_GMT_inner.csv'); plot3(x,y,z,'LineWidth',linewidth,'color','r'); hold on;
[x,y,z] = extract_crater_data('Crater_dataset_GMT_outer.csv'); plot3(x,y,z,'LineWidth',linewidth,'color','r'); hold on;
[x,y,z] = extract_geo_data('F1_data.txt'); plot3(x,y,z,'LineWidth',linewidth,'color','r'); hold on;
[x,y,z] = extract_geo_data('F2_data.txt'); plot3(x,y,z,'LineWidth',linewidth,'color','r'); hold on;
[x,y,z] = extract_geo_data('F3_data.txt'); plot3(x,y,z,'LineWidth',linewidth,'color','r'); hold on;
[x,y,z] = extract_geo_data('F3_data_2.txt'); plot3(x,y,z,'LineWidth',linewidth,'color','r'); hold on;
[x,y,z] = extract_geo_data('F3_data_3.txt'); plot3(x,y,z,'LineWidth',linewidth,'color','r'); hold on;
[x,y,z] = extract_geo_data('F3_data_4.txt'); plot3(x,y,z,'LineWidth',linewidth,'color','r'); hold on;
[x,y,z] = extract_geo_data('F4_data.txt'); plot3(x,y,z,'LineWidth',linewidth,'color','r'); hold on;

%axis([-15 90 -11 89 -30 0]);
%axis([0 70 0 70 -30 0]);
axis([xmin xmax ymin ymax zmin zmax]);
end

function [x,y,z] = extract_crater_data(infile)   
    fid=fopen(infile,'r'); [data,count]=fscanf(fid,'%g %g %g',[3,inf]); 
    fclose(fid);
    data=data';
    x = (data(:,1) + 70.6)*75.778;
    y = (data(:,2) - 47.2)*111.1743;
    z = zeros(count/3,1);  
end
        
function [x,y,z] = extract_geo_data(infile)   
    fid=fopen(infile,'r'); [data,count]=fscanf(fid,'%g %g',[2,inf]); 
    fclose(fid);
    data=data';
    x = (data(:,1) + 70.6)*75.778;
    y = (data(:,2) - 47.2)*111.1743;
    z = zeros(count/2,1);
end

function plot_box(xmin, xmax, ymin,ymax,zmin,zmax)
%axis([-15 90 -11 89 -30 0]);

hold on;
linewidth = 2;

x = [xmin xmax xmax xmin xmin];
y = [ymin ymin ymax ymax ymin];
plot(x,y ,'r','LineWidth',linewidth); hold on;

x = [xmin xmax xmax xmin xmin];
y = ymin*ones(1,5);
z = [zmax zmax zmin zmin zmax];
plot3(x,y ,z,'r','LineWidth',linewidth); hold on;

x = [xmin xmax xmax xmin xmin];
y = ymax*ones(1,5);
z = [zmax zmax zmin zmin zmax];
plot3(x,y ,z,'r','LineWidth',linewidth); hold on;

x = xmin*ones(1,5);
y = [ymin ymax ymax ymin ymin];
z = [zmax zmax zmin zmin zmax];
plot3(x,y ,z,'r','LineWidth',linewidth); hold on;

x = xmax*ones(1,5);
y = [ymin ymax ymax ymin ymin];
z = [zmax zmax zmin zmin zmax];
plot3(x,y ,z,'r','LineWidth',linewidth); hold on;

end

