clear all; close all; clc

model_no = 9;  real_fault_index = [1 2 3 11 12 13 14 15 16 17 18];%8
model_no = 14; real_fault_index = [1 2 4 5 8 9 11 12 13 14 16 17];      
model_no = 18; real_fault_index = [2 3 4 5 8 9 11 12 14 15 16];         
model_no = 25; real_fault_index = [1 2 3 6 8 11 13 15 16 17 18];        
model_no = 28; real_fault_index = [1 2 3 5 6 8  9 11 14 15 16 17];      
model_no = 29; real_fault_index = [1 2 3 4 8 10 11 13 14 16 17];        
model_no = 58; real_fault_index = [1 3 4 5 6 10 12 13 15 16 17 18];  %7 
% model_no = 61; real_fault_index = [1 2 3 4 6 7 8 9 10 11 12 ]; %16       
% model_no = 72; real_fault_index = [1 3 4 5 6 9 12 13 15 16 17];         
% model_no = 96; real_fault_index = [1 3 4 5  9 12 13 15 16 17 14  ];   %8 

linewidth = 1.5;
picname=['Fault Model ' num2str(model_no)];
Fig1 = figure('Name',picname,'Position', get(0, 'Screensize'),'color','w');

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

[x,y,z] = extract_crater_data('Crater_dataset_GMT_inner.csv'); plot(x,y,'LineWidth',linewidth,'color','r'); hold on;
text(59,25,0,{'N'},'color','b','FontWeight','bold','FontSize',20); hold on;
plot([(xmax-10) (xmax-10) (xmax-15) (xmax-10) (xmax-5)],...
    [5 20 15 20 15]+ymin,'LineWidth',linewidth,'color','r'); hold on;


for i = 1: length(real_fault_index)
    fault_no = real_fault_index(i);
    
x = mean(xv(fault_no,:));
y = mean(yv(fault_no,:));
z = mean(zv(fault_no,:));

if Dip(fault_no)>90
    Dip(fault_no)=180-Dip(fault_no);

    Strike(fault_no)=wrapTo360(Strike(fault_no)+180);
end
        
%[x y  L(fault_no) Strike(fault_no) Dip(fault_no)]
 
x1 = x - (0.5* L(fault_no)* sind(Strike(fault_no)));
x2 = x + (0.5* L(fault_no)* sind(Strike(fault_no)));

y1 = y - (0.5* L(fault_no)* cosd(Strike(fault_no)));
y2 = y + (0.5* L(fault_no)* cosd(Strike(fault_no)));
 

% [a,b]=sort(zv(fault_no,:));
% x1 = xv(fault_no,b(3));
% x2 = xv(fault_no,b(4));
% y1 = yv(fault_no,b(3));
% y2 = yv(fault_no,b(4));

plot([x1 x2],[y1 y2],'LineWidth',linewidth); axis equal; hold on;
text(x,y,num2str(fault_no),'FontSize',16); hold on;

end

xlim([0 65])
ylim([10 75])
grid on;
title(picname);
set(gca, 'fontsize', 18);

savefig(Fig1,['fig_fault_OADC_Model_No_' num2str(model_no) '.fig'])
F1    = getframe(Fig1);
imwrite(F1.cdata, ['fig_fault_OADC_Model_No_' num2str(model_no) '.png'], 'png')








function [x,y,z] = extract_crater_data(infile)   
    fid=fopen(infile,'r'); [data,count]=fscanf(fid,'%g %g %g',[3,inf]); 
    fclose(fid);
    data=data';
    x = (data(:,1) + 70.6)*75.778;
    y = (data(:,2) - 47.2)*111.1743;
    z = zeros(count/3,1);  
end