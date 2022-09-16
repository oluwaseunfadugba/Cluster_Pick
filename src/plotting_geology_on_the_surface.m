clear all; close all; clc;

Test = open('Simul.1.OADC.Final Model 4_editted.fig'); figure(Test); hold on
h = findobj(gcf,'Type', 'Line'); h.MarkerSize = 5;
 
set(0,'DefaultLegendAutoUpdate','off')
plot_box() ;
plot3([80 80 75 80 85],[0 20 15 20 15],[0 0 0 0 0],'LineWidth',3,'color','r'); hold on;

xt = [-10 0 10 30 0];
yt = [5 5 5 5 70];
zt = [0 0 0 0 0];
str = {'GNF','SLF','CHF','SSF','Impact Structure'};
text(xt,yt,zt,str,'color','b','FontWeight','bold','FontSize',16)

text(78.1,25,0,{'N'},'color','b','FontWeight','bold','FontSize',25)

[x,y,z] = extract_crater_data('Crater_dataset_GMT_outer.csv'); plot3(x,y,z,'LineWidth',3,'color','r'); hold on;
[x,y,z] = extract_geo_data('F1_data.txt'); plot3(x,y,z,'LineWidth',3,'color','r'); hold on;
[x,y,z] = extract_geo_data('F2_data.txt'); plot3(x,y,z,'LineWidth',3,'color','r'); hold on;
[x,y,z] = extract_geo_data('F3_data.txt'); plot3(x,y,z,'LineWidth',3,'color','r'); hold on;
[x,y,z] = extract_geo_data('F3_data_2.txt'); plot3(x,y,z,'LineWidth',3,'color','r'); hold on;
[x,y,z] = extract_geo_data('F3_data_3.txt'); plot3(x,y,z,'LineWidth',3,'color','r'); hold on;
[x,y,z] = extract_geo_data('F3_data_4.txt'); plot3(x,y,z,'LineWidth',3,'color','r'); hold on;
[x,y,z] = extract_geo_data('F4_data.txt'); plot3(x,y,z,'LineWidth',3,'color','r'); hold on;

axis([-15 90 -11 89 -30 0]);
title('Final Model')

%% Set up recording parameters (optional), and record
OptionZ.FrameRate=15;OptionZ.Duration=20;OptionZ.Periodic=true;
%CaptureFigVid([-48,16;-77,38; -20,10;-110,10;-190,80;-290,10;-380,10;-77,38], 'CSZ',OptionZ)

vw = [0,90;
-58,56;
0,90;
-15,60;
-28,28;
-30,25;
-50,30;
-59,60;
-77,31;
-30,25;
-62,31;
 -84,42;
 -44,26];

vw = [0,90;
      155,36;
      45,36;
      -83,60;
      -33,12;
      162,37;
      90,90;
      -38,30]; 
vw = [0,90;
    -7,75;
  -18,49;
  -41,36
  -46,30
  -46,25
  -52,24
  -64,31
  -90,50
  -121,57
  -151,62
  -166,60
  -177,53
  -200,36
  153,27
  135,52
  114,43
  90,52
  41,66
  10,66
  -7,48
  -17,30
  -27,23
  -42,30];
  
  
 vw = [0,90;
    -7,75;
  -18,49;
  -41,36
  -46,30
  -46,25
  -52,24
  -64,31
  -90,50
  -121,57
  -151,62
  -166,60
  -177,53
  -200,36
  -207,27
  -225,52
  -246,43
  -270,52
  -319,66
  -350,66
  -367,48
  -377,30
  -387,23
  -402,30];

CaptureFigVid(vw, 'CSZ3',OptionZ)

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

function plot_box()
hold on;

x = [-15 90 90 -15 -15];
y = [-11 -11 89 89 -11];
plot(x,y ,'r','LineWidth',3); hold on;

x = [-15 90 90 -15 -15];
y = -11*ones(1,5);
z = [0 0 -30 -30 0];
plot3(x,y ,z,'r','LineWidth',3); hold on;

x = [-15 90 90 -15 -15];
y = 89*ones(1,5);
z = [0 0 -30 -30 0];
plot3(x,y ,z,'r','LineWidth',3); hold on;

x = -15*ones(1,5);
y = [-11 89 89 -11 -11];
z = [0 0 -30 -30 0];
plot3(x,y ,z,'r','LineWidth',3); hold on;

x = 90*ones(1,5);
y = [-11 89 89 -11 -11];
z = [0 0 -30 -30 0];
plot3(x,y ,z,'r','LineWidth',3); hold on;

end
