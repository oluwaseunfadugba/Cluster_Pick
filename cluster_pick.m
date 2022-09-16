function varargout = cluster_pick(varargin)
% CLUSTER_PICK MATLAB code for cluster_pick.fig
%      CLUSTER_PICK, by itself, creates a new CLUSTER_PICK or raises the existing
%      singleton*.
%
%      H = CLUSTER_PICK returns the handle to a new CLUSTER_PICK or the handle to
%      the existing singleton*.
%
%      CLUSTER_PICK('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CLUSTER_PICK.M with the given input arguments.
%
%      CLUSTER_PICK('Property','Value',...) creates a new CLUSTER_PICK or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cluster_pick_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cluster_pick_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help cluster_pick

% Last Modified by GUIDE v2.5 15-Sep-2022 12:51:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cluster_pick_OpeningFcn, ...
                   'gui_OutputFcn',  @cluster_pick_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before cluster_pick is made visible.
function cluster_pick_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to cluster_pick (see VARARGIN)

global Cluster_number isauto_OADC ismanual_OADC

% Choose default command line output for cluster_pick
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes cluster_pick wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% Add the path to the source codes
addpath('src')

% Reset the plot_axes and uitable
h1=handles.density_plot_axes;
cla(h1,'reset');

h_table= handles.uitable1;
cla(h_table,'reset');
picks = [num2cell(false(1,1)) num2cell(zeros(1,4))];
set(h_table, 'Data', picks);

% Set default parameters for the sliders
set(handles.Psi_slider,'Value',0);
set(handles.Psi_edit,'String','0');
set(handles.Elevation_slider,'Value',90);
set(handles.Elevation_edit,'String','90');
set(handles.Azimuth_slider,'Value',0);
set(handles.Azimuth_edit,'String','0');

Cluster_number = [];
isauto_OADC = 0;
ismanual_OADC = 0;

% temporary remove the edit view controls
set(handles.Psi_edit,'visible','off');
set(handles.Psi_slider,'visible','off');
set(handles.Elevation_edit,'visible','off');
set(handles.Elevation_slider,'visible','off');
set(handles.Azimuth_edit,'visible','off');
set(handles.Azimuth_slider,'visible','off');
set(handles.text18,'visible','off');
set(handles.text3,'visible','off');
set(handles.text6,'visible','off');

set(handles.is_2D_simul_edit,'Enable','off') 

rotate3d on


% --- Outputs from this function are returned to the command line.
function varargout = cluster_pick_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function Psi_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Psi_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Psi_edit as text
%        str2double(get(hObject,'String')) returns contents of Psi_edit as a double

psi_ed=get(handles.Psi_edit,'String');
fprintf('psi_ed= %g\n',str2double(psi_ed));
set(handles.Psi_slider,'Value',str2double(psi_ed));


% --- Executes on slider movement.
function Psi_slider_Callback(hObject, eventdata, handles)
% hObject    handle to Psi_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

psi_angle=get(handles.Psi_slider,'Value');
fprintf('elevation= %g\n',psi_angle);
set(handles.Psi_edit,'String',num2str(psi_angle));




function Elevation_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Elevation_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Elevation_edit as text
%        str2double(get(hObject,'String')) returns contents of Elevation_edit as a double

elevation_ed=get(handles.Elevation_edit,'String');
fprintf('elevation_ed= %g\n',str2double(elevation_ed));
set(handles.Elevation_slider,'Value',str2double(elevation_ed));


function Azimuth_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Azimuth_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Azimuth_edit as text
%        str2double(get(hObject,'String')) returns contents of Azimuth_edit as a double

azimuth_ed=get(handles.Azimuth_edit,'String');
fprintf('azimuth_ed= %g\n',str2double(azimuth_ed));
set(handles.Azimuth_slider,'Value',str2double(azimuth_ed));



% --- Executes on slider movement.
function Elevation_slider_Callback(hObject, eventdata, handles)
% hObject    handle to Elevation_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

elevation=get(handles.Elevation_slider,'Value');
fprintf('elevation= %g\n',elevation);
set(handles.Elevation_edit,'String',num2str(elevation));



% --- Executes on slider movement.
function Azimuth_slider_Callback(hObject, eventdata, handles)
% hObject    handle to Azimuth_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
azimuth=get(handles.Azimuth_slider,'Value');
fprintf('azimuth= %g\n',azimuth);
set(handles.Azimuth_edit,'String',num2str(azimuth));


function read_file_menu_Callback(hObject, eventdata, handles)
% hObject    handle to read_file_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global data1 xc yc zc vec_plane xb_old yb_old zb_old xs ys zs N Nc 
global xclus yclus zclus Cluster_number N_cluster L W Strike Dip xv yv zv lambda3
global x_orig y_orig D infile

%  Initialize all working space
%initialize_all;
D=[1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0];
x_orig=0.0;
y_orig=0.0;

h1=handles.density_plot_axes;
cla(h1,'reset');

%  Read a seismicity file consisting of (x,y,z), depth is negative
%   Call uigetfile menu

[filename,dirname]=uigetfile('*.*');
if filename == 0; return;end 

infile=strcat(dirname,filename);
fid=fopen(infile,'r');

[data,count]=fscanf(fid,'%g %g %g',[3,inf]);

fclose(fid);

data1=data;
data=data';
[N,~]=size(data);

%  hypocentral locations, N=number of hypocenters
xs(1:N)=data(1:N,1);
ys(1:N)=data(1:N,2);
zs(1:N)=data(1:N,3);

%  Now plot the density plot on the GUI plot axes
h1=handles.density_plot_axes;
cla(h1,'reset')
set(h1,'NextPlot','add');

psi_angle=get(handles.Psi_slider,'Value');
elevation=get(handles.Elevation_slider,'Value');
azimuth=get(handles.Azimuth_slider,'Value');
zfactor=1;%get(handles.zoom_slider,'Value');

view_flag=0;

plot_density2(h1,psi_angle,elevation,azimuth,zfactor,view_flag);

view(3);
set(h1,'NextPlot','replace');


% --- Executes on button press in translate_pushbutton.
function translate_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to translate_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global PLOT_Latest x_poly y_poly x_orig y_orig n_orig DI D

%  Pick a new (x,y) relative origin for data

h1=handles.density_plot_axes;
set(h1,'Nextplot','add');

%  Pick a new origin point
[x_test,y_test]=ginput(1);

%  convert into original coordinates using the inverse rotation matrix
xm_old = [x_test ; y_test ; 0.0];
xm_new=DI*xm_old;

x_orig=xm_new(1) + x_orig;
y_orig=xm_new(2) + y_orig;

n_orig=n_orig+1;

%  Now plot the density plot on the GUI plot axes
cla(h1,'reset')
set(h1,'NextPlot','add');

psi_angle=get(handles.Psi_slider,'Value');%str2num(get(handles.rotate_about_x_edit,'String'));
elevation=get(handles.Elevation_slider,'Value');
azimuth=get(handles.Azimuth_slider,'Value');
zfactor=1;%get(handles.zoom_slider,'Value');

view_flag=0;

plot_density2(h1,psi_angle,elevation,azimuth,zfactor,view_flag);

set(h1,'NextPlot','replace');
rotate3d on
return

% --- Executes on button press in FM_file_browse.
function FM_file_browse_Callback(hObject, eventdata, handles)
% hObject    handle to FM_file_browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global FM_file

if get(handles.FM_file_ques,'Value') == 1
    
    %   Call uigetfile menu
    [filename,dirname]=uigetfile('*.*');
    if filename == 0; return;end 

    FM_file = strcat(dirname,filename);
else
    FM_file = '';
end


function plot_density_menu_Callback(hObject, eventdata, handles)
% hObject    handle to plot_density_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%  Plot the current density plot window in a separate figure
global h1 xs ys zs

h1=handles.density_plot_axes;
viewpt = h1.View;
cla(h1,'reset');
set(h1,'NextPlot','add');
    
plot3(h1,xs,ys,zs,'o','MarkerEdgeColor','k','MarkerFaceColor',[0 0.75 0.75]);

axis equal;
grid on;
%title('Input hypocenters');
xlabel('X km');
ylabel('Y km');
zlabel('Z km');
set(h1, 'fontsize', 18);
grid on;
view(3);

xlim([min(xs)-1 max(xs)+1])
ylim([min(ys)-1 max(ys)+1])
zlim([min(zs)-1 max(zs)+1])

view(h1,viewpt(1),viewpt(2));
rotate3d on

% --- Executes on button press in refresh_pushbutton.
function refresh_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to refresh_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%  refresh seismicity plot
global PLOT_Latest x_poly y_poly x_orig y_orig n_orig DI D
global xs ys zs xv yv zv 
global Kfaults Cluster_number picks isauto_OADC
 
h_table= handles.uitable1;
cla(h_table,'reset');

D=[ 1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0];
x_orig=0.0;
y_orig=0.0;

% Reset the plot handles
set(handles.Psi_slider,'Value',0);
set(handles.Psi_edit,'String','0');
set(handles.Elevation_slider,'Value',90);
set(handles.Elevation_edit,'String','90');
set(handles.Azimuth_slider,'Value',0);
set(handles.Azimuth_edit,'String','0');

%  Now plot the density plot on the GUI plot axes
h1=handles.density_plot_axes;
cla(h1,'reset')
set(h1,'NextPlot','add');

plot3(h1,xs,ys,zs,'o','MarkerEdgeColor','k','MarkerFaceColor',[0 0.75 0.75]);

if isauto_OADC == 1
    Cluster_number = Kfaults;
end 

if Cluster_number >= 1 

    flt_no = 1:Cluster_number;
    fColumn = num2cell(true(Cluster_number,1));

    picks(:,1) = fColumn;

    % Display fault parameters on uitable
    h_table= handles.uitable1;
    set(h_table, 'Data', picks);

    for k=flt_no
        fill3(h1,xv(k,1:4),yv(k,1:4),zv(k,1:4),'w','FaceAlpha',0.7,'FaceColor',[0.5 0.5 0.5]);
    end
end


hold off;

axis equal;
grid on;
%title('Input hypocenters');
xlabel('X km');
ylabel('Y km');
zlabel('Z km');
set(h1, 'fontsize', 18);
grid on;
view(3);
viewpt = h1.View;

xlim([min(xs)-1 max(xs)+1])
ylim([min(ys)-1 max(ys)+1])
zlim([min(zs)-1 max(zs)+1])

view(h1,viewpt(1),viewpt(2));
rotate3d on
set(h1,'NextPlot','replace');

% --- Executes on button press in pick_pushbutton.
function pick_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to pick_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%  Pick a polygon to delimit a fault seismicity cluster, subset the
%  seismicity, calculate a rectangular fault plane model using principal
%  components analysis

global data1 xc yc zc vec_plane xb_old yb_old zb_old xs ys zs N Nc
global PLOT_Latest x_poly y_poly x_orig y_orig n_orig DI D
global data2 xd yd zd
global xclus yclus zclus Cluster_number N_cluster L W Strike Dip xv yv zv lambda3
global L_man W_man Strike_man Dip_man xv_man yv_man zv_man vec_plane_man lambda3_man
global h1 h_table picks isauto_OADC ismanual_OADC Kfaults


% Display Warning to delete the results of manual picks
if isauto_OADC == 1
    opts.Interpreter = 'tex';
    opts.Default = 'Yes';
    reset_prog = questdlg('\bf \fontsize{15} Do you want to switch to manual pick? All auto picks will be lost!', ...
        'WARNING!', 'Yes','No thank you','Cancel',opts);

    if strcmp(reset_prog,'Yes')
        isauto_OADC = 0;
        
        h_table= handles.uitable1;
        cla(h_table,'reset');
        picks = [num2cell(false(1,1)) num2cell(zeros(1,4))];
 
        set(h_table, 'Data', picks);
        cla(h_table,'reset');
        
        %  Now plot the density plot on the GUI plot axes
        h1=handles.density_plot_axes;
        cla(h1,'reset')
        set(h1,'NextPlot','add');

        plot3(h1,xs,ys,zs,'o','MarkerEdgeColor','k','MarkerFaceColor',[0 0.75 0.75]);
        view(2);
        axis equal;
        grid on;
        %title('Input hypocenters');
        xlabel('X km');
        ylabel('Y km');
        zlabel('Z km');
        set(h1, 'fontsize', 18);
        grid on;
        
        xlim([min(xs)-1 max(xs)+1])
        ylim([min(ys)-1 max(ys)+1])
        zlim([min(zs)-1 max(zs)+1])

    end
end

% Proceed to manual cluster picking
if isauto_OADC == 0

    ismanual_OADC = 1;
    
    h1=handles.density_plot_axes;
    view(2);

    %  Pick the polygon
    [x_poly,y_poly]=getline(h1,'closed');

    %  subset the seismicity within the polygon
    % determine if a hypocenter is inside or outside the polygon
    IN=inpolygon(xd,yd,x_poly,y_poly);

    % subset the seismicity in this polygon
    % each transformed event is in the same order as the input data
    k=0;
    for kk=1:N
        if IN(kk) == 1
            k=k+1;
            xn(k)=xs(kk);
            yn(k)=ys(kk);
            zn(k)=zs(kk);
        end   
    end

    if isempty(Cluster_number)
        Cluster_number = 0;
        L_man = []; W_man = [];
        Strike_man = []; Dip_man = [];
        xv_man = []; yv_man = []; zv_man = [];
        vec_plane_man = []; lambda3_man = [];
    end

    %  save the clustered seismicity
    if k >= 0
        Cluster_number=Cluster_number+1;
        N_cluster(Cluster_number)=k;
        xclus(Cluster_number,1:N_cluster(Cluster_number))=xn(1:N_cluster(Cluster_number));
        yclus(Cluster_number,1:N_cluster(Cluster_number))=yn(1:N_cluster(Cluster_number));
        zclus(Cluster_number,1:N_cluster(Cluster_number))=zn(1:N_cluster(Cluster_number));  
    end

    %[Ln,Wn,Striken,Dipn,xvn,yvn,zvn,vec_planen,sum_lambda3qn] = New_recalcfault(Nt,xt,yt,zt,Kfaults)
    [Ln,Wn,Striken,Dipn,xvn,yvn,zvn,vec_planen,lambda3n] = ...
        New_recalcfault(N_cluster(Cluster_number),xclus(Cluster_number,1:N_cluster(Cluster_number)),...
        yclus(Cluster_number,1:N_cluster(Cluster_number)),zclus(Cluster_number,1:N_cluster(Cluster_number)),1);

    L_man = [L_man; Ln];
    W_man = [W_man; Wn];
    Strike_man = [Strike_man; Striken];
    Dip_man = [Dip_man; Dipn];
    xv_man = [xv_man; xvn];
    yv_man = [yv_man; yvn];
    zv_man = [zv_man; zvn];
    vec_plane_man = [vec_plane_man; vec_planen];
    lambda3_man = [lambda3_man; lambda3n];

    flt_no = 1:Cluster_number;
    plot_faults(h1,xs,ys,zs,xv_man,yv_man,zv_man,flt_no)

    view(2);

    % Displaying the faults on the uitable
    flt_details_man = [Strike_man  Dip_man L_man W_man mean(xv_man(:,:),2) mean(yv_man(:,:),2) ...
            mean(zv_man(:,:),2) xv_man(:,:) yv_man(:,:) zv_man(:,:)];

    fColumn_man = num2cell(true(Cluster_number,1));

    
    Kfaults = Cluster_number;
    xv = xv_man;
    yv = yv_man;
    zv = zv_man;
    
    L = L_man;
    W = W_man;
    Strike = Strike_man;
    Dip = Dip_man;
    lambda3 = lambda3_man;
    
    picks = [fColumn_man num2cell(flt_details_man(:,1:4))];

    h_table= handles.uitable1;
    set(h_table, 'Data', picks);
    rotate3d on
    set(h_table,'CellEditCallback',@onCellSelected);

end

function plot_faults_menu_Callback(hObject, eventdata, handles)
% hObject    handle to plot_faults_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%  Plot all seismicity and faults in a separate figure

global data1 xc yc zc vec_plane xb_old yb_old zb_old xs ys zs N Nc
global PLOT_Latest x_poly y_poly x_orig y_orig n_orig DI D
global data2 xd yd zd h1 isauto_OADC xt yt zt Nt picks
global xclus yclus zclus Cluster_number Kfaults N_cluster L W Strike Dip xv yv zv lambda3

% % plot a rendition of the data and the planes that fit the data

h1=handles.density_plot_axes;
viewpt = h1.View;
cla(h1,'reset');
set(h1,'NextPlot','add');
    
if isauto_OADC == 1
    Cluster_number = Kfaults;
    xclus = xt;
    yclus = yt;
    zclus = zt;
    N_cluster = Nt;
end 

plot3(h1,xs,ys,zs,'o','MarkerEdgeColor','k','MarkerFaceColor',[0 0.75 0.75]);

flt_no = 1:Cluster_number;
flt_no = flt_no(cell2mat(picks(:,1)));

for k=flt_no
    fill3(h1,xv(k,1:4),yv(k,1:4),zv(k,1:4),'w','FaceAlpha',0.7,'FaceColor',[0.5 0.5 0.5]);
end

hold off;

axis equal;
grid on;
%title('Input hypocenters');
xlabel('X km');
ylabel('Y km');
zlabel('Z km');
set(h1, 'fontsize', 18);
grid on;
view(3);

xlim([min(xs)-1 max(xs)+1])
ylim([min(ys)-1 max(ys)+1])
zlim([min(zs)-1 max(zs)+1])
rotate3d on
view(h1,viewpt(1),viewpt(2));
 
function plot_faults_only_menu_Callback(hObject, eventdata, handles)
% hObject    handle to plot_faults_only_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%  Plot all faults in a separate figure

global data1 xc yc zc vec_plane xb_old yb_old zb_old xs ys zs N Nc
global PLOT_Latest x_poly y_poly x_orig y_orig n_orig DI D
global data2 xd yd zd h1 isauto_OADC xt yt zt Nt picks
global xclus yclus zclus Cluster_number Kfaults N_cluster L W Strike Dip xv yv zv lambda3

h1=handles.density_plot_axes;
viewpt = h1.View;
cla(h1,'reset');
set(h1,'NextPlot','add');
    
if isauto_OADC == 1
    Cluster_number = Kfaults;
    xclus = xt;
    yclus = yt;
    zclus = zt;
    N_cluster = Nt;
end 

flt_no = 1:Cluster_number;
flt_no = flt_no(cell2mat(picks(:,1)));

for k=flt_no
    fill3(h1,xv(k,1:4),yv(k,1:4),zv(k,1:4),'w','FaceAlpha',0.7,'FaceColor',[0.5 0.5 0.5]);
end

hold off;

axis equal;
grid on;
%title('Input hypocenters');
xlabel('X km');
ylabel('Y km');
zlabel('Z km');
set(h1, 'fontsize', 18);
grid on;
view(3);

xlim([min(xs)-1 max(xs)+1])
ylim([min(ys)-1 max(ys)+1])
zlim([min(zs)-1 max(zs)+1])

rotate3d on

view(h1,viewpt(1),viewpt(2));

function plot_clusters_faults_menu_Callback(hObject, eventdata, handles)
% hObject    handle to plot_clusters_faults_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%  Plot cluster seismicity with faults in a separate figure

global data1 xc yc zc vec_plane xb_old yb_old zb_old xs ys zs N Nc
global PLOT_Latest x_poly y_poly x_orig y_orig n_orig DI D
global data2 xd yd zd h1 isauto_OADC  xt yt zt picks Nt
global xclus yclus zclus Cluster_number Kfaults N_cluster L W Strike Dip xv yv zv lambda3

h1=handles.density_plot_axes;
viewpt = h1.View;
cla(h1,'reset');
set(h1,'NextPlot','add');
  
if isauto_OADC == 1
    Cluster_number = Kfaults;
    xclus = xt;
    yclus = yt;
    zclus = zt;
    N_cluster = Nt;
end 

flt_no = 1:Cluster_number;
flt_no = flt_no(cell2mat(picks(:,1)));
           
for k=flt_no
   
    kc=N_cluster(k);
    plot3(h1,xclus(k,1:kc),yclus(k,1:kc),zclus(k,1:kc),'o','MarkerEdgeColor','k','MarkerFaceColor','k');
    fill3(h1,xv(k,1:4),yv(k,1:4),zv(k,1:4),'w','FaceAlpha',0.7,'FaceColor',[0.5 0.5 0.5]);
end

hold off;

axis equal;
grid on;
%title('Input hypocenters');
xlabel('X km');
ylabel('Y km');
zlabel('Z km');
set(h1, 'fontsize', 18);
grid on;
view(3);

xlim([min(xs)-1 max(xs)+1])
ylim([min(ys)-1 max(ys)+1])
zlim([min(zs)-1 max(zs)+1])
rotate3d on
view(h1,viewpt(1),viewpt(2));

function write_file_menu_Callback(hObject, eventdata, handles)
% hObject    handle to write_file_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%  Write the clustered seismicity and fault plane models to a file

global xclus yclus zclus Cluster_number Kfaults N_cluster L W Strike Dip xv yv zv lambda3
global isauto_OADC  xt yt zt Nt picks
global infile arg_simul_tag

global i arg_kmin arg_kmax arg_dip_threshold arg_clus_mineqs 
global arg_is3D_simul arg_err_av arg_FM_file
    
    
%  uiputfile dialog

[filename,dirname]=uiputfile('*.*');
if filename == 0; return;end
    
faultfile=strcat(dirname,filename);

fid=fopen(faultfile,'w');

fprintf(fid,'%s\n','Cluster Model File');
fprintf(fid,'%s%s\n','Data: ', infile);
fprintf(fid,'%s%s\n','Simultag: ', arg_simul_tag);
fprintf(fid,'\n');

if isauto_OADC == 1
    fprintf(fid,'%s\n','Automatic Fault Detection');
    fprintf(fid,'\n');
    fprintf(fid,'%s\n','Important User Settings');
    fprintf(fid,'%s%s\n','Rnd_seed: ', num2str(i));
    fprintf(fid,'%s%s\n','kmin: ', num2str(arg_kmin));
    fprintf(fid,'%s%s\n','kmax: ', num2str(arg_kmax));
    fprintf(fid,'%s%s\n','Dip threshold: ', num2str(arg_dip_threshold));
    fprintf(fid,'%s%s\n','Min eqs per cluster: ', num2str(arg_clus_mineqs));
    fprintf(fid,'%s%s\n','is_3D_simul: ', num2str(arg_is3D_simul));
    fprintf(fid,'%s%s\n','err_av: ', num2str(arg_err_av));
    fprintf(fid,'%s%s\n','FM_file: ', arg_FM_file);
    fprintf(fid,'\n');
else
    fprintf(fid,'%s\n','Manual Fault Picking:');
    fprintf(fid,'\n');
end


if isauto_OADC == 1
    Cluster_number = Kfaults;
    xclus = xt;
    yclus = yt;
    zclus = zt;
    N_cluster = Nt;
end 

flt_no = 1:Cluster_number;
flt_no = flt_no(cell2mat(picks(:,1)));

fprintf(fid,'%s%s\n','No of Faults: ', num2str(Cluster_number));

if Cluster_number ~= length(flt_no)
    fprintf(fid,'%s%s\n','But printing only faults: ', num2str(flt_no));
end

fprintf(fid,'\n');
fprintf(fid,'%s\n','Legend:');
fprintf(fid,'%s\n','For each fault');
fprintf(fid,'%s\n','           -  Fault number');
fprintf(fid,'%s\n','           -  Length Width Strike Dip Thickness');
fprintf(fid,'%s\n','           -  x coord of the fault plane corners');
fprintf(fid,'%s\n','           -  y coord of the fault plane corners');
fprintf(fid,'%s\n','           -  z coord of the fault plane corners');
fprintf(fid,'%s\n','           -  Number of hypocenters in the cluster');
fprintf(fid,'%s\n','           -  xt yt zt - coordinate of the hypocenters in the cluster');
      
fprintf(fid,'\n');
fprintf(fid,'\n');

flt_no = 1:Cluster_number;
flt_no = flt_no(cell2mat(picks(:,1)));

for k=flt_no
    
    %  cluster index
    fprintf(fid,'%i \n',k);
    
    %  fault geometry
    fprintf(fid,'%12.5g %12.5g %12.5g %12.5g %12.5g \n',[L(k) W(k) Strike(k) Dip(k) lambda3(k)]);
    
    %  rectangular fault vertices
    fprintf(fid,'%12.5g %12.5g %12.5g %12.5g \n',xv(k,1:4));
    fprintf(fid,'%12.5g %12.5g %12.5g %12.5g \n',yv(k,1:4));
    fprintf(fid,'%12.5g %12.5g %12.5g %12.5g \n',zv(k,1:4));
    
    %  cluster seismicity
    n=N_cluster(k);
    fprintf(fid,'%i \n',n);
    
    for kk=1:n
        fprintf(fid,'%12.5g %12.5g %12.5g \n',[xclus(k,kk) yclus(k,kk) zclus(k,kk)]);
    end
    
end

fclose(fid);
  


% --------------------------------------------------------------------
function Generate_synthetic_dataset_Callback(hObject, eventdata, handles)
% hObject    handle to Generate_synthetic_dataset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Generate_Synthetic_Dataset


% --- Executes on button press in Fault_Clustering_push_button.
function Fault_Clustering_push_button_Callback(hObject, eventdata, handles)
% hObject    handle to Fault_Clustering_push_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global h_table picks h1 Kfaults FM_file infile arg_simul_tag
global xs ys zs xv yv zv isauto_OADC ismanual_OADC xt yt zt Nt lambda3
global Strike Dip L W 

global i arg_kmin arg_kmax arg_dip_threshold arg_clus_mineqs 
global arg_is3D_simul arg_err_av arg_FM_file
    
% Display Warning to delete the results of manual picks
if ismanual_OADC == 1
    opts.Interpreter = 'tex';
    opts.Default = 'Yes';
    reset_prog = questdlg('\bf \fontsize{15} Do you want to auto pick? All manual picks will be lost!', ...
        'WARNING!', 'Yes','No thank you','Cancel',opts);

    if strcmp(reset_prog,'Yes')
        ismanual_OADC = 0;
    end
end
 
if ismanual_OADC == 0
    
    isauto_OADC = 1;
    
    % Default inputs
    arg_N_loop = 1;  
    arg_PLOT_FLAG0 = 0;   % =0, no plots at all
    arg_PLOT_FLAG1 = 0;   % =0, no intermediate loop plots of data and planes
    arg_comb_coplanar = 1; % =0, Do not check or combine coplananr faults
    arg_plot_avg_FM = 0; % =1, Plot intermediate average FM focalsphere

    % User inputs
    arg_hypo_infile = infile; %'testdata.txt';%'Simul_hypos.txt';
    arg_FM_file = FM_file; %'';%'All_write_mechanisms.txt'; 
    
    arg_dip_threshold = str2double(get(handles.Dip_th_edit,'String'));
    arg_clus_mineqs = str2double(get(handles.mineqs_clus_edit,'String')); 
    arg_is3D_simul = get(handles.is_3D_simul_edit,'Value');
    i = str2double(get(handles.md_seed_edit,'String'));
    arg_kmin = str2double(get(handles.kmin_edit,'String'));
    arg_kmax = str2double(get(handles.kmax_edit,'String'));
    arg_err_av = str2double(get(handles.max_th_edit,'String'));
    arg_simul_tag = get(handles.Simultag_edit,'String');

    if arg_is3D_simul == 1
    %     % Run OADC_3D

        OADC_3D_now_12_4_20(i,arg_kmin,arg_kmax,arg_err_av,...
            arg_hypo_infile,arg_clus_mineqs,arg_N_loop,arg_simul_tag,...
            arg_dip_threshold,arg_FM_file,arg_PLOT_FLAG0,arg_PLOT_FLAG1,...
            arg_comb_coplanar, arg_plot_avg_FM);

        % ------- remove OADC intermediate figures --------
        Figures = findobj( 'Type', 'Figure' , '-not' , 'Tag' , get( handles.output , 'Tag' ) );
        NFigures = length( Figures );

        for nFigures = 1 : NFigures
          close(Figures(nFigures));
        end
        %--------------------------------

        if arg_comb_coplanar == 0
            data=load([pwd '/' arg_simul_tag '_results/' arg_simul_tag '.saved_variables.mat']);
        else
            data=load([pwd '/' arg_simul_tag '_results/' arg_simul_tag '.saved_variables_coplanar.mat']);
        end

        Kfaults=data.Kfaults;
        xs=data.xs; ys=data.ys; zs=data.zs;
        xv=data.xv(1:Kfaults,:); yv=data.yv(1:Kfaults,:); zv=data.zv(1:Kfaults,:);
        xt=data.xt; yt=data.yt; zt=data.zt;
        Nt=data.Nt(1:Kfaults);
        Strike=data.Strike(1:Kfaults); Dip=data.Dip(1:Kfaults);
        L=data.L(1:Kfaults); W=data.W(1:Kfaults); 
        lambda3=data.lambda3(1:Kfaults);

        % generate fault details matrix for display
        flt_details = [Strike'  Dip' L' W' mean(xv(:,:),2) mean(yv(:,:),2) ...
            mean(zv(:,:),2) xv(:,:) yv(:,:) zv(:,:)];

        h1=handles.density_plot_axes;
        cla(h1,'reset');
        set(h1,'NextPlot','add');

        flt_no = 1:Kfaults;
        plot_faults(h1,xs,ys,zs,xv,yv,zv,flt_no)
        view(3);
        
        fColumn = num2cell(true(Kfaults,1));

        picks = [fColumn num2cell(flt_details(:,1:4))];

        % Display fault parameters on uitable
        h_table= handles.uitable1;
        set(h_table, 'Data', picks);
        rotate3d on
        set(h_table,'CellEditCallback',@onCellSelected);

    else
        % Run OADC_2D

        
        
    end
end

% --------------------------------------------------------------------
function Cluster_seismicity_dataset_Callback(hObject, eventdata, handles)
% hObject    handle to Cluster_seismicity_dataset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Cluster_Seismicity_Dataset






% ---------------- TO DO EDIT BELOW THIS LINE -----------------------------

function Dip_th_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Dip_th_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Dip_th_edit as text
%        str2double(get(hObject,'String')) returns contents of Dip_th_edit as a double


% --- Executes during object creation, after setting all properties.
function Dip_th_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Dip_th_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mineqs_clus_edit_Callback(hObject, eventdata, handles)
% hObject    handle to mineqs_clus_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mineqs_clus_edit as text
%        str2double(get(hObject,'String')) returns contents of mineqs_clus_edit as a double


% --- Executes during object creation, after setting all properties.
function mineqs_clus_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mineqs_clus_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function write_clusterredfile_menu_Callback(hObject, eventdata, handles)
% hObject    handle to write_clusterredfile_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function file_menu_Callback(hObject, eventdata, handles)
% hObject    handle to file_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function plot_menu_Callback(hObject, eventdata, handles)
% hObject    handle to plot_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function Simultag_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Simultag_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
% Hints: get(hObject,'String') returns contents of Simultag_edit as text
%        str2double(get(hObject,'String')) returns contents of Simultag_edit as a double

% --- Executes during object creation, after setting all properties.
function Simultag_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Simultag_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
%
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in is_2D_simul_edit.
function is_2D_simul_edit_Callback(hObject, eventdata, handles)
% hObject    handle to is_2D_simul_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
% Hint: get(hObject,'Value') returns toggle state of is_2D_simul_edit

% --- Executes on button press in is_3D_simul_edit.
function is_3D_simul_edit_Callback(hObject, eventdata, handles)
% hObject    handle to is_3D_simul_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
% Hint: get(hObject,'Value') returns toggle state of is_3D_simul_edit

% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
% Hint: get(hObject,'Value') returns toggle state of radiobutton3

% --- Executes on button press in radiobutton4.
function radiobutton4_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of radiobutton4

function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double

% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double

% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double

% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function kmin_edit_Callback(hObject, eventdata, handles)
% hObject    handle to kmin_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
% Hints: get(hObject,'String') returns contents of kmin_edit as text
%        str2double(get(hObject,'String')) returns contents of kmin_edit as a double

% --- Executes during object creation, after setting all properties.
function kmin_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kmin_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function kmax_edit_Callback(hObject, eventdata, handles)
% hObject    handle to kmax_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
% Hints: get(hObject,'String') returns contents of kmax_edit as text
%        str2double(get(hObject,'String')) returns contents of kmax_edit as a double

% --- Executes during object creation, after setting all properties.
function kmax_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kmax_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function max_th_edit_Callback(hObject, eventdata, handles)
% hObject    handle to max_th_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
% Hints: get(hObject,'String') returns contents of max_th_edit as text
%        str2double(get(hObject,'String')) returns contents of max_th_edit as a double

% --- Executes during object creation, after setting all properties.
function max_th_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to max_th_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function md_seed_edit_Callback(hObject, eventdata, handles)
% hObject    handle to md_seed_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
% Hints: get(hObject,'String') returns contents of md_seed_edit as text
%        str2double(get(hObject,'String')) returns contents of md_seed_edit as a double

% --- Executes during object creation, after setting all properties.
function md_seed_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to md_seed_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in FM_file_ques.
function FM_file_ques_Callback(hObject, eventdata, handles)
% hObject    handle to FM_file_ques (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
% Hint: get(hObject,'Value') returns toggle state of FM_file_ques



% --- Executes on button press in FM_file_ques2.
function FM_file_ques2_Callback(hObject, eventdata, handles)
% hObject    handle to FM_file_ques2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
% Hint: get(hObject,'Value') returns toggle state of FM_file_ques2


% --- Executes during object creation, after setting all properties.
function Azimuth_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Azimuth_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function Elevation_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Elevation_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes during object creation, after setting all properties.
function Azimuth_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Azimuth_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes during object creation, after setting all properties.
function Psi_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Psi_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function Psi_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Psi_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes during object creation, after setting all properties.
function Elevation_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Elevation_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
