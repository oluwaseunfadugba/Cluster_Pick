function varargout = Cluster_Seismicity_Dataset(varargin)
% CLUSTER_SEISMICITY_DATASET MATLAB code for Cluster_Seismicity_Dataset.fig
%      CLUSTER_SEISMICITY_DATASET, by itself, creates a new CLUSTER_SEISMICITY_DATASET or raises the existing
%      singleton*.
%
%      H = CLUSTER_SEISMICITY_DATASET returns the handle to a new CLUSTER_SEISMICITY_DATASET or the handle to
%      the existing singleton*.
%
%      CLUSTER_SEISMICITY_DATASET('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CLUSTER_SEISMICITY_DATASET.M with the given input arguments.
%
%      CLUSTER_SEISMICITY_DATASET('Property','Value',...) creates a new CLUSTER_SEISMICITY_DATASET or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Cluster_Seismicity_Dataset_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Cluster_Seismicity_Dataset_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Cluster_Seismicity_Dataset

% Last Modified by GUIDE v2.5 15-Sep-2022 12:48:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Cluster_Seismicity_Dataset_OpeningFcn, ...
                   'gui_OutputFcn',  @Cluster_Seismicity_Dataset_OutputFcn, ...
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


% --- Executes just before Cluster_Seismicity_Dataset is made visible.
function Cluster_Seismicity_Dataset_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Cluster_Seismicity_Dataset (see VARARGIN)

global perc_value is_declus is_3D_declus is_2D_declus
global xs_decl ys_decl zs_decl 

% Choose default command line output for Cluster_Seismicity_Dataset
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Setting default values
perc_value=5.0;
set(handles.perc_slider,'Value',perc_value);
set(handles.perc_edit,'String',num2str(perc_value));

is_declus = 0;

cla(handles.axes1);view(handles.axes1,2);
cla(handles.axes2);view(handles.axes1,2);
cla(handles.axes3);view(handles.axes1,2);
cla(handles.axes4);view(handles.axes1,2);

set(handles.is_3D_declus,'Value',1); is_3D_declus = 1;
set(handles.is_2D_declus,'Value',0); is_2D_declus = 0;

xs_decl = [];
ys_decl = [];
zs_decl = [];

% UIWAIT makes Cluster_Seismicity_Dataset wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Cluster_Seismicity_Dataset_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --------------------------------------------------------------------
% --- Executes on slider movement.
function perc_slider_Callback(hObject, eventdata, handles)
% hObject    handle to perc_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

global is_declus 
global is_3D_declus is_2D_declus

perc_value=get(handles.perc_slider,'Value');
set(handles.perc_edit,'String',num2str(round(perc_value,1)));

if is_declus == 1
    if is_3D_declus == 1
        perform_clustering_edit();
    elseif is_2D_declus == 1
        perform_clustering_edit_2D();
    end
end

% --------------------------------------------------------------------
function perc_edit_Callback(hObject, eventdata, handles)
% hObject    handle to perc_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of perc_edit as text
%        str2double(get(hObject,'String')) returns contents of perc_edit as a double

global is_declus 
global is_3D_declus is_2D_declus

perc_value=str2double(get(handles.perc_edit,'String'));
set(handles.perc_slider,'Value',round(perc_value,1));

if is_declus == 1
    if is_3D_declus == 1
        perform_clustering_edit();
    elseif is_2D_declus == 1
        perform_clustering_edit_2D();
    end
end

% --------------------------------------------------------------------
% --- Executes on button press in default_perc.
function default_perc_Callback(hObject, eventdata, handles)
% hObject    handle to default_perc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global is_declus 
global is_3D_declus is_2D_declus

perc_value=5.0;
set(handles.perc_slider,'Value',perc_value);
set(handles.perc_edit,'String',num2str(perc_value));

if is_declus == 1
    if is_3D_declus == 1
        perform_clustering_edit();
    elseif is_2D_declus == 1
        perform_clustering_edit_2D();
    end
end

% --------------------------------------------------------------------
function Read_Seis_Data_Callback(hObject, eventdata, handles)
% hObject    handle to Read_Seis_Data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global xs_decl ys_decl zs_decl h1_decl
global infilename_decl
global is_3D_declus is_2D_declus

h1_decl=handles.axes1;
cla(h1_decl,'reset');

%  Read a seismicity file consisting of (x,y,z), depth is negative
%   Call uigetfile menu
[infilename_decl,dirname]=uigetfile('*.*');
if infilename_decl == 0; return;end 

infile_decl=strcat(dirname,infilename_decl);
fid=fopen(infile_decl,'r');

if is_3D_declus == 1

    [data,~]=fscanf(fid,'%g %g %g',[3,inf]);

    fclose(fid);

    data=data';
    [N,~]=size(data);

    %  hypocentral locations, N=number of hypocenters
    xs_decl(1:N)=data(1:N,1);
    ys_decl(1:N)=data(1:N,2);
    zs_decl(1:N)=data(1:N,3);

    %  Now plot the density plot on the GUI plot axes
    h1_decl=handles.axes1;
    cla(h1_decl,'reset')
    set(h1_decl,'NextPlot','add');

    scatter3(h1_decl,xs_decl,ys_decl,zs_decl,'filled','MarkerEdgeColor','k',...
            'MarkerFaceColor',[0 .75 .75])
    axis(h1_decl,'equal');
    grid(h1_decl,'on');
    xlabel(h1_decl,'X km');
    ylabel(h1_decl,'Y km');
    zlabel(h1_decl,'Z km');

    xx = round((max(xs_decl)-min(xs_decl))/5);
    yy = round((max(ys_decl)-min(ys_decl))/5);
    zz = round((max(zs_decl)-min(zs_decl))/5);
    
    xlim(h1_decl,[min(xs_decl)-xx max(xs_decl)+xx])
    ylim(h1_decl,[min(ys_decl)-yy max(ys_decl)+yy])
    zlim(h1_decl,[min(zs_decl)-zz max(zs_decl)+zz])

    view(h1_decl,3);hold on;
    
elseif is_2D_declus == 1
    [data,~]=fscanf(fid,'%g %g',[2,inf]);

    fclose(fid);

    data=data';
    [N,~]=size(data);

    %  hypocentral locations, N=number of hypocenters
    xs_decl(1:N)=data(1:N,1);
    ys_decl(1:N)=data(1:N,2);

    %  Now plot the density plot on the GUI plot axes
    h1_decl=handles.axes1;
    cla(h1_decl,'reset')
    set(h1_decl,'NextPlot','add');

    scatter(h1_decl,xs_decl,ys_decl,'filled','MarkerEdgeColor','k',...
            'MarkerFaceColor',[0 .75 .75])
    axis(h1_decl,'equal');
    grid(h1_decl,'on');
    xlabel(h1_decl,'X km');
    ylabel(h1_decl,'Y km');

    xx = round((max(xs_decl)-min(xs_decl))/5);
    yy = round((max(ys_decl)-min(ys_decl))/5);
    
    xlim(h1_decl,[min(xs_decl)-xx max(xs_decl)+xx])
    ylim(h1_decl,[min(ys_decl)-yy max(ys_decl)+yy])

    view(h1_decl,2);hold on;
end

set(h1_decl,'NextPlot','add');

% --------------------------------------------------------------------
function Pick_Random_Dataset_Polygon_Callback(hObject, eventdata, handles)
% hObject    handle to Pick_Random_Dataset_Polygon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global xs_decl ys_decl zs_decl 
global xs_now ys_now zs_now
global xs_sort ys_sort zs_sort vol_sort
global xr yr zr
global nkh cdf_rnd vr
global nks vs cdf_cat
global h1_decl h2_decl h3_decl h4_decl h_perc_edit
global is_declus 
global is_3D_declus is_2D_declus

h1_decl=handles.axes1;
view(h1_decl,2);

%  Pick the polygon
[x_poly_decl,y_poly_decl]=getline(h1_decl,'closed'); hold on;

cla(h1_decl,'reset')
set(h1_decl,'NextPlot','add');

if is_3D_declus == 1
    scatter3(h1_decl,xs_decl,ys_decl,zs_decl,'filled','MarkerEdgeColor','k',...
        'MarkerFaceColor',[0 .75 .75]); hold on
elseif is_2D_declus == 1
    scatter(h1_decl,xs_decl,ys_decl,'filled','MarkerEdgeColor','k',...
        'MarkerFaceColor',[0 .75 .75]); hold on
end   
    
plot(h1_decl,x_poly_decl,y_poly_decl,'r-'); hold on;

view(h1_decl,2);hold on;

axis(h1_decl,'equal');
grid(h1_decl,'on');
xlabel(h1_decl,'X km');
ylabel(h1_decl,'Y km');

if is_3D_declus == 1
    zlabel(h1_decl,'Z km');
end
% 
xx = round((max(xs_decl)-min(xs_decl))/5);
yy = round((max(ys_decl)-min(ys_decl))/5);
    
xlim(h1_decl,[min(xs_decl)-xx max(xs_decl)+xx])
ylim(h1_decl,[min(ys_decl)-yy max(ys_decl)+yy])
grid(h1_decl,'on');
view(h1_decl,2)

% Generating ramdomized catalog
if is_3D_declus == 1
    rnd_in_z = 1; % Randomized the hypocenters in depth as well
    [xr,yr,zr]=randomcat(xs_decl,ys_decl,zs_decl,x_poly_decl,y_poly_decl,rnd_in_z);

    % Remove earthquakes outside the polygon
    IN=inpolygon(xs_decl,ys_decl,x_poly_decl,y_poly_decl);
    xs_now = xs_decl(IN);
    ys_now = ys_decl(IN);
    zs_now = zs_decl(IN);

    % determine tetrahedra volumes
    [xs_sort, ys_sort, zs_sort, vol_sort] = calc_tetvol(xs_now,ys_now,zs_now);
    [~, ~, ~, volrand_sort] = calc_tetvol(xr',yr',zr');

    %  Compute the normalized cumulative probability density function of the
    %  random and original catalog tetrahedra volumes
    [cdf_rnd,vr]=ecdf(volrand_sort); nkh=length(vr);
    [cdf_cat,vs]=ecdf(vol_sort); nks=length(vs);

elseif is_2D_declus == 1 
    [xr,yr]=randomcat_2D(xs_decl,ys_decl,x_poly_decl,y_poly_decl);

    % Remove earthquakes outside the polygon
    IN=inpolygon(xs_decl,ys_decl,x_poly_decl,y_poly_decl);
    xs_now = xs_decl(IN);
    ys_now = ys_decl(IN);

    % determine tetrahedra volumes
    [xs_sort, ys_sort, vol_sort] = calc_tetvol_2D(xs_now,ys_now);
    [~, ~, volrand_sort] = calc_tetvol_2D(xr',yr');

    %  Compute the normalized cumulative probability density function of the
    %  random and original catalog tetrahedra volumes
    [cdf_rnd,vr]=ecdf(volrand_sort); nkh=length(vr);
    [cdf_cat,vs]=ecdf(vol_sort); nks=length(vs);
    
end   
% 
% get necessary handles
h1_decl=handles.axes1;
h2_decl=handles.axes2;
h3_decl=handles.axes3;
h4_decl=handles.axes4;
h_perc_edit=handles.perc_edit;

% perform declustering

if is_3D_declus == 1
    perform_clustering_edit();
elseif is_2D_declus == 1
    perform_clustering_edit_2D();
end

is_declus = 1;

% --------------------------------------------------------------------
function Write_Clusterred_Data_Callback(hObject, eventdata, handles)
% hObject    handle to Write_Clusterred_Data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global xs_clu ys_clu zs_clu
global infilename_decl h_perc_edit
global V05 NV05_cat
global xs_diffuse ys_diffuse zs_diffuse
global is_3D_declus is_2D_declus

%  uiputfile dialog

[outfilename_decl,dirname]=uiputfile('*.*');
if outfilename_decl == 0; return;end
    
outputfile_decl=strcat(dirname,outfilename_decl);

fid=fopen(outputfile_decl,'w');
fprintf(fid,'%s\n','Decluster Result File');
fprintf(fid,'%s%s\n','Data: ', infilename_decl);
fprintf(fid,'%s%s\n','Simultag: ', outfilename_decl(1:end-4));
fprintf(fid,'\n');

perc = get(h_perc_edit,'String');

if is_3D_declus == 1
    fprintf(fid,'%s%s%s%s\n','Volume at ',perc,' percent probability (km cubed) = ',num2str(V05)');
    fprintf(fid,'%s%s\n','Probability at target volume (e.g., NV05) = ',num2str(NV05_cat)');
    fprintf(fid,'\n');
    fprintf(fid,'%s%s\n','No of clusterred events: ', num2str(length(xs_clu)));
    fprintf(fid,'\n');
    fprintf(fid,'%s\n','x_coord     y_coord     z_coord');
    fprintf(fid,'\n');
    for kk=1:length(xs_clu)
        fprintf(fid,'%12.5f %12.5f %12.5f\n',[xs_clu(kk) ys_clu(kk) zs_clu(kk)]);
    end
    fprintf(fid,'\n');
    fprintf(fid,'\n');

    % Writing diffuse events
    fprintf(fid,'%s%s\n','No of diffuse events: ', num2str(length(xs_diffuse)));
    fprintf(fid,'\n');
    fprintf(fid,'%s\n','x_coord     y_coord     z_coord');
    fprintf(fid,'\n');
    for kk=1:length(xs_diffuse)
        fprintf(fid,'%12.5f %12.5f %12.5f\n',[xs_diffuse(kk) ys_diffuse(kk) zs_diffuse(kk)]);
    end

elseif is_2D_declus == 1
    
    fprintf(fid,'%s%s%s%s\n','Area at ',perc,' percent probability (km squared) = ',num2str(V05)');
    fprintf(fid,'%s%s\n','Probability at target area (e.g., NA05) = ',num2str(NV05_cat)');
    fprintf(fid,'\n');
    fprintf(fid,'%s%s\n','No of clusterred events: ', num2str(length(xs_clu)));
    fprintf(fid,'\n');
    fprintf(fid,'%s\n','x_coord     y_coord ');
    fprintf(fid,'\n');
    for kk=1:length(xs_clu)
        fprintf(fid,'%12.5f %12.5f %12.5f\n',[xs_clu(kk) ys_clu(kk)]);
    end
    fprintf(fid,'\n');
    fprintf(fid,'\n');

    % Writing diffuse events
    fprintf(fid,'%s%s\n','No of diffuse events: ', num2str(length(xs_diffuse)));
    fprintf(fid,'\n');
    fprintf(fid,'%s\n','x_coord     y_coord ');
    fprintf(fid,'\n');
    for kk=1:length(xs_diffuse)
        fprintf(fid,'%12.5f %12.5f %12.5f\n',[xs_diffuse(kk) ys_diffuse(kk)]);
    end
    
end

% Printing figure to file
F1 = getframe(handles.axes1); imwrite(F1.cdata, [outfilename_decl(1:end-4) '_Input_Data.png'], 'png')
F2 = getframe(handles.axes2); imwrite(F2.cdata, [outfilename_decl(1:end-4) '_Clustered_Data.png'], 'png')
F3 = getframe(handles.axes3); imwrite(F3.cdata, [outfilename_decl(1:end-4) '_CDF.png'], 'png')
F4 = getframe(handles.axes4); imwrite(F4.cdata, [outfilename_decl(1:end-4) '_Input_ClusteredData.png'], 'png')

savefig(gcf,[outfilename_decl(1:end-4) '.fig'])

% Move all figures to a folder name with the output filename
save([outfilename_decl(1:end-4) '_workspace'])

eval(sprintf('%s%s%s','! mkdir ',outfilename_decl(1:end-4),'_results'))
eval(sprintf('%s%s%s %s%s','! mv ',outfilename_decl(1:end-4), '*',outfilename_decl(1:end-4),'_results'))


% --- Executes on button press in is_2D_declus.
function is_2D_declus_Callback(hObject, eventdata, handles)
% hObject    handle to is_2D_declus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of is_2D_declus

global is_3D_declus is_2D_declus xs_decl ys_decl zs_decl

if is_3D_declus == 1 && length(xs_decl)>1
 
    % Display Warning to delete the results of manual picks
    opts.Interpreter = 'tex';
    opts.Default = 'Yes';
    reset_prog = questdlg('\bf \fontsize{15} Do you want to switch to 2D? All 3D data will be lost!', ...
        'WARNING!', 'Yes','No thank you','Cancel',opts);

    if strcmp(reset_prog,'Yes')
       	is_3D_declus = 0;
        is_2D_declus = get(hObject,'Value');

        xs_decl = [];
        ys_decl = [];
        zs_decl = [];

        cla(handles.axes1);view(handles.axes1,2);
        cla(handles.axes2);view(handles.axes1,2);
        cla(handles.axes3);view(handles.axes1,2);
        cla(handles.axes4);view(handles.axes1,2);

    end
    
else
    set(handles.is_3D_declus,'Value',0); is_3D_declus = 0;
    set(handles.is_2D_declus,'Value',1); is_2D_declus = 1;

    cla(handles.axes1);view(handles.axes1,2);
    cla(handles.axes2);view(handles.axes1,2);
    cla(handles.axes3);view(handles.axes1,2);
    cla(handles.axes4);view(handles.axes1,2);
end
  

% --- Executes on button press in is_3D_declus.
function is_3D_declus_Callback(hObject, eventdata, handles)
% hObject    handle to is_3D_declus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of is_3D_declus

global is_3D_declus is_2D_declus xs_decl ys_decl zs_decl

if is_2D_declus == 1 && length(xs_decl)>1
 
    % Display Warning to delete the results of manual picks
    opts.Interpreter = 'tex';
    opts.Default = 'Yes';
    reset_prog = questdlg('\bf \fontsize{15} Do you want to switch to 3D? All 2D data will be lost!', ...
        'WARNING!', 'Yes','No thank you','Cancel',opts);

    if strcmp(reset_prog,'Yes')
        is_3D_declus = get(hObject,'Value');
        is_2D_declus = 0;

        xs_decl = [];
        ys_decl = [];
        zs_decl = [];
        
        cla(handles.axes1);view(handles.axes1,2);
        cla(handles.axes2);view(handles.axes1,2);
        cla(handles.axes3);view(handles.axes1,2);
        cla(handles.axes4);view(handles.axes1,2);
    end
    
else
    set(handles.is_3D_declus,'Value',1); is_3D_declus = 1;
    set(handles.is_2D_declus,'Value',0); is_2D_declus = 0;
    
    cla(handles.axes1);view(handles.axes1,2);
    cla(handles.axes2);view(handles.axes1,2);
    cla(handles.axes3);view(handles.axes1,2);
    cla(handles.axes4);view(handles.axes1,2);

end
  

% -------------------------------------------------------------------------
% --------------- DO NOT EDIT BELOW THIS LINE -----------------------------
% -------------------------------------------------------------------------

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function perc_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to perc_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function perc_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to perc_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
