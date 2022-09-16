function varargout = Generate_Synthetic_Dataset(varargin)
% GENERATE_SYNTHETIC_DATASET MATLAB code for Generate_Synthetic_Dataset.fig
%      GENERATE_SYNTHETIC_DATASET, by itself, creates a new GENERATE_SYNTHETIC_DATASET or raises the existing
%      singleton*.
%
%      H = GENERATE_SYNTHETIC_DATASET returns the handle to a new GENERATE_SYNTHETIC_DATASET or the handle to
%      the existing singleton*.
%
%      GENERATE_SYNTHETIC_DATASET('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GENERATE_SYNTHETIC_DATASET.M with the given input arguments.
%
%      GENERATE_SYNTHETIC_DATASET('Property','Value',...) creates a new GENERATE_SYNTHETIC_DATASET or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Generate_Synthetic_Dataset_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Generate_Synthetic_Dataset_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Generate_Synthetic_Dataset

% Last Modified by GUIDE v2.5 13-Jun-2021 23:01:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Generate_Synthetic_Dataset_OpeningFcn, ...
                   'gui_OutputFcn',  @Generate_Synthetic_Dataset_OutputFcn, ...
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


% --- Executes just before Generate_Synthetic_Dataset is made visible.
function Generate_Synthetic_Dataset_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Generate_Synthetic_Dataset (see VARARGIN)

global h1 xs_3D ys_3D zs_3D tag_3D picks_3D is_2D is_3D picks_2D

% Choose default command line output for Generate_Synthetic_Dataset
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Add the path to the source codes
addpath('src')

% % Reset the plot_axes and uitable
h1=handles.plot_seis;
cla(h1,'reset');
% set(h1,'YTickLabel',[]);
% set(h1,'YTick',[]);

h_table= handles.uitable1;
cla(h_table,'reset');
 
% Set default parameters for the sliders
set(handles.is_3D,'Value',1);
set(handles.is_2D,'Value',0);

is_2D = 0;
is_3D = 1;

xs_3D=[]; ys_3D=[]; zs_3D=[]; tag_3D = [];

picks_3D=[num2cell(false(1,1)) {0} {0} {0} {0} {0} {0} {0} {0} {0}];
picks_2D=[num2cell(false(1,1)) {0} {0} {0} {0} {0} {0}];

h_table= handles.uitable1;
set(h_table, 'Data', picks_3D);
set(h_table, 'ColumnName', {' ', 'Strike', 'Dip', 'L', 'W', 'Err', 'nhypos', 'xb', 'yb', 'zb'});
set(h_table,'CellEditCallback',@onCellSelected_gen_synthetic);

% UIWAIT makes Generate_Synthetic_Dataset wait for user response (see UIRESUME)
% uiwait(handles.figure1);



% --- Outputs from this function are returned to the command line.
function varargout = Generate_Synthetic_Dataset_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in is_3D.
function is_3D_Callback(hObject, eventdata, handles)
% hObject    handle to is_3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global is_3D is_2D picks_3D 

h1=handles.plot_seis;
cla(h1,'reset');

% Hint: get(hObject,'Value') returns toggle state of is_3D
is_3D = get(hObject,'Value');
is_2D = 0;

h_table= handles.uitable1;
set(h_table, 'ColumnName', {' ', 'Strike', 'Dip', 'L', 'W', 'Err', 'nhypos', 'xb', 'yb', 'zb'});

% Determine synthetic datasets and plot them.
create_syn_and_plot_datasets

set(h_table, 'Data', picks_3D);
set(h_table,'CellEditCallback',@onCellSelected_gen_synthetic);



% --- Executes on button press in is_2D.
function is_2D_Callback(hObject, eventdata, handles)
% hObject    handle to is_2D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global is_2D  is_3D picks_2D

h1=handles.plot_seis;
cla(h1,'reset');

% Hint: get(hObject,'Value') returns toggle state of is_2D
is_2D = get(hObject,'Value');
is_3D = 0;

h_table= handles.uitable1;
set(h_table, 'ColumnName', {' ', 'Strike', 'L', 'Err', 'nhypos', 'xb', 'yb'});

% Determine synthetic datasets and plot them.
create_syn_and_plot_datasets

set(h_table, 'Data', picks_2D);
set(h_table,'CellEditCallback',@onCellSelected_gen_synthetic);

    
   
% --- Executes on button press in Add_push.
function Add_push_Callback(hObject, eventdata, handles)
% hObject    handle to Add_push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global picks_3D picks_2D Kfaults_3D Kfaults_2D is_3D is_2D 

data=get(handles.uitable1, 'data');
data(end+1,:)={0}; %if data is a cell or  
set(handles.uitable1, 'data', data);

if is_3D == 1
    Kfaults_3D = Kfaults_3D+1;

    picks_3D = [picks_3D; [num2cell(false(1,1)) {0} {0} {0} {0} {0} {0} {0} {0} {0}]];

    h_table= handles.uitable1;
    set(h_table, 'Data', picks_3D);
    set(h_table,'CellEditCallback',@onCellSelected_gen_synthetic);

elseif is_2D == 1
    Kfaults_2D = Kfaults_2D+1;

    picks_2D = [picks_2D; [num2cell(false(1,1)) {0} {0} {0} {0} {0} {0}]];

    h_table= handles.uitable1;
    set(h_table, 'Data', picks_2D);
    set(h_table,'CellEditCallback',@onCellSelected_gen_synthetic);

end    
    
    
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


% --- Executes on button press in Export_Push.
function Export_Push_Callback(hObject, eventdata, handles)
% hObject    handle to Export_Push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global picks_3D xs_3D ys_3D zs_3D tag_3D is_3D
global picks_2D xs_2D ys_2D  tag_2D is_2D 

%  uiputfile dialog
[filename,dirname]=uiputfile('*.txt');
if filename == 0; return;end

faultfile = strcat(dirname,filename);
fid = fopen(faultfile,'w');

%  Write cluster seismicity to file
if is_3D == 1
    Kfaults_3D = length(picks_3D(:,1));
    flt_no_3D = 1:Kfaults_3D;
    flt_no_3D = flt_no_3D(cell2mat(picks_3D(:,1)));

    for k=flt_no_3D
        xs_n_3D = xs_3D(k==tag_3D);
        ys_n_3D = ys_3D(k==tag_3D);
        zs_n_3D = zs_3D(k==tag_3D);

        for kk=1:length(xs_n_3D)
            fprintf(fid,'%12.5g %12.5g %12.5g \n',[xs_n_3D(kk) ys_n_3D(kk) zs_n_3D(kk)]);
        end
    end

elseif is_2D==1
    Kfaults_2D = length(picks_2D(:,1));
    flt_no_2D = 1:Kfaults_2D;
    flt_no_2D = flt_no_2D(cell2mat(picks_2D(:,1)));

    for k=flt_no_2D
        xs_n_2D = xs_2D(k==tag_2D);
        ys_n_2D = ys_2D(k==tag_2D);

        for kk=1:length(xs_n_2D)
            fprintf(fid,'%12.5g %12.5g \n',[xs_n_2D(kk) ys_n_2D(kk)]);
        end
    end
end

fclose(fid);
 


% --- Executes on button press in Default_inputs_push.
function Default_inputs_push_Callback(hObject, eventdata, handles)
% hObject    handle to Default_inputs_push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global h_table   
global picks_3D xs_3D ys_3D zs_3D xv_3D yv_3D zv_3D tag_3D is_3D Kfaults_3D
global picks_2D xs_2D ys_2D zs_2D xv_2D yv_2D zv_2D tag_2D is_2D Kfaults_2D

if is_3D == 1
    % Displaying the faults on the uitable
    Strike_3D = [45 45 135]';
    Dip_3D = [90 90 45]';
    L_3D = [100 100 100]';
    W_3D = [50 50 50]';
    Err_3D = [1 1 1]';
    nhypos_3D = [200 200 200]';
    xb_3D = [ 0 0 0]';
    yb_3D = [20 -20 0]';
    zb_3D = [0 0 0]';

    Kfaults_3D = length(Strike_3D);
    flt_details_man = [Strike_3D  Dip_3D L_3D W_3D Err_3D nhypos_3D xb_3D yb_3D zb_3D];
    fColumn_man = num2cell(true(Kfaults_3D,1));

    picks_3D = [fColumn_man num2cell(flt_details_man(:,1:9))];

    % Determine synthetic datasets and plot them.
    create_syn_and_plot_datasets

    h_table= handles.uitable1;
    set(h_table, 'Data', picks_3D);
    set(h_table,'CellEditCallback',@onCellSelected_gen_synthetic);

elseif is_2D == 1
      
    % Displaying the faults on the uitable
    Strike_2D = [45 45 135]';
    L_2D = [100 100 100]';
    Err_2D = [1 1 1]';
    nhypos_2D = [200 200 200]';
    xb_2D = [ 0 0 0]';
    yb_2D = [20 -20 0]';

    Kfaults_2D = length(Strike_2D);
    flt_details_man = [Strike_2D  L_2D Err_2D nhypos_2D xb_2D yb_2D];
    fColumn_man = num2cell(true(Kfaults_2D,1));

    picks_2D = [fColumn_man num2cell(flt_details_man(:,1:6))];

    % Determine synthetic datasets and plot them.
    create_syn_and_plot_datasets

    h_table= handles.uitable1;
    set(h_table, 'Data', picks_2D);
    set(h_table,'CellEditCallback',@onCellSelected_gen_synthetic);

end


% --- Executes on button press in Generate_Datasets_push.
function Generate_Datasets_push_Callback(hObject, eventdata, handles)
% hObject    handle to Generate_Datasets_push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global h_table h1 
global picks_3D xs_3D ys_3D zs_3D tag_3D Kfaults_3D
global picks_2D xs_2D ys_2D  tag_2D Kfaults_2D

viewpt = h1.View;
create_syn_and_plot_datasets
view(h1,viewpt(1),viewpt(2));
