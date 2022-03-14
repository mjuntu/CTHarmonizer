function varargout = kernelmatch(varargin)
% KERNELMATCH MATLAB code for kernelmatch.fig
%      KERNELMATCH, by itself, creates a new KERNELMATCH or raises the existing
%      singleton*.
%
%      H = KERNELMATCH returns the handle to a new KERNELMATCH or the handle to
%      the existing singleton*.
%
%      KERNELMATCH('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in KERNELMATCH.M with the given input arguments.
%
%      KERNELMATCH('Property','Value',...) creates a new KERNELMATCH or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before kernelmatch_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to kernelmatch_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help kernelmatch

% Last Modified by GUIDE v2.5 14-May-2021 14:19:55

% Nina Hänninen & Mikael Juntunen, 2021

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @kernelmatch_OpeningFcn, ...
                   'gui_OutputFcn',  @kernelmatch_OutputFcn, ...
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


% --- Executes just before kernelmatch is made visible.
function kernelmatch_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to kernelmatch (see VARARGIN)

% Choose default command line output for kernelmatch
handles.output = hObject;

%load('DATA/results.mat');
% results_load = matfile('DATA/results.mat');
% varname = who(results_load);
% results = varname{1,1};
%fieldnames(results_load.(results)) 

% for ii = 1:length(fieldnames(results))
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load in the result.mat file
handles.function_path = GetExecutableFolder();
% Load in the the result file containing the NPSs, MTFs etc
if isdeployed
    % User is running an executable in standalone mode. 
    % Identify the DATA directory, load in the results file and test that
    % the directory is otherwise ok (contains the slice images)
    %% 1. Identify the DATA directory
    waitfor(msgbox('Cannot locate the DATA directory. Please identify the path.')); % For macos, the below title text in uigetdir does not work, so we have this additional messagebox  
    handles.data_path = uigetdir('','Cannot locate the DATA directory. Please identify the path.');
    results = handles.data_path;
    %% 2. Try to load in the results file and example slice
    try
        load(sprintf('%s/results.mat',handles.data_path));
        handles.results = results;
        % Input the list of CT models:
        list_models = fieldnames(handles.results);
        % Remove the data path from the list of variables
        idx = find(cellfun(@isempty,strfind(list_models,'data_path')));
        list_models = list_models(idx);
        handles.list_models = list_models;
        
        %
        % Check if the example slices do actually exist in the data path
        %
        files = fullfile(sprintf('%s/Slices/%s/%s/%s/Slice_1.dcm',handles.data_path,handles.list_models{1},results.(list_models{1}).Studies{1},results.(list_models{1}).Series_name{1}));
        if ~exist(files) % No data found --> no valid DATA directory
            waitfor(warndlg('No valid DATA directory given, please check.','Warning'));    
        end
    catch
        waitfor(warndlg('No valid DATA directory given, please check.','Warning'));    
    end
else
    % User is running an m-file from the MATLAB integrated development environment (regular MATLAB).
    load(sprintf('%s/DATA/results.mat',handles.function_path));
    handles.results = results;

    % Input the list of CT models:
    list_models = fieldnames(handles.results);
    % Remove the data path from the list of variables
    idx = find(cellfun(@isempty,strfind(list_models,'data_path')));
    list_models = list_models(idx);
    handles.list_models = list_models;
    %
    % Check if the example slices do actually exist in the data path
    %
    try
        handles.data_path = results.data_path;
        % Check that the path actually exists (makes sure that no old paths
        % exist)
        files = fullfile(sprintf('%s/Slices/%s/%s/%s/Slice_1.dcm',handles.data_path,handles.list_models{1},results.(list_models{1}).Studies{1},results.(list_models{1}).Series_name{1}));

        if ~exist(files) % No data found --> no valid DATA directory
            waitfor(msgbox('Cannot locate the DATA directory. Please identify the path.')); % For macos, the below title text in uigetdir does not work, so we have this additional messagebox  
            handles.data_path = uigetdir('','Cannot locate the DATA directory. Please identify the path.');
            results.data_path = handles.data_path;
        end
    catch
        waitfor(msgbox('Cannot locate the DATA directory. Please identify the path.')); % For macos, the below title text in uigetdir does not work, so we have this additional messagebox  
        handles.data_path = uigetdir('','Cannot locate the DATA directory. Please identify the path.');
        results.data_path = handles.data_path;
    end
    % Check that the given folder was valid. If still no valid directory, give
    % a warning
    try
        files = fullfile(sprintf('%s/Slices/%s/%s/%s/Slice_1.dcm',handles.data_path,handles.list_models{1},results.(list_models{1}).Studies{1},results.(list_models{1}).Series_name{1}));

        if ~exist(files) % No data found --> no valid DATA directory
            results = rmfield(results,'data_path');
            waitfor(warndlg('No valid DATA directory given, please check.','Warning'));    
        end
        % Update the results.mat file
        if isdeployed 
            waitfor(msgbox('Going to data save.'));    
            % User is running an executable in standalone mode. 
            save('results.mat', 'results', '-append');
        else
            % User is running an m-file from the MATLAB integrated development environment (regular MATLAB).
            save(sprintf('%s/results.mat',handles.data_path),'results','-append');
        end 
    catch
        results = rmfield(results,'data_path');
        waitfor(warndlg('No valid DATA directory given, please check.','Warning'));    
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update the results.mat file
if ~isdeployed % For deployed version, data cannot be saved in this manneer
    % User is running an m-file from the MATLAB integrated development environment (regular MATLAB).
    save(sprintf('%s/results.mat',handles.data_path),'results','-append');
end 

list_models = fieldnames(handles.results);
% Remove the data path from the list of variables
idx = find(cellfun(@isempty,strfind(list_models,'data_path')));
list_models = list_models(idx);
handles.list_models = list_models;

% Empty first option added
list_models = [{'Choose...'}; list_models];
set(handles.targetmodel, 'String', list_models);
%set(handles.targetmodel, 'String', fieldnames(results));
set(handles.matchmodel, 'String', list_models);

% Set default window parameters
    WL_value_default = 40;
    WW_value_default = 500;
    
    handles.WL_value_default = WL_value_default;
    handles.WW_value_default = WW_value_default;
    handles.WL_value = WL_value_default;
    handles.WW_value = WW_value_default;
    
    % default image
    handles.fig_nro = 1;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes kernelmatch wait for user response (see UIRESUME)
% uiwait(handles.figure1);






% --- Outputs from this function are returned to the command line.
function varargout = kernelmatch_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in matchmodel.
function matchmodel_Callback(hObject, eventdata, handles)
% hObject    handle to matchmodel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns matchmodel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from matchmodel
CTmodel = get(hObject,'Value');
CTmodel_contents = cellstr(get(hObject,'String'));
CTmodel_name = CTmodel_contents{get(hObject,'Value')};
set(handles.text5,'String',CTmodel_name);
handles.target_CTmodel = CTmodel_name;

if ~any(strcmp(handles.list_models,CTmodel_name))
    axes(handles.axes2);
    cla
else

% ACTUAL MATCHING

[matching_kernel, handles] = match_kernel(handles,CTmodel_name);
set(handles.text2,'String',matching_kernel);

end


% Update handles structure
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function matchmodel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to matchmodel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in manualmatch.
function manualmatch_Callback(hObject, eventdata, handles)
% hObject    handle to manualmatch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% close current window
%closereq;


% --- Executes on button press in advanced.
function advanced_Callback(hObject, eventdata, handles)
% hObject    handle to advanced (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if isfield(handles,'matching_kernel')
    if ~isempty(handles.matching_kernel)
fontsize = 12;
fontname = 'Times new roman';

%close all force;
%mkdir(figLoc);
figure


matched_idx = handles.idxs_to_menu_sort;
m = handles.m_sort;
to_be_matched = handles.matching_model_struct;
reference = handles.model_struct;

% target kernel index
idx_ref = handles.target_image_indx(1);

% matching kernel index
ii_ref = handles.match_rank;


    % Find the reference images and matched kernel images
    %     Iref = reference.ims{idx_ref}; % OLD METHOD
    % Reference image
    Iref = dicomread(sprintf('%s/Slices/%s/%s/%s/Slice_1.dcm',handles.data_path,handles.chosen_target_model,reference.Studies{idx_ref},reference.Series_name{idx_ref}));
    Iref = cat(3,Iref,dicomread(sprintf('%s/Slices/%s/%s/%s/Slice_2.dcm',handles.data_path,handles.chosen_target_model,reference.Studies{idx_ref},reference.Series_name{idx_ref})));
    info = dicominfo(sprintf('%s/Slices/%s/%s/%s/Slice_1.dcm',handles.data_path,handles.chosen_target_model,reference.Studies{idx_ref},reference.Series_name{idx_ref}));
    Iref = double(Iref)*info.RescaleSlope + info.RescaleIntercept;
    % Matched image 
    Imatched = dicomread(sprintf('%s/Slices/%s/%s/%s/Slice_1.dcm',handles.data_path,handles.target_CTmodel,to_be_matched.Studies{matched_idx(ii_ref)},to_be_matched.Series_name{matched_idx(ii_ref)}));
    Imatched = cat(3,Imatched,dicomread(sprintf('%s/Slices/%s/%s/%s/Slice_2.dcm',handles.data_path,handles.target_CTmodel,to_be_matched.Studies{matched_idx(ii_ref)},to_be_matched.Series_name{matched_idx(ii_ref)})));
    info = dicominfo(sprintf('%s/Slices/%s/%s/%s/Slice_1.dcm',handles.data_path,handles.target_CTmodel,to_be_matched.Studies{matched_idx(ii_ref)},to_be_matched.Series_name{matched_idx(ii_ref)}));
    Imatched = double(Imatched)*info.RescaleSlope + info.RescaleIntercept;
    LIMS = [-200,200];
    LIMS2 = [-200,200];
%     figure;
    ha = tight_subplot(3,2,[.01 .01],[.1 .1],[.2 .2]);
    % Display first slice with homogeneous region
    axes(ha(1));imshow(Iref(:,:,1),LIMS);title(sprintf('%s\n%s',reference.Model,reference.Recon_names{idx_ref}),'interpreter','none','fontsize',fontsize,'fontname',fontname);
    axes(ha(2));imshow(Imatched(:,:,1),LIMS);
    title(sprintf('%s\n%s\n m = %4.3f',to_be_matched.Model,to_be_matched.Recon_names{matched_idx(ii_ref)},m(ii_ref)),'interpreter','none','fontsize',fontsize,'fontname',fontname);
    % Display second slice with contrast targets
    axes(ha(3));imshow(Iref(:,:,2),LIMS2);%title(sprintf('%s\n%s %s (index = %d)',devices(4).name,reference.Recon_info{ii_ref},reference.Kernel_names{ii_ref},ii_ref),'interpreter','none');
    axes(ha(4));imshow(Imatched(:,:,2),LIMS2);%title(sprintf('%s\n%s (index = %d)',devices(5).name,to_be_matched.Recon_info{matched_idx(ii_ref)},matched_idx(ii_ref)));
   
    axes(ha(5));
    plot(reference.MTF(idx_ref).freq,reference.MTF(idx_ref).MTF,'k');hold on;
    plot(to_be_matched.MTF(matched_idx(ii_ref)).freq,to_be_matched.MTF(matched_idx(ii_ref)).MTF,'k--');hold off;
    set(gca,'ytick',[]);
    xlabel('Spatial frequency (1/mm)','fontsize',fontsize,'fontname',fontname);
    ylabel('MTF','fontsize',fontsize,'fontname',fontname);
    legend(reference.Model,to_be_matched.Model,'interpreter','none','fontsize',fontsize-5,'fontname',fontname)
    axes(ha(6));
    plot(reference.NPS(idx_ref).freq,reference.NPS(idx_ref).NPS/max(reference.NPS(idx_ref).NPS),'k');hold on;
    plot(to_be_matched.NPS(matched_idx(ii_ref)).freq,to_be_matched.NPS(matched_idx(ii_ref)).NPS/max(to_be_matched.NPS(matched_idx(ii_ref)).NPS),'k--');hold off;
    set(gca,'ytick',[]);
    xlabel('Spatial frequency (1/mm)','fontsize',fontsize,'fontname',fontname);
    ylabel('nNPS','fontsize',fontsize,'fontname',fontname);
    legend(reference.Model,to_be_matched.Model,'interpreter','none','fontsize',fontsize-5,'fontname',fontname)
    %drawnow;
    %print(sprintf('%s/%d.png',figLoc,ii_ref),'-dpng');
    %close gcf;
    end
end



% --- Executes on slider movement.
function WL_slider_Callback(hObject, eventdata, handles)
% Adjust Window Length with slider

% hObject    handle to WL_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


%       % CT levels:
%    the upper grey level (x) is calculated via WL + (WW � 2)
%    the lower grey level (y) is calculated via WL - (WW � 2)

%I = handles.Datainput;
% handles.WL_value = get(handles.WL_slider,'Value');
% handles.WW_value = get(handles.WW_slider,'Value');

set(handles.error_text,'String','')
WL_value = get(hObject,'Value');
handles.WL_value = WL_value;

WW_value = get(handles.WW_slider,'Value');
axes(handles.axes1);
caxis([(WL_value - (WW_value/2)) (WL_value + (WW_value/2))]);
axes(handles.axes2);
caxis([(WL_value - (WW_value/2)) (WL_value + (WW_value/2))]);
set(handles.WL_edit,'String',WL_value);

guidata(hObject, handles);




% --- Executes during object creation, after setting all properties.
function WL_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WL_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function WW_slider_Callback(hObject, eventdata, handles)
% Adjust Window Width with slider

% hObject    handle to WW_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

%I = handles.Datainput;
set(handles.error_text,'String','')
WW_value = get(hObject,'Value');
handles.WW_value = WW_value;

WL_value = get(handles.WL_slider,'Value');
axes(handles.axes1);
caxis([(WL_value - (WW_value/2)) (WL_value + (WW_value/2))]);
axes(handles.axes2);
caxis([(WL_value - (WW_value/2)) (WL_value + (WW_value/2))]);
set(handles.WW_edit,'String',WW_value);

guidata(hObject, handles);




% --- Executes during object creation, after setting all properties.
function WW_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WW_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function WL_edit_Callback(hObject, eventdata, handles)
% Adjust Window Length by inputting value in the field

% hObject    handle to WL_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of WL_edit as text
%        str2double(get(hObject,'String')) returns contents of WL_edit as a double

WL_value = get(hObject,'String');
WL_max = get(handles.WL_slider,'Max');
WL_min = get(handles.WL_slider,'Min');
try
WL_value = str2double(WL_value);

if WL_value <= WL_max && WL_value >= WL_min
    set(handles.error_text,'String','')
    handles.WL_value = WL_value;
    
    WW_value = get(handles.WW_slider,'Value');
    axes(handles.axes1);
    caxis([(WL_value - (WW_value/2)) (WL_value + (WW_value/2))]);
    axes(handles.axes2);
    caxis([(WL_value - (WW_value/2)) (WL_value + (WW_value/2))]);
    set(handles.WL_slider,'Value',WL_value);
else
    set(handles.error_text,'String',['Please input numeric W-L value between ', num2str(WL_min), ' and ', num2str(WL_max)]);
end
catch
    set(handles.error_text,'String',['Please input numeric W-L value between ', num2str(WL_min), ' and ', num2str(WL_max)]);
end

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function WL_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WL_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function WW_edit_Callback(hObject, eventdata, handles)
% Adjust Window Width by inputting value in the field

% hObject    handle to WW_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of WW_edit as text
%        str2double(get(hObject,'String')) returns contents of WW_edit as a double

WW_value = get(hObject,'String');
WW_max = get(handles.WW_slider,'Max');
WW_min = get(handles.WW_slider,'Min');
try
WW_value = str2double(WW_value);

if WW_value <= WW_max && WW_value >= WW_min
    set(handles.error_text,'String','')
    handles.WW_value = WW_value;
    
    WL_value = get(handles.WL_slider,'Value');
    axes(handles.axes1);
    caxis([(WL_value - (WW_value/2)) (WL_value + (WW_value/2))]);
    axes(handles.axes2);
    caxis([(WL_value - (WW_value/2)) (WL_value + (WW_value/2))]);
    set(handles.WW_slider,'Value',WW_value);
else
    set(handles.error_text,'String',['Please input numeric W-W value between ', num2str(WW_min), ' and ', num2str(WW_max)]);
end
catch
    set(handles.error_text,'String',['Please input numeric W-W value between ', num2str(WW_min), ' and ', num2str(WW_max)]);
end

guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function WW_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WW_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Reset_Button.
function Reset_Button_Callback(hObject, eventdata, handles)
% Reset to default window width and level

% hObject    handle to Reset_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% DatainputHU = handles.DatainputHU;
% axes(handles.axes1);
% imshow(DatainputHU);

try
%Datainput_info = handles.Datainput_info;
% set slider values
WL_value = handles.WL_value_default;
WW_value = handles.WW_value_default;
handles.WL_value = WL_value;
handles.WW_value = WW_value;
set(handles.WL_slider,'Value',WL_value);
set(handles.WL_edit,'String',WL_value);
set(handles.WW_slider,'Value',WW_value);
set(handles.WW_edit,'String',WW_value);
axes(handles.axes1);
caxis([(WL_value - (WW_value/2)) (WL_value + (WW_value/2))]);
axes(handles.axes2);
caxis([(WL_value - (WW_value/2)) (WL_value + (WW_value/2))]);

guidata(hObject, handles);
catch
    set(handles.error_text,'String','Please input target data')
end


% --- Executes on selection change in targetmodel.
function targetmodel_Callback(hObject, eventdata, handles)
% hObject    handle to targetmodel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns targetmodel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from targetmodel


contents = cellstr(get(hObject,'String'));
chosen_target_model = contents{get(hObject,'Value')};
handles.chosen_target_model = chosen_target_model;

results = handles.results;

fig_nro = handles.fig_nro;

if isfield(results,chosen_target_model)
    % Check which CT model is inputted
    model_struct = results.(chosen_target_model);
    handles.model_struct = model_struct;
    
    % Find slice thickness values for chosen model and input them into corresponding popupmenu
    slice_thick_values = unique(model_struct.Slice_thickness);
    
    % Empty possible previous inputs
    set(handles.slicethicknessmenu, 'Value', 1);
    
    % Set the new slice thickness values
    %slice_thick_values = [{'Choose...'}; slice_thick_values];
    set(handles.slicethicknessmenu, 'String', slice_thick_values);
    slice_thickness = slice_thick_values(1);
    handles.slice_thickness = slice_thickness;
    
    
    % Find slice CTDIvol values for chosen model and input them into corresponding popupmenu
    CTDIvol_values = unique(round(model_struct.CTDIvol));
    
    % Empty possible previous inputs
    set(handles.CTDIvolmenu, 'Value', 1);
    
    % Set the new values
    set(handles.CTDIvolmenu, 'String', CTDIvol_values);
    CTDIvol = CTDIvol_values(1);
    handles.CTDIvol = CTDIvol;
    
    % Find Body part values for chosen model and input them into corresponding popupmenu
    bodypart_values = unique(model_struct.Body_part);
    
    % Empty possible previous inputs
    set(handles.bodypartmenu, 'Value', 1);
    
    % Set the new values
    set(handles.bodypartmenu, 'String', bodypart_values);
    body_part = bodypart_values{1};
    handles.body_part = body_part;
    
    
    % Find kernel values for chosen parameters, and input them into corresponding popupmenu
    [handles, target_kernels] = show_target_kernel_list(handles);
    
    kerneli = target_kernels(1);
    handles.kerneli = kerneli;
    set(handles.chosen_kernel, 'String', kerneli);
    
    target_image_indx = find_image(model_struct, slice_thickness, CTDIvol, body_part, kerneli);
    handles.target_image_indx = target_image_indx;
    
    handles = show_imag_axes1(model_struct, target_image_indx, handles, fig_nro);
    
    if isfield(handles, 'target_CTmodel') % If matching model is defined
        if any(strcmp(handles.list_models,handles.target_CTmodel))
            % ACTUAL MATCHING
            
            CTmodel_name = handles.target_CTmodel;
            [matching_kernel, handles] = match_kernel(handles,CTmodel_name);
            set(handles.text2,'String',matching_kernel);
        end
    end
    
else
    % Set default texts
    set(handles.slicethicknessmenu, 'String','Slice Thickness');
    set(handles.CTDIvolmenu, 'String', 'CTDIvol');
    set(handles.bodypartmenu, 'String', 'Body Part');
    set(handles.kernelmenu, 'String', 'Kernel');
    axes(handles.axes1);
    cla
    axes(handles.axes2);
    cla
end

guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function targetmodel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to targetmodel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in slicethicknessmenu.
function slicethicknessmenu_Callback(hObject, eventdata, handles)
% hObject    handle to slicethicknessmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns slicethicknessmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from slicethicknessmenu

contents = cellstr(get(hObject,'String'));
slice_thickness = contents{get(hObject,'Value')};
slice_thickness = str2double(slice_thickness);

if ~isnan(slice_thickness)
    handles.slice_thickness = slice_thickness;
    
    [handles] = show_target_kernel_list(handles);
    
    handles = show_kernel_and_match(handles);
    
      
    
else
%     set(handles.CTDIvolmenu, 'Value', 1);
%     set(handles.CTDIvolmenu, 'String', 'CTDIvolmenu');
end

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function slicethicknessmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slicethicknessmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in CTDIvolmenu.
function CTDIvolmenu_Callback(hObject, eventdata, handles)
% hObject    handle to CTDIvolmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns CTDIvolmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from CTDIvolmenu

contents = cellstr(get(hObject,'String'));
CTDIvol = contents{get(hObject,'Value')};
CTDIvol = str2double(CTDIvol);

if ~isnan(CTDIvol)
    handles.CTDIvol = CTDIvol;
    
    [handles] = show_target_kernel_list(handles);
    handles = show_kernel_and_match(handles);
    
    
end


guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function CTDIvolmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CTDIvolmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in bodypartmenu.
function bodypartmenu_Callback(hObject, eventdata, handles)
% hObject    handle to bodypartmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns bodypartmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from bodypartmenu

contents = cellstr(get(hObject,'String'));
body_part = contents{get(hObject,'Value')};

if isfield(handles,'model_struct')
    handles.body_part = body_part;
    
    %model_struct = handles.model_struct;
    
    [handles] = show_target_kernel_list(handles);
    
    
    handles = show_kernel_and_match(handles);
    
end

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function bodypartmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bodypartmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in kernelmenu.
function kernelmenu_Callback(hObject, eventdata, handles)
% hObject    handle to kernelmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns kernelmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from kernelmenu

contents = cellstr(get(hObject,'String'));
kerneli = contents{get(hObject,'Value')};


    handles.kerneli = kerneli;
    set(handles.chosen_kernel, 'String', kerneli);
    
    handles = show_kernel_and_match(handles);
    
    
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function kernelmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kernelmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% Find CT image from the 'results' structure based on given variables
function imag = find_image(model_struct, slice, ctdi, bodypart, kernel)
    
slice_indxs = model_struct.Slice_thickness == slice;
ctdi_indxs = round(model_struct.CTDIvol) == ctdi;
bodypart_indxs = strcmp(model_struct.Body_part,bodypart);
kernel_indxs = strcmp(model_struct.Recon_names,kernel);

the_indxs = slice_indxs & ctdi_indxs & bodypart_indxs & kernel_indxs;
imag = find(the_indxs);


% Plot image on axes1 (Target) based on image index
% fig_nro defines which of the two images is shown, either 1 or 2
function handles = show_imag_axes1(model_struct, target_image_indx, handles,fig_nro)


if isempty(target_image_indx)   % indx is empty
    set(handles.error_text,'String','Such combination does not exist, please choose some other paramaters')
    axes(handles.axes1);
    cla
    
else
    if length(target_image_indx) > 1    % multiple hits
        
        %set(handles.error_text,'String','Multiple images found, showing only first one')
    else
        set(handles.error_text,'String','')
    end
    
    % Original method with results structure containing the image data
%     target_image = model_struct.ims{1,target_image_indx(1)};
    % New method 
    target_image = dicomread(sprintf('%s/Slices/%s/%s/%s/Slice_1.dcm',handles.data_path,handles.chosen_target_model,model_struct.Studies{target_image_indx(1)},model_struct.Series_name{target_image_indx(1)}));
    target_image = cat(3,target_image,dicomread(sprintf('%s/Slices/%s/%s/%s/Slice_2.dcm',handles.data_path,handles.chosen_target_model,model_struct.Studies{target_image_indx(1)},model_struct.Series_name{target_image_indx(1)})));
    info = dicominfo(sprintf('%s/Slices/%s/%s/%s/Slice_1.dcm',handles.data_path,handles.chosen_target_model,model_struct.Studies{target_image_indx(1)},model_struct.Series_name{target_image_indx(1)}));
    target_image = double(target_image)*info.RescaleSlope + info.RescaleIntercept;
    
    axes(handles.axes1);
    imshow(target_image(:,:,fig_nro));
    

       
    WL_value = handles.WL_value;
    WW_value = handles.WW_value;
    set(handles.WL_slider,'Value',WL_value);
    set(handles.WL_edit,'String',WL_value);
    set(handles.WW_slider,'Value',WW_value);
    set(handles.WW_edit,'String',WW_value);
    handles.WL_value = WL_value;
    handles.WW_value = WW_value;
    caxis([(WL_value - (WW_value/2)) (WL_value + (WW_value/2))]);
 
    
end


% Plot image on axes1 (Target) based on image index
% fig_nro defines which of the two images is shown, either 1 or 2
function handles = show_imag_axes2(model_struct, target_image_indx, handles, fig_nro)


if isempty(target_image_indx)   % indx is empty
    set(handles.error_text,'String','No image index defined')
    axes(handles.axes1);
    cla
    
else
    if length(target_image_indx) > 1    % multiple hits
        
        set(handles.error_text,'String','Multiple images found, showing only first one')
    else
        set(handles.error_text,'String','')
    end
    % OLD method to obtain the target image
    %target_image = model_struct.ims{1,target_image_indx(1)};
    % New method
    target_image = dicomread(sprintf('%s/Slices/%s/%s/%s/Slice_1.dcm',handles.data_path,handles.target_CTmodel,model_struct.Studies{target_image_indx(1)},model_struct.Series_name{target_image_indx(1)}));
    target_image = cat(3,target_image,dicomread(sprintf('%s/Slices/%s/%s/%s/Slice_2.dcm',handles.data_path,handles.target_CTmodel,model_struct.Studies{target_image_indx(1)},model_struct.Series_name{target_image_indx(1)})));
    info = dicominfo(sprintf('%s/Slices/%s/%s/%s/Slice_1.dcm',handles.data_path,handles.target_CTmodel,model_struct.Studies{target_image_indx(1)},model_struct.Series_name{target_image_indx(1)}));
    target_image = double(target_image)*info.RescaleSlope + info.RescaleIntercept;
    
    
    axes(handles.axes2);
    imshow(target_image(:,:,fig_nro));
    

       
    WL_value = handles.WL_value;
    WW_value = handles.WW_value;
    caxis([(WL_value - (WW_value/2)) (WL_value + (WW_value/2))]);
 
    
end


% --- Executes on button press in fignro_button.
function fignro_button_Callback(hObject, eventdata, handles)
% Change the phantom image

% hObject    handle to fignro_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fignro_button

fign = get(hObject,'Value');

if fign
    fig_nro = 2;
    handles.fig_nro = fig_nro;
    set(handles.fignro_button,'String','Image 1')
else
    fig_nro = 1;
    handles.fig_nro = fig_nro;
    set(handles.fignro_button,'String','Image 2')
end

   
    model_struct = handles.model_struct;
    slice_thickness = handles.slice_thickness;
    CTDIvol = handles.CTDIvol;
    body_part = handles.body_part;
    kerneli = handles.kerneli;
    target_image_indx = handles.target_image_indx;
    
    handles = show_imag_axes1(model_struct, target_image_indx, handles, fig_nro);
    
    if isfield(handles,'matching_kernel')
    if ~isempty(handles.matching_kernel)
        matching_model_struct = handles.matching_model_struct;
        matching_kernel = handles.matching_kernel;
        handles = show_imag_axes2(matching_model_struct, matching_kernel, handles, fig_nro);
    end
    end

guidata(hObject, handles);


% Find matching kernel from the 'results' structure based on given variables
function [matching_kernel, handles] = match_kernel(handles,matched_model)

% Values from Winslow et al
% ad = 10;cd1 = 0.5;
% at = 150;ct = -0.05;
% ar = 150;cr = -0.05;
% am = 100;cm = -0.2;
% Values used in Juntunen et al.
ad = 10;cd1 = 0.5;
at = 30;ct = -0.2;
ar = 30;cr = -0.2;
am = 100;cm = -0.2;
    
results = handles.results;
% model_struct, slice, ctdi, bodypart, kernel

if isfield(results,matched_model)
    % Check which CT model is inputted
    matching_model_struct = results.(matched_model);
    handles.matching_model_struct = matching_model_struct;
else
    % Some error
end

if isfield(handles,'target_image_indx')
if ~isempty(handles.target_image_indx)

     slice_thickness = handles.slice_thickness;
     CTDIvol = handles.CTDIvol;
     body_part = handles.body_part;
%     kerneli = handles.kerneli;

% Initialize the matching process between Siemens Flash and DRIVE
width_goal = 512;
reference = handles.model_struct;

% Limit the analysis to measurements with the correct bow-tie filter and
% slice thickness
% idxs_ref = find(reference.Width == width_goal & ~cellfun(@isempty,strfind(reference.Body_part,'MEDIUM FILTER')) & reference.Slice_thickness > 2);% & ~cellfun(@isempty,strfind(reference.Recon_info,'FBP')));
% idxs_ref = find(reference.Width == width_goal & ~cellfun(@isempty,strfind(reference.Body_part,'ABDOMEN')) & reference.Slice_thickness < 2);% & ~cellfun(@isempty,strfind(reference.Recon_info,'FBP')));
idx_ref = handles.target_image_indx(1);

% Data to be matched
to_be_matched = matching_model_struct;

% Please NOTE: The name of the matched body part changes depending on the
% vendor --> the commented part below works for GE devices
%idxs_to_be_matched = find(to_be_matched.Width == width_goal & ~cellfun(@isempty,strfind(to_be_matched.Body_part,'ABDOMEN')) & to_be_matched.Slice_thickness < 2);% & ~cellfun(@isempty,strfind(to_be_matched.Recon_info,'B')));
%idxs_to_be_matched = find(to_be_matched.Width == width_goal & ~cellfun(@isempty,strfind(to_be_matched.Body_part,'MEDIUM FILTER')) & to_be_matched.Slice_thickness < 2);% & ~cellfun(@isempty,strfind(reference.Recon_info,'FBP')));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADD HERE any other matching body-part values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(strfind(body_part,'ABDOMEN')) |  ~isempty(strfind(body_part,'BODY FILTER')) | ~isempty(strfind(body_part,'MEDIUM FILTER'))
    idxs_to_be_matched = find(to_be_matched.Width == width_goal & ( ~cellfun(@isempty,( strfind(to_be_matched.Body_part,'ABDOMEN'))) | ~cellfun(@isempty,( strfind(to_be_matched.Body_part,'MEDIUM FILTER')))  | ~cellfun(@isempty,( strfind(to_be_matched.Body_part,'BODY FILTER'))) ) & to_be_matched.Slice_thickness < slice_thickness+1 & to_be_matched.Slice_thickness > slice_thickness-1 ...
        & to_be_matched.CTDIvol < CTDIvol+5 & to_be_matched.CTDIvol > CTDIvol-5);% & ~cellfun(@isempty,strfind(to_be_matched.Recon_info,'B')));
elseif ~isempty(strfind(body_part,'BRAIN')) | ~isempty(strfind(body_part,'HEAD'))
    idxs_to_be_matched = find(to_be_matched.Width == width_goal & ( ~cellfun(@isempty,( strfind(to_be_matched.Body_part,'BRAIN'))) | ~cellfun(@isempty,( strfind(to_be_matched.Body_part,'HEAD'))) ) & to_be_matched.Slice_thickness < slice_thickness+1 & to_be_matched.Slice_thickness > slice_thickness-1 ...
        & to_be_matched.CTDIvol < CTDIvol+5 & to_be_matched.CTDIvol > CTDIvol-5);
end
% if isempty(idxs_to_be_matched) && strcmp(body_part,'MEDIUM FILTER')
%     idxs_to_be_matched = find(to_be_matched.Width == width_goal & ( ~cellfun(@isempty,( strfind(to_be_matched.Body_part,'ABDOMEN'))) | ~cellfun(@isempty,( strfind(to_be_matched.Body_part,'MEDIUM FILTER')))  | ~cellfun(@isempty,( strfind(to_be_matched.Body_part,'BODY FILTER'))) ) & to_be_matched.Slice_thickness < slice_thickness+1 & to_be_matched.Slice_thickness > slice_thickness-1 ...
%         & to_be_matched.CTDIvol < CTDIvol+5 & to_be_matched.CTDIvol > CTDIvol-5);% & ~cellfun(@isempty,strfind(to_be_matched.Recon_info,'B')));
% elseif isempty(idxs_to_be_matched) && strcmp(body_part,'ABDOMEN')
%     idxs_to_be_matched = find(to_be_matched.Width == width_goal & ~cellfun(@isempty,strfind(to_be_matched.Body_part,'MEDIUM FILTER')) & to_be_matched.Slice_thickness < slice_thickness+1 & to_be_matched.Slice_thickness > slice_thickness-1 ...
%         & to_be_matched.CTDIvol < CTDIvol+5 & to_be_matched.CTDIvol > CTDIvol-5);% & ~cellfun(@isempty,strfind(to_be_matched.Recon_info,'B')));
% elseif isempty(idxs_to_be_matched) && strcmp(body_part,'BRAIN')
%     idxs_to_be_matched = find(to_be_matched.Width == width_goal & ~cellfun(@isempty,strfind(to_be_matched.Body_part,'HEAD')) & to_be_matched.Slice_thickness < slice_thickness+1 & to_be_matched.Slice_thickness > slice_thickness-1 ...
%         & to_be_matched.CTDIvol < CTDIvol+5 & to_be_matched.CTDIvol > CTDIvol-5);% & ~cellfun(@isempty,strfind(to_be_matched.Recon_info,'B')));
% elseif isempty(idxs_to_be_matched) && strcmp(body_part,'HEAD')
%     idxs_to_be_matched = find(to_be_matched.Width == width_goal & ~cellfun(@isempty,strfind(to_be_matched.Body_part,'BRAIN')) & to_be_matched.Slice_thickness < slice_thickness+1 & to_be_matched.Slice_thickness > slice_thickness-1 ...
%         & to_be_matched.CTDIvol < CTDIvol+5 & to_be_matched.CTDIvol > CTDIvol-5);% & ~cellfun(@isempty,strfind(to_be_matched.Recon_info,'B')));
% end


% Perform the actual matching process
m = zeros(numel(reference),numel(to_be_matched)); % Matching parameter. The higher the value, the better the match
Ds = m;
Rs = m;
Ts = m;
Ms = m;

idxs_ref = idx_ref;
ii_kernel_ref = 1;
    for ii_kernel_matched = 1:numel(idxs_to_be_matched)
        % Skip the matching if the width 
        % Difference in dose
        d = (reference.CTDIvol(idxs_ref(ii_kernel_ref)) - to_be_matched.CTDIvol(idxs_to_be_matched(ii_kernel_matched)))/reference.CTDIvol(idxs_ref(ii_kernel_ref));
        
        %% OLD method based on the weighted average frequency and
        % MTF10/MTF50
        % Difference in noise texture. 
%         t = (reference.FR_NPS(idxs_ref(ii_kernel_ref)) - to_be_matched.FR_NPS(idxs_to_be_matched(ii_kernel_matched)));
%         % RMSE
%         %t = norm(reference.NPS(idxs_ref(ii_kernel_ref)).NPS/max(reference.NPS(idxs_ref(ii_kernel_ref)).NPS) - to_be_matched.NPS(idxs_to_be_matched(ii_kernel_matched)).NPS/max(to_be_matched.NPS(idxs_to_be_matched(ii_kernel_matched)).NPS),2)/norm(reference.NPS(idxs_ref(ii_kernel_ref)).NPS/max(reference.NPS(idxs_ref(ii_kernel_ref)).NPS),2);
% 
%         % Difference in spatial resolution. 
%         f_to_be_matched = (to_be_matched.f_50(idxs_to_be_matched(ii_kernel_matched)) + to_be_matched.f_10(idxs_to_be_matched(ii_kernel_matched)))/2;
%         f_ref = (reference.f_50(idxs_ref(ii_kernel_ref)) + reference.f_10(idxs_ref(ii_kernel_ref)))/2;
% 
%         r = (f_to_be_matched-f_ref);
        %% NEW method based on RMSE
        % NPS
        NPSref = interp1(reference.NPS(idxs_ref(ii_kernel_ref)).freq,reference.NPS(idxs_ref(ii_kernel_ref)).NPS/max(reference.NPS(idxs_ref(ii_kernel_ref)).NPS),(0:0.05:max([max(reference.NPS(idxs_ref(ii_kernel_ref)).freq),max(to_be_matched.NPS(idxs_to_be_matched(ii_kernel_matched)).freq)])));
        NPS_to_be_matched = interp1(to_be_matched.NPS(idxs_to_be_matched(ii_kernel_matched)).freq,to_be_matched.NPS(idxs_to_be_matched(ii_kernel_matched)).NPS/max(to_be_matched.NPS(idxs_to_be_matched(ii_kernel_matched)).NPS),(0:0.05:max([max(reference.NPS(idxs_ref(ii_kernel_ref)).freq),max(to_be_matched.NPS(idxs_to_be_matched(ii_kernel_matched)).freq)])));
        NPSref(isnan(NPSref)) = 0;NPS_to_be_matched(isnan(NPS_to_be_matched)) = 0;
        % Difference in noise texture
        t = norm(NPSref-NPS_to_be_matched,2)/norm(NPSref,2); 

        % MTF
        MTFref = interp1(reference.MTF(idxs_ref(ii_kernel_ref)).freq,reference.MTF(idxs_ref(ii_kernel_ref)).MTF/max(reference.MTF(idxs_ref(ii_kernel_ref)).MTF),(0:0.05:max([max(reference.MTF(idxs_ref(ii_kernel_ref)).freq),max(to_be_matched.MTF(idxs_to_be_matched(ii_kernel_matched)).freq)])));
        if ~isempty(to_be_matched.MTF(idxs_to_be_matched(ii_kernel_matched)).freq)
            MTF_to_be_matched = interp1(to_be_matched.MTF(idxs_to_be_matched(ii_kernel_matched)).freq,to_be_matched.MTF(idxs_to_be_matched(ii_kernel_matched)).MTF/max(to_be_matched.MTF(idxs_to_be_matched(ii_kernel_matched)).MTF),(0:0.05:max([max(reference.MTF(idxs_ref(ii_kernel_ref)).freq),max(to_be_matched.MTF(idxs_to_be_matched(ii_kernel_matched)).freq)])));
        
            MTFref(isnan(MTFref)) = 0;MTF_to_be_matched(isnan(MTF_to_be_matched)) = 0;

            r = norm(MTFref-MTF_to_be_matched,2)/norm(MTFref,2);
        else
            r = 1e3;
        end
        % Noise magnitude term
        % mval = (reference.NPS(idxs_ref(ii_kernel_ref)).NPS_SD - to_be_matched.NPS(idxs_to_be_matched(ii_kernel_matched)).NPS_SD);

        %% Evaluate the matching function value
        D = 1./(1 + exp(-ad*(d+cd1))); % Dose term
        T =  (1./(1 + exp(-at*(t - ct))) - 1./(1 + exp(-at*(t + ct))) ); % Noise texture
%         R = 1./(1 + exp(-ar*(r - cr))); % Resolution term from winslow
        % Modified resolution term using a sigmoid function
        R =  (1./(1 + exp(-ar*(r - cr))) - 1./(1 + exp(-ar*(r + cr))) ); 
        % Include noise magnitude term
        %M = (1./(1 + exp(-am*(mval - cm))) - 1./(1 + exp(-am*(mval + cm))) );
        %% Store results
        Ds(ii_kernel_ref,ii_kernel_matched) = D;
        Rs(ii_kernel_ref,ii_kernel_matched) = R;
        Ts(ii_kernel_ref,ii_kernel_matched) = T;
        %Ms(ii_kernel_ref,ii_kernel_matched) = M;
        m(ii_kernel_ref,ii_kernel_matched) = D*R*T;
    end
    
    
    % Identify matching kernel for each target FBP kernel
matched_idx = [];
if sum(m) ~= 0
    for ii_row = 1:size(m,1)
        idxs = find(m(ii_row,:) == max(m(ii_row,:)));
        matched_idx(ii_row) = idxs(1);
    end
    
    m_value = m(matched_idx);
    matching_kernel = idxs_to_be_matched(matched_idx);
    handles.idxs_to_be_matched = idxs_to_be_matched;
    handles.matching_kernel = matching_kernel;
    handles.m_vector = m;
    handles.m_value = m_value;
    
    fig_nro = handles.fig_nro;
    
    handles = show_imag_axes2(to_be_matched, matching_kernel, handles, fig_nro);
    
    set(handles.text18,'String',to_be_matched.Recon_names(matching_kernel));
    set(handles.text19,'String',to_be_matched.Body_part(matching_kernel));
    set(handles.text_m,'String',m_value);
    set(handles.text24,'String',to_be_matched.Slice_thickness(matching_kernel));
    set(handles.text27,'String',round(to_be_matched.CTDIvol(matching_kernel)));
    
    idxs_to_menu = idxs_to_be_matched(m > 0.1);
    m2 = m(m > 0.1);
    
    [m_sort morder] = sort(m2,'descend');
    idxs_to_menu_sort = idxs_to_menu(morder);
    if ~isempty(idxs_to_menu_sort)
        %set(handles.matchres_menu,'String','')
        set(handles.matchres_menu, 'String', to_be_matched.Recon_names(idxs_to_menu_sort));
    else
        set(handles.error_text,'String','No good match found')
        set(handles.matchres_menu, 'String', 'Bad result');
    end
    
    handles.idxs_to_menu_sort = idxs_to_menu_sort;
    handles.m_sort = m_sort;
    handles.match_rank = 1;
      
   
else
    matching_kernel = [];   % leave empty if no match found
    handles.matching_kernel = matching_kernel;
    set(handles.error_text,'String','No matching kernel found')
    set(handles.matchres_menu, 'String', [{'No match...'}]);
    cla(handles.axes2)
    set(handles.text18,'String','')
    set(handles.text19,'String','')
    set(handles.text_m,'String','')
    set(handles.text24,'String','')
    set(handles.text27,'String','')
end
else
    matching_kernel = [];   % leave empty if no target image found
    set(handles.text18,'String','')
    set(handles.text19,'String','')
    set(handles.text_m,'String','')
    set(handles.text24,'String','')
    set(handles.text27,'String','')
    set(handles.matchres_menu, 'String', [{'No match...'}]);
end
else
matching_kernel = [];   % leave empty if no target image found
end


% Find the inputted kernel and its match
function handles = show_kernel_and_match(handles)

if isfield(handles,'model_struct')
    fig_nro = handles.fig_nro;
    
    slice_thickness = handles.slice_thickness;
    model_struct = handles.model_struct;
    CTDIvol = handles.CTDIvol;
    body_part = handles.body_part;
    kerneli = handles.kerneli;
    target_image_indx = find_image(model_struct, slice_thickness, CTDIvol, body_part, kerneli);
    handles.target_image_indx = target_image_indx;
    
    
    if isempty(target_image_indx)   % indx is empty
        set(handles.error_text,'String','No such combination found, please choose some othet parameters')
        axes(handles.axes1);
        cla
        axes(handles.axes2);
        cla
        set(handles.text18,'String','')
        set(handles.text19,'String','')
        set(handles.text_m,'String','')
        set(handles.text24,'String','')
        set(handles.text27,'String','')
        
    else
        handles = show_imag_axes1(model_struct, target_image_indx, handles, fig_nro);
        
        if isfield(handles, 'target_CTmodel') % If matching model is defined
            if any(strcmp(handles.list_models,handles.target_CTmodel))
                % ACTUAL MATCHING
                
                CTmodel_name = handles.target_CTmodel;
                [matching_kernel, handles] = match_kernel(handles,CTmodel_name);
                set(handles.text2,'String',matching_kernel);
            end
        end
        
    end
end



function [handles, target_kernels] = show_target_kernel_list(handles)

 model_struct = handles.model_struct;
 body_part = handles.body_part;
 CTDIvol = handles.CTDIvol;
 slice_thickness = handles.slice_thickness;
    
    % Find kernel values for matching bodypart, and input them into corresponding popupmenu
    idxs_for_kernel = find(~cellfun(@isempty,strfind(model_struct.Body_part,body_part)));
    idxs_for_kernel2 = find(round(model_struct.CTDIvol) == CTDIvol);
    idxs_for_kernel3 = find(model_struct.Slice_thickness == slice_thickness);
    
    idintersect1 = intersect(idxs_for_kernel,idxs_for_kernel2);
    idintersect = intersect(idintersect1,idxs_for_kernel3);
      
    kernel_values = unique(model_struct.Recon_names(idintersect));
    
    % Empty possible previous inputs
    set(handles.kernelmenu, 'Value', 1);
    
    % Set the new values
    set(handles.kernelmenu, 'String', kernel_values);
    
    if isfield(handles,'kerneli')
        kerneli = handles.kerneli;
        % Check if current kernel value is not found and input default if not
        if ~any(strcmp(kernel_values,kerneli))
            kerneli = kernel_values(1);
            handles.kerneli = kerneli;
        end
        set(handles.chosen_kernel, 'String', kerneli);
    end
   
    
    target_kernels = kernel_values;


% --- Executes on selection change in matchres_menu.
function matchres_menu_Callback(hObject, eventdata, handles)
% show matching results

% hObject    handle to matchres_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns matchres_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from matchres_menu

match_rank = get(hObject,'Value');
match_rank_contents = cellstr(get(hObject,'String'));
match_rank_name = match_rank_contents{get(hObject,'Value')};
handles.match_rank = match_rank;

%set(handles.text5,'String',CTmodel_name);
%handles.target_CTmodel = CTmodel_name;
to_be_matched = handles.matching_model_struct;

if ~any(strcmp(to_be_matched.Recon_names,match_rank_name))
    %axes(handles.axes2);
    %cla
else
    fig_nro = handles.fig_nro;
    idxs_to_menu_sort = handles.idxs_to_menu_sort;
    m_sort = handles.m_sort;
    matching_kernel = idxs_to_menu_sort(match_rank);
    handles.matching_kernel = matching_kernel;
    
    handles = show_imag_axes2(to_be_matched, matching_kernel, handles, fig_nro);
    
    set(handles.text18,'String',to_be_matched.Recon_names(matching_kernel));
    set(handles.text19,'String',to_be_matched.Body_part(matching_kernel));
    set(handles.text_m,'String',m_sort(match_rank));
    set(handles.text24,'String',to_be_matched.Slice_thickness(matching_kernel));
    set(handles.text27,'String',round(to_be_matched.CTDIvol(matching_kernel)));
    
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function matchres_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to matchres_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Utility functions
function executableFolder = GetExecutableFolder() 
try
    if isdeployed 
        % User is running an executable in standalone mode. 
        executableFolder = ctfroot;
% 			fprintf(1, '\nIn function GetExecutableFolder(), currentWorkingDirectory = %s\n', executableFolder);
    else
        % User is running an m-file from the MATLAB integrated development environment (regular MATLAB).
        executableFolder = pwd; 
    end 
catch ME
    errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
        ME.stack(1).name, ME.stack(1).line, ME.message);
    uiwait(warndlg(errorMessage));
end
return;
