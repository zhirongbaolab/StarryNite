function varargout = lineage_launcher_v2(varargin)
%LINEAGE_LAUNCHER_V2 M-file for lineage_launcher_v2.fig
%      LINEAGE_LAUNCHER_V2, by itself, creates a new LINEAGE_LAUNCHER_V2 or raises the existing
%      singleton*.
%
%      H = LINEAGE_LAUNCHER_V2 returns the handle to a new LINEAGE_LAUNCHER_V2 or the handle to
%      the existing singleton*.
%
%      LINEAGE_LAUNCHER_V2('Property','Value',...) creates a new LINEAGE_LAUNCHER_V2 using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to lineage_launcher_v2_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      LINEAGE_LAUNCHER_V2('CALLBACK') and LINEAGE_LAUNCHER_V2('CALLBACK',hObject,...) call the
%      local function named CALLBACK in LINEAGE_LAUNCHER_V2.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to runbutton (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help lineage_launcher_v2

% Last Modified by GUIDE v2.5 07-Jan-2019 14:35:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @lineage_launcher_v2_OpeningFcn, ...
                   'gui_OutputFcn',  @lineage_launcher_v2_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before lineage_launcher_v2 is made visible.
function lineage_launcher_v2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no outputname args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line outputname for lineage_launcher_v2
handles.outputname = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes lineage_launcher_v2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);
handles.ROIcount=0;
handles.ROIs={};
guidata(hObject,handles)

% --- Outputs from this function are returned to the command line.
function varargout = lineage_launcher_v2_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning outputname args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line outputname from handles structure
varargout{1} = handles.outputname;


% --- Executes on button press in runbutton.
function runbutton_Callback(hObject, eventdata, handles)
% hObject    handle to runbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

suffix=get(handles.embryosuffix,'String');
parameterfilenames={};
imagefilename=get(handles.imagefilename,'String');

embryodir=handles.embryodir;
wrongslashes=findstr(embryodir,'\');
embryodir(wrongslashes)='/';
outputdir=get(handles.outputdir,'String');
if(iscell(outputdir))
    outputdir=outputdir{1};
end
wrongslashes=findstr(outputdir,'\');
outputdir(wrongslashes)='/';
if ~isempty(outputdir)
    outputdir=[outputdir,'/'];
end

if isempty(outputdir)
    %if none specified use base image location pluus image base name
    outputdirectory=[embryodir,handles.embryodirname,'/'];
else
    %if specified use specified plust image base name
    outputdirectory=[outputdir,handles.embryodirname,'/'];
end
%make this directory
mkdir(outputdirectory);

for i=1:length(handles.ROIs)
%for each embryo ROI
%make copy of param file with emb name
%append configuration matching gui to it
%first one gets sliceoutput=true
%close it
%paramfilename=[outputdirectory,'/paramfile_',handles.embryoname,suffix,num2str(i,'%04d'),'.txt'];


paramfilename=[outputdirectory,'paramfile_',handles.embryodirname,suffix,num2str(i,'%04d'),'.txt'];

%paramfilename=[get(handles.outputdir,'String'),'\paramfile_',num2str(i,'%04d'),'_',suffix,'.txt'];
wrongslashes=findstr(paramfilename,'\');
paramfilename(wrongslashes)='/';

parameterfilenames{i}=paramfilename;
copyfile(get(handles.parameterfilename,'String'),paramfilename);

file=fopen(paramfilename,'a');

if(get(handles.splitimage,'Value'))
    fprintf(file,'splitstack=true;\n\r ');
%    fprintf(file,'newscope=true;\n\r');
%    fprintf(file,'SIMPLETIFF=false;\n\r');
     
else
    %{
    if(get(handles.matlabim,'Value'))
        fprintf(file,'newscope=false;\n\r');
        fprintf(file,'SIMPLETIFF=false;\n\r');
        fprintf(file,'MATLAB_STACK=true;\n\r');
        
    else
    %}
        %if(get(handles.simpletiff,'Value'))
       fprintf(file,'splitstack=false;\n\r ');

 %           fprintf(file,'newscope=false;\n\r');
           % fprintf(file,'MATLAB_STACK=false;\n\r');
 %           fprintf(file,'SIMPLETIFF=true;\n\r');

        %end
    %end 
end

fprintf(file,'%Parameter overwrites generated by ROI interface:\n\r');
if(get(handles.splitimage,'Value'))
    if(get(handles.green,'Value'))
        fprintf(file,'rednuclei=false;\n\r ');
       % fprintf(file,'LSM_channel=1;\n\r');
    else
        fprintf(file,'rednuclei=true;\n\r ');
       % fprintf(file,'LSM_channel=2;\n\r');
    end
end

if(get(handles.flipimagesbox,'Value'))
     fprintf(file,'flipstack=true;\n\r ');

else
     fprintf(file,'flipstack=false;\n\r ');

end

fprintf(file,['start_time=',get(handles.starttime,'String'),';\n\r ']);
fprintf(file,['end_time=',get(handles.endtime,'String'),';\n\r ']);

if(i==1&&get(handles.makeslices,'Value'))
    fprintf(file,'outputSlice=true;\n\r ');
else
    fprintf(file,'outputSlice=false;\n\r ');
end

%output ROI i
fprintf(file,'ROI=true;\n\r ');
points=round(handles.ROIs{i}.getPosition());
fprintf(file,['ROIxmin=',num2str(max(1,min(points(:,1)))),';\n\r ']);
fprintf(file,['ROIxmax=',num2str(max(points(:,1))),';\n\r ']);
fprintf(file,['ROIymin=',num2str(max(1,min(points(:,2)))),';\n\r ']);
fprintf(file,['ROIymax=',num2str(max(points(:,2))),';\n\r ']);


%build up a string for the polygonal ROI to put in parameter file
bigstring='ROIpoints=[';
for p=1:size(points,1);
    bigstring=[bigstring,num2str(points(p,1)),' ',num2str(points(p,2)),' ; '];
end
bigstring=[bigstring,']'];
fprintf(file,[bigstring,';\n\r ']);%save point array in parameter file

fclose(file);
end

for i=1:length(handles.ROIs)
%for each embryo ROI
%call main detction with images 
%call tracking
%put everything in right final location
%this is done with call to central driver script
%{
embryodir=handles.embryodir;
wrongslashes=findstr(embryodir,'\');
embryodir(wrongslashes)='/';
%}
outputdir=get(handles.outputdir,'String');
if(iscell(outputdir))
    outputdir=outputdir{1};
end
wrongslashes=findstr(outputdir,'\');
outputdir(wrongslashes)='/';
outputdir=[outputdir,'/'];
points=round(handles.ROIs{i}.getPosition());
tic
%%detect_track_driver_allmatlab(parameterfilenames{i},embryodir,handles.embryoname,[suffix,num2str(i)],outputdirectory,points,false);%,lineageparameterfile);
%this is the good one
detect_track_driver_allmatlab_v2(parameterfilenames{i},imagefilename,[suffix,num2str(i)],outputdirectory,points,false);

toc
end
'all embryos completed'

% --- Executes on button press in reset.
function reset_Callback(hObject, eventdata, handles)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when selected object is changed in colorchannel.
function colorchannel_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in colorchannel 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
%if changed update to one of default param files for red or green 
%decided no longer needed
%if(get(handles.green,'Value'))
%set(handles.parameterfilename,'String','L:\bin\starryniteII\matlab-parameters-file-newscope_middle_weak_thresholds_journalV_besthack_nodata.txt');
%else
%    set(handles.parameterfilename,'String','./');
%end


function embryosuffix_Callback(hObject, eventdata, handles)
% hObject    handle to embryosuffix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of embryosuffix as text
%        str2double(get(hObject,'String')) returns contents of embryosuffix as a double



% --- Executes during object creation, after setting all properties.
function embryosuffix_CreateFcn(hObject, eventdata, handles)
% hObject    handle to embryosuffix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function parameterfilename_Callback(hObject, eventdata, handles)
% hObject    handle to parameterfilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of parameterfilename as text
%        str2double(get(hObject,'String')) returns contents of parameterfilename as a double



% --- Executes during object creation, after setting all properties.
function parameterfilename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to parameterfilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in parameterbutton.
function parameterbutton_Callback(hObject, eventdata, handles)
% hObject    handle to parameterbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName,FilterIndex] = uigetfile('*.*');
set(handles.parameterfilename,'String',[PathName,FileName]);



function imagefilename_Callback(hObject, eventdata, handles)
% hObject    handle to imagefilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of imagefilename as text
%        str2double(get(hObject,'String')) returns contents of imagefilename as a double


 
% --- Executes during object creation, after setting all properties.
function imagefilename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imagefilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in imagebutton.
function imagebutton_Callback(hObject, eventdata, handles)
% hObject    handle to imagebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName,FilterIndex] = uigetfile('l:/*.*');

set(handles.imagefilename,'String',[PathName,FileName]);
ind=max(strfind(FileName,'_'));
if (isempty(ind))
    ind=regexp(FileName,'\d' ,'once');
end

ind2=max(findstr(FileName,'.'));
framenum=str2num(FileName(ind+2:ind2-1));
ind3=(findstr(PathName,filesep));

ind3=ind3(length(ind3)-1);%take second to last occurance of dir delimater
handles.embryodirname=PathName(ind3+1:length(PathName)-1);%end dir name
handles.embryodir=PathName;%vs full path

%{
%looks for the 


set(handles.imagefilename,'String',[PathName,FileName(1:ind-1)]);
suffix=FileName(ind2+1:length(FileName));
 handles.lsm=false;

if (strcmp(suffix,'mat'))
    load([PathName,FileName]);
    X=stack;
    clear stack;
    handles.matlabformat=true;
    im=max(X,[],3);
else
    if(strcmp(suffix,'lsm'))
        handles.lsm=true;
         handles.matlabformat=false;
         set(handles.imagefilename,'String',[PathName,FileName]);
         X=(loadCellStackLSMtime([PathName,FileName],1,1,inf));
        im=max(X,[],3);
    else
  %}      

if(get(handles.splitimage,'Value'))
    if(get(handles.green,'Value'))
        X=im2double(((loadCellStackMetamorph([PathName,FileName(1:ind-1)],framenum,2,30,[0,0,0,0],false))));
    else
        X=im2double(((loadCellStackMetamorph([PathName,FileName(1:ind-1)],framenum,1,30,[0,0,0,0],false))));
    end
    im=max(X,[],3);
    im=fliplr(im);
else
    X=loadSimpleStackTiff([PathName,FileName]);
    %simple tiff or slice
    if(size(X,3)>1)
        im=max(X,[],3);
    else
        im=X;
    end
end
%{
if (get(handles.simpletiff,'Value'))
            X=loadSimpleStackTiff([PathName,FileName]);
            'simple tiff'
            im=max(X,[],3);
        else
            if(get(handles.green,'Value'))
                X=im2double(((loadCellStackMetamorph([PathName,FileName(1:ind-1)],framenum,2,30,[0,0,0,0],false))));
            else
                X=im2double(((loadCellStackMetamorph([PathName,FileName(1:ind-1)],framenum,1,30,[0,0,0,0],false))));
            end
            im=max(X,[],3);
            im=fliplr(im);
        end
    end
end
%}
minvalr=prctile(reshape(im,[1,numel(im)]),25);
maxvalr=prctile(reshape(im,[1,numel(im)]),99); 
im(im>maxvalr)=maxvalr;
im(im<minvalr)=minvalr;

hold on;
i=imagesc(im);
axis tight;
set(handles.axes1,'YDir','reverse');
colormap(gray);
set(i,'HitTest','off');
%if(isfield(handles,'lsm')&&handles.lsm)
%    handles.embryoname=FileName;
%else
%handles.embryoname=FileName(1:ind-1);
%end
guidata(hObject,handles)


% --- Executes on mouse press over axes background.
function axes1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.ROIcount=handles.ROIcount+1;
handles.ROIs{handles.ROIcount}=impoly(handles.axes1);
guidata(hObject,handles)



function endtime_Callback(hObject, eventdata, handles)
% hObject    handle to endtime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of endtime as text
%        str2double(get(hObject,'String')) returns contents of endtime as a double


% --- Executes during object creation, after setting all properties.
function endtime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to endtime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function starttime_Callback(hObject, eventdata, handles)
% hObject    handle to starttime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of starttime as text
%        str2double(get(hObject,'String')) returns contents of starttime as a double


% --- Executes during object creation, after setting all properties.
function starttime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to starttime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in outputdirbutton.
function outputdirbutton_Callback(hObject, eventdata, handles)
% hObject    handle to outputdirbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path = uigetdir('L:\');
set(handles.outputdir,'String',path);


function outputdir_Callback(hObject, eventdata, handles)
% hObject    handle to outputdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of outputdir as text
%        str2double(get(hObject,'String')) returns contents of outputdir as a double


% --- Executes during object creation, after setting all properties.
function outputdir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to outputdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in makeslices.
function makeslices_Callback(hObject, eventdata, handles)
% hObject    handle to makeslices (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of makeslices


% --- Executes on button press in clearROIbutton.
function clearROIbutton_Callback(hObject, eventdata, handles)
% hObject    handle to clearROIbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%clear all existing ROIS
for i=1:length (handles.ROIs)
    handles.ROIs{i}.delete;
end
handles.ROIcount=0;
handles.ROIs={};
guidata(hObject,handles)

% --- Executes on button press in flipimagesbox.
function flipimagesbox_Callback(hObject, eventdata, handles)
% hObject    handle to flipimagesbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of flipimagesbox


% --- Executes on button press in splitimage.
function splitimage_Callback(hObject, eventdata, handles)
% hObject    handle to splitimage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of splitimage
if(get(handles.splitimage,'Value'))
    set(handles.flipimagesbox,'Value',true);
else
     set(handles.flipimagesbox,'Value',false);
end
