function varargout = trainint_launcher(varargin)
% TRAININT_LAUNCHER MATLAB code for trainint_launcher.fig
%      TRAININT_LAUNCHER, by itself, creates a new TRAININT_LAUNCHER or raises the existing
%      singleton*.
%
%      H = TRAININT_LAUNCHER returns the handle to a new TRAININT_LAUNCHER or the handle to
%      the existing singleton*.
%
%      TRAININT_LAUNCHER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRAININT_LAUNCHER.M with the given input arguments.
%
%      TRAININT_LAUNCHER('Property','Value',...) creates a new TRAININT_LAUNCHER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before trainint_launcher_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to trainint_launcher_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help trainint_launcher

% Last Modified by GUIDE v2.5 17-Apr-2013 18:10:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @trainint_launcher_OpeningFcn, ...
                   'gui_OutputFcn',  @trainint_launcher_OutputFcn, ...
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


% --- Executes just before trainint_launcher is made visible.
function trainint_launcher_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to trainint_launcher (see VARARGIN)

% Choose default command line output for trainint_launcher
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes trainint_launcher wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = trainint_launcher_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
