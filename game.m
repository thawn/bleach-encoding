function varargout = game(varargin)
% game MATLAB code for game.fig
%      game, by itself, creates a new game or raises the existing
%      singleton*.
%
%      H = game returns the handle to a new game or the handle to
%      the existing singleton*.
%
%      game('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in game.M with the given input arguments.
%
%      game('Property','Value',...) creates a new game or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before game_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to game_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help game

% Last Modified by GUIDE v2.5 08-Jun-2018 16:26:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @game_OpeningFcn, ...
                   'gui_OutputFcn',  @game_OutputFcn, ...
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


% --- Executes just before game is made visible.
function game_OpeningFcn(hObject, eventdata, handles, varargin )
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to game (see VARARGIN)

% Choose default command line output for game
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes game wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = game_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5


% --- Executes on button press in checkbox6.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox6


% --- Executes on button press in checkbox7.
function checkbox7_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox7


% --- Executes on button press in checkbox8.
function checkbox8_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox8


% --- Executes on button press in checkbox9.
function checkbox9_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox9


% --- Executes on button press in checkbox10.
function checkbox10_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox10


% --- Executes on button press in checkbox11.
function checkbox11_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox11


% --- Executes on button press in checkbox12.
function checkbox12_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox12


% --- Executes on button press in checkbox13.
function checkbox13_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox13


% --- Executes on button press in checkbox14.
function checkbox14_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox14


% --- Executes on button press in checkbox15.
function checkbox15_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox15


% --- Executes on button press in checkbox16.
function checkbox16_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox16


% --- Executes on button press in checkbox17.
function checkbox17_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox17


% --- Executes on button press in checkbox18.
function checkbox18_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox18


% --- Executes on button press in checkbox19.
function checkbox19_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox19


% --- Executes on button press in checkbox20.
function checkbox20_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox20


% --- Executes on button press in checkbox21.
function checkbox21_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox21


% --- Executes on button press in checkbox22.
function checkbox22_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox22


% --- Executes on button press in checkbox23.
function checkbox23_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox23


% --- Executes on button press in checkbox24.
function checkbox24_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox24


% --- Executes on button press in checkbox25.
function checkbox25_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox25


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val1 = get(handles.checkbox1,'value');
val2 = get(handles.checkbox2,'value');
val3 = get(handles.checkbox3,'value');
val4 = get(handles.checkbox4,'value');
val5 = get(handles.checkbox5,'value');
val6 = get(handles.checkbox6,'value');
val7 = get(handles.checkbox7,'value');
val8 = get(handles.checkbox8,'value');
val9 = get(handles.checkbox9,'value');
val10 = get(handles.checkbox10,'value');
val11 = get(handles.checkbox11,'value');
val12 = get(handles.checkbox12,'value');
val13 = get(handles.checkbox13,'value');
val14 = get(handles.checkbox14,'value');
val15 = get(handles.checkbox15,'value');
val16 = get(handles.checkbox16,'value');
val17 = get(handles.checkbox17,'value');
val18 = get(handles.checkbox18,'value');
val19 = get(handles.checkbox19,'value');
val20 = get(handles.checkbox20,'value');
val21 = get(handles.checkbox21,'value');
val22 = get(handles.checkbox22,'value');
val23 = get(handles.checkbox23,'value');
val24 = get(handles.checkbox24,'value');
val25 = get(handles.checkbox25,'value');
Pattern1 = evalin( 'base','Pattern1');
Pattern2 = evalin( 'base','Pattern2');
Pattern3 = evalin( 'base','Pattern3');
Pattern4 = evalin( 'base','Pattern4');
Pattern5 = evalin( 'base','Pattern5');
Pattern6 = evalin( 'base','Pattern6');
Pattern7 = evalin( 'base','Pattern7');
Pattern8 = evalin( 'base','Pattern8');
Pattern9 = evalin( 'base','Pattern9');
Pattern10 = evalin( 'base','Pattern10');
Pattern11 = evalin( 'base','Pattern11');
Pattern12 = evalin( 'base','Pattern12');
Pattern13 = evalin( 'base','Pattern13');
Pattern14 = evalin( 'base','Pattern14');
Pattern15 = evalin( 'base','Pattern15');
Pattern16 = evalin( 'base','Pattern16');
Pattern17 = evalin( 'base','Pattern17');
Pattern18 = evalin( 'base','Pattern18');
Pattern19 = evalin( 'base','Pattern19');
Pattern20 = evalin( 'base','Pattern20');
Pattern21 = evalin( 'base','Pattern21');
Pattern22 = evalin( 'base','Pattern22');
Pattern23 = evalin( 'base','Pattern23');
Pattern24 = evalin( 'base','Pattern24');
Pattern25 = evalin( 'base','Pattern25');
output = [];
output1 = [];
meanintensity = [];
times =10;
for i = 1:times
    BES = BleachEncodingSim('MicrotubuleLength',10000);
    if val1 == 1
       BES.bleach(Pattern1, 0.05);
    end
    if val2 == 1
       BES.bleach(Pattern2, 0.05);
    end
    if val3 == 1
       BES.bleach(Pattern3, 0.05);
    end
    if val4 == 1
       BES.bleach(Pattern4, 0.05);
    end
    if val5 == 1
       BES.bleach(Pattern5, 0.05);
    end
    if val6 == 1
       BES.bleach(Pattern6, 0.05);
    end
    if val7 == 1
       BES.bleach(Pattern7, 0.05);
    end
    if val8 == 1
       BES.bleach(Pattern8, 0.05);
    end
    if val9 == 1
       BES.bleach(Pattern9, 0.05);
    end
    if val10 == 1
       BES.bleach(Pattern10, 0.05); 
    end
    if val11 == 1
       BES.bleach(Pattern11, 0.05); %0.05
    end
    if val12 == 1
       BES.bleach(Pattern12, 0.05); %0.05
    end
    if val13 == 1
       BES.bleach(Pattern13, 0.05); %0.05
    end
    if val14 == 1
       BES.bleach(Pattern14, 0.05); %0.05
    end
    if val15 == 1
       BES.bleach(Pattern15, 0.05); %0.05
    end
    if val16 == 1
       BES.bleach(Pattern16, 0.05); %0.06
    end
    if val17 == 1
       BES.bleach(Pattern17, 0.07); %0.07
    end
    if val18 == 1
       BES.bleach(Pattern18, 0.06); %0.07
    end
    if val19 == 1
       BES.bleach(Pattern19, 0.03); 
    end
    if val20 == 1
       BES.bleach(Pattern20, 0.07); 
    end
    if val21 == 1
       BES.bleach(Pattern21, 0.09); 
    end
    if val22 == 1
       BES.bleach(Pattern22, 0.08); 
    end
    if val23 == 1
       BES.bleach(Pattern23, 0.07); 
    end
    if val24 == 1
       BES.bleach(Pattern24, 0.06); 
    end
    if val25 == 1
       BES.bleach(Pattern25, 0.05);
    end
    [output,meanintensity] = BES.plotFFTAnalysis(i,output,meanintensity);
end
fig1 = figure('Position',[200 150 1000 600]);
ax1 = axes(fig1);
if times == 1
   plot(output,'LineStyle','none','Marker','*');
else
   boxplot(output);
end 
grid(ax1,'on');
set(ax1,'XAxisLocation','top','XTick',...
    [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30]);
annotation('textbox',[0.7 0.8 0.2 0.1],'VerticalAlignment','middle','String',...
{sprintf('mean intensity: %10.0f', mean(meanintensity))});
fig2 = figure('Position',[400 150 1000 600]);
ax2 = axes(fig2);
BES.plotLineScan;

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.checkbox1,'value',0);
set(handles.checkbox2,'value',0);
set(handles.checkbox3,'value',0);
set(handles.checkbox4,'value',0);
set(handles.checkbox5,'value',0);
set(handles.checkbox6,'value',0);
set(handles.checkbox7,'value',0);
set(handles.checkbox8,'value',0);
set(handles.checkbox9,'value',0);
set(handles.checkbox10,'value',0);
set(handles.checkbox11,'value',0);
set(handles.checkbox12,'value',0);
set(handles.checkbox13,'value',0);
set(handles.checkbox14,'value',0);
set(handles.checkbox15,'value',0);
set(handles.checkbox16,'value',0);
set(handles.checkbox17,'value',0);
set(handles.checkbox18,'value',0);
set(handles.checkbox19,'value',0);
set(handles.checkbox20,'value',0);
set(handles.checkbox21,'value',0);
set(handles.checkbox22,'value',0);
set(handles.checkbox23,'value',0);
set(handles.checkbox24,'value',0);
set(handles.checkbox25,'value',0);
