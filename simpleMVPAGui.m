%% Building blocks for a GUI to help run plmvpaLite without scripting
% for reference: https://blogs.mathworks.com/pick/2007/12/28/matlab-basics-guis-without-guide/

function simpleMVPAGui

h.fig = figure('position', [1100 30 360 60])

%% run classification
h.buttonMVPA = uicontrol('style', 'pushbutton', ...
    'position', [15 10 100 40] , ...
    'string', 'Classification');

%add callback
set(h.buttonMVPA, 'callback', {@addButtonMVPAtype, h})

%% run RSA
h.buttonRSA = uicontrol('style', 'pushbutton', ...
    'position', [245 10 100 40] , ...
    'string', 'RSA');

%add callback
set(h.buttonRSA, 'callback', {@addButtonRSAtype, h})




function h = addButtonMVPAtype(hObject, eventdata, h)

h.buttonMVPAtype = uicontrol('style', 'pushbutton', ...
    'position', [130 10 100 40] , ...
    'string', 'Type'); 
set(h.buttonMVPAtype, 'callback', {@removeButton, h});
set(h.buttonMVPA, 'enable', 'off'); %turn off option to click on MVPA

function h = removeButton(hObject, eventdata, h)

delete(h.buttonMVPAtype)
h = rmfield(h, 'buttonMVPAtype');
set(h.buttonMVPA, 'enable', 'on');