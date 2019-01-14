%% Building blocks for a GUI to help run plmvpaLite without scripting
% for reference: https://blogs.mathworks.com/pick/2007/12/28/matlab-basics-guis-without-guide/

function simpleMVPAGui

h.fig = figure('position', [1100 30 360 60])

%% run classification
h.buttonMVPA = uicontrol('style', 'pushbutton', ...
    'position', [15 10 100 40] , ...
    'string', 'Classification');

%add callback
set(h.buttonMVPA, 'callback', {@addButton, h})

%% run RSA
h.buttonRSA = uicontrol('style', 'pushbutton', ...
    'position', [245 10 100 40] , ...
    'string', 'RSA');

%add callback
set(h.buttonRSA, 'callback', {@addButton, h})




function h = addButton(hObject, eventdata, h)

h.buttonTwo = uicontrol('style', 'pushbutton', ...
    'position', [130 10 100 40] , ...
    'string', 'Remove button'); 
set(h.buttonTwo, 'callback', {@removeButton, h});
set(h.buttonOne, 'enable', 'off');

function h = removeButton(hObject, eventdata, h)

delete(h.buttonTwo)
h = rmfield(h, 'buttonTwo');
set(h.buttonOne, 'enable', 'on');