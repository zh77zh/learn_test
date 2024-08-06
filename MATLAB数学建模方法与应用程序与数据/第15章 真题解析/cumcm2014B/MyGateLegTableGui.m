function MyGateLegTableGui

OldHandle = findobj( 'Type', 'figure', 'Tag', 'MyGateLegTable' ) ;
if ishandle( OldHandle )
    close( OldHandle ) ;
end
fig = figure ;
init_MyGateLegTable(fig) ; 

%--------------------------------------------------------------------------
function init_MyGateLegTable(fig)

set(fig,'units','normalized','position',[0.0344,0.0924,0.9275,0.8359],...
    'menubar','none','name','折叠桌模拟制作系统',...
    'numbertitle','off','color',[0.925 0.914 0.847],'tag','MyGateLegTable');
axes('pos',[0.063 0.34 0.6 0.618],'tag','axes1');
axis([-25,25,0,50]);
rwm = uicontextmenu;
uimenu(rwm,'label','画网格','callback','grid on');
uimenu(rwm,'label','清除网格','callback','grid off');
line_width=uimenu(rwm,'label','设置线宽');
uimenu(line_width,'label','0.5','callback',...
    ['hand_line=findall(gcf,''type'',''line'');',...
    'set(hand_line,''linewidth'',0.5);']);
uimenu(line_width,'label','1','callback',...
    ['hand_line=findall(gcf,''type'',''line'');',...
    'set(hand_line,''linewidth'',1);']);
uimenu(line_width,'label','2','callback',...
    ['hand_line=findall(gcf,''type'',''line'');',...
    'set(hand_line,''linewidth'',2);']);
uimenu(line_width,'label','3','callback',...
    ['hand_line=findall(gcf,''type'',''line'');',...
    'set(hand_line,''linewidth'',3);']);
uimenu(line_width,'label','4','callback',...
    ['hand_line=findall(gcf,''type'',''line'');',...
    'set(hand_line,''linewidth'',4);']);
uimenu(line_width,'label','6','callback',...
    ['hand_line=findall(gcf,''type'',''line'');',...
    'set(hand_line,''linewidth'',6);']);
uimenu(line_width,'label','8','callback',...
    ['hand_line=findall(gcf,''type'',''line'');',...
    'set(hand_line,''linewidth'',8);']);
uimenu(line_width,'label','10','callback',...
    ['hand_line=findall(gcf,''type'',''line'');',...
    'set(hand_line,''linewidth'',10);']);
uimenu(rwm,'label','复制当前图像','callback',@CopyFig);
uimenu(rwm,'label','保存当前图像','callback',@SaveFig);
uimenu(rwm,'label','打印当前图像','callback',@PrintFig);
set(gca,'uicontextmenu',rwm);

panel0 = uipanel(fig,'units','normalized',...
    'pos',[0.038 0.002 0.346 0.265],'title','请输入相关参数',...
    'fontsize',10,'fontweight','bold','bac',[0.925 0.914 0.847]);
uicontrol(panel0,'style','text','units','normalized',...
    'pos',[0.03 0.76 0.2 0.15],'string','桌半宽 R = ',...
    'fontsize',12,'fontweight','bold','fontunits','normalized',...
    'HorizontalAlignment','left','bac',[0.925 0.914 0.847])
uicontrol(panel0,'style','edit','units','normalized','pos',[0.26 0.72 0.18 0.22],...
    'fontsize',12,'fontunits','normalized','string','25',...
    'tag','edit1','backgroundcolor',[1 1 1],'callback',@ResizeAxesRange);
uicontrol(panel0,'style','text','units','normalized',...
    'pos',[0.03 0.43 0.2 0.15],'string','桌高 H = ',...
    'fontsize',12,'fontweight','bold','fontunits','normalized',...
    'HorizontalAlignment','left','bac',[0.925 0.914 0.847])
uicontrol(panel0,'style','edit','units','normalized','pos',[0.26 0.39 0.18 0.22],...
    'fontsize',12,'fontunits','normalized','string','53',...
    'tag','edit2','backgroundcolor',[1 1 1]);
uicontrol(panel0,'style','text','units','normalized',...
    'pos',[0.03 0.07 0.2 0.15],'string','桌厚 t = ',...
    'fontsize',12,'fontweight','bold','fontunits','normalized',...
    'HorizontalAlignment','left','bac',[0.925 0.914 0.847])
uicontrol(panel0,'style','edit','units','normalized','pos',[0.26 0.03 0.18 0.22],...
    'fontsize',12,'fontunits','normalized','string','3',...
    'tag','edit3','backgroundcolor',[1 1 1]);

uicontrol(panel0,'style','text','units','normalized',...
    'pos',[0.52 0.76 0.2 0.15],'string','木条宽 w = ',...
    'fontsize',12,'fontweight','bold','fontunits','normalized',...
    'HorizontalAlignment','left','bac',[0.925 0.914 0.847])
uicontrol(panel0,'style','edit','units','normalized','pos',[0.77 0.72 0.18 0.22],...
    'fontsize',12,'fontunits','normalized','string','2.5',...
    'tag','edit4','backgroundcolor',[1 1 1]');
uicontrol(panel0,'style','text','units','normalized',...
    'pos',[0.52 0.43 0.2 0.15],'string','板半长 L = ',...
    'fontsize',12,'fontweight','bold','fontunits','normalized',...
    'HorizontalAlignment','left','bac',[0.925 0.914 0.847])
uicontrol(panel0,'style','edit','units','normalized','pos',[0.77 0.39 0.18 0.22],...
    'fontsize',12,'fontunits','normalized','string','',...
    'tag','edit5','backgroundcolor',[1 1 1]);
uicontrol(panel0,'style','text','units','normalized',...
    'pos',[0.52 0.07 0.2 0.15],'string','轴位置 b = ',...
    'fontsize',12,'fontweight','bold','fontunits','normalized',...
    'HorizontalAlignment','left','bac',[0.925 0.914 0.847])
uicontrol(panel0,'style','edit','units','normalized','pos',[0.77 0.03 0.18 0.22],...
    'fontsize',12,'fontunits','normalized','string','',...
    'tag','edit6','backgroundcolor',[1 1 1]);

panel1 = uipanel(fig,'units','normalized',...
    'pos',[0.39 0.002 0.28 0.265],'title','请输入边缘函数',...
    'fontsize',10,'fontweight','bold','bac',[0.925 0.914 0.847]);
uicontrol(panel1,'style','text','units','normalized',...
    'pos',[0.05 0.78 0.4 0.15],'string','f(x) = @(x)',...
    'fontsize',12,'fontweight','bold','fontunits','normalized',...
    'HorizontalAlignment','left','bac',[0.925 0.914 0.847])
uicontrol(panel1,'style','edit','units','normalized','pos',[0.06 0.46 0.88 0.25],...
    'fontsize',12,'fontunits','normalized','string','sqrt(R^2-x^2)',...
    'tag','edit7','backgroundcolor',[1 1 1],'callback',@PlotEdgeCurve);

uicontrol(panel1,'style','push','units','normalized',...
    'pos',[0.04 0.08 0.28 0.23],'string','手  绘',...
    'fontsize',12,'fontweight','bold',...
    'fontunits','normalized','callback',@SetHandDrawnFun)
uicontrol(panel1,'style','push','units','normalized',...
    'pos',[0.36 0.08 0.28 0.23],'string','参数优化',...
    'fontsize',12,'fontweight','bold',...
    'fontunits','normalized','callback',@CalculateOptimalParameters)
uicontrol(panel1,'style','push','units','normalized',...
    'pos',[0.68 0.08 0.28 0.23],'string','效果展示',...
    'fontsize',12,'fontweight','bold',...
    'fontunits','normalized','callback',@DynamicGraphDisplay)

panel2 = uipanel(fig,'units','normalized',...
    'pos',[0.704 0.002 0.286 0.974],'title','输出各木条开槽相关结果',...
    'fontsize',10,'fontweight','bold','bac',[0.925 0.914 0.847]);
uicontrol(panel2,'style','push','units','normalized',...
    'pos',[0.02 0.01 0.25 0.06],'string','旋转图形',...
    'fontsize',12,'fontweight','bold',...
    'fontunits','normalized','callback','rotate3d on')
uicontrol(panel2,'style','push','units','normalized',...
    'pos',[0.37 0.01 0.25 0.06],'string','关闭旋转',...
    'fontsize',12,'fontweight','bold',...
    'fontunits','normalized','callback','rotate3d off')
uicontrol(panel2,'style','push','units','normalized',...
    'pos',[0.73 0.01 0.25 0.06],'string','清除图形',...
    'fontsize',12,'fontweight','bold',...
    'fontunits','normalized','callback',@ClearAxes)

uitable(panel2,'units','normalized',...
    'position',[0.0168 0.08 0.9642 0.9115],...
    'data',cell(40,5),...
    'ColumnName',{'开槽起点','开槽终点','开槽长度'},...
    'ColumnEditable',true,...
    'tag','datatable','fontsize',10);

handles = guihandles(gcf);
guidata(gcf,handles);
%--------------------------------------------------------------------------
function ResizeAxesRange(~, ~)
handles = guihandles(gcf);
R = str2double(get(handles.edit1,'string'));
if isempty(R) || R == 0 || isnan(R)
    R = 25;
end
axis([-R,R,0,2*R]);

function ClearAxes(~, ~)
cla; 
view(2);
axis normal;
handles = guihandles(gcf);
R = str2double(get(handles.edit1,'string'));
if isempty(R) || R == 0 || isnan(R)
    R = 25;
end
axis([0,R,0,2*R]);
setappdata(gcf,'xy_point',[]);
setappdata(gcf,'Lb',[]);

function PlotEdgeCurve(~, ~)
handles = guihandles(gcf);
R = str2double(get(handles.edit1,'string'));
if isempty(R) || R == 0 || isnan(R)
    R = 25;
end
FunStr = get(handles.edit7,'string');
if isempty(FunStr)
    FunStr = 'sqrt(R^2-x.^2)';
end
FunStr = vectorize(FunStr);
EdgeFun = eval(['@(x)',FunStr]);
x = linspace(-R,R,100);
y = EdgeFun(x);
plot(x,y,'r','linewidth',2);

function SetHandDrawnFun(~, ~)
handles = guidata(gcf);
set(handles.axes1,'ButtonDownFcn',{@PlotEdgeCurveByHand 'start'});
set(handles.edit7,'string','');
R = str2double(get(handles.edit1,'string'));
if isempty(R) || R == 0 || isnan(R)
    R = 25;
end
cla
view(2);
axis normal;
axis([0,R,0,2*R]);

%%-------------------------------------------------------------------------
%鼠标手绘
%%-------------------------------------------------------------------------
function PlotEdgeCurveByHand(~, ~, action)

switch action
    case 'start'        
        point=get(gca,'CurrentPoint'); 
        xy_point=getappdata(gcf,'xy_point');
        xy_point=[xy_point;point(1,1:2)];
        setappdata(gcf,'xy_point',xy_point);
        line(point(1,1),point(1,2),'clipping','on',...
            'linewidth',3,'color','r');  % 'erasemode','background',
        data_point=point(1,1:2);        
        setappdata(gcf,'data_point',data_point);
        set(gcbf,'WindowButtonMotionFcn',{@PlotEdgeCurveByHand 'move'});
        set(gcbf,'WindowButtonUpFcn',{@PlotEdgeCurveByHand 'stop'});        
    case 'move'
        data_point=getappdata(gcf,'data_point');
        point=get(gca,'CurrentPoint'); 
        xy_point=getappdata(gcf,'xy_point');
        xy_point=[xy_point;point(1,1:2)];
        setappdata(gcf,'xy_point',xy_point);
        line([data_point(1),point(1,1)],[data_point(2),point(1,2)],...
            'clipping','on',...  % 'erasemode','background',
            'linewidth',3,'color','r'); 
        data_point=point(1,1:2);
        setappdata(gcf,'data_point',data_point);
    case 'stop'
        set(gcbf,'WindowButtonMotionFcn','');
        set(gcbf,'WindowButtonUpFcn','');
end

function CalculateOptimalParameters(~, ~)
handles = guidata(gcf);
xy_point = getappdata(gcf,'xy_point');
R = str2double(get(handles.edit1,'string'));
if isempty(R) || R == 0 || isnan(R)
    R = 25;
end
H = str2double(get(handles.edit2,'string'));
if isempty(H) || H == 0 || isnan(H)
    H = 53;
end
t = str2double(get(handles.edit3,'string'));
if isempty(t) || t == 0 || isnan(t)
    t = 3;
end
w = str2double(get(handles.edit4,'string'));
if isempty(w) || w == 0 || isnan(w)
    w = 2.5;
end
if ~isempty(xy_point)
    EdgeFun = @(x)MyCsapi(xy_point(:,1),xy_point(:,2),abs(x));
else
    FunStr = get(handles.edit7,'string');
    if isempty(FunStr)
        FunStr = 'sqrt(R^2-x.^2)';
    end
    FunStr = vectorize(FunStr);
    EdgeFun = eval(['@(x)',FunStr]);
end
H = H-t;
Problem = MakeObjFun(R,H,w/2,t,EdgeFun);
Lb = fmincon(Problem);
set(handles.edit5,'string',num2str(Lb(1)));
set(handles.edit6,'string',num2str(Lb(2)));

function DynamicGraphDisplay(~, ~)
handles = guidata(gcf);
R = str2double(get(handles.edit1,'string'));
if isempty(R) || R == 0 || isnan(R)
    R = 25;
end
H = str2double(get(handles.edit2,'string'));
if isempty(H) || H == 0 || isnan(H)
    H = 53;
end
t = str2double(get(handles.edit3,'string'));
if isempty(t) || t == 0 || isnan(t)
    t = 3;
end
w = str2double(get(handles.edit4,'string'));
if isempty(w) || w == 0 || isnan(w)
    w = 2.5;
end

L = str2double(get(handles.edit5,'string'));
b = str2double(get(handles.edit6,'string'));
if isempty(L) || isempty(b) || isnan(L) || isnan(b)
    warndlg('您需要指定参数L和b的值，或者通过参数优化求L和b的值！！！','警告','replace');
    return;
end

xy_point = getappdata(gcf,'xy_point');
if ~isempty(xy_point)
    EdgeFun = @(x)MyCsapi(xy_point(:,1),xy_point(:,2),abs(x));
else
    FunStr = get(handles.edit7,'string');
    if isempty(FunStr)
        FunStr = 'sqrt(R^2-x.^2)';
    end
    FunStr = vectorize(FunStr);
    EdgeFun = eval(['@(x)',FunStr]);
end
axes(handles.axes1);
cla;
opengl software
H = H-t;
Result = PlotGateLegTable(R,L,H,w/2,t,EdgeFun,b);
set(handles.datatable,'data',Result);

function CopyFig(~, ~)
% 复制图片子函数
handles = guidata(gcf);
new_fig = figure('units','normalized',...
    'position',[0.29429,0.31901,0.409956,0.546875],...
    'Renderer','painters',...
    'visible','off');
new_axes = copyobj(handles.axes1,new_fig);
set(new_axes,'position',[0.13 0.11 0.775 0.815]);
hgexport(new_fig, '-clipboard');
delete(new_fig);

function SaveFig(~, ~)
% 保存图片子函数
handles = guidata(gcf);
geshi = {'*.emf','Enhanced metafile (*.emf)';
    '*.jpg','JPEG image (*.jpg)';
    '*.*','All Files (*.*)'};
[FileName, FilePath] = uiputfile(geshi,'保存图像文件','untitled.emf');
new_fig = figure('visible','off');
new_axes = copyobj(handles.axes1,new_fig);
set(new_axes,'position',[0.1 0.1 0.8 0.8]);
if ~isequal([FileName,FilePath],[0,0])
    FileFullName=[FilePath FileName];
    saveas(gca,FileFullName);
end
delete(new_fig);

function PrintFig(~, ~)
% 打印图片子函数
handles = guidata(gcf);
new_fig = figure('visible','on');
new_axes = copyobj(handles.axes1,new_fig);
%set(new_axes,'position',[0.1 0.15 0.85 0.8]);
set(new_axes,'position',[0.13 0.11 0.775 0.815]);

function yi = MyCsapi(x,y,xi)
yi = csapi(x,y,xi);
[xmin,id] = min(x);
yi(xi<xmin) = y(id);
[xmax,id] = max(x);
yi(xi>xmax) = y(id);