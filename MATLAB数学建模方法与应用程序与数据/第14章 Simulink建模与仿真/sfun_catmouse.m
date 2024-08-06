function [sys,x0,str,ts] = sfun_catmouse(t,x,u,flag,a,b,c)
% 猫追老鼠的Simulink仿真
% 绘制演示动画的S-Function
% 2019.5.4修改，xiezhh（谢中华）

switch flag
    case 0
        [sys,x0,str,ts] = mdlInitializeSizes(a,b,c);
    case 3
        sys = mdlOutputs(t,x,u,c);
    case 9
        sys = mdlTerminate;
    case { 1, 2, 4 }
        sys = [];
    otherwise
        error('Unhandled flag');
end

function [sys,x0,str,ts] = mdlInitializeSizes(a,b,c)
% 初始化
sizes = simsizes;
sizes.NumContStates  = 0;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 0;
sizes.NumInputs      = -1;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;
sys = simsizes(sizes);
str = [];
x0 = [];
ts  = [-1 0];
if a < b
    return
end
maxy = 0.85*b*c/(a-b);
maxy = 15;
% 图形初始化
OldHandle = findobj('Type','figure','Tag','catmouse') ;
if ishandle( OldHandle )
    close( OldHandle );
end
fig = figure('units','normalized',...
    'position',[0.25,0.2,0.5,0.6],...
    'name','猫追老鼠的Simulink仿真',...
    'numbertitle','off',...
    'color',[0.8 0.8 0.8],...
    'tag','catmouse');
ax = axes('parent',fig,...
    'position',[0.03 0.1 0.8 0.8]);
hpoint1 = line(0,0,'Color',[0 0 1],...
    'Marker','.',...
    'MarkerSize',40,...
    'parent',ax);
hpoint2 = line(c,0,'MarkerFaceColor',[0 1 0],...
    'Marker','h',...
    'MarkerSize',15,...
    'parent',ax);
hline = line(0,0,'Color',[1 0 0],...
    'linewidth',2,...
    'parent',ax);
line([c,c],[0,maxy],'LineWidth',2);
hcat = text(-0.8,0,'猫','FontSize',12);
hmouse = text(c+0.3,0,'鼠','FontSize',12);
uicontrol(fig,'style','text',...
    'units','normalized',...
    'position',[0.81 0.2 0.15 0.05],...
    'string','时间:','fontsize',13,...
    'fontweight','bold',...
    'backgroundcolor',[0.8 0.8 0.8],...
    'HorizontalAlignment','left');
hedit = uicontrol(fig,'style','edit',...
    'units','normalized',...
    'position',[0.81 0.1 0.15 0.08],...
    'fontsize',13,'string','0 秒',...
    'backgroundcolor',[1 1 1]);
axis equal 
axis([0,c+1,0,maxy])
title('猫追老鼠的动画演示','FontSize',15,...
    'FontWeight','Bold')
xlabel('X');
ylabel('Y');
setappdata(fig,'handles',...
    [hpoint1,hpoint2,hline,hcat,hmouse,hedit]);
setappdata(fig,'xdata',0);
setappdata(fig,'ydata',0);
drawnow;


function sys = mdlOutputs(t,~,u,c)
sys = [];
fig = findobj('Type','figure',...
    'Tag','catmouse');
if ~isempty(fig)
    h = getappdata(fig,'handles');
    hpoint1 = h(1);
    hpoint2 = h(2);
    hline = h(3);
    hc = h(4);
    hm = h(5);
    hedit = h(6);
    xdata = [getappdata(fig,'xdata');u(1)];
    ydata = [getappdata(fig,'ydata');u(3)];
    set(hpoint1,'xdata',u(1),'ydata',u(3));
    set(hpoint2,'xdata',c,'ydata',u(2));
    set(hline,'xdata',xdata,'ydata',ydata);
    set(hc,'position',[u(1)-0.8,u(3),0]);
    set(hm,'position',[c+0.3,u(2),0]);
    set(hedit,'string',[num2str(t),' 秒']);
    setappdata(fig,'xdata',xdata);
    setappdata(fig,'ydata',ydata);
    drawnow;
%     f = getframe(gcf);  
%     imind = frame2im(f);
%     [imind,cm] = rgb2ind(imind,256);
%     if u == 0    
%        imwrite(imind,cm,'猫追老鼠的动画演示.gif','gif', 'Loopcount',inf,'DelayTime',0.01);
%     else
%        imwrite(imind,cm,'猫追老鼠的动画演示.gif','gif','WriteMode','append','DelayTime',0.01);
%     end
else
    return;
end

function sys = mdlTerminate
sys = [];