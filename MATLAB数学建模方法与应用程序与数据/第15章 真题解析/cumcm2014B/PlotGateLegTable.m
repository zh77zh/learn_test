function Result = PlotGateLegTable(R,L,H,w,t,EdgeFun,b)
% 绘制折叠桌折叠的动画效果图，并计算各木条开槽长度
% 坐标原点在桌子的下表面
% R：桌子半径
% L：平板半长
% H：桌子展开后高度（桌子的厚度不计算在内）
% w：木条半宽
% t：桌子的厚度
% EdgeFun：桌子边缘函数
% b：钢筋位置参数（最外侧木条外端点到钢筋的距离）
% Result：开槽起点，开槽终点，开槽长度
xi = -(R-w):2*w:R-w;    % 木条中心线x坐标
yi = real(EdgeFun(xi)); % 桌面边缘y坐标
a = yi(1);              % 桌面上最外侧木条的一半长
if nargin < 7
    b = (L-a)/2;        % 默认的b参数
end
% 若没有输出，则绘制折叠桌折叠过程动态效果图
if ~nargout
    figure('name','折叠桌折叠过程动态效果图','numbertitle','off');
end

n = numel(xi);          % 木条个数
fac = [1 2 3 4;2 6 7 3;4 3 7 8;1 5 8 4;1 2 6 5;5 6 7 8];% 构成木条各个面的顶点编号
verti1 = cell(1,n);  % 桌面上各木条的顶点坐标
verti2 = verti1;     % 左侧桌腿的顶点坐标
verti3 = verti1;     % 右侧桌腿的顶点坐标
h1 = zeros(1,n);
h2 = h1;
for i = 1:n
    vx = xi(i) + w*[-1 1 1 -1 -1 1 1 -1]'; % 顶点x坐标
    vz = t*[-1 -1 -1 -1 0 0 0 0]';         % 顶点z坐标
    vy1 = yi(i)*[1 1 -1 -1 1 1 -1 -1]';    % 顶点y坐标
    vy2 = [L,L,yi(i),yi(i),L,L,yi(i),yi(i)]'; % 顶点y坐标
    verti1{i} = [vx,vy1,vz];  % 桌面上各木条的顶点坐标
    patch('faces',fac,'vertices',verti1{i},'FaceColor',[1,1,0]); % 画桌面木条
    verti2{i} = [vx,vy2,vz];  % 左侧桌腿的顶点坐标
    % 画左侧桌腿
    h1(i) = patch('faces',fac,'vertices',verti2{i},'FaceColor',[0,0,1]);
    verti3{i} = [vx,-vy2,vz]; % 右侧桌腿的顶点坐标
    % 画右侧桌腿
    h2(i) = patch('faces',fac,'vertices',verti3{i},'FaceColor',[0,0,1]);
end
view([-60,24]);  % 设置视点位置
xlabel('x');ylabel('y');zlabel('z');  % 坐标轴标签
axis equal; axis([-R,R,-L,L,-5,H+5]); % 设置坐标轴显示属性
set(gca,'ZDir','reverse','color','none'); % 设置z轴正向向下
hold on;   % 图形保持
% 画桌脚边缘线
hline = plot3(xi,L*ones(size(xi)),zeros(size(xi)),'r','linewidth',3);

hei = linspace(0,H,24);  % 定义桌高向量
hei(hei>(L-a)) = [];     % 去除不合理桌高
m = numel(hei);          % 桌高个数
xTop = zeros(m,n);
yTop = xTop;
zTop = xTop;
D = xTop;
if ~nargout, f = getframe(gca); Frame = repmat(f.cdata,[1,1,1,m]);end
% 通过循环更新木条位置
for j = 1:m      % 对桌高进行循环
    for k = 1:n  % 对木条进行循环
        [CosBetai,di] = SubFcn(yi(k),hei(j));  % 计算夹角余弦及di
        SinBetai = real(sqrt(1-CosBetai^2));   % 夹角正弦 
        Oi = repmat([xi(k),yi(k),0],[8,1]);    % 左侧铰链结合点坐标
        Rot = [1,0,0;0,CosBetai,SinBetai;0,-SinBetai,CosBetai]; % 旋转矩阵
        vert_Rot = (verti2{k}-Oi)*Rot+Oi;  % 求左侧桌腿旋转后顶点坐标
        set(h1(k),'vertices',vert_Rot);    % 更新左侧桌腿的顶点坐标
        xTop(j,k) = xi(k);  % 桌脚边缘线x坐标
        yTop(j,k) = yi(k)+(L-yi(k))*CosBetai; % 桌脚边缘线y坐标
        zTop(j,k) = (L-yi(k))*SinBetai;  % 桌脚边缘线z坐标
        D(j,k) = di;
        % 更新桌脚边缘线的三维坐标
        set(hline,'XData',xTop(j,:),'YData',yTop(j,:),'ZData',zTop(j,:))
        Oi = repmat([xi(k),-yi(k),0],[8,1]); % 右侧铰链结合点坐标
        Rot = [1,0,0;0,CosBetai,-SinBetai;0,SinBetai,CosBetai]; % 旋转矩阵
        vert_Rot = (verti3{k}-Oi)*Rot+Oi;  % 求右侧桌腿旋转后顶点坐标
        set(h2(k),'vertices',vert_Rot);  % 更新右侧桌腿的顶点坐标
    end
    if j == 1
        % 绘制桌脚边缘线滑过的曲面
        hsurf = mesh(xTop([j;j],:),yTop([j;j],:),zTop([j;j],:),'FaceAlpha',0);
    else
        % 更新桌脚边缘线滑过的曲面
        set(hsurf,'XData',xTop(1:j,:),'YData',yTop(1:j,:),'ZData',zTop(1:j,:));
    end
    if ~nargout, f = getframe(gca); Frame(:,:,:,j) = f.cdata;end
    pause(0.3);
    drawnow;
end

% 平展状态下各木条开槽起点和终点的y坐标
yNotch = [(L-b)*ones(n,1),(D(end,:)+yi)'];
Result1 = [yNotch,abs(diff(yNotch,1,2))];  % 开槽起点，开槽终点，开槽长度
if nargout
    Result = Result1;  % 输出开槽起点，开槽终点，开槽长度结果
else
    % 绘制各木条开槽位置示意图
    figure('name','各木条开槽位置示意图','numbertitle','off');
    for i = 1:n
        vert = verti1{i}(5:8,:);
        % 绘制桌面木条投影
        patch('faces',1:4,'vertices',vert,'FaceColor',[0.5,0.5,0.5]);
        vert = verti2{i}(5:8,:);
        % 绘制左侧桌腿投影
        patch('faces',1:4,'vertices',vert,'FaceColor',[1,1,1]);
        vert = verti3{i}(5:8,:);
        % 绘制右侧桌腿投影
        patch('faces',1:4,'vertices',vert,'FaceColor',[1,1,1]);
    end
    view(2);
    hold on;
    plot(xi,yi,'r','linewidth',2);  % 绘制桌面边缘线
    plot(xi,-yi,'r','linewidth',2); % 绘制桌面边缘线
    plot(xi,[yNotch,-yNotch],'r','linewidth',2); % 绘制开槽位置线
    xlabel('x'); ylabel('y');  % 坐标轴标签
    axis equal;  axis([-R,R,-L,L]);  % 设置坐标轴显示属性
    
    % 绘制折叠桌动态变化过程的示意图（抓取9幅静态图片）
    figure('name','折叠桌动态变化过程的示意图');
    id = round(linspace(1,m,9));
    montage(Frame(:,:,:,id),'Size',[3,3]);
    % 将计算结果写入excel文件
    Result2 = [{'开槽起点','开槽终点','开槽长度'};num2cell(Result1)];
    xlswrite('各木条开槽起点和终点的y坐标.xlsx',Result2);
end

%-------------------------------------------------------------
%  子函数
%-------------------------------------------------------------
function [CosBeta,d] = SubFcn(y,z)
    % 求各木条旋转角度的余弦值;
    CosTheta = sqrt((L-a)^2-z^2)/(L-a);
    d = sqrt((y-a)^2+(L-a-b)^2-2*(y-a)*(L-a-b)*CosTheta);
    if y == a
        CosBeta = CosTheta;
    else
        CosBeta = -((y-a)^2+d^2-(L-a-b)^2)/(2*(y-a)*d);
    end
end
end