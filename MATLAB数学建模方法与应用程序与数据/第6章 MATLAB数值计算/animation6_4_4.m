%% examp6.4-4
% 用有限元法求解波动方程,生成gif格式动画文件
model = createpde(1);  % 创建包含一个方程的微分方程模型
geometryFromEdges(model,@squareg);  % 创建正方形求解区域

applyBoundaryCondition(model,'Edge',[2,4],'g',0);  % 左右边值条件
applyBoundaryCondition(model,'Edge',[1,3],'g',0);  % 上下边值条件
Me = generateMesh(model,'Hmax',0.08,'GeometricOrder','linear');  % 划分网格

% 方程参数
specifyCoefficients(model,'m',1,'d',0,'c',0.01,'a',0,'f',@framp)

% 初值条件
u0 = 0;
ut0 = 0;
setInitialConditions(model,u0,ut0)

tlist = linspace(0,20,61);  % 定义时间向量
u1 = solvepde(model,tlist);  % 方程求解

% 结果可视化
[X,Y] = meshgrid(linspace(-1,1,60)); % 对x,y轴作矩形网格划分
Uxy = interpolateSolution(u1,X,Y,1:numel(tlist)); % 把求解结果插值到矩形网格点
U = reshape(Uxy(:,1),size(X));
figure;
h = surf(X,Y,U);  % 绘制第一个时刻的插值曲面
axis([-1.1,1.1,-1.1,1.1,-0.8,0.8]);
view(-30,70);  % 设置视点位置
colormap(jet);  % 设置颜色矩阵
shading interp; % 插值染色
%light('pos',[0.6,-0.6,20]);
camlight;  % 加入光源
lighting gouraud;  % 设置光照模式
xlabel('x');ylabel('y'),zlabel('u');

filename = '水波扩散动画.gif';
Frame = getframe(gcf);
IM = frame2im(Frame);
[IM,map] = rgb2ind(IM,256);
imwrite(IM,map,filename,'gif', 'Loopcount',inf,'DelayTime',0.1);

% 水波扩散的动态展示
for i = 1:numel(tlist)
    U = reshape(Uxy(:,i),size(X));
    set(h,'ZData',U);  % 更新坐标
    Frame = getframe(gcf);
    IM = frame2im(Frame);
    [IM,map] = rgb2ind(IM,256);
    imwrite(IM,map,filename,'gif','WriteMode','append','DelayTime',0.1);
    pause(0.1);
end