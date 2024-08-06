%% examp6.4-5
% 用有限元法求解膜振动方程,生成gif格式动画文件
model2 = createpde(1);
geometryFromEdges(model2,@squareg);   % 创建正方形求解区域

% 设置第左右边界的边值条件（Dirichlet 边值条件）
applyBoundaryCondition(model2,'Edge',[2,4],'u',0);
% 设置上下边界的边值条件（Neumann 边值条件）
applyBoundaryCondition(model2,'Edge',[1,3],'g',0);
Me = generateMesh(model2,'Hmax',0.1,'GeometricOrder','linear');  % 划分网格

% 方程参数
specifyCoefficients(model2,'m',1,'d',0,'c',1,'a',0,'f',0);
% 初始条件
u0fun = @(location)atan(cos(pi/2*location.x));
ut0fun = @(location)3*sin(pi*location.x).*exp(cos(pi*location.y));
setInitialConditions(model2,u0fun,ut0fun);

tlist = linspace(0,6,41);                      % 定义时间向量
u2 = solvepde(model2,tlist);                   % 方程求解

XY = u2.Mesh.Nodes';
Tri = u2.Mesh.Elements';
UXY = u2.NodalSolution;
figure;
h = trisurf(Tri,XY(:,1),XY(:,2),UXY(:,1));     % 绘制三角网面图
axis([-1.1,1.1,-1.1,1.1,-3,3]);                % 设置坐标轴范围
xlabel('x');ylabel('y');zlabel('u(x,y)');

filename = '膜振动动画.gif';
Frame = getframe(gcf);
IM = frame2im(Frame);
[IM,map] = rgb2ind(IM,256);
imwrite(IM,map,filename,'gif', 'Loopcount',inf,'DelayTime',0.1);

for i = 1:numel(tlist)
    set(h(1),'Vertices',[XY(:,1),XY(:,2),UXY(:,i)]);  % 更新坐标
    Frame = getframe(gcf);
    IM = frame2im(Frame);
    [IM,map] = rgb2ind(IM,256);
    imwrite(IM,map,filename,'gif','WriteMode','append','DelayTime',0.1);
    pause(0.1);
end