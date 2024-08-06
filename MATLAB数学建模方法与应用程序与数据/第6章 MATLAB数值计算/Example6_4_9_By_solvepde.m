%% 调用parabolic函数求解偏微分方程组（examp6.4-9）
% 1. 创建包含两个方程的微分方程模型
model6 = createpde(2);

% 2. 创建正方形求解区域
% 其中3为矩形编号，4是边数，矩形区域顶点x坐标为[0,1,1,0],y坐标为[0,0,1,1]
R1 = [3,4,0,1,1,0,0,0,1,1]'; 
geom = R1;
ns = char('R1')';             % 求解区域名称
sf = 'R1';                    % 构造求解区域的公式
g = decsg(geom,sf,ns);        % 创建几何体
geometryFromEdges(model6,g);  % 创建自定义求解区域
figure;
pdegplot(model6,'EdgeLabels','on');   % 绘制求解区域，并显示边界标签 
axis equal;
axis([-0.1,1.1,-0.1,1.1]);

% 3. 设置边值条件（Dirichlet 边值条件）
% u1的右边界对应的边值条件
applyBoundaryCondition(model6,'Edge',2,'u',1,'EquationIndex',1);
% u2的左边界对应的边值条件
applyBoundaryCondition(model6,'Edge',4,'u',0,'EquationIndex',2);

% 4. 划分网格
Me = generateMesh(model6,'Hmax',0.1,'GeometricOrder','linear'); 

% 5. 设置方程参数
m = 0;
d = 1;
c = [1/80;0;0;1/91;0;0];
a = 0;
specifyCoefficients(model6,'m',m,'d',d,'c',c,'a',a,'f',@ffun);

% 6. 设置初始条件
u0 = [1;0];
setInitialConditions(model6,u0);

% 7. 方程求解
tlist = linspace(0,2,31);  % 定义时间向量
u6 = solvepde(model6,tlist);  % 方程求解

% 8. 结果可视化
xi = Me.Nodes(1,:);  % 网格点x坐标
[xi,id] = sort(xi);
U = u6.NodalSolution;
U1 = squeeze(U(id,1,:))';
U2 = squeeze(U(id,2,:))';
figure;
surf(xi,tlist,U1);
shading interp;
title('u1(x,t)'); xlabel('Distance x'); ylabel('Time t');

figure;
surf(xi,tlist,U2);
shading interp;
title('u2(x,t)'); xlabel('Distance x'); ylabel('Time t');