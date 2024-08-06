%--------------------------------------------------------------------------
%  第6章  MATLAB数值计算
%--------------------------------------------------------------------------
% CopyRight：xiezhh


%% examp6.1-1
h = 0.01;
x = 0:h:2*pi;
y = sin(x);

dy_dx1 = diff(y)./diff(x);
dy_dx2 = gradient(y,h);

figure;
plot(x,y);
hold on
plot(x(1:end-1),dy_dx1,'k:');
plot(x,dy_dx2,'r--');
legend('y = sin(x)','导函数曲线（diff）','导函数曲线（gradient）');
xlabel('x'); ylabel('正弦曲线及导函数曲线')

%% examp6.1-2
t1 = linspace(0,2*pi,60);
x1 = cos(t1); y1 = sin(t1);
s1 = abs(trapz(x1,y1))

t2 = linspace(0,2*pi,200);
x2 = cos(t2); y2 = sin(t2);
s2 = abs(trapz(x2,y2))

t3 = linspace(0,2*pi,2000);
x3 = cos(t3); y3 = sin(t3);
s3 = abs(trapz(x3,y3))

%% examp6.1-3
fun1 = @(x)exp(-x.^2);
s1 = integral(fun1,0,1)

fun2 = @(x,y)x.*sqrt(10-y.^2);
yfun1 = @(x)x-2;
yfun2 = @(x)2-sin(x);
s2 = integral2(fun2,-1,2,yfun1,yfun2)

fun3 = @(x,y,z)x.*y.*z;
yfun1 = @(x)x;
yfun2 = @(x)2*x;
zfun1 = @(x,y)x.*y;
zfun2 = @(x,y)2*x.*y;
s3 = integral3(fun3,1,2,yfun1,yfun2,zfun1,zfun2)

%% examp6.1-4
fxy = @(x,y)exp(-x.^2)./(x.^2+y.^2);
fy1 = @(y)integral(@(x)fxy(x,y),-1,1)^2;
fy2 = @(y)arrayfun(@(t)fy1(t),y);
fun1 = @(y)2*y.*exp(-y.^2).*fy2(y);
s1 = integral(fun1,0.2,1)

fun2 = @(x1,x2,x3,x4)exp(x1.*x2.*x3.*x4);
f_x1 = @(x1)integral3(@(x2,x3,x4)fun2(x1,x2,x3,x4),0,1,0,1,0,1);
f_x1 = @(x1)arrayfun(@(t)f_x1(t),x1);
s2 = integral(f_x1,0,1)

%% examp6.2-1
A = [1 -1 2 -3;-2 2 1 1;-1 1 8 -8];
b = [3;-1;6];
x = A\b

%% examp6.2-2
p = [2, -3, 5, -10];
x = roots(p)

%% examp6.2-3
fun = @(x)-x.*sin(5*exp(1-x.^2));
figure;
fplot(fun,[-1 1]);
grid on;
[x,fval] = fzero(fun,[0.2,0.4])
hold on;
plot(x,fval,'ko');
xlabel('x');
ylabel('$$ y = -xsin(5e^{1-x^2}) $$','Interpreter','latex');

%% examp6.2-4
fun = @(X)[X(1) - X(2) - exp(-X(1)); -X(1) + 2*X(2) - exp(-X(2))];
x0 = [1,1];
options = optimset('Display','iter'); %显示迭代过程
[x,fval] = fsolve(fun,x0,options)

%% examp6.2-5
xyt = [500    3300    21    9
       300     200    19   29
       800    1600    14   51
      1400    2200    13   17
      1700     700    11   46
      2300    2800    14   47
      2500    1900    10   14
      2900     900    11   46
      3200    3100    17   57
      3400     100    16   49];
x = xyt(:,1);
y = xyt(:,2);
Minutes = xyt(:,3);
Seconds = xyt(:,4);
T = Minutes + Seconds/60; 
modelfun = @(b) sqrt((x-b(1)).^2+(y-b(2)).^2+b(3).^2)/(60*b(4))+b(5)-T;
b0 = [1000 100 10 1 1];
options = optimoptions('fsolve','Display','none',...
    'Algorithm','Levenberg-Marquardt');
[Bval,Fval] = fsolve(modelfun,b0,options)

%% examp6.3-1
a = 14;b = 8;c = 10;
f = @(t,x)sqrt((c-x(1))^2+(b*t-x(2))^2); 
fun = @(t,x)[a*(c-x(1))/f(t,x);a*(b*t-x(2))/f(t,x)];
tspan = linspace(0,1.06,100);
x0 = [0;0];
[t,x] = ode45(fun,tspan,x0);

figure;
hpoint1 = line(0,0,'Color',[0 0 1],'Marker',...
    '.','MarkerSize',40);
hpoint2 = line(c,0,'MarkerFaceColor',[0 1 0],...
    'Marker','p','MarkerSize',15);
hline = line(0,0,'Color',[1 0 0],'linewidth',2);
line([c c],[0 c],'LineWidth',2);
hcat = text(-0.8,0,'猫','FontSize',12);
hmouse = text(c+0.3,0,'鼠','FontSize',12);
xlabel('X'); ylabel('Y');
axis([0 c+1 0 9.5]);

for i = 1:size(x,1)
    ymouse = t(i)*b;
    set(hpoint1,'xdata',x(i,1),'ydata',x(i,2));
    set(hpoint2,'xdata',c,'ydata',ymouse);
    set(hline,'xdata',x(1:i,1),'ydata',x(1:i,2));
    set(hcat,'Position',[x(i,1)-0.8,x(i,2),0]);
    set(hmouse,'Position',[c+0.3,ymouse,0]);
    pause(0.1);
end

%% examp6.3-2
fun = @(t,y,mu)[y(2);mu*(1-y(1)^2)*y(2)-y(1)];
tspan = [0,30];%时间区间
y0 = [1 0];
ColorOrder = {'r','b','k'};
LineStyle = { '-','--',':'};
figure; ha1 = axes; hold on;
figure; ha2 = axes; hold on;
for mu = 1:3
    [t,y] = ode45(fun,tspan,y0,[],mu);
    plot(ha1,t,y(:,1),'color',ColorOrder{mu},'LineStyle',LineStyle{mu});
    plot(ha2,y(:,1),y(:,2),'color',ColorOrder{mu},'LineStyle',LineStyle{mu});
end
xlabel(ha1,'t'); ylabel(ha1,'x(t)');
legend(ha1,'\mu = 1','\mu = 2','\mu = 3');
hold off
xlabel(ha2,'位移'); ylabel(ha2,'速度');
legend(ha2,'\mu = 1','\mu = 2','\mu = 3');
hold off

%% examp6.3-3
fun = @(t,y,dy)[dy(1)-y(2);
                dy(2)*sin(y(4))+dy(4)^2+2*y(1)*y(3)-y(1)*dy(2)*y(4);
                dy(3)-y(4);
                y(1)*dy(2)*dy(4)+cos(dy(4))-3*y(2)*y(3)];

t0 = 0;         % 自变量的初值
y0 = [1;0;0;1]; % 状态变量初值向量y0
% fix_y0用来指定初值向量y0的元素是否可以改变。1表示对应元素不能改变，0为可以改变
fix_y0 = [1;1;1;1]; % 本例中y0的值都给出了，因此都不能改变，所有fix_y0全为1
dy0 = [0;3;1;0];    % 猜测一下一阶导数dy的初值dy0;
% 由于本例中一阶导数dy的初值dy0是猜测的，都可以改变，因此fix_dy0 全部为0
fix_dy0 = [0;0;0;0];
% 调用decic函数来决定y和dy的初值
[y02,dy02] = decic(fun,t0,y0,fix_y0,dy0,fix_dy0);

%求解微分方程
[t,y] = ode15i(fun,[0,5],y02,dy02); % y02和dy02由decic输出
% 结果图示
figure;
plot(t,y(:,1),'k-','linewidth',2);
hold on
plot(t,y(:,2),'k--','linewidth',2);
plot(t,y(:,3),'k-.','linewidth',2);
plot(t,y(:,4),'k:','linewidth',2);
% 图例,位置自动选择最佳位置
L = legend('y_1(t)','y_2(t)','y_3(t)','y_4(t)','Location','best');
set(L,'fontname','Times New Roman');
xlabel('t');ylabel('y(t)');

%% examp6.3-4
lags = [1,3];       % 延迟常数向量
history = [0,0,1];  % 小于初值时的历史函数
tspan = [0,8];      % 时间区间
% 方法一：调用dde23函数求解
sol = dde23(@ddefun,lags,history,tspan); 
% % 方法二：调用ddesd函数求解
% sol = ddesd(@ddefun,lags,history,tspan); 

% 画图呈现结果
figure;
plot(sol.x,sol.y(1,:),'k-','linewidth',2);
hold on
plot(sol.x,sol.y(2,:),'k-.','linewidth',2);
plot(sol.x,sol.y(3,:),'k-*','linewidth',1);
hold off
% 图例,位置自动选择最佳位置
L = legend('y_1(t)','y_2(t)','y_3(t)','Location','best');
set(L,'fontname','Times New Roman');   % 设置图例字体
xlabel('t');ylabel('y(t)');            % 添加坐标轴标签

%% examp6.3-6
% 微分方程组所对应的匿名函数
BvpOdeFun  = @(t,y)[y(2)
                    2*y(2)*cos(t)-y(1)*sin(4*t)-cos(3*t)];
% 边界条件所对应的匿名函数。
% 边界条件为 y1(0) = 1, y1(4) = 2，这里0,4分别对应y的下边界和上边界。
% 这里ylow(1)表示y1(0)，yup(1)表示y1(4)，类似的y2(0)和y2(4)分别用ylow(2)和yup(2)表示
BvpBcFun = @(ylow,yup)[ylow(1)-1; yup(1)-2];

T = linspace(0,4,10); % 为调用bvpinit生成初始化网格作准备
% 对状态变量y作出初始假设，由于y1(0) = 1,y1(4) = 2，可选取一个满足上述条件的函数
% y1(t) = 1+t/4来作为对y1(t)的初始假设，从而其导数1/4作为对y2(t)的初始假设
BvpYinit = @(t)[ 1+t/4; 1/4 ];
solinit = bvpinit(T,BvpYinit); % 调用bvpinit函数生成初始解

sol = bvp4c(BvpOdeFun,BvpBcFun,solinit); % 调用bvp4c求解,也可以换成bvp5c
tint = linspace(0,4,100);
Stint = deval(sol,tint); % 根据得到的sol利用deval函数求出[0,4]区间内更多其他的解

% 画图呈现结果
figure;
plot(tint,Stint(1,:),'k-','linewidth',2);
hold on
plot(tint,Stint(2,:),'k:','linewidth',2);
% 图例,位置自动选择最佳位置
L = legend('y_1(t)','y_2(t)','Location','best');
set(L,'fontname','Times New Roman');   % 设置图例字体
xlabel('t');ylabel('y(t)');            % 添加坐标轴标签

%% examp6.4-1
u1 = ones(1,49);
%  根据差分方程构造目标函数（方程组）
objfun = @(u)([u(2:end,:);u1]+[u1;u(1:end-1,:)]+...
       [u(:,2:end),u1']+[0*u1',u(:,1:end-1)])/4-u;
U0 = rand(49);
[Uin,Error] = fsolve(objfun,U0);  % 求解内点温度
U = zeros(size(U0)+2);
U(:,end) = 1;
U(1,:) = 1;
U(end,:) = 1;
U(2:end-1,2:end-1) = Uin;
[X,Y] = meshgrid(linspace(0,1,51));
figure;
surf(X,Y,U);
xlabel('X'); ylabel('Y'); zlabel('U(X,Y)');

%% examp6.4-2
U = zeros(100);  % 初值矩阵
t = (1:100)/100; x = t;  % t和x的划分向量
U(1,:) = sin(t);  % 下边界条件
U(end,:) = cos(t);  % 上边界条件
U(:,1) = x;  % 初值条件
b2 = 0.001; dx = 0.01;dt = 0.01;r = b2*dt/dx^2;  % 参数
% 差分方程求解
for j = 1:99
       U(2:99,j+1) = (1-2*r)*U(2:99,j)+r*(U(3:100,j)+U(1:98,j));
end
[T,X] = meshgrid(t);  % 网格矩阵
figure;
surf(T,X,U);  % 绘制面图
xlabel('T');  ylabel('X');  zlabel('U(T,X)');  % 坐标轴标签

%% examp6.4-3
u = zeros(301);                 % 定义零矩阵
dt = 1/300; c = 0.03;           % 参数
x = linspace(0,1,301); t = x;  % t和x的划分向量
u(:,1) = x.*(1-x)/10;           % 初值
u(1,:) = sin(t);                % 边值
v = sin(2*pi*x);                % 初速度
% 计算u(i,2)
u(2:300,2) = (1-c)*u(2:300,1) + ...
       1/2*c*(u(3:301,1)+u(1:299,1)) + v(2:300)'*dt;
figure;
h = plot(x,u(:,1));             % 绘制初始弦曲线
axis([-0.1,1.1,-1,1]);          % 设置坐标轴范围
xlabel('x');ylabel('U(x)');     % 坐标轴标签
% 用有限差分法求解方程，并动态展示求解结果
for j = 3:301
       u(2:300,j) = 2*(1-c)*u(2:300,j-1)+c*(u(3:301,j-1)+...
           u(1:299,j-1))-u(2:300,j-2);
       set(h,'YData',u(:,j));      % 更新弦上各点位移
       pause(0.1);                 % 暂停0.1秒
end
%text(0.8,0.8,['t = ',num2str((j-1)/300,'%3.2f')])

%% examp6.4-4
% 用有限元法求解波动方程
model = createpde(1);  % 创建包含一个方程的微分方程模型
geometryFromEdges(model,@squareg);  % 创建正方形求解区域
figure;
pdegplot(model,'EdgeLabels','on');  % 绘制求解区域
axis equal;
axis([-1.1,1.1,-1.1,1.1]);

applyBoundaryCondition(model,'Edge',[2,4],'g',0);  % 左右边值条件
applyBoundaryCondition(model,'Edge',[1,3],'g',0);  % 上下边值条件
Me = generateMesh(model,'Hmax',0.08,'GeometricOrder','linear');  % 划分网格
figure;
pdeplot(model);     % 显示网格图
axis equal;
axis([-1.1,1.1,-1.1,1.1]);

% 方程参数
specifyCoefficients(model,'m',1,'d',0,'c',0.01,'a',0,'f',@framp0);

% 初值条件
u0 = 0;
ut0 = 0;
setInitialConditions(model,u0,ut0);

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
% 水波扩散的动态展示
for i = 1:numel(tlist)
    U = reshape(Uxy(:,i),size(X));
    set(h,'ZData',U);  % 更新坐标
    pause(0.1);
end

%% examp6.4-5
model2 = createpde(1);
geometryFromEdges(model2,@squareg);   % 创建正方形求解区域
figure;
pdegplot(model2,'EdgeLabels','on');   % 绘制求解区域，并显示边界标签
axis equal;
axis([-1.1,1.1,-1.1,1.1]);

% 设置第左右边界的边值条件（Dirichlet 边值条件）
applyBoundaryCondition(model2,'Edge',[2,4],'u',0);
% 设置上下边界的边值条件（Neumann 边值条件）
applyBoundaryCondition(model2,'Edge',[1,3],'g',0);
Me = generateMesh(model2,'Hmax',0.1,'GeometricOrder','linear');  % 划分网格

% 方程参数
specifyCoefficients(model2,'m',1,'d',0,'c',1,'a',0,'f',0);
% 初值条件
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
for i = 1:numel(tlist)
    set(h(1),'Vertices',[XY(:,1),XY(:,2),UXY(:,i)]);  % 更新坐标
    pause(0.1);
end

%% examp6.4-6
% 1. 创建包含一个方程的微分方程模型
model3 = createpde(1);

% 2. 创建圆形求解区域
geometryFromEdges(model3,@circleg);
figure;
pdegplot(model3,'EdgeLabels','on');   % 绘制求解区域，并显示边界标签 
axis equal;
axis([-1.1,1.1,-1.1,1.1]);

% 3. 设置边值条件（Dirichlet 边值条件）
NumEdges = model3.Geometry.NumEdges;  % 求解区域边界数
applyBoundaryCondition(model3,'Edge',1:NumEdges,'u',0);

% 4. 划分网格
generateMesh(model3,'Hmax',0.02,'GeometricOrder','linear'); 

% 5. 设置方程参数
specifyCoefficients(model3,'m',0,'d',1,'c',1,'a',0,'f',0);

% 6. 设置初值条件
u0fun = @(location)sqrt(location.x.^2 +location.y.^2) <= 0.4;
setInitialConditions(model3,u0fun);

% 7. 方程求解
tlist = linspace(0,0.1,21);  % 定义时间向量
u3 = solvepde(model3,tlist);  % 方程求解

% 8. 结果可视化
U = u3.NodalSolution;
figure;
umax = max(U(:));  % 最大温度
umin = min(U(:));  % 最小温度
% 热扩散的动态展示
for i = 1:numel(tlist)
    pdeplot(model3,'XYData',U(:,i)); % 绘制温度分布图
    caxis([umin, umax]);             % 设置坐标系颜色范围
    axis equal;
    axis([-1.1,1.1,-1.1,1.1]);
    xlabel('x');ylabel('y');
    pause(0.1);
end

%% examp6.4-7
% 1. 创建包含一个方程的微分方程模型
model4 = createpde;

% 2. 从外部文件导入求解区域的几何描述
importGeometry(model4,'Block.stl'); 
figure;
h = pdegplot(model4,'FaceLabels','on'); % 绘制求解区域
h(1).FaceAlpha = 0.5; % 设置透明度值为0.5

% 3. 设置边值条件
% x = 0, x = 100, z = 0, z = 50所对应的边值条件
applyBoundaryCondition(model4,'Face',1:4,'u',0);
% y = 0所对应的边值条件
applyBoundaryCondition(model4,'Face',6,'g',-1);
% y = 20所对应的边值条件
applyBoundaryCondition(model4,'Face',5,'g',1);

% 4. 划分网格
generateMesh(model4);

% 5. 设置方程参数
f = @(location,state)log(1 + location.x + location.y./(1+location.z));
specifyCoefficients(model4,'m',0,'d',0,'c',1,'a',0,'f',f);

% 6. 方程求解
u4 = solvepde(model4);

% 7. 结果可视化
p = u4.Mesh.Nodes;  % 三角网顶点坐标
% 对x,y,z坐标轴进行网格划分
xi = linspace(min(p(1,:)),max(p(1,:)),60);
yi = linspace(min(p(2,:)),max(p(2,:)),60);
zi = linspace(min(p(3,:)),max(p(3,:)),60);
[X,Y,Z] = meshgrid(xi,yi,zi);  % 生成矩形网格数据
% 把求解结果插值到矩形网格点
V = interpolateSolution(u4,X,Y,Z);
% 把插值结果转为三维数组
V = reshape(V,size(X));

figure;
% 绘制切片图
slice(X,Y,Z,V,xi,yi,zi);
shading interp;  % 插值染色
alpha(0.03);  % 设置透明度
xlabel('x');
ylabel('y');
zlabel('z');
colorbar;  % 添加颜色条
axis equal;

%% examp6.4-8
% 1. 创建包含一个方程的微分方程模型
model5 = createpde;

% 2. 创建圆形求解区域
geometryFromEdges(model5,@circleg);

% 3. 设置边值条件
boundaryfun = @(location,state)location.x.^2;  % 定义边界函数
NumEdges = model5.Geometry.NumEdges;
applyBoundaryCondition(model5,'Edge',1:NumEdges,...
    'u',boundaryfun,'Vectorized','on');

% 4. 划分网格
generateMesh(model5,'Hmax',0.1,'GeometricOrder','linear');

% 5. 设置方程参数
c = @(location,state)1./sqrt(1 + state.ux.^2 + state.uy.^2);
specifyCoefficients(model5,'m',0,'d',0,'c',c,'a',0,'f',0);

% 6. 方程求解
u5 = solvepde(model5);

% 7. 结果可视化
figure;
U = u5.NodalSolution;
pdeplot(model5,'xydata',U,'zdata',U);

%% examp6.4-9
m = 0;                                     % 方程中的m参数
x = linspace(0,1,30);                      % 定义x向量
t = linspace(0,2,30);                      % 定义t向量
sol = pdepe(m,@pdefun,@pdeic,@pdebc,x,t);  % 方程求解
u1 = sol(:,:,1); u2 = sol(:,:,2);          % 分别提取u1和u2的结果
% 结果可视化
figure;
surf(x,t,u1);                              % 绘制u1关于x和t的三维曲面
title('u1(x,t)'); xlabel('Distance x'); ylabel('Time t')
figure;
surf(x,t,u2);                              % 绘制u2关于x和t的三维曲面
title('u2(x,t)'); xlabel('Distance x'); ylabel('Time t')

%% examp6.5-1
syms a h1 h2 h3 h4 x y
r1 = sqrt(h1^2 + x^2 + y^2);
r2 = sqrt(h2^2 + (a-x)^2 + y^2);
r3 = sqrt(h3^2 + (a-x)^2 + (a-y)^2);
r4 = sqrt(h4^2 + x^2 + (a-y)^2);
Ixy = h1/r1^3 + h2/r2^3 + h3/r3^3 + h4/r4^3;
Ixy = subs(Ixy,{a,h1,h2,h3,h4},{20,7,9,8,10});
figure;
fsurf(Ixy,[0,20,0,20],'ShowContours','on')
xlabel('x');ylabel('y');zlabel('I(x,y)')
colorbar;
view(19,18);

dx = diff(Ixy,x);
dy = diff(Ixy,y);
dxfun = matlabFunction(dx);
dyfun = matlabFunction(dy);
Ixyfun = matlabFunction(Ixy);
eq = @(x)1000*[dxfun(x(1),x(2));dyfun(x(1),x(2))];
Xg = fsolve(eq,[10,10]);
X1 = fsolve(eq,[0,0]);
X2 = fsolve(eq,[20,0]);
X3 = fsolve(eq,[20,20]);
X4 = fsolve(eq,[0,20]);
X = [Xg;X1;X2;X3;X4;0,0;20,0;20,20;0,20];
Ival = Ixyfun(X(:,1),X(:,2));
id = (1:9)';
xi = X(:,1);
yi = X(:,2);
VarNames = {'序号','x坐标','y坐标','照度值I'};
Result = table(id,xi,yi,Ival,'VariableNames',VarNames)

figure;
fsurf(Ixy,[0,20,0,20],'ShowContours','on')
xlabel('x');ylabel('y');zlabel('I(x,y)')
colorbar;
view(19,18);
hold on
stem3(X(1:5,1),X(1:5,2),Ival(1:5),'ro','filled');

%% examp6.5-2
% 指数增长模型
syms x(t) lambda x0
x(t) = dsolve(diff(x) == lambda*x, x(0) == x0)

% ------------SI模型---------------
% 符号解（解析解）
syms I(t) lambda I0
I(t) = dsolve(diff(I) == lambda*I*(1-I),I(0) == I0);
I = simplify(I)

% 数值解
ts = 0:0.5:15;                             % 时间
I0 = 0.02;                      % I(t)的初值
lambda = 1;  % 每个患病者每天有效接触的平均人数，称为日接触率
SIFun = @(t,I)lambda * I .* (1-I);
[t,It] = ode45(SIFun,ts,I0);           % 调用ode45函数求解模型
% 结果可视化
figure
St = 1-It;                           % 易感者
plot(t,St,'--b',t,It,'r','LineWidth',2); % 绘制S(t)和I(t)曲线
grid on                                % 显示网格线
xlabel('时间/天');       % x轴标签
legend('易感者 S(t)', '患病者 I(t)','Location','best')  % 添加图例

% ------------SIS模型---------------
% 符号解（解析解）
syms I(t) lambda mu I0
% 当lambda ≠ mu时
I1(t) = dsolve(diff(I) == lambda*I*(1-I)-mu*I,I(0) == I0);
I1 = simplify(I1,'Steps',326)
% 当lambda = mu时
I2(t) = dsolve(diff(I) == -lambda*I^2,I(0) == I0)

% 数值解
ts = 0:0.5:15;                             % 定义时间向量
I0 = 0.02;                      % I(t)的初值
lambda = 1;  % 每个患病者每天有效接触的平均人数，称为日接触率
mu = 0.3;    % 日治愈率
SISFun = @(t,I)lambda*I.*(1-I)-mu*I;
[t,It] = ode45(SISFun,ts,I0);           % 调用ode45函数求解模型
% 结果可视化
figure
St = 1-It;                           % 易感者
plot(t,St,'--b',t,It,'r','LineWidth',2); % 绘制S(t)和I(t)曲线
grid on                                % 显示网格线
xlabel('时间/天');       % x轴标签
legend('易感者 S(t)', '患病者 I(t)','Location','best')  % 添加图例

% ------------SIR模型---------------
ts = 0:0.5:25;                                       % 定义时间向量
S0 = 0.98;                                           % S(t)的初值
I0 = 0.02;                                           % I(t)的初值
R0 = 0;                                              % R(t)的初值
lambda = 1;                                          % 日接触率
mu = 0.3;                                            % 日治愈率
SIRFun = @(t,SIR)[-lambda*SIR(1).*SIR(2);...
    lambda*SIR(1).*SIR(2)-mu*SIR(2);mu*SIR(2)];      % 定义模型对应的匿名函数
[t,SIR] = ode45(SIRFun,ts,[S0,I0,R0]);              % 调用ode45函数求解模型
% 结果可视化
figure                                               % 新建图窗
St = SIR(:,1);                                       % 易感者比例
It = SIR(:,2);                                       % 患病者比例
Rt = SIR(:,3);                                       % 移出者比例
plot(t,St,'--b',t,It,'r',t,Rt,':k','LineWidth',2); % 绘制S(t)、I(t)和R(t)曲线
grid on                                              % 显示网格线
xlabel('时间/天');                                   % 添加x轴标签
legend('易感者 S(t)', '患病者 I(t)',...
    '康复者 R(t)','Location','best')                 % 添加图例

%% examp6.5-3
m = 2;                                     % 方程中的m参数
r = linspace(0,9,41)/100;                  % 定义距离向量r
t = linspace(0,10,41)*3600;                % 定义时间向量t
T = pdepe(m,@pdefun2,@pdeic2,@pdebc2,r,t);   % 模型求解

% 结果可视化
figure;
surf(r*100,t/3600,T);                     % 绘制温度T关于距离r和时间t的三维曲面
colorbar;
zlim([min(T(:))-2,max(T(:))+2]);
xlabel('距离 r （cm）');
ylabel('时间 t （h）')
zlabel('温度 T （℃）'); 
view(116,33);
text(8,-0.5,34.5,'开始冷藏','Rotation',40)
text(0,5,26,'瓜心')
text(7,2,10.5,'瓜皮')

% 绘制西瓜等效切面温度扩散动画
theta = linspace(0,2*pi,61);
X = 100*r'*cos(theta);
Y = 100*r'*sin(theta);
Ti = T(1,:)'*ones(size(theta));
xmin = min(X(:))-1;
xmax = max(X(:))+1;
ymin = min(Y(:))-1;
ymax = max(Y(:))+1;
X_bac = linspace(xmin,xmax,50);
Y_bac = linspace(ymin,ymax,50);
[X_bac,Y_bac] = meshgrid(X_bac,Y_bac);
figure;
% 绘制背景温度面
Tinf = 6;     % 终极温度 [℃]，即冰箱冷藏室温度
surf(X_bac,Y_bac,Tinf*ones(size(X_bac)),...
    'FaceColor','interp','EdgeColor','k');
hold on;
% 绘制零时刻西瓜等效温度切面
hs = surf(X,Y,Ti,'FaceColor','interp','EdgeColor','k');
axis equal;
axis([xmin,xmax,ymin,ymax])
xlabel('X');
ylabel('Y');
caxis([Tinf,max(T(:))]);
colorbar;
view(2);
ht = title(['冷藏 ',num2str(t(1)/3600,'%2.1f'),' 小时']);

% 通过循环绘制各时刻西瓜等效温度切面
for i = 2:numel(t)
    Ti = T(i,:)'*ones(size(theta));
    set(hs,'ZData',Ti);
    set(ht,'String',['冷藏',num2str(t(i)/3600,'%2.1f'),'小时']);
    pause(0.2);
end