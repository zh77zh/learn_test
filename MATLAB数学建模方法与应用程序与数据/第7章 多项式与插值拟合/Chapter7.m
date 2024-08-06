%--------------------------------------------------------------------------
%  第7章  多项式与插值拟合
%--------------------------------------------------------------------------

%% examp7.1-1 多项式拟合
%--------------------散点图------------------
[Data,Textdata] = xlsread('食品零售价格.xls');
x = Data(:,1);
y = Data(:,3);
timestr = Textdata(3:end,2);
figure;
plot(x,y,'k.','Markersize',15);
set(gca,'XTick',1:2:numel(x),'XTickLabel',timestr(1:2:end));
set(gca,'XTickLabelRotation',-30);
xlabel('时间(x)');
ylabel('食品零售价格分类指数');
%-------------------4阶多项式拟合--------------------
[p4,S4] = polyfit(x,y,4)
r = poly2sym(p4);
r = vpa(r,5)
 
%--------------------更高阶多项式拟合---------------------
[p5,S5] = polyfit(x,y,5);
S5.normr
[p6,S6] = polyfit(x,y,6);
S6.normr
[p7,S7] = polyfit(x,y,7);
S7.normr
[p8,S8] = polyfit(x,y,8);
S8.normr
[p9,S9] = polyfit(x,y,9);
S9.normr

%-------------------拟合效果图----------------------
figure;
plot(x,y,'k.','Markersize',15);
set(gca,'XTick',1:2:numel(x),'XTickLabel',timestr(1:2:end));
set(gca,'XTickLabelRotation',-90);
xlabel('时间(x)');
ylabel('食品零售价格分类指数');
hold on;
yd4 = polyval(p4,x);
yd6 = polyval(p6,x);
yd8 = polyval(p8,x);
yd9 = polyval(p9,x);
plot(x,yd4,'r:+');
plot(x,yd6,'g--s');
plot(x,yd8,'b-.d');
plot(x,yd9,'m-p');
legend('原始散点','4次多项式拟合','6次多项式拟合','8次多项式拟合','9次多项式拟合')

%% examp7.4-1 一维插值
x0 = [0,3,5,7,9,11,12,13,14,15];
y01 = [0,1.8,2.2,2.7,3.0,3.1,2.9,2.5,2.0,1.6];
y02 = [0,1.2,1.7,2.0,2.1,2.0,1.8,1.2,1.0,1.6];
x = 0:0.1:15;
ysp1 = interp1(x0,y01,x,'spline');
ysp2 = interp1(x0,y02,x,'spline');
figure;
plot([x0,x0],[y01,y02],'o');
hold on;
plot(x,ysp1,'r',x,ysp2,'r');
xlabel('X')
ylabel('Y')
legend('插值节点','三次样条插值','location','northwest') 


%% examp7.4-2 一维插值
fun = @(x)sin(pi*x/2).*(x>=-1&x<1) + x.*exp(1-x.^2).*(x>=1 | x<-1);
%%----------------区间[0,1]上的三次样条插值------------------
x01 = linspace(0,1,6);
y01 = fun(x01); 
x1 = linspace(0,1,20);
pp1 = csape(x01,[1,y01,0],'complete');
y1 = fnval(pp1,x1); 
%%----------------区间[1,3]上的三次样条插值------------------
x02 = linspace(1,3,8);
y02 = fun(x02);   
x2 = linspace(1,3,30); 
pp2 = csape(x02,[0,y02,0.01],[1,2]);
y2 = fnval(pp2,x2);
%%-----------------------绘图---------------------
figure;
plot([x01,x02],[y01,y02],'ko');
hold on;
plot([x1,x2],fun([x1,x2]),'k','linewidth',2);
plot([x1,x2],[y1,y2],'--','linewidth',2);
xlabel('X');
ylabel('Y = f(x)');
legend('插值节点','原函数图像','三次样条插值');

%% examp7.4-3 二维网格节点插值
x = 100:100:500;
y = 100:100:400;
[X,Y] = meshgrid(x,y);
Z = [450  478  624  697  636
        420  478  630  712  698
        400  412  598  674  680
        310  334  552  626  662];
xd = 100:20:500;
yd = 100:20:400;
[Xd1,Yd1] = meshgrid(xd,yd);
[Xd2,Yd2] = ndgrid(xd,yd);

figure;  % 新建图形窗口
% -------------- 调用interp2函数作三次样条插值-------------------
Zd1 = interp2(X,Y,Z,Xd1,Yd1,'spline');
subplot(1,2,1);
surf(Xd1,Yd1,Zd1);
xlabel('X'); ylabel('Y'); zlabel('Z'); title('interp2')

% ---------调用griddedInterpolant函数作三次样条插值--------------
F = griddedInterpolant({x,y},Z','spline');
Zd2 = F(Xd2,Yd2);
subplot(1,2,2);
surf(Xd2,Yd2,Zd2);
xlabel('X'); ylabel('Y'); zlabel('Z'); title('griddedInterpolant')

%% examp7.4-4 二维散乱节点插值
xyz = xlsread('cumcm2011A.xls',1,'B4:D322');
Cd = xlsread('cumcm2011A.xls',2,'C4:C322');
x = xyz(:,1);
y = xyz(:,2);
z = xyz(:,3);
xd = linspace(min(x),max(x),60);
yd = linspace(min(y),max(y),60);
[Xd,Yd] = meshgrid(xd,yd);
% ------------调用griddata函数作散乱节点插值---------------
Zd1 = griddata(x,y,z,Xd,Yd);
Cd1 = griddata(x,y,Cd,Xd,Yd);
figure;
subplot(1,2,1);
surf(Xd,Yd,Zd1,Cd1);
shading interp;
xlabel('X'); ylabel('Y'); zlabel('Z'); title('griddata');
colorbar;

% ------------调用scatteredInterpolant函数作散乱节点插值------------
F1 = scatteredInterpolant(x,y,z,'linear','none');
Zd2 = F1(Xd,Yd);  % 计算插值点处的海拔高度
F2 = scatteredInterpolant(x,y,Cd,'linear','none');
Cd2 = F2(Xd,Yd);
subplot(1,2,2);
surf(Xd,Yd,Zd2,Cd2);
shading interp;
xlabel('X'); ylabel('Y'); zlabel('Z');title('scatteredInterpolant'); 
colorbar;

%% examp7.4-5 高维插值
data = xlsread('温度场.xlsx');
x = data(:,1);
y = data(:,2);
z = data(:,3);
v = data(:,4);
xi = linspace(min(x),max(x),60);
yi = linspace(min(y),max(y),60);
zi = linspace(min(z),max(z),60);
[Xd,Yd,Zd] = meshgrid(xi,yi,zi);
F = scatteredInterpolant(x,y,z,v);
Vd = F(Xd,Yd,Zd);

id = [1:15:60,60];  % 设置切片位置下标
figure;
% 绘制切片图
slice(Xd,Yd,Zd,Vd,xi(id),yi(id),zi(id));
shading interp;  % 插值染色
alpha(0.5);  % 设置透明度
xlabel('x');ylabel('y');zlabel('z')
axis equal
colorbar;  % 添加颜色条

figure;
MyIsosurface(Xd,Yd,Zd,Vd,800);

%% examp7.5-1 声呐定位
%------------------------------
% 1. 根据附件1中数据求目标物体所在的大致区域
%------------------------------
AEd = xlsread('声呐扫描数据.xlsx',1);
[x1,y1,z1] = sph2cart(AEd(:,1),AEd(:,2),AEd(:,3));
xyz_Range = minmax([x1,y1,z1]');
T1 = table(xyz_Range(:,1),xyz_Range(:,2),...
    'VariableNames',{'最小值','最大值'},...
    'RowNames',{'x','y','z'})
V = range(x1)*range(y1)*range(z1)

%------------------------------
% 2. 根据附件2中数据求目标物体的球心坐标与球半径
%------------------------------
AEd = xlsread('声呐扫描数据.xlsx',2);
[x2,y2,z2] = sph2cart(AEd(:,1),AEd(:,2),AEd(:,3));
A = [x2.^2+y2.^2+z2.^2, -2*x2, -2*y2, -2*z2];
B = ones(size(x2));
beta = A\B;
xyz0 = beta(2:4)/beta(1);
r = sqrt(1/beta(1) + xyz0'*xyz0);
T2 = table(xyz0(1),xyz0(2),xyz0(3),r,...
    'VariableNames',{'x0','y0','z0','r'})


%------------------------------
% 3. 根据附件3中数据求目标物体的运行轨迹
%------------------------------
AEd = xlsread('声呐扫描数据.xlsx',3);
[x3,y3,z3] = sph2cart(AEd(:,2),AEd(:,3),AEd(:,4));
figure;
plot3(x3,y3,z3,'.');
grid on;
axis equal;
xlabel('x(m)'); ylabel('y(m)'); zlabel('z(m)');

t = AEd(:,1);
D = zeros(10,3);
for i = 1:10
    id = (t == i);
    xi = x3(id); yi = y3(id); zi = z3(id);
    A = [xi.^2+yi.^2+zi.^2, -2*xi, -2*yi, -2*zi];
    B = ones(size(xi));
    beta = A\B;
    D(i,:) = beta(2:4)/beta(1);
end
t = (1:10)';
T3 = table(t,D(:,1),D(:,2),D(:,3),...
    'VariableNames',{'t','x','y','z'})

p1 = polyfit(T3.t, T3.x, 1);
p2 = polyfit(T3.t, T3.z, 2);
str1 = char(vpa(poly2sym(p1,sym('t')),5));
str1 = ['x = ',str1];
str2 = char(vpa(poly2sym(p2,sym('t')),5));
str2 = ['z = ',str2];
tnew = linspace(0,11,30);
xnew = polyval(p1,tnew);
znew = polyval(p2,tnew);
figure;
subplot(1,2,1);
plot(T3.t, T3.x,'ko');
hold on;
plot(tnew, xnew,'b');
text(1,-70,str1);
grid on;
xlabel('t(min)'); ylabel('x(m)');

subplot(1,2,2);
plot(T3.t, T3.z,'ko');
hold on;
plot(tnew, znew,'b');
text(1,130,str2);
grid on;
xlabel('t(min)'); ylabel('z(m)');

% 预测目标物体下一时刻的位置
tnext = 11;
xnext = polyval(p1,tnext);
ynext = mean(T3.y);
znext = polyval(p2,tnext);
T = table(tnext,xnext,ynext,znext,...
    'VariableNames',{'t','x','y','z'})

%% examp7.5-2 国土面积测量问题
data = xlsread('河南省省界线经纬度坐标.xlsx');
lat = data(:,1);
lon = data(:,2);
L = sqrt(diff(lon).^2 + diff(lat).^2);
L = cumsum([0;L]);
Li = linspace(0,max(L),1000);
lati = spline(L,lat,Li);
loni = spline(L,lon,Li);
R = 6371;
lat0 = 0;
lon0 = 0;
lati(end+1) = lati(1);
loni(end+1) = loni(1);
[phi,w] = distance('gc',lat0,lon0,lati,loni);
phi = phi*pi/180; w = w*pi/180;
ds = (1-cos((phi(1:end-1)+phi(2:end))/2)).*diff(w);
SD1 = R^2*abs(sum(ds))        % 方法1
SD2 = areaint(lati,loni,R)    % 方法2
