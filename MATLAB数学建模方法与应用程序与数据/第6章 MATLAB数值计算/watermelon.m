function watermelon
% 绘制一个西瓜
%   August 4  2017, editted by ZhongHua-Xie, Tianjin University of Science and Technology.

% ----绘制西瓜面----
% 西瓜相关参数
a = 10;
b = 10;
c = 15;
% 通过球面坐标定义西瓜外表面网格坐标
Pj = linspace(-pi/2,pi/2,60);
Tj = linspace(0,2*pi,60);
[P,T] = meshgrid(Pj,Tj);
X = a*cos(P).*cos(T);
Y = b*cos(P).*sin(T);
Z = c*sin(P);
id = (Z > 0) & (T > 3*pi/2 & T < 2*pi);
Z(id) = NaN;  % 镂空处理
C = sin(T*10 + rand(size(T)));  % 西瓜外表面颜色索引矩阵
% 自定义颜色矩阵
ColorMat = linspace(0,0.2,64)';
ColorMat = [ColorMat,-3.5*ColorMat + 1,ColorMat];
surf(X,Y,Z,C,'FaceColor','interp','EdgeColor','none')  % 绘制西瓜外表面
colormap(ColorMat)  % 设置西瓜外表面颜色
hold on;
view(31,22)
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');

% ----绘制瓜秧（三维螺线）----
Z2 = linspace(0,pi,100);
X2 = Z2.*cos(5*Z2)/2;
Y2 = Z2.*sin(5*Z2)/2;
plot3(X2,Y2,Z2 + 15,'Color',[0.2,0.7,0.2],'LineWidth',3);

% ----绘制西瓜外表面镂空部分对应的内切面----
[P,R] = meshgrid(linspace(0,pi/2,30),linspace(1,0,30));
X = a*cos(P).*R*cos(3*pi/2);
Y = b*cos(P).*R*sin(3*pi/2);
Z = c*sin(P).*R;
surf(X,Y,Z,'FaceColor','r','EdgeColor','none','FaceAlpha',1)
rng(3)
id = randi(900,5,1);
plot3(X(id),Y(id),Z(id),'k.','MarkerSize',12)  % 绘制西瓜籽

X = a*cos(P).*R*cos(2*pi);
Y = b*cos(P).*R*sin(2*pi);
Z = c*sin(P).*R;
surf(X,Y,Z,'FaceColor','r','EdgeColor','none','FaceAlpha',1)
plot3(X(id),Y(id),Z(id),'k.','MarkerSize',12) % 绘制西瓜籽

[T,R] = meshgrid(linspace(3*pi/2,2*pi,30),linspace(1,0,30));
X = a*cos(T).*R;
Y = b*sin(T).*R;
Z = 0*R;
surf(X,Y,Z,'FaceColor','r','EdgeColor','none','FaceAlpha',1)
plot3(X(id),Y(id),Z(id),'k.','MarkerSize',12) % 绘制西瓜籽

% ----绘制坐标轴线----
quiver3(0,0,0,1,0,0,a,'Color',[0.5,0.5,0.2],'ShowArrowHead','off')
quiver3(0,0,0,0,-1,0,b,'Color',[0.5,0.5,0.2],'ShowArrowHead','off')
quiver3(0,0,0,0,0,1,c,'Color',[0.5,0.5,0.2],'ShowArrowHead','off')

caxis([-1,1])
% camlight
% lighting phong