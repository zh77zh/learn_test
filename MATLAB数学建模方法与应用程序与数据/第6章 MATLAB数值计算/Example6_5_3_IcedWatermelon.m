function Example6_5_3_IcedWatermelon
% 把西瓜放冰箱里，多长时间才能冻透？
% 求解冰镇西瓜过程中的温度扩散问题（热传导方程）
%       1   ∂T     1  ∂        ∂T
%      ---.---- = ---.---( r^2.---)
%       a   ∂t    r^2 ∂r       ∂r
%
%   初始条件：
%            t = 0, 0 ≤ r ≤ R, T = T0
%   边值条件：
%            t > 0, r = 0, ∂T/∂r = 0
%            t > 0, r = R, h(T - Tinf) + k.∂T/∂r = 0
%
%   1. 假设西瓜的外形是一个半径为R的完美的球体
%   2. 假设西瓜没有瓜皮，没有瓜籽，瓜瓤也是完全均匀的，其密度、导热能力、比热容
%      全是均匀的，一切物性参数都不随温度变化而变化
%   3. 冰箱是一个恒温为Tinf ℃的环境
%
%   August 4  2017, editted by ZhongHua-Xie, Tianjin University of Science and Technology.

% 相关参数
h = 5;        % 西瓜与静止冷空气的对流换热系数 [w/(m^2.K)]
k = 0.48;     % 西瓜的导热系数 [w/(m.K)]
p = 918;      % 西瓜的密度 [kg/m^3]
Cp = 3990;    % 热容 [J/(kg.K)]
a = k/(p*Cp); % 西瓜的热扩散系数 [m^2/s]
T0 = 32;      % 初始温度 [℃]
Tinf = 6;     % 终极温度 [℃]，即冰箱冷藏室温度

m = 2;                                     % 方程中的m参数
r = linspace(0,9,41)/100;                  % 定义距离向量r
t = linspace(0,10,41)*3600;                % 定义时间向量t
T = pdepe(m,@pdefun,@pdeic,@pdebc,r,t);    % 方程求解

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
num = size(T,1);
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
surf(X_bac,Y_bac,Tinf*ones(size(X_bac)),...
    'FaceColor','interp','EdgeColor','none');
hold on;
% 绘制零时刻西瓜等效温度切面
hs = surf(X,Y,Ti,'FaceColor','interp','EdgeColor','none');
axis equal;
axis([xmin,xmax,ymin,ymax])
xlabel('X');
ylabel('Y');
caxis([Tinf,max(T(:))]);
colorbar;
view(2);
ht = title(['冷藏 ',num2str(t(1)/3600,'%2.1f'),' 小时']);

filename = '冰镇西瓜温度扩散动画.gif';
f = getframe(gcf);
IM = frame2im(f);
[IM,map] = rgb2ind(IM,256);
imwrite(IM,map,filename,'gif', 'Loopcount',inf,'DelayTime',0.2);

% 通过循环绘制各时刻西瓜等效温度切面
for i = 2:num
    Ti = T(i,:)'*ones(size(theta));
    set(hs,'ZData',Ti);
    set(ht,'String',['冷藏 ',num2str(t(i)/3600,'%2.1f'),' 小时']);
    pause(0.2);
    
    f = getframe(gcf);
    IM = frame2im(f);
    [IM,map] = rgb2ind(IM,256);
    imwrite(IM,map,filename,'gif','WriteMode','append','DelayTime',0.2);   
end

    function [c,f,s] = pdefun(r,t,T,dT)
        % 偏微分方程函数
        c = 1/a;
        f = dT;
        s = 0;
    end

    function u0 = pdeic(r)
        % 初值条件函数
        u0 = T0;
    end

    function [pa,qa,pb,qb] = pdebc(ra,Ta,rb,Tb,t)
        % 边值条件函数
        pa = 0;
        qa = 1;
        pb = h*(Tb - Tinf);
        qb = k;
    end

end
