m = 2;                                     % 方程中的m参数
r = linspace(0,9,41)/100;                  % 定义距离向量r
t = linspace(0,10,41)*3600;                % 定义时间向量t
T = pdepe(m,@pdefun2,@pdeic2,@pdebc2,r,t);   % 模型求解

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
Tinf = 6;     % 终极温度 [℃]，即冰箱冷藏室温度
k = 1;
for i = [1 13 25 41]
    subplot(2,2,k)
    % 绘制背景温度面
    surf(X_bac,Y_bac,Tinf*ones(size(X_bac)),...
        'FaceColor','interp','EdgeColor','none');
    hold on;
    % 绘制零时刻西瓜等效温度切面
    Ti = T(i,:)'*ones(size(theta));
    hs = surf(X,Y,Ti,'FaceColor','interp','EdgeColor','none');
    axis equal;
    axis([xmin,xmax,ymin,ymax])
    xlabel('X');
    ylabel('Y');
    caxis([Tinf,max(T(:))]);
    colorbar;
    view(2);
    ht = title(['冷藏 ',num2str(t(i)/3600,'%2.1f'),' 小时']);
    k = k+1;
end