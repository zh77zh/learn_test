% 图6.4-3程序
u = zeros(301);                 % 定义零矩阵
dt = 1/300; c = 0.03;           % 参数
x = linspace(0,1,301); t = x';  % t和x的划分向量
u(:,1) = x.*(1-x)/10;           % 初值
u(1,:) = sin(t);                % 边值
v = sin(2*pi*x);                % 初速度
% 计算u(i,2)
u(2:300,2) = (1-c)*u(2:300,1) + ...
       1/2*c*(u(1:299,1)+u(3:301,1)) + v(2:300)'*dt;

% 用有限差分法求解方程
for j = 3:301
       u(2:300,j) = 2*(1-c)*u(2:300,j-1)+c*(u(3:301,j-1)+...
           u(1:299,j-1))-u(2:300,j-2);
end

tk = 0:0.2:1;
for k = 1:numel(tk)
    subplot(2,3,k)
    id = (t == tk(k));
    plot(x,u(:,id));                % 绘制初始弦曲线
    axis([-0.1,1.1,-1,1]);          % 设置坐标轴范围
    xlabel('x');ylabel('U(x)');     % 坐标轴标签
    text(0.8,0.8,['t = ',num2str(tk(k))])
end