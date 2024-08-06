function f = framp(p,t,u,time)
% p    --- 三角网顶点坐标
% t    --- 三角网顶点编号
% u    --- u(x,y,t)
% time --- 时间
n = size(t,2);  % 网格个数
f = zeros(1,n); % 定义0向量
it1 = t(1,:);   % 每一个网格上第1个顶点编号
it2 = t(2,:);   % 每一个网格上第2个顶点编号
it3 = t(3,:);   % 每一个网格上第3个顶点编号
xpts = (p(1,it1)+p(1,it2)+p(1,it3))/3;  % 每一个网格中心点x坐标
ypts = (p(2,it1)+p(2,it2)+p(2,it3))/3;  % 每一个网格中心点y坐标
[~,id] = min(xpts.^2 + ypts.^2);  % 寻找位于求解区域的中心的网格
f(id) = 5*sin(1.7*time);  % 函数赋值
end