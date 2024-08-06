function MyGateLegTable1
% 绘制桌高53 cm，桌面直径50 cm的圆形折叠桌的折叠效果图

R = 25;
L = 120/2;
w = 2.5/2;
t = 3;
H = 53-t;
EdgeFun = @(x)sqrt(R^2-x.^2);
PlotGateLegTable(R,L,H,w,t,EdgeFun);
end