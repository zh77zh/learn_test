function MyGateLegTable2
% 绘制桌高70 cm，桌面直径80 cm的圆形折叠桌的折叠效果图

R = 40;
w = 2.5/2;
t = 3;
H = 70-t;
EdgeFun = @(x)sqrt(R^2-x.^2);
Problem = MakeObjFun(R,H,w,t,EdgeFun);
[Lb,fval] = fmincon(Problem)
L = Lb(1);
b = Lb(2);

PlotGateLegTable(R,L,H,w,t,EdgeFun,b);
end