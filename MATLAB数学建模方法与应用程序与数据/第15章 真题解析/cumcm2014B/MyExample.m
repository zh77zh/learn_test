%% 用户自定义参数及桌面边缘线（唇形桌面）
R = 40;
w = 2.5/2; 
t = 3;
H = 83-t;
EdgeFun = @(x)20-15*sin(6*exp(-abs(x/25)));
Problem = MakeObjFun(R,H,w,t,EdgeFun);
[Lb,fval] = fmincon(Problem)
L = Lb(1);
b = Lb(2);
PlotGateLegTable(R,L,H,w,t,EdgeFun,b);