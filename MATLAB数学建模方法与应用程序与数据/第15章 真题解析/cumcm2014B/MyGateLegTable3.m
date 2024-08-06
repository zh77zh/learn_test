function MyGateLegTable3
% 绘制桌高60 cm，桌面直径60 cm的八边形折叠桌的折叠效果图

R = 30;
w = 2.5/2;
t = 3;
H = 60-t;
k = tand(15);
EdgeFun = @(x)(30*(1+k)+x).*(x>=-30 & x<-30*k) + 30*(x>=-30*k & x<30*k) + ...
    (30*(1+k)-x).*(x>=30*k & x<=30);
Problem = MakeObjFun(R,H,w,t,EdgeFun);
[Lb,fval] = fmincon(Problem)
L = Lb(1);
b = Lb(2);

PlotGateLegTable(R,L,H,w,t,EdgeFun,b);
end