function MyGateLegTable2
% ��������70 cm������ֱ��80 cm��Բ���۵������۵�Ч��ͼ

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