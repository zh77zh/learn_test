function MyGateLegTable1
% ��������53 cm������ֱ��50 cm��Բ���۵������۵�Ч��ͼ

R = 25;
L = 120/2;
w = 2.5/2;
t = 3;
H = 53-t;
EdgeFun = @(x)sqrt(R^2-x.^2);
PlotGateLegTable(R,L,H,w,t,EdgeFun);
end