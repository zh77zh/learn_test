% 平板折叠桌平展状态示意图
% 绘制折叠桌折叠的动画效果图，并计算各木条开槽长度
% R：桌子半径
% L：平板半长
% h：桌子展开后高度（桌子的厚度不计算在内）
% w：木条半宽
% EdgeFun:桌子外边缘函数
% t：桌子的厚度

R = 25;
L = 60;
h = 50;
w = 2.5/2;
t = 3;
EdgeFun = @(x)sqrt(R^2-x.^2);
xi = -(R-w):2*w:R-w;
yi = EdgeFun(xi);
a = yi(1);
b = (L-a)/2;
n = numel(xi);
fac = [1 2 3 4;2 6 7 3;4 3 7 8;1 5 8 4;1 2 6 5;5 6 7 8]; 

for i = 1:n
    verti1 = [xi(i)-w,yi(i),-t;xi(i)+w,yi(i),-t;xi(i)+w,-yi(i),-t;xi(i)-w,-yi(i),-t;...
        xi(i)-w,yi(i),0;xi(i)+w,yi(i),0;xi(i)+w,-yi(i),0;xi(i)-w,-yi(i),0];
    patch('faces',fac,'vertices',verti1,'FaceColor',[1,1,0]);
    verti2{i} = [xi(i)-w,L,-t;xi(i)+w,L,-t;xi(i)+w,yi(i),-t;xi(i)-w,yi(i),-t;...
        xi(i)-w,L,0;xi(i)+w,L,0;xi(i)+w,yi(i),0;xi(i)-w,yi(i),0];
    h1(i) = patch('faces',fac,'vertices',verti2{i},'FaceColor',[1,0,0]);
    verti3{i} = [xi(i)-w,-L,-t;xi(i)+w,-L,-t;xi(i)+w,-yi(i),-t;xi(i)-w,-yi(i),-t;...
        xi(i)-w,-L,0;xi(i)+w,-L,0;xi(i)+w,-yi(i),0;xi(i)-w,-yi(i),0];
    h2(i) = patch('faces',fac,'vertices',verti3{i},'FaceColor',[1,0,0]);
end
axis equal;
view(90,90);
axis([-R,R,-L,L]);
axis off;
hold on;
Bet = linspace(0,2*pi,100);
plot(R*cos(Bet),R*sin(Bet),'b');

CosTheta = sqrt((L-a)^2-h^2)/(L-a);
di = sqrt((yi-a).^2+(L-a-b)^2-2*(yi-a)*(L-a-b)*CosTheta);   % 各木条的开槽长度
w = w/2;
% 绘制各木条开槽位置示意图
for i = 1:n
    verti1 = [xi(i)-w,yi(i)+di(i),-t;xi(i)+w,yi(i)+di(i),-t;xi(i)+w,L-b,-t;xi(i)-w,L-b,-t;...
        xi(i)-w,yi(i)+di(i),0;xi(i)+w,yi(i)+di(i),0;xi(i)+w,L-b,0;xi(i)-w,L-b,0];
    patch('faces',fac,'vertices',verti1,'FaceColor',[1,1,1]);
    verti2 = [xi(i)-w,-yi(i)-di(i),-t;xi(i)+w,-yi(i)-di(i),-t;xi(i)+w,b-L,-t;xi(i)-w,b-L,-t;...
        xi(i)-w,-yi(i)-di(i),0;xi(i)+w,-yi(i)-di(i),0;xi(i)+w,b-L,0;xi(i)-w,b-L,0];
    patch('faces',fac,'vertices',verti2,'FaceColor',[1,1,1]);
end
