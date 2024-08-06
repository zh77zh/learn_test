function Problem = MakeObjFun(R,H,w,t,EdgeFun)
% 计算长方形平板材料和圆形折叠桌的最优设计加工参数
% R：桌子半径
% L：平板半长
% H：桌子折叠后高度（桌子的厚度不计算在内）
% w：木条半宽
% t：桌子的厚度
% EdgeFun：桌子边缘函数
% Problem：优化问题
%    目标函数：min 开槽长度  &  max 转动惯量
%    约束条件：
%          a + H < L < 5R
%          a + b < L < 5R
%          0 < Theta < pi/2
%          0 < Beta  < pi/2
%          (a+R)/2 < yA1 < R    % yA1是最外侧木条外端点 A 的 y 坐标
%          zAi < zA1            % zAi是第i根木条外端点 A 的 z 坐标，

xi = -(R-w):2*w:R-w;    % 木条中心线x坐标
yi = real(EdgeFun(xi)); % 桌面边缘y坐标
a = yi(1);              % 桌面上最短木条一半长
ymaxi = max(yi);        % yi最大值
n = numel(xi);          % 木条数目
h = H;                  % 桌高（不考虑桌子厚度）
Problem.objective = @ObjFcn;    % 设置目标函数
Problem.nonlcon  = @NlinConFcn; % 设置非线性约束函数
Problem.x0 = [H+a,H+a-ymaxi];   % 设置初值
Problem.Aineq = [-1,0;1,0;-1,1];% 线性不等式约束系数矩阵
Problem.bineq = -[a+h;-5*R;ymaxi+t];  % 线性不等式约束常数向量
Problem.Aeq = [];               % 线性等式约束系数矩阵
Problem.beq = [];               % 线性等式约束常数向量
Problem.lb = [max([R,h+a]);0];  % 变量下界
Problem.ub = [];                % 变量上界
Problem.solver = 'fmincon';     % 设置求解器函数
Problem.options = '';           % 设置优化参数

function fval = ObjFcn(Lb)
% 目标函数：min 开槽长度/转动惯量
% Lb = [L,b]
% L：平板半长
% b：钢筋位置参数

L = Lb(1);
b = Lb(2);
J = 0;  % 转动惯量初始值
yz0 = [a+sqrt((L-a)^2-h^2),h];  % 转轴的yz坐标
D = zeros(1,n);
for k = 1:n
    J = J + sum((yz0-[0,-t/2]).^2)*2*yi(k);  % 桌面上各木条的转动惯量
    % 平展状态下y轴正向一侧桌腿两端坐标
    verti2 = [xi(k),L,-t/2;xi(k),yi(k),-t/2]; 
    % 计算第k根桌腿与桌面夹角余弦及铰链到钢筋的距离
    [CosBetai,d] = SubFcn(yi(k),h,L,b);
    SinBetai = sqrt(1-CosBetai^2);  % 夹角正弦
    Ai = repmat([xi(k),yi(k),0],[2,1]);   % 旋转点（铰链处）坐标
    Rot = [1,0,0;0,CosBetai,SinBetai;0,-SinBetai,CosBetai]; % 旋转矩阵
    % 计算第k根桌腿旋转后两端点坐标
    vert_Rot = (verti2-Ai)*Rot+Ai;
    xyzLeft = mean(vert_Rot);  % 第k根桌腿重心坐标
    J = J + sum((yz0-xyzLeft(2:3)).^2)*(L-yi(k)); % y轴正向一侧桌腿的转动惯量
    D(k) = d;
    % 以下计算y轴负向一侧桌腿的转动惯量
    verti2 = [xi(k),-L,-t/2;xi(k),-yi(k),-t/2]; 
    Ai = repmat([xi(k),-yi(k),0],[2,1]);
    Rot = [1,0,0;0,CosBetai,-SinBetai;0,SinBetai,CosBetai];
    vert_Rot = (verti2-Ai)*Rot+Ai;
    xyzRight = mean(vert_Rot);
    J = J + sum((yz0-xyzRight(2:3)).^2)*(L-yi(k));
end

Li = abs(D-(L-b-yi));   % 各木条的开槽长度
Ltotal = sum(Li);       % 开槽总长度
fval = Ltotal/J;        % 目标函数值：开槽总长度/转动惯量
end

%-------------------------------------------------------------
%  子函数
%-------------------------------------------------------------
function [CosBeta,d] = SubFcn(y,z,L,b)
    % 求各木条旋转角度的负余弦值;
    CosTheta = sqrt((L-a)^2-z^2)/(L-a);
    d = sqrt((y-a)^2+(L-a-b)^2-2*(y-a)*(L-a-b)*CosTheta);
    if y == a
        CosBeta = CosTheta;
    else
        CosBeta = -((y-a)^2+d^2-(L-a-b)^2)/(2*(y-a)*d);
    end
end

function [C,Ceq] = NlinConFcn(Lb)
    % 非线性约束函数
    L = Lb(1);  % 桌子半长
    b = Lb(2);  % 钢筋位置参数
    [CosBeta,di] = SubFcn(EdgeFun(0),h,L,b);  % 旋转角度的负余弦值
    Beta = acos(-CosBeta); % 旋转角度
    SinTheta = h/(L-a);  % theta角的正弦值
    Theta = asin(SinTheta); % theta角
    yA1 = a + (L-a)*cos(Theta);  % 最外侧木条外端点 A 的 y 坐标(a+R)/2 < yA1 < R+t
    Ceq = di-(L-a-b);  % 非线性等式约束（等腰三角形结构）
    C = [-Beta;Beta-pi/2;-Theta;Theta-pi/2;yA1-R;(a+R)/2-yA1]; % 非线性不等式约束
    zLeftbottom = zeros(1,n);
    % 通过循环计算各桌脚的z坐标
    for k = 1:n
        [CosBetai,di] = SubFcn(yi(k),h,L,b);
        SinBetai = sqrt(1-CosBetai^2);
        zLeftbottom(k) = (L-yi(k))*SinBetai;
    end
    % 非线性不等式约束
    C = [C;max(zLeftbottom(2:end-1))-zLeftbottom(1);di'+yi'-L];
end
end