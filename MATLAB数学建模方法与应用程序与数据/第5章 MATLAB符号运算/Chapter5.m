%--------------------------------------------------------------------------
%  第5章  MATLAB符号计算
%--------------------------------------------------------------------------
% CopyRight：xiezhh


%% examp5.1-1
a = sym('6.01');         % 定义符号常数
b = sym('b','real');     % 定义实数域上的符号变量
A = [1, 2; 3, 4];        % 定义数值矩阵
B = sym(A);              % 把数值矩阵转为符号矩阵
C = sym('c%d%d',[3,4])

syms  x  y               % 同时定义多个复数域上的符号变量
syms  z  positive        % 定义正实数域上的符号变量
syms  f(x,y)             % 定义符号函数
f(x,y) = x + y^2;        % 指定符号函数表达式
c = f(1, 2) 

zv = solve(z^2 == 1, z)  % 求方程z^2 = 1的解（只有大于零的解）

syms  z           % 撤销对符号变量取值域的限定，将其恢复为复数域上的符号变量
zv = solve(z^2 == 1, z)

%% examp5.1-2
syms x                     % 定义符号变量
assume(x>0 & x<5);         % 对符号变量的取值域进行限定，0<x<5
assumeAlso(x,'integer');   % 对符号变量的取值域增加别的限定，x取整数
assumptions(x)             % 查看符号变量取值域的限定

result = solve(x^2>12)     % 求解不等式

syms  x

%% examp5.1-3
syms a b c x y z          % 定义多个符号变量
f1 = a*x^2+b*x-c;         % 创建符号表达式f1
f2 = sin(x)*cos(y);       % 创建符号表达式f2
f3 = (x+y)/z;             % 创建符号表达式f3
f4 = [x+1, x^2; x^3, x^4] % 创建符号表达式矩阵f4

f5 = f4'                  % 符号表达式矩阵的共轭转置（'）

f6 = f4.'

%% examp5.1-4
syms x y                % 定义符号变量
f1 = abs(x) >= 0        % 创建符号表达式

f2 = x^2 + y^2 == 1     % 创建符号表达式

f3 = ~(y - sqrt(x) > 0) % 创建符号表达式

f4 = x > 0 | y < -1     % 创建符号表达式

f5 = x > 0 & y < -1     % 创建符号表达式

%% examp5.1-5
syms x
f = abs(x) >= 0;                % 创建符号表达式
result1 = isAlways(f)           % 判断不等式|x|>=0是否成立

result2 = isequaln(abs(x), x)   % 判断|x|是否等于x

assume(x>0);                    % 限定x>0
result3 = isequaln(abs(x), x)   % 重新判断|x|是否等于x

syms x                             % 撤销对符号变量取值域的限定

%% examp5.1-6
syms x y
f = factor(x^3-y^3)

fa = factor(sym('12345678901234567890'))

%% examp5.1-7
syms x y
f = (x+y)*(x^2+y^2+1);
collect(f,y)

%% examp5.1-8
syms x y a b
f = [cos(x+y); (a+b)*exp((a-b)^2)];
expand(f)

%% examp5.1-9
syms x
f1 = sqrt(4/x^2+4/x+1);
g1 = simplify(f1)                                % 按默认设置进行化简

g2 = simplify(f1,'IgnoreAnalyticConstraints',1)  % 忽略分析约束进行化简

pretty(g2)                                       % 把符号表达式显示为数学公式形式

f2 = cos(3*acos(x));
g3 = simplify(f2, 'Steps', 4)                    % 进行4步化简

%% examp5.1-10
syms f(x) g(x)
f(x) = exp(x);
g(x) = sin(x);
y1(x) = f(g(x))
y2(x) = compose(f,g)
y = y1(pi)

%% examp5.1-11
syms f(x)
f(x) = exp(sin(x));
g(x) = finverse(f(x))

%% examp5.1-12
syms f(x)                      % 定义符号函数
f(x) = log(sym(5.2))*exp(x);   % 指定符号函数表达式
y = f(3)                       % 计算符号函数在x = 3处的函数值

y1 = double(y)                 % 把符号数转为双精度数

y2 = vpa(y,10)                 % 以10位有效数字形式显示符号数

x = 3;                         % 指定x的值
y3 = eval(f)                   % 执行MATLAB运算，得到函数值

%% examp5.1-13
syms a b x y
f = a*sin(x)+b;               % 定义符号表达式
f1 = subs(f,sin(x),log(y))  % 符号项替换

% 变量替换方式一
f2 = subs(f1,[a,b],[2,5])     % 同时替换变量a，b的值

% 变量替换方式二
f3 = subs(f1,{a,b},{2,5})     % 同时替换变量a，b的值

%% examp5.1-14
syms a b x                                    % 定义符号变量
f = a*sin(x)+b;                               % 定义符号表达式
y = subs(f, {a,b,x}, {2, 5, 1:3})             % 同时替换多个符号变量的值
y = double(y)                                 % 将计算结果转为双精度值

%% examp5.1-15
syms a b x
f(x) = symfun(a*sin(x)+b, x);     % 把符号表达式转为符号函数
y = f(1:3)

%% examp5.1-16
syms a b c d x
f = a*(x+b)^c+d;                          % 定义符号表达式
g = subs(f,{a,b,c,d},{2,-1,sym(1/2),3});  % 同时替换多个变量
FunFromSym1 = matlabFunction(g)           % 将符号表达式转为匿名函数

y = FunFromSym1(10)                       % 调用匿名函数计算函数值

% 将符号表达式转为M文件函数FunFromSym2.m
matlabFunction(g,'file',[pwd,'\FunFromSym2.m'],...
    'vars',{'x'},'outputs',{'y'});
y = FunFromSym2(10)                       % 调用M函数计算函数值

%% examp5.1-17
syms f(x)               % 定义符号函数
f(x) = 1/log(abs(x));   % 指定符号函数表达式
fplot(f,[-6,6]);       % 绘制函数图形
xlabel('x');
ylabel('$$ f(x) = \frac{1}{ln|x|} $$','Interpreter','Latex');

%% examp5.1-18
syms f(x,y)
f(x,y) = x*y/(x^2+y^2);
fsurf(f,[-1,1,-1,1])
xlabel('x');ylabel('y');zlabel('z')
hold on
plot3(0,0,0,'r.','MarkerSize',20)
title('$$f(x,y)=\frac{xy}{x^2+y^2}$$','Interpreter','latex')
view(127,34)

%% examp5.2-1
syms n a b c k x y
xn = (-1)^n/(n+1)^2;
L1 = limit(xn,n,inf)

f1 = sin(a*x)/(a*x);
L2 = limit(f1,x,0,'left')

f2 = (1-2/x)^(k*x);
L3 = limit(f2,x,inf)

f3 = a/(1+x^2+y^2);
L4 = limit(limit(f3,x,b),y,c)

%% examp5.2-2
syms x y
f(x) = sin(x)^2;
df = diff(f,x)
df_1 = df(1)

ddf = diff(f,x,2)

F(x,y) = cos(x+sin(y))-sin(y);
dy = -diff(F,x)/diff(F,y)

%% examp5.2-3
syms x1 x2
f = [x1+x2;x2*log(x1)];
v = [x1;x2];
jac = jacobian(f,v)

%% examp5.2-4
syms x
f = exp(x); 
g = taylor(f, x, 0, 'Order', 6)

%% examp5.2-5
syms k
f1 = (k-2)/2^k;
s1 = symsum(f1,k,3,inf)

f2 = [1/(2*k+1)^2,  (-1)^k/3^k];
s2 = symsum(f2,k,1,inf)

%% examp5.2-6
syms x y z a
F = int(x*log(a*x),x)

f1 = sqrt(1-x^2);
s1 = int(f1,x,-1,1)

f2 = exp(-x^2/2);
s2 = int(f2,x,-inf,inf)

f3 = (x+y)/z;
s3 = int(int(int(f3,z,x*y,2*x*y),y,x,2*x),x,1,2)
s4 = double(s3) 

%% examp5.3-1
syms x
Result1 = solve(x^3 - 2*x^2 + 4*x == 8, x)

Result2 = solve(sin(x) + cos(2*x) == 1, x)

[Result3,params,conditions] = solve(sin(x) + cos(2*x) == 1, x, 'ReturnConditions',true)

Result4 = solve(x + x*exp(x) == 10, x)

%% examp5.3-2
syms x y
[X,Y] = solve([1/x^3 + 1/y^3 == 28, 1/x + 1/y == 4], [x,y]) 

%% examp5.3-3
syms y(x)
Y = dsolve(diff(y,2) == x+y)

%% examp5.3-4
syms y(t)
Y = dsolve(diff(y) == 1 + y^2, y(0) == 1)

Y = dsolve(diff(y) == 1 + y^2, y(0) == 1, 'IgnoreAnalyticConstraints', false)

%% examp5.3-5
syms y(x)
Y = dsolve(x*diff(y,2)-3*diff(y) == x^2, [y(1) == 0, y(5) == 0])

h = fplot(Y,[-1,6]);
set(h,'color','k','LineWidth',2,'LineStyle','--');
hold on;
plot([1 5],[0,0],'p','color','r','markersize',12); %画微分方程的两个边值点
text(1,1,'y(1) = 0'); %图上标注边值条件
text(4,1,'y(5) = 0');
title('');
hold off;

%% examp5.3-6
syms y(x)
Dy = diff(y);
eq = x*diff(y,2)-3*Dy == x^2;
Y = dsolve(eq, [y(1) == 0, Dy(5) == 1])

%% examp5.3-7
syms x(t)  y(t)
[X, Y] = dsolve(diff(x) == y, diff(y) == -x)

%% examp5.4-1
syms x(t) y(t) lambda mu % 定义符号函数
eq1 = [diff(x) == -lambda*x,diff(y) == lambda*x-mu*y];  % 定义微分方程
cond = [x(0) == 1100,y(0) == 0];                        % 定义初值条件
[x(t),y(t)] = dsolve(eq1,cond)                          % 求解微分方程

lambda = solve(x(5) == x(0)/2)                   % 求解血液对药物的吸收率系数
syms y2(t) mu tau a
y2(t) = dsolve(diff(y2) == -mu*y2, y2(tau) == a)
mu = solve(y2(tau+6) == y2(tau)/2,mu)  % 求解血液对药物的排除率

x(t) = subs(x,'lambda',lambda);
y(t) = subs(y,{'lambda','mu'},{lambda,mu});

y2h = double(subs(y,2))           % t = 2h 时血液中药量
tB = double(solve(y(t) == 200))   % 严重中毒时刻
tC = double(solve(y(t) == 400))   % 致命中毒时刻
tD = double(solve(diff(y) == 0))  % 不加施救情形下，血药量达到峰值时刻
ymax = double(subs(y,tD))         % 最大血药量

% 结果可视化
figure
fplot(x,[0,25],'b--')  % 胃肠道中药量变化曲线
hold on
fplot(y,[0,25],'r')  % 血液中药量变化曲线
grid on
set(gca,'XMinorGrid','on','YMinorGrid','on')
xlabel('时间 t（h）')
ylabel('药量（mg）')
plot([2,tB,tC,tD],[y2h,200,400,ymax],'b*')
text([2,tB,tC,tD]+0.5,[y2h,200,400,ymax],...
    {'A','B','C','D'})
legend('x(t)','y(t)')

% ---------------施救方案：口服活性炭-------------------
syms x(t) z(t)
lambda = log(2)/5;
mu = 2*log(2)/6;
eq2 = [diff(x) == -lambda*x,diff(z) == lambda*x-mu*z];  % 微分方程
cond2 = [x(0) == 1100,z(2) == 236.5588];  % 初值条件
[x(t),z(t)] = dsolve(eq2,cond2);
z = vpa(z,5)

% 结果可视化
figure
fplot(x,[0,25],'b--')  % 胃肠道中药量变化曲线
hold on
fplot(y,[0,25],'r-.')  % 血液中药量变化曲线（不施救）
fplot(z,[2,25],'k')  % 血液中药量变化曲线（施救）
grid on
set(gca,'XMinorGrid','on','YMinorGrid','on')
xlabel('时间 t（h）')
ylabel('药量（mg）')
legend('x(t)','y(t)','z(t)')

t4 = double(solve(diff(z) == 0))  % 施救情形下，血药量达到峰值时刻
zmax = double(subs(z,t4))         % 最大血药量


% ---------------施救方案：体外血液透析-------------------
syms x(t) z(t)                                         % 定义符号函数
lambda = log(2)/5;                                     % 定义lambda
mu = 6*log(2)/6;                                       % 定义mu
eq2 = [diff(x) == -lambda*x,diff(z) == lambda*x-mu*z]; % 定义微分方程
cond2 = [x(0) == 1100,z(2) == 236.5588];               % 定义初值条件
[x(t),z(t)] = dsolve(eq2,cond2);                       % 求解微分方程
z = vpa(z,5)                                           % 以5位有效数字形式显示

% 结果可视化
figure                                                 % 创建figure窗口
fplot(x,[0,25],'b--')                                  % 胃肠道中药量变化曲线
hold on                                              % 开启图形保持
fplot(y,[0,25],'r-.')                                % 血液中药量变化曲线（不施救）
fplot(z,[2,25],'k')                                  % 血液中药量变化曲线（施救）
grid on                                              % 添加主网格线
set(gca,'XMinorGrid','on','YMinorGrid','on')         % 添加次网格线
xlabel('时间 t（h）');  ylabel('药量（mg）');           % 添加坐标轴标签
legend('x(t)','y(t)','z(t)')                         % 添加图例
