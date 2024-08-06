%--------------------------------------------------------------------------
%  ��5��  MATLAB���ż���
%--------------------------------------------------------------------------
% CopyRight��xiezhh


%% examp5.1-1
a = sym('6.01');         % ������ų���
b = sym('b','real');     % ����ʵ�����ϵķ��ű���
A = [1, 2; 3, 4];        % ������ֵ����
B = sym(A);              % ����ֵ����תΪ���ž���
C = sym('c%d%d',[3,4])

syms  x  y               % ͬʱ�������������ϵķ��ű���
syms  z  positive        % ������ʵ�����ϵķ��ű���
syms  f(x,y)             % ������ź���
f(x,y) = x + y^2;        % ָ�����ź������ʽ
c = f(1, 2) 

zv = solve(z^2 == 1, z)  % �󷽳�z^2 = 1�Ľ⣨ֻ�д�����Ľ⣩

syms  z           % �����Է��ű���ȡֵ����޶�������ָ�Ϊ�������ϵķ��ű���
zv = solve(z^2 == 1, z)

%% examp5.1-2
syms x                     % ������ű���
assume(x>0 & x<5);         % �Է��ű�����ȡֵ������޶���0<x<5
assumeAlso(x,'integer');   % �Է��ű�����ȡֵ�����ӱ���޶���xȡ����
assumptions(x)             % �鿴���ű���ȡֵ����޶�

result = solve(x^2>12)     % ��ⲻ��ʽ

syms  x

%% examp5.1-3
syms a b c x y z          % ���������ű���
f1 = a*x^2+b*x-c;         % �������ű��ʽf1
f2 = sin(x)*cos(y);       % �������ű��ʽf2
f3 = (x+y)/z;             % �������ű��ʽf3
f4 = [x+1, x^2; x^3, x^4] % �������ű��ʽ����f4

f5 = f4'                  % ���ű��ʽ����Ĺ���ת�ã�'��

f6 = f4.'

%% examp5.1-4
syms x y                % ������ű���
f1 = abs(x) >= 0        % �������ű��ʽ

f2 = x^2 + y^2 == 1     % �������ű��ʽ

f3 = ~(y - sqrt(x) > 0) % �������ű��ʽ

f4 = x > 0 | y < -1     % �������ű��ʽ

f5 = x > 0 & y < -1     % �������ű��ʽ

%% examp5.1-5
syms x
f = abs(x) >= 0;                % �������ű��ʽ
result1 = isAlways(f)           % �жϲ���ʽ|x|>=0�Ƿ����

result2 = isequaln(abs(x), x)   % �ж�|x|�Ƿ����x

assume(x>0);                    % �޶�x>0
result3 = isequaln(abs(x), x)   % �����ж�|x|�Ƿ����x

syms x                             % �����Է��ű���ȡֵ����޶�

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
g1 = simplify(f1)                                % ��Ĭ�����ý��л���

g2 = simplify(f1,'IgnoreAnalyticConstraints',1)  % ���Է���Լ�����л���

pretty(g2)                                       % �ѷ��ű��ʽ��ʾΪ��ѧ��ʽ��ʽ

f2 = cos(3*acos(x));
g3 = simplify(f2, 'Steps', 4)                    % ����4������

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
syms f(x)                      % ������ź���
f(x) = log(sym(5.2))*exp(x);   % ָ�����ź������ʽ
y = f(3)                       % ������ź�����x = 3���ĺ���ֵ

y1 = double(y)                 % �ѷ�����תΪ˫������

y2 = vpa(y,10)                 % ��10λ��Ч������ʽ��ʾ������

x = 3;                         % ָ��x��ֵ
y3 = eval(f)                   % ִ��MATLAB���㣬�õ�����ֵ

%% examp5.1-13
syms a b x y
f = a*sin(x)+b;               % ������ű��ʽ
f1 = subs(f,sin(x),log(y))  % �������滻

% �����滻��ʽһ
f2 = subs(f1,[a,b],[2,5])     % ͬʱ�滻����a��b��ֵ

% �����滻��ʽ��
f3 = subs(f1,{a,b},{2,5})     % ͬʱ�滻����a��b��ֵ

%% examp5.1-14
syms a b x                                    % ������ű���
f = a*sin(x)+b;                               % ������ű��ʽ
y = subs(f, {a,b,x}, {2, 5, 1:3})             % ͬʱ�滻������ű�����ֵ
y = double(y)                                 % ��������תΪ˫����ֵ

%% examp5.1-15
syms a b x
f(x) = symfun(a*sin(x)+b, x);     % �ѷ��ű��ʽתΪ���ź���
y = f(1:3)

%% examp5.1-16
syms a b c d x
f = a*(x+b)^c+d;                          % ������ű��ʽ
g = subs(f,{a,b,c,d},{2,-1,sym(1/2),3});  % ͬʱ�滻�������
FunFromSym1 = matlabFunction(g)           % �����ű��ʽתΪ��������

y = FunFromSym1(10)                       % ���������������㺯��ֵ

% �����ű��ʽתΪM�ļ�����FunFromSym2.m
matlabFunction(g,'file',[pwd,'\FunFromSym2.m'],...
    'vars',{'x'},'outputs',{'y'});
y = FunFromSym2(10)                       % ����M�������㺯��ֵ

%% examp5.1-17
syms f(x)               % ������ź���
f(x) = 1/log(abs(x));   % ָ�����ź������ʽ
fplot(f,[-6,6]);       % ���ƺ���ͼ��
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
plot([1 5],[0,0],'p','color','r','markersize',12); %��΢�ַ��̵�������ֵ��
text(1,1,'y(1) = 0'); %ͼ�ϱ�ע��ֵ����
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
syms x(t) y(t) lambda mu % ������ź���
eq1 = [diff(x) == -lambda*x,diff(y) == lambda*x-mu*y];  % ����΢�ַ���
cond = [x(0) == 1100,y(0) == 0];                        % �����ֵ����
[x(t),y(t)] = dsolve(eq1,cond)                          % ���΢�ַ���

lambda = solve(x(5) == x(0)/2)                   % ���ѪҺ��ҩ���������ϵ��
syms y2(t) mu tau a
y2(t) = dsolve(diff(y2) == -mu*y2, y2(tau) == a)
mu = solve(y2(tau+6) == y2(tau)/2,mu)  % ���ѪҺ��ҩ����ų���

x(t) = subs(x,'lambda',lambda);
y(t) = subs(y,{'lambda','mu'},{lambda,mu});

y2h = double(subs(y,2))           % t = 2h ʱѪҺ��ҩ��
tB = double(solve(y(t) == 200))   % �����ж�ʱ��
tC = double(solve(y(t) == 400))   % �����ж�ʱ��
tD = double(solve(diff(y) == 0))  % ����ʩ�������£�Ѫҩ���ﵽ��ֵʱ��
ymax = double(subs(y,tD))         % ���Ѫҩ��

% ������ӻ�
figure
fplot(x,[0,25],'b--')  % θ������ҩ���仯����
hold on
fplot(y,[0,25],'r')  % ѪҺ��ҩ���仯����
grid on
set(gca,'XMinorGrid','on','YMinorGrid','on')
xlabel('ʱ�� t��h��')
ylabel('ҩ����mg��')
plot([2,tB,tC,tD],[y2h,200,400,ymax],'b*')
text([2,tB,tC,tD]+0.5,[y2h,200,400,ymax],...
    {'A','B','C','D'})
legend('x(t)','y(t)')

% ---------------ʩ�ȷ������ڷ�����̿-------------------
syms x(t) z(t)
lambda = log(2)/5;
mu = 2*log(2)/6;
eq2 = [diff(x) == -lambda*x,diff(z) == lambda*x-mu*z];  % ΢�ַ���
cond2 = [x(0) == 1100,z(2) == 236.5588];  % ��ֵ����
[x(t),z(t)] = dsolve(eq2,cond2);
z = vpa(z,5)

% ������ӻ�
figure
fplot(x,[0,25],'b--')  % θ������ҩ���仯����
hold on
fplot(y,[0,25],'r-.')  % ѪҺ��ҩ���仯���ߣ���ʩ�ȣ�
fplot(z,[2,25],'k')  % ѪҺ��ҩ���仯���ߣ�ʩ�ȣ�
grid on
set(gca,'XMinorGrid','on','YMinorGrid','on')
xlabel('ʱ�� t��h��')
ylabel('ҩ����mg��')
legend('x(t)','y(t)','z(t)')

t4 = double(solve(diff(z) == 0))  % ʩ�������£�Ѫҩ���ﵽ��ֵʱ��
zmax = double(subs(z,t4))         % ���Ѫҩ��


% ---------------ʩ�ȷ���������ѪҺ͸��-------------------
syms x(t) z(t)                                         % ������ź���
lambda = log(2)/5;                                     % ����lambda
mu = 6*log(2)/6;                                       % ����mu
eq2 = [diff(x) == -lambda*x,diff(z) == lambda*x-mu*z]; % ����΢�ַ���
cond2 = [x(0) == 1100,z(2) == 236.5588];               % �����ֵ����
[x(t),z(t)] = dsolve(eq2,cond2);                       % ���΢�ַ���
z = vpa(z,5)                                           % ��5λ��Ч������ʽ��ʾ

% ������ӻ�
figure                                                 % ����figure����
fplot(x,[0,25],'b--')                                  % θ������ҩ���仯����
hold on                                              % ����ͼ�α���
fplot(y,[0,25],'r-.')                                % ѪҺ��ҩ���仯���ߣ���ʩ�ȣ�
fplot(z,[2,25],'k')                                  % ѪҺ��ҩ���仯���ߣ�ʩ�ȣ�
grid on                                              % �����������
set(gca,'XMinorGrid','on','YMinorGrid','on')         % ��Ӵ�������
xlabel('ʱ�� t��h��');  ylabel('ҩ����mg��');           % ����������ǩ
legend('x(t)','y(t)','z(t)')                         % ���ͼ��
