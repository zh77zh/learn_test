%--------------------------------------------------------------------------
%  第2章   MATLAB程序设计
%--------------------------------------------------------------------------
% CopyRight：xiezhh


%% examp2.1-1
A = input('请输入三角形的三条边：');
if A(1) + A(2) > A(3) & A(1) + A(3) > A(2) & A(2) + A(3) > A(1)
    p = (A(1) + A(2) + A(3)) / 2;
    s = sqrt(p*(p - A(1))*(p - A(2))*(p - A(3)));
    disp(['该三角形面积为：' num2str(s)]);
else
    disp('不能构成一个三角形。')
end

%% examp2.1-2
num = input('请输入一个数：');
switch num
    case -1
        disp('I am a teacher.');
    case 0
        disp('I am a student.');
    case 1
        disp('You are a teacher.');
    otherwise
        disp('You are a student.');
end

%% examp2.1-3
% 程序1：for循环
y = 0;
for i = 1:inf 
    y = y + i^2;
    if  y > 2000
        break;
    end
end
n = i - 1 
y = y - i^2

% 程序2：while循环
y = 0;
i = 0;
while  y <= 2000
    i = i + 1;
    y = y + i^2;
end
n = i-1 
y = y-i^2

%% examp2.1-4
A = [1,2,3;4,5,6]; B = [7,8,9;10,11,12];
try
    X = A*B
catch ME
    disp(ME.message);    % 显示出错原因
end

%% examp2.2-2
fun1 = @(x,y) cos(x).*sin(y)
x = [0,1,2];
y = [-1,0,1];
z = fun1(x,y)

%% examp2.2-3
fun2 = @(x)(x>=-1 & x<0).*sin(pi*x.^2)+(x>=0).*exp(1-x);
fun2([-0.5,0,0.5]) 