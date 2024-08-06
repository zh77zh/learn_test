%--------------------------------------------------------------------------
%  第1章   MATLAB数组运算
%--------------------------------------------------------------------------
% CopyRight：xiezhh

%% examp1.2-1
x = 1
y = 1+2+sqrt(9)
z = 'Hellow World !!!'

(7189+(1021-913)*80)/64^0.5

pi
pi = 1
clear pi
pi

pi = 1;
clear = 2;
clear pi

iskeyword('for')
iskeyword('xiezhh')

whos

%% examp1.3-1
x = [1  -1.65  2.2  -3.1];
y1 = abs(x)
y2 = sin(x)
y3 = round(x)
y4 = floor(x)
y5 = ceil(x)
y6 = min(x)
y7 = mean(x)
y8 = range(x)
y9 = sign(x)

%% examp1.4-1
x = [1,0,2,-3   5]

%% examp1.4-2
y = [5;2;0]

%% examp1.4-3
x = 1:10 
y = 1:2:10

%% examp1.4-4
x = linspace(1, 10, 10)

%% examp1.4-5
X = []

%% examp1.4-6
A = [1, 2, 3;4  5  6;7  8, 9]

%% examp1.4-7
x = A(:)'

%% examp1.4-8
y = 1:18 ;                % 定义长度为18的向量y
B = reshape(y, [3, 6])

%% examp1.4-9
x1 = 1:3;      % 定义一个向量x1
x2 = 4:6;      % 定义一个向量x2
C = [x1; x2]

%% examp1.4-10
A = [1, 2;3, 4];      % 定义2*2的矩阵A
B = repmat(A,[2,3])

%% examp1.4-11
x = ['abc'; 'def'; 'ghi'] 
size(x)

%% examp1.4-12
x = 2i+5
y = [1  2  3; 4  5  6]*i+7
a = [1  2; 3  4];
b = [5  6; 7  8];
c = complex(a,b)

%% examp1.4-13
syms a b c d
x = [a  b; c  d]
y = [1  2  3; 4  5  6];  
y = sym(y)
z = sym('a%d%d',[2,3])

%% examp1.4-14
A = zeros(3)
B = ones(3,5)
C = eye(3,5)
D = diag([1 2 3])
E = diag(D)
F = rand(3)
G = magic(3)

%% examp1.4-15
% 定义一个2行，2列，2页的3维数组
x(1:2, 1:2, 1)=[1  2; 3  4];
x(1:2, 1:2, 2)=[5  6; 7  8];

%% examp1.4-16
A1 = [1  2; 3  4];
A2 = [5  6; 7  8];
A = cat(3, A1, A2)

%% examp1.4-17
x = reshape(1:12, [2, 2, 3]) 

%% examp1.4-18
x = repmat([1  2; 3  4], [1 1 2])

%% examp1.4-19
x = [1  2  3; 4  5  6; 7  8  9];
y1 = x(1, 2)
y2 = x(2:3, 1:2)
y3 = x(1, :)
y4 = x(:, 1:2)
x(2,:) = []

%% examp1.4-20
x = [1  2  3; 4  5  6; 7  8  9];  % 定义一个3行3列的矩阵x
y5 = x(3:6)

%% examp1.4-21
A = rand(3)
x = A(A>0.5)

%% examp1.4-22
c1 = {[1  2; 3  4], 'xiezhh', 10; [5  6  7], ['abc';'def'], 'I LOVE MATLAB'}

%% examp1.4-23
c2 = cell(2,4)
c2{2, 3} = [1  2  3]

%% examp1.4-24
c = {[1  2], 'xie', 'xiezhh'; 'MATLAB', [3  4; 5  6], 'I LOVE MATLAB'}
c(2, 2)
c{2, 2}

c = {[1  2],  'xiezhh'; 'MATLAB', [3  4; 5  6]};
celldisp(c)

%% examp1.4-25
% 通过直接赋值方式定义一个1行2列的结构体数组
struct1(1).name = 'xiezhh';
struct1(2).name = 'heping';
struct1(1).age = 31;
struct1(2).age = 22;
struct1 

%% examp1.4-26
struct2 = struct('name', {'xiezhh', 'heping'}, 'age',{31, 22})
struct2(1).name

%% examp1.4-27
Name = {'Smith';'Johnson';'Williams';'Jones';'Brown'}; 
Age = [38;43;38;40;49]; 
Height = [71;69;64;67;64];
Weight = [176;163;131;133;119];
BP = [124 93; 109 77; 125 83; 117 75; 122 80];
D = dataset({Age,'Age'},{Height,'Height'},{Weight,'Weight'},...
    {BP,'BloodPressure'},'ObsNames',Name)
x = D(1,:)
y = double(x)
H = D.Height

%% examp1.4-28
Name = {'Smith';'Johnson';'Williams';'Jones';'Brown'};
Age = [38;43;38;40;49];
Height = [71;69;64;67;64];
Weight = [176;163;131;133;119]; 
BloodPressure = [124 93; 109 77; 125 83; 117 75; 122 80]; 
T = table(Age,Height,Weight,BloodPressure,... 
    'RowNames',Name)
T1 = T(4,:)
T2 = T(:,2:3)
H = T.Height

%% examp1.4-29
A1 = rand(60,50);
B1 = mat2cell(A1, [10 20 30], [25 25])
C1 = cell2mat(B1);
isequal(A1,C1)
A2 = [1 2 3 4;5 6 7 8;9 10 11 12];
B2 = num2cell(A2)
C = {'Heping', 'Tianjin', 22;  'Xiezhh', 'Xingyang', 31}
fields = {'Name', 'Address', 'Age'};
S = cell2struct(C, fields, 2)
CS = struct2cell(S)
isequal(C,CS')
x = [1;2;3;4;5];
x = cellstr(num2str(x));
y = strcat('xiezhh', x, '.txt')
TS = struct2table(S)

%% examp1.5-1
A = [1  2; 3  4];
B = [5  6; 7  8];
C = A+B
D = A-B

%% examp1.5-2
A = [1  2  3; 4  5  6];
B = [1  1  1  1; 2  2  2  2; 3  3  3  3];
C = A*B 
D = [1  1  1; 2  2  2];
E = A.*D 

%% examp1.5-3
A = [2  3  8; 1  -2  -4; -5  3  1];
b = [-5; 3; 2];
x = A\b
B = A;
C = A./B

%% examp1.5-4
A = [1  2; 3  4];
B = A ^ 2
C = A .^ 2
D = A .^ A

%% examp1.5-5
A = [1  2; 3  4];
B = [2  2; 2  2];
C1 = A > B 
C2 = A ~= B
C3 = A >=2

%% examp1.5-6
A = [0  0  1  2];
B = [0  -2  0  1];
C1 = A | B
C2 = A & B
C3 = ~ A
C4 = xor(A, B)
x = 5;
y = 0;
x || y
x && y

%% examp1.5-7
A = [1  2  3; 4  5  6; 7  8  9]
B = A'

%% examp1.5-8
A = [1  2  3; 4  5  6; 7  8  9];
B1 = flipud(A)
B2 = fliplr(A)
B3 = rot90(A)

%% examp1.5-9
A = [1  2; 3  4];
d1 = det(A)
syms a b c d
B = [a  b; c  d];
d2 = det(B)

%% examp1.5-10
A = [1  2; 3  4];
Ai = inv(A)
syms a b c d
B = [a  b; c  d];
Bi = inv(B) 
C = [1  2  3; 4  5  6];
Cpi = pinv(C)
D = C * Cpi * C

%% examp1.5-11
A = [5  0  4; 3  1  6; 0  2  3];
d = eig(A)
[V, D] = eig(A)
[Vs, Ds] = eig(sym(A))

%% examp1.5-12
A = [1  2  3; 4  5  6; 7  8  9];
t = trace(A)
r = rank(A)