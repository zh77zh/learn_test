function [Ls,Lq,Ws,Wq,P] = QueuingSystem(lambda,mu,model,VarT)
% 排队模型（M/M/c/N/m）求解函数，M表示负指数分布，c表示服务台的个数， 
% N表示排队系统的容量限制，m表示顾客源数目。
%
%   [Ls,Lq,Ws,Wq,P] = QueuingSystem(lambda,mu,model)
%      lambda： 平均到达率
%      mu： 平均服务率
%      model：排队模型，三个元素的向量[c,N,m]
%      VarT：服务时间的方差，此参数仅用于M/G/1/inf/inf模型
%      Ls：平均队长
%      Lq：平均等待队长
%      Ws：平均逗留时间
%      Wq：平均等待时间
%      P：稳态概率，P = [p0,p1,p2,……]
% Example：
%      [Ls,Lq,Ws,Wq,P] = QueuingSystem(0.394,1.89,[1,inf,inf])
%      [Ls,Lq,Ws,Wq,P] = QueuingSystem(5.5/48,1/8,[1,inf,inf],16)

if nargin == 2
    model = [1,inf,inf];
end

c = model(1);  % 服务台个数
N = model(2);  % 系统容量
m = model(3);  % 顾客源数目
if isinf(m)
    if isinf(N)
        if nargin < 4
            [Ls,Lq,Ws,Wq,P] = M_M_C(lambda,mu,c);
        else
            [Ls,Lq,Ws,Wq,P] = M_G_1(lambda,mu,VarT);
        end
    else
        [Ls,Lq,Ws,Wq,P] = M_M_C_N(lambda,mu,c,N);
    end
else
    [Ls,Lq,Ws,Wq,P] = M_M_C_m_m(lambda,mu,c,m);
end
P = P(P >= 0.0001);

function [Ls,Lq,Ws,Wq,P] = M_M_C(lambda,mu,c)
% 求解M/M/C/inf/inf排队模型
rho = lambda/(c*mu);
if rho >= 1
    warning('服务强度大于1，队列无限长！')
    [Ls,Lq,Ws,Wq,P] = deal([]);
    return
end
k = 0:c-1;
p0 = (sum((c*rho).^k./factorial(k))+(c*rho)^c/(1-rho)/factorial(c))^(-1);
pn = zeros(100,1);
for i = 1:100
    if i <= c
        pn(i) = (c*rho)^i*p0/factorial(i);
    else
        pn(i) = c^c*rho^i*p0/factorial(c);
    end
end
Lq = (c*rho)^c*rho*p0/(1-rho)^2/factorial(c);
Ls = Lq+c*rho;
Ws = Ls/lambda;
Wq = Lq/lambda;
P = [p0;pn];

function [Ls,Lq,Ws,Wq,P] = M_M_C_N(lambda,mu,c,N)
% 求解M/M/C/N/inf排队模型
rho = lambda/(c*mu);
k = 0:c-1;
if rho == 1
    p0 = (sum(c.^k./factorial(k)) + c^c*(N-c+1)/factorial(c))^(-1);
else
    p0 = sum((c*rho).^k./factorial(k));
    p0 = (p0 + c^c*(rho^c-rho^(N+1))/(1-rho)/factorial(c))^(-1);
end
pn = zeros(N,1);
for i = 1:N
    if i <= c
        pn(i) = (c*rho)^i*p0/factorial(i);
    else
        pn(i) = c^c*rho^i*p0/factorial(c);
    end
end
if rho == 1
    Lq = N*(N-1)/(2*(N+1));
else
    Lq = (c*rho)^c*rho*p0*(1-rho^(N-c)-(N-c)*(1-rho)*rho^(N-c));
    Lq = Lq/(1-rho)^2/factorial(c);
end
Ls = Lq+c*rho*(1-pn(N));
lambda_e = lambda*(1-pn(N));
Ws = Ls/lambda_e;
Wq = Lq/lambda_e;
P = [p0;pn];

function [Ls,Lq,Ws,Wq,P] = M_M_C_m_m(lambda,mu,c,m)
% 求解M/M/C/m/m排队模型
rho = m*lambda/(c*mu);
k1 = 0:c;
k2 = c+1:m;
p0 = sum((c*rho/m).^k1./factorial(k1)./factorial(m-k1));
p0 = p0+c^c*sum((rho/m).^k2./factorial(m-k2))/factorial(c);
p0 = p0^(-1)/factorial(m);
pn = zeros(m,1);
for i = 1:m
    if i <= c
        pn(i) = (c*rho/m)^i*p0*nchoosek(m,i);
    else
        pn(i) = (c*rho/m)^i*p0*nchoosek(m,i)*factorial(i)*c^(c-i);
        pn(i) = pn(i)/factorial(c);
    end
end
Ls = sum((1:m)'.*pn);
lambda_e = lambda*(m-Ls);
Lq = Ls - lambda_e/mu;
Ws = Ls/lambda_e;
Wq = Lq/lambda_e;
P = [p0;pn];

function [Ls,Lq,Ws,Wq,P] = M_G_1(lambda,mu,VarT)
% 求解M/G/1/inf/inf排队模型
% G表示服务时间服从一般分布
rho = lambda/mu;
P = 1-rho;
Lq = (rho^2 + lambda^2*VarT)/(2*(1-rho));
Ls = Lq + rho;
Ws = Ls/lambda;
Wq = Lq/lambda;