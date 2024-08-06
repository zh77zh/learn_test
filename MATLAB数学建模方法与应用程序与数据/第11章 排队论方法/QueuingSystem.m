function [Ls,Lq,Ws,Wq,P] = QueuingSystem(lambda,mu,model,VarT)
% �Ŷ�ģ�ͣ�M/M/c/N/m����⺯����M��ʾ��ָ���ֲ���c��ʾ����̨�ĸ����� 
% N��ʾ�Ŷ�ϵͳ���������ƣ�m��ʾ�˿�Դ��Ŀ��
%
%   [Ls,Lq,Ws,Wq,P] = QueuingSystem(lambda,mu,model)
%      lambda�� ƽ��������
%      mu�� ƽ��������
%      model���Ŷ�ģ�ͣ�����Ԫ�ص�����[c,N,m]
%      VarT������ʱ��ķ���˲���������M/G/1/inf/infģ��
%      Ls��ƽ���ӳ�
%      Lq��ƽ���ȴ��ӳ�
%      Ws��ƽ������ʱ��
%      Wq��ƽ���ȴ�ʱ��
%      P����̬���ʣ�P = [p0,p1,p2,����]
% Example��
%      [Ls,Lq,Ws,Wq,P] = QueuingSystem(0.394,1.89,[1,inf,inf])
%      [Ls,Lq,Ws,Wq,P] = QueuingSystem(5.5/48,1/8,[1,inf,inf],16)

if nargin == 2
    model = [1,inf,inf];
end

c = model(1);  % ����̨����
N = model(2);  % ϵͳ����
m = model(3);  % �˿�Դ��Ŀ
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
% ���M/M/C/inf/inf�Ŷ�ģ��
rho = lambda/(c*mu);
if rho >= 1
    warning('����ǿ�ȴ���1���������޳���')
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
% ���M/M/C/N/inf�Ŷ�ģ��
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
% ���M/M/C/m/m�Ŷ�ģ��
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
% ���M/G/1/inf/inf�Ŷ�ģ��
% G��ʾ����ʱ�����һ��ֲ�
rho = lambda/mu;
P = 1-rho;
Lq = (rho^2 + lambda^2*VarT)/(2*(1-rho));
Ls = Lq + rho;
Ws = Ls/lambda;
Wq = Lq/lambda;