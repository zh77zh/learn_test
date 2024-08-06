function [Lmax,CI,CR,W] = AHP(A)
% 层次分析法
% [Lmax,CI,CR,W] = AHP(A)
%    A：判断矩阵（两两比较矩阵）
%    Lmax：A的最大特征值
%    CI：一致性指标
%    CR：一致性比率指标
%    W：权重向量
%
% Example:
%    A = [1 2 3;1/2 1 4;1/3 1/4 1];
%    [Lmax,CI,CR,W] = AHP(A)

Lmax = []; CI = []; CR = []; W = [];
[m,n] = size(A);
if m ~= n
    warning('判断矩阵应为方阵')
    return
end
id = tril(true(m),-1);
aij = A(id);
B = A';
bij = B(id);
if any(abs(aij.*bij-1)>1e-10)
    warning('判断矩阵应为正互反矩阵')
    return
end
% 内置的随机一致性指标
RI_Vec = [0,0,0.52,0.89,1.12,1.26,1.36,1.41,...
    1.46,1.49,1.52,1.54,1.56,1.58,1.59,1.60];
% 若判断矩阵的阶数小于17，则用内置的随机一致性指标，否则重新计算
if m < 17
    RI = RI_Vec(m);
else
    RI = MyRI(m);  % 调用自编函数计算随机一致性指标
end
[V,L] = eig(A,'vector');  % 求判断矩阵的特征值与特征向量
L = real(L);              % 特征值取实部
[Lmax,id] = max(L);       % 最大特征值
CI = (Lmax-m)/(m-1);      % 计算一致性指标
CR = CI/RI;               % 计算一致性比率指标
W = V(:,id);              % 最大特征值对应的特征向量
W = W/sum(W);             % 特征向量归一化，即权重向量
if CR < 0.1
    disp('通过一致性检验')
else
    disp('未通过一致性检验')
end
