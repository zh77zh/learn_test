function RI = MyRI(n)
% 通过1000次随机模拟计算n阶正互反矩阵的随机一致性指标
% 
if n == 1
    RI = 0;
else
    A = eye(n);
    id = tril(true(n),-1);
    population = [1:9,1./(2:9)];
    k = n*(n-1)/2;
    m = 1000;
    L = zeros(1,m);
    for i = 1:m
        aij = randsample(population,k,1);
        A(id) = aij;
        B = A';
        B(id) = 1./aij;
        L(i) = max(real(eig(B)));
    end
    RI = (mean(L)-n)/(n-1);
end