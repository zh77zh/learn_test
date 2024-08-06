function [c,f,s] = pdefun(x,t,u,du)
% 微分方程的描述函数
c = [1;1];
f = [1/80;1/91].*du;
y = u(1) - u(2);
F = exp(y)-exp(-y);
s = [-F; F];
end