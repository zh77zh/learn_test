function [c,f,s] = pdefun(x,t,u,du)
% ΢�ַ��̵���������
c = [1;1];
f = [1/80;1/91].*du;
y = u(1) - u(2);
F = exp(y)-exp(-y);
s = [-F; F];
end