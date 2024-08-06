function [pa,qa,pb,qb] = pdebc(xa,ua,xb,ub,t)
% 边值条件描述函数
pa = [0;ua(2)];
qa = [1;0];
pb = [ub(1)-1;0];
qb = [0;1];
end