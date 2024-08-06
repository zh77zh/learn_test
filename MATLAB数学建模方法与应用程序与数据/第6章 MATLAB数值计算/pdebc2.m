function [pa,qa,pb,qb] = pdebc2(ra,Ta,rb,Tb,t)
% 冰镇西瓜模型的边值条件描述函数
% ra,rb  --- 下边界和上边界条件中的r值
% Ta,Tb  --- 下边界和上边界条件中的T值
% t      --- 时间
% pa,pb  --- 下边界和上边界条件中的p函数值
% qa,qb  --- 下边界和上边界条件中的q函数值

h = 5;        % 西瓜与静止冷空气的对流换热系数 [w/(m^2.K)]
Tinf = 6;     % 终极温度 [℃]，即冰箱冷藏室温度
k = 0.48;     % 西瓜的导热系数 [w/(m.K)]
pa = 0;
qa = 1;
pb = h*(Tb - Tinf);
qb = k;
end