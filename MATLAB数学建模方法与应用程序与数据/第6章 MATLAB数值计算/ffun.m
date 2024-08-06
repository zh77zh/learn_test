function f = ffun(location,state)
% 例6.4-9的f函数
% location  --- 包含如下字段的结构体变量
%    location.x --- x坐标
%    location.y --- y坐标
%    location.z --- z坐标
% state  --- 包含如下字段的结构体变量
%    state.u    --- 状态变量u
%    state.ux   --- u对x的一阶偏导
%    state.uy   --- u对y的一阶偏导
%    state.uz   --- u对z的一阶偏导
%    state.time --- 时间变量

n = numel(location.x);   % 网格点个数
u1 = state.u(1);         % u1
u2 = state.u(2);         % u2
f = [exp(u2-u1)-exp(u1-u2);exp(u1-u2)-exp(u2-u1)].*ones(1,n);
end