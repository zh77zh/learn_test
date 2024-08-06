function f = framp0(location, state)
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

f = 5*sin(1.7*state.time).*(location.x.^2 + location.y.^2 <= 0.001);
end