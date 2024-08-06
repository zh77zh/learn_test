function f = ObjFun_BaseStation(x,R)
%------------------------------------
%   手机基站优化的目标函数
%------------------------------------
%   x：基站圆心坐标向量
%   R：基站圆半径

C = reshape(x,2,[]);
N = numel(R);
f = 0;

% 通过循环计算所有圆相交的面积总和
for i = 1:N-1
    for j = i+1:N
        % 计算两圆相交面积
        aij = AreaIntersect(C(:,i),R(i),C(:,j),R(j));
        if isequal(C(:,i),C(:,j))
            % 如果圆心重合，加上一个大的惩罚项
            f = f + aij + 1e8;
        else
            % 若圆心不重合，加上一个惩罚项
            f = f + aij + 1/40*1/sum((C(:,i)-C(:,j)).^2);
        end
    end
end
end

function s = AreaIntersect(C1,r1,C2,r2)
%------------------------------------
%   计算两圆相交面积的子函数
%------------------------------------
% C1和C2为两圆的圆心坐标（列向量），r1和r2为两圆半径

d = norm(C1-C2);  % 计算圆心距离
if d >= r1+r2
    % 两圆相离
    s = 0;
elseif d >= 0 && d <= abs(r1-r2)
    % 两圆相套
    s = pi*min(r1,r2)^2;
else
    % 两圆相交
    s = r1^2*MiddleTerm(r1,r2,d) + r2^2*MiddleTerm(r2,r1,d);
end
end

function f = MiddleTerm(ri,rj,dij)
% ----计算中间项的子函数----
f = (ri^2+dij^2-rj^2)/(2*ri*dij);
f = acos(f) - f*sqrt(1-f^2);
end