function V = OilVolumnFun2(Alp,Bet,h,r,A,B,C)
% 第二问：有变位时油量与变位角度及油位高度的关系函数，针对实际储油罐（两端为球冠的圆柱体）
%   Alp ---- 纵向倾斜角度（单位为°）
%   Bet ---- 横向偏转角度（单位为°）
%   h   ---- 油位计显示的油位高度（单位：m）
%   r   ---- 油罐截面圆半径（单位：m）
%   A   ---- 球冠高度（单位：m）
%   B   ---- 探针位置参数B（单位：m）
%   C   ---- 探针位置参数C（单位：m）
%   V   ---- 油罐内油量（单位：m?）

R = A/2 + r^2/(2*A);
a1 = (R-A)*cosd(Alp) - r*sind(Alp);
b1 = (R-A)*sind(Alp) + r*cosd(Alp);
a2 = (A-R)*cosd(Alp) - r*sind(Alp);
b2 = (A-R)*sind(Alp) + r*cosd(Alp);

V = zeros(size(h));
for i = 1:numel(h)
    H = r+(h(i)-r)*cosd(Bet);
    if H >= 0 && H <= 2*r-B*tand(Alp)
        Vleft = 2*integral2(@(y,x)sqrt(R^2-(x-a1).^2-(y-b1).^2),0,H*cosd(Alp)+B*sind(Alp),@(y)a1-sqrt(R^2-(y-b1).^2),@(y)-y*tand(Alp));
    elseif H > 2*r-B*tand(Alp) && H <= 2*r
        Vleft = pi*A^2*(A/6+r^2/(2*A));
    end
    if H >= 0 && H <= C*tand(Alp)
        Vright = 0;
    elseif H > C*tand(Alp) && H <= 2*r
        Vright = 2*integral2(@(y,x)sqrt(R^2-(x-a2).^2-(y-b2).^2),0,H*cosd(Alp)-C*sind(Alp),@(y)-y*tand(Alp),@(y)a2+sqrt(R^2-(y-b2).^2));
    end
    Vmid = OilVolumnFun1(Alp,H,r,r,B,C);
    V(i) = Vleft + Vmid + Vright;
end
V = V(:);