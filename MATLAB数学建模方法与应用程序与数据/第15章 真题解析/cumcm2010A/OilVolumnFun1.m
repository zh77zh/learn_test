function V = OilVolumnFun1(Alp,h,a,b,B,C)
% 第一问：有变位时油量与倾斜角度及油位高度的关系函数，针对小椭圆型储油罐（两端平头的椭圆柱体）
%   Alp ---- 纵向倾斜角度（单位为°）
%   h   ---- 油位高度（单位：m）
%   a   ---- 油罐截面椭圆长半轴（单位：m）
%   b   ---- 油罐截面椭圆短半轴（单位：m）
%   B   ---- 探针位置参数B（单位：m）
%   C   ---- 探针位置参数C（单位：m）
%   V   ---- 油罐内油量（单位：m?）

tcot = cotd(Alp);
ttan = tand(Alp);

V = zeros(size(h));
for i = 1:numel(h)
    if h(i) >=0 && h(i) <= C*ttan
        V(i) = 2*a*integral2(@(x,y)sqrt(1-y.^2/b^2),-B,h(i)*tcot,-b,@(x)h(i)-x*ttan-b);
    elseif h(i) > C*ttan && h(i) <= 2*b-B*ttan
        V(i) = 2*a*integral2(@(x,y)sqrt(1-y.^2/b^2),-B,C,-b,@(x)(h(i)-x*ttan)-b);
    elseif h(i) > 2*b-B*ttan && h(i) <= 2*b
        V(i) = pi*a*b*(B-(2*b-h(i))*tcot) + 2*a*integral2(@(x,y)sqrt(1-y.^2/b^2),(h(i)-2*b)*tcot,C,-b,@(x)h(i)-x*ttan-b);
    end
end
V = V(:);