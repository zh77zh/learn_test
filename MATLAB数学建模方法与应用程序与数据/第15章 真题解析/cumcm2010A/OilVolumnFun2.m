function V = OilVolumnFun2(Alp,Bet,h,r,A,B,C)
% �ڶ��ʣ��б�λʱ�������λ�Ƕȼ���λ�߶ȵĹ�ϵ���������ʵ�ʴ��͹ޣ�����Ϊ��ڵ�Բ���壩
%   Alp ---- ������б�Ƕȣ���λΪ�㣩
%   Bet ---- ����ƫת�Ƕȣ���λΪ�㣩
%   h   ---- ��λ����ʾ����λ�߶ȣ���λ��m��
%   r   ---- �͹޽���Բ�뾶����λ��m��
%   A   ---- ��ڸ߶ȣ���λ��m��
%   B   ---- ̽��λ�ò���B����λ��m��
%   C   ---- ̽��λ�ò���C����λ��m��
%   V   ---- �͹�����������λ��m?��

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