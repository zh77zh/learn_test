function V = OilVolumnFun1(Alp,h,a,b,B,C)
% ��һ�ʣ��б�λʱ��������б�Ƕȼ���λ�߶ȵĹ�ϵ���������С��Բ�ʹ��͹ޣ�����ƽͷ����Բ���壩
%   Alp ---- ������б�Ƕȣ���λΪ�㣩
%   h   ---- ��λ�߶ȣ���λ��m��
%   a   ---- �͹޽�����Բ�����ᣨ��λ��m��
%   b   ---- �͹޽�����Բ�̰��ᣨ��λ��m��
%   B   ---- ̽��λ�ò���B����λ��m��
%   C   ---- ̽��λ�ò���C����λ��m��
%   V   ---- �͹�����������λ��m?��

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