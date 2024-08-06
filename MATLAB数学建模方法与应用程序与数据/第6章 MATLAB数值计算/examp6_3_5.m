function examp6_3_5
    function v = ddex3hist(t)
        %��ʷ����
        v = [ log(t); 1./t];
    end
    function d = ddex3delay(t,y)
        %�ӳٺ���
        d = exp(1 - y(2));
    end
    function dydt = ddex3de(t,y,Z)
        %�ӳ�΢�ַ��̺���.����ֻ��һ���ӳ�����ZΪ1�е�������y2(exp(1-y2(t)))�ӳ�
        %��exp(1-y2(t))����y2���ǵڶ���״̬���������y2(exp(1-y2(t)))��Z(2)����ʾ
        dydt = [ y(2); -Z(2)*y(2)^2*exp(1 - y(2))];
    end
t0 = 0.1;
tfinal = 5;
tspan = [t0, tfinal];%���ʱ�䷶Χ
sol = ddesd(@ddex3de,@ddex3delay,@ddex3hist,tspan);
%׼ȷ��
texact = linspace(t0,tfinal);
yexact = ddex3hist(texact);
%���»�ͼ���ֽ��
figure;
plot(sol.x,sol.y(1,:),'o','markersize',7);
hold on
plot(sol.x,sol.y(2,:),'*','markersize',7);
plot(texact,yexact(1,:),'k-','linewidth',2);
plot(texact,yexact(2,:),'k:','linewidth',2);
%ͼ��,λ���Զ�ѡ�����λ��
L = legend('{\ity}_1,ddesd','{\ity}_2,ddesd','{\ity}_1,������',...
    '{\ity}_2,������','Location','best');
set(L,'fontname','Times New Roman');
hold off
xlabel('\fontname{����}ʱ��t','fontsize',16);
ylabel('\fontname{����}y�Ľ�','fontsize',16);
title('ddesd���ͽ�����Ա�ͼ');
end