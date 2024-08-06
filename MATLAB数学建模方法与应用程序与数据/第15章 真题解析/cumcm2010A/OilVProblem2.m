%*************************************************
%    �ڶ��ʣ�ʵ�ʴ��͹ޣ�����Ϊ��ڵ�Բ���壩
%*************************************************
%% ���������б�ǶȺͺ���ƫת�Ƕ�
A = 1;
B = 2;
C = 6;
r = 1.5;
data = xlsread('����A����2��ʵ�ʲɼ����ݱ�.xls',1);
h = data(:,5)/1000;
u = data(1:end,4)/1000;
h_top = h(1:301);
u_top = u(2:301);
Vfun = @(Alp,Bet)OilVolumnFun2(Alp,Bet,h_top,r,A,B,C);
DV = @(Alp,Bet)-diff(Vfun(Alp,Bet));
ObjFun = @(Angle)sum((DV(Angle(1),Angle(2))-u_top).^2);
[Alp_Bet,fval] = fminunc(ObjFun,[2,2])

%% �������Ч��ͼ
hnew = linspace(0,3,60);
Vnew = 1000*OilVolumnFun2(Alp_Bet(1),Alp_Bet(2),hnew,r,A,B,C);
V = data(:,6);
[hs,id] = sort(h);
figure;
plot(hs,V(id),'k',hnew,Vnew,'k--')
xlabel('��λ�߶�h��m��');
ylabel('������V��L��');
legend('��ʾ����','��������','Location','NorthWest')

%% ģ�������ļ���
h_bottom = h(303:end);
u_bottom = 1000*u(304:end);
Alp = Alp_Bet(1);
Bet = Alp_Bet(2);
Vhat = 1000*OilVolumnFun2(Alp,Bet,h_bottom,r,A,B,C);
uhat = -diff(Vhat);
R = corr(u_bottom,uhat)
RMSE = sqrt(sum((u_bottom - uhat).^2)/numel(uhat))
figure
plot(u_bottom,uhat,'+')
hline = refline([1,0]);
set(hline,'Color','r')
xlabel('ʵ�ʳ�������L��');
ylabel('���۳�������L��');
title(['���ϵ�� R = ',num2str(R)])

%% ģ�Ͳ����������ȷ���
u_all = 1000*u([2:301,304:end]);
RmseFun = @(x,y)sqrt(sum((x-y).^2)/numel(x));
Alpi = linspace(0,5,30);
Betj = [2,4,6];
RMSE = zeros(30,3);
for i = 1:30
    for j = 1:3
        % �������ǰ��������
        V_top = 1000*OilVolumnFun2(Alpi(i),Betj(j),h_top,r,A,B,C);
        % ������ͺ���������
        V_bottom = 1000*OilVolumnFun2(Alpi(i),Betj(j),h_bottom,r,A,B,C);
        uij = -[diff(V_top);diff(V_bottom)];  % ���۳�����
        RMSE(i,j) = RmseFun(u_all,uij);       % ���������
    end
end
figure
plot(Alpi,RMSE(:,1),'r',Alpi,RMSE(:,2),'g--',Alpi,RMSE(:,3),'b:')
xlabel('������б��\alpha');
ylabel('��������L��');
legend('\beta = 2��','\beta = 4��','\beta = 6��','Location','NorthWest')

Beti = linspace(1,8,30);
Alpj = [2.0,2.1,2.2];
RMSE = zeros(30,3);
for i = 1:30
    for j = 1:3
        % �������ǰ��������
        V_top = 1000*OilVolumnFun2(Alpj(j),Beti(i),h_top,r,A,B,C);
        % ������ͺ���������
        V_bottom = 1000*OilVolumnFun2(Alpj(j),Beti(i),h_bottom,r,A,B,C);
        uij = -[diff(V_top);diff(V_bottom)];  % ���۳�����
        RMSE(i,j) = RmseFun(u_all,uij);       % ���������
    end
end
figure
plot(Beti,RMSE(:,1),'r',Beti,RMSE(:,2),'g--',Beti,RMSE(:,3),'b:')
xlabel('����ƫת��\beta');
ylabel('��������L��');
legend('\alpha = 2.0��','\alpha = 2.1��','\alpha = 2.2��','Location','NorthWest')