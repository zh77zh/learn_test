%*************************************************
%    ��һ�ʣ�С��Բ�ʹ��͹ޣ�����ƽͷ����Բ���壩
%*************************************************
%% ���ǹ����ޱ�λ����
Alp = 0;
a = 0.89;
b = 0.6;
B = 0.4;
C = 2.05;
data = xlsread('����A����1��ʵ��ɼ����ݱ�.xls',1);
h = data(:,4)/1000;  % ��λ�߶�
Vreal = data(:,3)+262;   % ������ʵ��ֵ��
Vcal = 1000*OilVolumnFun1(Alp,h,a,b,B,C);  % ����������ֵ��

figure;
plot(h,Vreal);
hold on;
plot(h,Vcal,'r--');
legend('ʵ��ֵ','����ֵ','Location','NorthWest');
title('�ޱ�λ�͹���λ�߶�����ݹ�ϵ����');
xlabel('��λ�߶�h��m��');
ylabel('������L��');

%% ���ǹ���������б����
Alp = 4.1;
a = 0.89;
b = 0.6;
B = 0.4;
C = 2.05;
data = xlsread('����A����1��ʵ��ɼ����ݱ�.xls',3);
h = data(:,4)/1000;  % ��λ�߶�
Vreal = data(:,3)+215;   % ������ʵ��ֵ��
Vcal = 1000*OilVolumnFun1(Alp,h,a,b,B,C);  % ����������ֵ��

figure;
plot(h,Vreal);
hold on;
plot(h,Vcal,'r--');
legend('ʵ��ֵ','����ֵ','Location','NorthWest');
title('�б�λ�͹���λ�߶�����ݹ�ϵ����');
xlabel('��λ�߶�h��m��');
ylabel('������L��');
