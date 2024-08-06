%*************************************************
%    第一问：小椭圆型储油罐（两端平头的椭圆柱体）
%*************************************************
%% 考虑罐体无变位情形
Alp = 0;
a = 0.89;
b = 0.6;
B = 0.4;
C = 2.05;
data = xlsread('问题A附件1：实验采集数据表.xls',1);
h = data(:,4)/1000;  % 油位高度
Vreal = data(:,3)+262;   % 油量（实测值）
Vcal = 1000*OilVolumnFun1(Alp,h,a,b,B,C);  % 油量（理论值）

figure;
plot(h,Vreal);
hold on;
plot(h,Vcal,'r--');
legend('实测值','理论值','Location','NorthWest');
title('无变位油罐油位高度与罐容关系曲线');
xlabel('油位高度h（m）');
ylabel('油量（L）');

%% 考虑罐体纵向倾斜情形
Alp = 4.1;
a = 0.89;
b = 0.6;
B = 0.4;
C = 2.05;
data = xlsread('问题A附件1：实验采集数据表.xls',3);
h = data(:,4)/1000;  % 油位高度
Vreal = data(:,3)+215;   % 油量（实测值）
Vcal = 1000*OilVolumnFun1(Alp,h,a,b,B,C);  % 油量（理论值）

figure;
plot(h,Vreal);
hold on;
plot(h,Vcal,'r--');
legend('实测值','理论值','Location','NorthWest');
title('有变位油罐油位高度与罐容关系曲线');
xlabel('油位高度h（m）');
ylabel('油量（L）');
