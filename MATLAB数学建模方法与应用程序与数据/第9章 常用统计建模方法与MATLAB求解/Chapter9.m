%--------------------------------------------------------------------------
%  ��9��  ����ͳ�ƽ�ģ������MATLAB���
%--------------------------------------------------------------------------
% CopyRight��xiezhh


%% examp9.1-1
T = readtable('���ɼ�.xls','ReadRowNames',false);
T.Properties.VariableNames = {'Class','StudentId','Height',...
    'Weight','Rank','VC','Score1','Score2','Score3'};
whichstats = {'mean','std','min','max'};
T1 = T(:,{'Class','Height'});
statarray = grpstats(T1,'Class',whichstats)
T2 = T(:,{'Class','Weight'});
statarray = grpstats(T2,'Class',whichstats)
T3 = T(:,{'Class','VC'});
statarray = grpstats(T3,'Class',whichstats)

T4 = T(:,{'Height','Weight','VC','Score1','Score2','Score3'});
T4 = table2array(T4);
id = any(isnan(T4),2);
T4(id,:) = [];
corrcoef(T4)

T5 = T.Rank;
tabulate(T5);

% ͳ��ͼ
% ����ͼ
VC = T.VC;                            % ��ȡ�λ�������
group = T.Class;                      % ��ȡ�༶����
figure;                               % �½�ͼ�δ���
boxplot(VC,group)                     % ��������ͼ
ylabel('�λ���')                       % y���ǩ

% Ƶ��ֱ��ͼ����ܶ�����
figure;                     % �½�ͼ�δ���
[f, xc] = ecdf(VC);         % ����ecdf��������xc���ľ���ֲ�����ֵf
ecdfhist(f, xc);            % ����Ƶ��ֱ��ͼ
xlabel('�λ���');           % ΪX��ӱ�ǩ
ylabel('f(x)');             % ΪY��ӱ�ǩ
hold on
[f2,xc2] = ksdensity(VC);   % ���ܶȹ���
plot(xc2,f2,'r')            % ���ƺ��ܶ�����
legend('Ƶ��ֱ��ͼ','���ܶ�����','Location','NorthEast');  % ���ͼ��

% ��̬����ͼ
figure;          % �½�ͼ�δ���
normplot(VC);    % ������̬����ͼ

% ����ɢ��ͼ
figure;
gscatter(T.Height,T.Weight,T.Class,'br','o*');  % ���Ʒ���ɢ��ͼ
xlabel('���'); 
ylabel('����')

% ����ɢ��ͼ����
figure;
data = [T.Height,T.Weight,T.VC];            % ��ȡ��ߡ����غͷλ�������
Group = T.Class;
Clr = 'br';
Sym = 'o*';
Siz = [6,6];
Leg = 'on';
Dispopt = 'stairs';
VarNames = {'���','����','�λ���'};           % ������Ԫ������
gplotmatrix (data,[],Group,Clr,Sym,Siz,Leg,Dispopt,VarNames);    % ���Ʒ���ɢ��ͼ����


% ��ά��ͼ
RankStr = T.Rank;
tab = tabulate(RankStr);
x = cell2mat(tab(:,3));
explode = [0,0,1,0,0];
labels = tab(:,1);
figure;
pie3(x,explode,labels);

%% examp9.2-1 ���÷ֲ���������
%+++++++++++++++++++++++++����normfit�������+++++++++++++++++++++++++++++++
x = [15.14  14.81  15.11  15.26  15.08  15.17  15.12  14.95  15.05  14.87];
[muhat,sigmahat,muci,sigmaci] = normfit(x,0.1)

%++++++++++++++++++++++++++++����mle�������++++++++++++++++++++++++++++++++
x = [15.14  14.81  15.11  15.26  15.08  15.17  15.12  14.95  15.05  14.87];
[mu_sigma,mu_sigma_ci] = mle(x,'distribution','norm','alpha',0.1)

%% examp9.2-2 �Զ���ֲ���������
x = [0.7917,0.8448,0.9802,0.8481,0.7627
        0.9013,0.9037,0.7399,0.7843,0.8424
        0.9842,0.7134,0.9959,0.6444,0.8362
        0.7651,0.9341,0.6515,0.7956,0.8733];
x = x(:);
PdfFun = @(x,theta) theta*x.^(theta-1).*(x>0 & x<1);
[phat,pci] = mle(x,'pdf',PdfFun,'start',1)

%% examp9.2-3 ���������
rand('seed',1);
randn('seed',1);
x = normrnd(35,5,600,1);
y = evrnd(20,2,400,1);
data = [x;y];
pdffun = @(t,mu1,sig1,mu2,sig2)0.6*normpdf(t,mu1,sig1)+0.4*evpdf(t,mu2,sig2);
[phat,pci] = mle(data,'pdf',pdffun,'start',[10,10,10,10],...
    'lowerbound',[-inf,0,-inf,0],'upperbound',[inf,inf,inf,inf])

%% examp9.3-1 �����׼����֪ʱ�ĵ�����̬�����ֵ�ļ���
%++++++++++++++++++++++++++++++++˫�����++++++++++++++++++++++++++++++++++
x = [97  102  105  112  99  103  102  94  100  95  105  98  102  100  103];
mu0 = 100;
Sigma = 2;
Alpha = 0.05;
[h,p,muci,zval] = ztest(x,mu0,Sigma,Alpha)

%++++++++++++++++++++++++++++++++�������++++++++++++++++++++++++++++++++++
x = [97  102  105  112  99  103  102  94  100  95  105  98  102  100  103];
mu0 = 100;
Sigma = 2;
Alpha = 0.05;
tail = 'right';
[h,p,muci,zval] = ztest(x,mu0,Sigma,Alpha,tail)

%% examp9.3-2 �����׼��δ֪ʱ�ĵ�����̬�����ֵ�ļ���
x = [49.4  50.5  50.7  51.7  49.8  47.9  49.2  51.4  48.9];
mu0 = 50;
Alpha = 0.05;
[h,p,muci,stats] = ttest(x,mu0,Alpha)

%% examp9.3-3 �����׼��δ֪ʱ��������̬�����ֵ�ıȽϼ��飨����������
x = [20.1,  20.0,  19.3,  20.6,  20.2,  19.9,  20.0,  19.9,  19.1,  19.9];
y = [18.6,  19.1,  20.0,  20.0,  20.0,  19.7,  19.9,  19.6,  20.2];
alpha = 0.05;
tail = 'both';
vartype = 'equal';
[h,p,muci,stats] = ttest2(x,y,alpha,tail,vartype)

%% examp9.3-4 �����׼��δ֪ʱ��������̬�����ֵ�ıȽϼ��飨���������
x = [80.3,68.6,72.2,71.5,72.3,70.1,74.6,73.0,58.7,78.6,85.6,78.0];
y = [74.0,71.2,66.3,65.3,66.0,61.6,68.8,72.6,65.7,72.6,77.1,71.5];
Alpha = 0.05;
tail = 'both';
[h,p,muci,stats] = ttest(x,y,Alpha,tail)

%% examp9.3-5 �����ֵδ֪ʱ�ĵ�����̬���巽��ļ���
x = [49.4  50.5  50.7  51.7  49.8  47.9  49.2  51.4  48.9];
var0 = 1.5;
alpha = 0.05;
tail = 'both';
[h,p,varci,stats] = vartest(x,var0,alpha,tail)

%% examp9.3-6 �����ֵδ֪ʱ��������̬���巽��ıȽϼ���
x = [20.1,  20.0,  19.3,  20.6,  20.2,  19.9,  20.0,  19.9,  19.1,  19.9];
y = [18.6,  19.1,  20.0,  20.0,  20.0,  19.7,  19.9,  19.6,  20.2];
alpha = 0.05;
tail = 'both';
[h,p,varci,stats] = vartest2(x,y,alpha,tail)

%% examp9.4-3  �γ̼���
x = xlsread('2012˫ɫ�򿪽�����.xls',1,'I2:I98');
[h,p,stats] = runstest(x,[],'method','approximate') 

%% examp9.4-4  ���ż���1
x = [-ones(69,1);zeros(23,1);ones(108,1)];
p = signtest(x)

%% examp9.4-5  ���ż���2
x = [80.3,68.6,72.2,71.5,72.3,70.1,74.6,73.0,58.7,78.6,85.6,78.0];
y = [74.0,71.2,66.3,65.3,66.0,61.6,68.8,72.6,65.7,72.6,77.1,71.5];
p = signtest(x,y)

%% examp9.4-6  Wilcoxon�����ȼ���
x = [20.21,19.95,20.15,20.07,19.91,19.99,20.08,20.16,...
        19.99,20.16,20.09,19.97,20.05,20.27,19.96,20.06];
[p,h,stats] = signrank(x,20)

%% examp9.4-7  Mann-Whitney�Ⱥͼ���
x = [133,112,102,129,121,161,142,88,115,127,96,125];
y = [71,119,101,83,107,134,92];
[p,h,stats] = ranksum(x,y,'method','approximate')

%% examp9.4-8  �ֲ���������
data = xlsread('���ɼ�.xls');
VC = data(:,4);

[mu1, sigma1] = normfit(VC)
[h,p,stats] = chi2gof(VC,'cdf',{'normcdf', mu1, sigma1})

[h,p,kstat,critval] = lillietest(VC)

parmhat = lognfit(VC);
mu2 = parmhat(1)
sigma2 = parmhat(2)
[h,p,stats] = chi2gof(VC,'cdf',{'logncdf', mu2, sigma2})

CDF = [VC, logncdf(VC, mu2, sigma2)];
[h,p,ksstat,cv] = kstest(VC,CDF)

%% examp9.4-9  ���������
[~,~,rawdata] = xlsread('��������.xls');  % ��ȡԭʼ����
% ��ȡ����״���ͳ������ݣ���������״��������Ϊ��һ�У�������������Ϊ�ڶ���
data = rawdata(2:end,[7,3]);
myCrossTab(data);

%% examp9.5-1  �������
[x,y] = xlsread('�ߵ���ѧ�ɼ�.xls');
score = x(:,5);
college = y(2:end,3);
college_id = x(:,4);
for i = 1:6
    scorei = score(college_id == i);
    [h,p] = lillietest(scorei);
    result(i,:) = p;
end
result

[p,stats] = vartestn(score,college)

[p,table,stats] = anova1(score,college)

[c,m,h,gnames] = multcompare(stats);
VarNames = {'ѧԺ1','ѧԺ2','��������','���ֵ��','��������','pֵ'};
T1 = [gnames(c(:,1:2)),num2cell(c(:,3:end))];
T1 = cell2table(T1,'VariableNames',VarNames)

T2 = [gnames,num2cell(m)];
T2 = cell2table(T2,'VariableNames',{'ѧԺ','ƽ���ɼ�','��ֵ��׼��'})

%% examp9.6-1  һԪ���Իع�
data = xlsread('�������������.xls');
x = data(:,5);
y = data(:,10);
figure;
plot(x, y, 'k.', 'Markersize', 15);
xlabel('����������(x)');
ylabel('����������(y)');
%����x��y���������ϵ������R
R = corrcoef(x, y)

mdl1 = fitlm(x,y)    % ģ�����

figure;
mdl1.plot;
xlabel('����������(x)');
ylabel('����������(y)');
title('');
legend('ԭʼɢ��','�ع�ֱ��','��������');

xnew = [0.035,0.04]';
ynew = mdl1.predict(xnew)

Res = mdl1.Residuals;
Res_Stu = Res.Studentized;
id = find(abs(Res_Stu)>2);  
mdl2 = fitlm(x,y, 'Exclude',id)

figure;
mdl2.plot;
xlabel('����������(x)');
ylabel('����������(y)');
title('');
legend('�޳��쳣���ݺ�ɢ��','�ع�ֱ��','��������');

%% examp9.6-2 һԪ�����Իع�
HeadData = xlsread('��ͯ­�Է������ָ��.xls'); 
x = HeadData(:, 4); 
y = HeadData(:, 9);
figure;
plot(x, y, 'k.');
xlabel('����(x)');
ylabel('ͷΧ(y)');

HeadCirFun = @(beta, x)beta(1)*exp(beta(2)./(x+beta(3)));
beta0 = [53,-0.2604,0.6276];
nlm1 = fitnlm(x,y,HeadCirFun,beta0)

xnew = linspace(0,16,50)';
ynew = nlm1.predict(xnew);
figure;
plot(x, y, 'k.');
hold on;
plot(xnew, ynew, 'linewidth', 3);
xlabel('����(x)');
ylabel('ͷΧ(y)');
legend('ԭʼ����ɢ��','�����Իع�����');

% ����������䴦ͷΧԤ��ֵ��Ԥ��ֵ��95%��������͹۲�ֵ��95%Ԥ������
x0 = 10;
[yp,ypci1] = nlm1.predict(x0,'Prediction','curve')
[~,ypci2] = nlm1.predict(x0,'Prediction','observation')

cftool(x,y)

%% examp9.6-3 ��Ԫ���Ժ͹������Իع�
data = xlsread('���������������.xls');
X = data(:,3:7);
y = data(:,2);
[R,P] = corrcoef([y,X])
VarNames = {'y','x1','x2','x3','x4','x5'};
matrixplot(R,'FigShap','h','FigSize','full','FigStyle' ,'Triu', ...
    'ColorBar','on','XVar', VarNames,'YVar',VarNames);

mmdl1 = fitlm(X,y)

figure;
subplot(1,2,1);
mmdl1.plotResiduals('histogram');
title('(a) �в�ֱ��ͼ');
xlabel('�в�r');ylabel('f(r)');
subplot(1,2,2);
mmdl1.plotResiduals('probability');
title('(b) �в���̬����ͼ');
xlabel('�в�');ylabel('����');

Res3 = mmdl1.Residuals;
Res_Stu3 = Res3.Studentized;
id3 = find(abs(Res_Stu3)>2);

Model = 'poly10101';
mmdl2 = fitlm(X,y,Model,'Exclude',id3)

Model = 'poly22222';
mmdl3 = fitlm(X,y,Model)

figure;
plot(y,'ko');
hold on
plot(mmdl1.predict(X),':');
plot(mmdl2.predict(X),'r-.');
plot(mmdl3.predict(X),'k');
legend('y��ԭʼɢ��','5Ԫ���Իع����',...
    '3Ԫ���Իع����','��ȫ���λع����');
xlabel('y�Ĺ۲����'); 
ylabel('y'); 

mmdl4 = LinearModel.stepwise(X,y, 'poly22222')
yfitted = mmdl4.Fitted;
figure;
plot(y,'ko');
hold on
plot(yfitted,':','linewidth',2);
legend('y��ԭʼɢ��','�𲽻ع����');
xlabel('y�Ĺ۲����');
ylabel('y');    

model = [0 0 0 0 0
         1 0 0 0 0
         0 1 0 0 0
         0 0 0 0 1
         2 0 0 0 0
         1 1 0 0 0
         0 1 1 0 0
         1 0 0 1 0
         0 0 0 2 0
         1 0 0 0 1
         0 1 0 0 1
         0 0 1 0 1
         0 0 0 0 2];
mmdl5 = fitlm(X,y,model)

%% examp9.6-4 ��Ԫ�����Իع�
modelfun = @(b,x)sqrt((x(:,1)-b(1)).^2+(x(:,2)-b(2)).^2+b(3).^2)/(60*b(4))+b(5);
xyt = [500    3300    21    9
       300     200    19   29
       800    1600    14   51
      1400    2200    13   17
      1700     700    11   46
      2300    2800    14   47
      2500    1900    10   14
      2900     900    11   46
      3200    3100    17   57
      3400     100    16   49];
xy = xyt(:,1:2); Minutes = xyt(:,3); Seconds = xyt(:,4);
T = Minutes + Seconds/60;

b0 = [1000 100 1 1 1];
mnlm = fitnlm(xy,T,modelfun,b0) 

%% examp9.7-2 Q�;������
%***************************��ȡ���ݣ������б�׼��***************************
[X,textdata] = xlsread('�ֵ��������˾�����֧��.xls');
obslabel = textdata(4:end,1);
X = zscore(X);

%******************************* �ֲ����� **********************************
y = pdist(X);
Z = linkage(y,'ward');
H = dendrogram(Z,0,'orientation','right','labels',obslabel);
set(H,'LineWidth',2,'Color','k');
xlabel('��׼�����루Ward������')

%******************************* �������� **********************************
eva = evalclusters(X,'linkage','silhouette','KList',[2:6])

%*********************����clusterdata��������һ������************************
id1 = clusterdata(X,'linkage','ward','maxclust',3);
obslabel(id1 == 1)
obslabel(id1 == 2)
obslabel(id1 == 3)

%******************************* K��ֵ���� **********************************
startdata = X(1:3,:);
id2 = kmeans(X,3,'Start',startdata);
obslabel(id2 == 1)
obslabel(id2 == 2)
obslabel(id2 == 3)

%% examp9.7-3 R�;������
%*************************��ȡ���ݣ���תΪ��������***************************
[X,textdata] = xlsread('ȫ����װ��׼.xls');
y = 1 - X(tril(true(size(X)),-1))';

%***********************����linkage��������ϵͳ������************************
Z = linkage(y,'average');

%****************************** ���ƾ�������ͼ *****************************
varlabel = textdata(2:end,1);
H = dendrogram(Z,0,'orientation','right','labels',varlabel);
set(H,'LineWidth',2,'Color','k');
xlabel('������루��ƽ������')


%% examp9.8-1 �б����
% 1. ��ȡ���ݲ��������ݱ�
T = readtable('�����.xls','ReadRowNames',true);
T.Properties.VariableNames = {'x1','x2','y'};
T_train = T(1:15,:);

% 2. ѵ��������
ResponseVarName = 'y';
Mdl = fitcdiscr(T_train,ResponseVarName)

% 3. ��ѵ���õķ���������Ԥ��
label = predict(Mdl,T)

% 4. ���Ʒ���ɢ��ͼ��������
x1 = T.x1;
x2 = T.x2;
x1i = linspace(min(x1),max(x1),60);
x2i = linspace(min(x2),max(x2),60);
[x1Grid,x2Grid] = meshgrid(x1i,x2i);
[~,Score] = predict(Mdl,[x1Grid(:),x2Grid(:)]);
figure;
contour(x1Grid,x2Grid,reshape(Score(:,1),size(x1Grid)),1,'LineColor','k');
hold on
gscatter(x1,x2,T.y,'rbk','*.p');

%% examp9.9-1 ���ɷַ���
% 1. ��ȡ����
[data,textdata] = xlsread('2016��������ƽ������.xls');
cityname = textdata(2:end,1);

% 2. ���ɷַ���
[coeff,score,latent,~,explained] = pca(data)
cumsum(explained)

% 3. ����ǰ�������ɷֵ÷ֵ�ɢ��ͼ
figure;
plot(score(:,1),score(:,2),'o');
hold on
text(score(:,1),score(:,2),cityname)
xlabel('��һ���ɷֵ÷֣�����ɷ֣�')
ylabel('�ڶ����ɷֵ÷֣����ȳɷ֣�')