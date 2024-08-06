%--------------------------------------------------------------------------
%  第9章  常用统计建模方法与MATLAB求解
%--------------------------------------------------------------------------
% CopyRight：xiezhh


%% examp9.1-1
T = readtable('体测成绩.xls','ReadRowNames',false);
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

% 统计图
% 箱线图
VC = T.VC;                            % 提取肺活量数据
group = T.Class;                      % 提取班级数据
figure;                               % 新建图形窗口
boxplot(VC,group)                     % 绘制箱线图
ylabel('肺活量')                       % y轴标签

% 频率直方图与核密度曲线
figure;                     % 新建图形窗口
[f, xc] = ecdf(VC);         % 调用ecdf函数计算xc处的经验分布函数值f
ecdfhist(f, xc);            % 绘制频率直方图
xlabel('肺活量');           % 为X轴加标签
ylabel('f(x)');             % 为Y轴加标签
hold on
[f2,xc2] = ksdensity(VC);   % 核密度估计
plot(xc2,f2,'r')            % 绘制核密度曲线
legend('频率直方图','核密度曲线','Location','NorthEast');  % 添加图例

% 正态概率图
figure;          % 新建图形窗口
normplot(VC);    % 绘制正态概率图

% 分组散点图
figure;
gscatter(T.Height,T.Weight,T.Class,'br','o*');  % 绘制分组散点图
xlabel('身高'); 
ylabel('体重')

% 分组散点图矩阵
figure;
data = [T.Height,T.Weight,T.VC];            % 提取身高、体重和肺活量数据
Group = T.Class;
Clr = 'br';
Sym = 'o*';
Siz = [6,6];
Leg = 'on';
Dispopt = 'stairs';
VarNames = {'身高','体重','肺活量'};           % 变量名元胞数组
gplotmatrix (data,[],Group,Clr,Sym,Siz,Leg,Dispopt,VarNames);    % 绘制分组散点图矩阵


% 三维饼图
RankStr = T.Rank;
tab = tabulate(RankStr);
x = cell2mat(tab(:,3));
explode = [0,0,1,0,0];
labels = tab(:,1);
figure;
pie3(x,explode,labels);

%% examp9.2-1 常用分布参数估计
%+++++++++++++++++++++++++调用normfit函数求解+++++++++++++++++++++++++++++++
x = [15.14  14.81  15.11  15.26  15.08  15.17  15.12  14.95  15.05  14.87];
[muhat,sigmahat,muci,sigmaci] = normfit(x,0.1)

%++++++++++++++++++++++++++++调用mle函数求解++++++++++++++++++++++++++++++++
x = [15.14  14.81  15.11  15.26  15.08  15.17  15.12  14.95  15.05  14.87];
[mu_sigma,mu_sigma_ci] = mle(x,'distribution','norm','alpha',0.1)

%% examp9.2-2 自定义分布参数估计
x = [0.7917,0.8448,0.9802,0.8481,0.7627
        0.9013,0.9037,0.7399,0.7843,0.8424
        0.9842,0.7134,0.9959,0.6444,0.8362
        0.7651,0.9341,0.6515,0.7956,0.8733];
x = x(:);
PdfFun = @(x,theta) theta*x.^(theta-1).*(x>0 & x<1);
[phat,pci] = mle(x,'pdf',PdfFun,'start',1)

%% examp9.2-3 多参数估计
rand('seed',1);
randn('seed',1);
x = normrnd(35,5,600,1);
y = evrnd(20,2,400,1);
data = [x;y];
pdffun = @(t,mu1,sig1,mu2,sig2)0.6*normpdf(t,mu1,sig1)+0.4*evpdf(t,mu2,sig2);
[phat,pci] = mle(data,'pdf',pdffun,'start',[10,10,10,10],...
    'lowerbound',[-inf,0,-inf,0],'upperbound',[inf,inf,inf,inf])

%% examp9.3-1 总体标准差已知时的单个正态总体均值的检验
%++++++++++++++++++++++++++++++++双侧检验++++++++++++++++++++++++++++++++++
x = [97  102  105  112  99  103  102  94  100  95  105  98  102  100  103];
mu0 = 100;
Sigma = 2;
Alpha = 0.05;
[h,p,muci,zval] = ztest(x,mu0,Sigma,Alpha)

%++++++++++++++++++++++++++++++++单侧检验++++++++++++++++++++++++++++++++++
x = [97  102  105  112  99  103  102  94  100  95  105  98  102  100  103];
mu0 = 100;
Sigma = 2;
Alpha = 0.05;
tail = 'right';
[h,p,muci,zval] = ztest(x,mu0,Sigma,Alpha,tail)

%% examp9.3-2 总体标准差未知时的单个正态总体均值的检验
x = [49.4  50.5  50.7  51.7  49.8  47.9  49.2  51.4  48.9];
mu0 = 50;
Alpha = 0.05;
[h,p,muci,stats] = ttest(x,mu0,Alpha)

%% examp9.3-3 总体标准差未知时的两个正态总体均值的比较检验（独立样本）
x = [20.1,  20.0,  19.3,  20.6,  20.2,  19.9,  20.0,  19.9,  19.1,  19.9];
y = [18.6,  19.1,  20.0,  20.0,  20.0,  19.7,  19.9,  19.6,  20.2];
alpha = 0.05;
tail = 'both';
vartype = 'equal';
[h,p,muci,stats] = ttest2(x,y,alpha,tail,vartype)

%% examp9.3-4 总体标准差未知时的两个正态总体均值的比较检验（配对样本）
x = [80.3,68.6,72.2,71.5,72.3,70.1,74.6,73.0,58.7,78.6,85.6,78.0];
y = [74.0,71.2,66.3,65.3,66.0,61.6,68.8,72.6,65.7,72.6,77.1,71.5];
Alpha = 0.05;
tail = 'both';
[h,p,muci,stats] = ttest(x,y,Alpha,tail)

%% examp9.3-5 总体均值未知时的单个正态总体方差的检验
x = [49.4  50.5  50.7  51.7  49.8  47.9  49.2  51.4  48.9];
var0 = 1.5;
alpha = 0.05;
tail = 'both';
[h,p,varci,stats] = vartest(x,var0,alpha,tail)

%% examp9.3-6 总体均值未知时的两个正态总体方差的比较检验
x = [20.1,  20.0,  19.3,  20.6,  20.2,  19.9,  20.0,  19.9,  19.1,  19.9];
y = [18.6,  19.1,  20.0,  20.0,  20.0,  19.7,  19.9,  19.6,  20.2];
alpha = 0.05;
tail = 'both';
[h,p,varci,stats] = vartest2(x,y,alpha,tail)

%% examp9.4-3  游程检验
x = xlsread('2012双色球开奖数据.xls',1,'I2:I98');
[h,p,stats] = runstest(x,[],'method','approximate') 

%% examp9.4-4  符号检验1
x = [-ones(69,1);zeros(23,1);ones(108,1)];
p = signtest(x)

%% examp9.4-5  符号检验2
x = [80.3,68.6,72.2,71.5,72.3,70.1,74.6,73.0,58.7,78.6,85.6,78.0];
y = [74.0,71.2,66.3,65.3,66.0,61.6,68.8,72.6,65.7,72.6,77.1,71.5];
p = signtest(x,y)

%% examp9.4-6  Wilcoxon符号秩检验
x = [20.21,19.95,20.15,20.07,19.91,19.99,20.08,20.16,...
        19.99,20.16,20.09,19.97,20.05,20.27,19.96,20.06];
[p,h,stats] = signrank(x,20)

%% examp9.4-7  Mann-Whitney秩和检验
x = [133,112,102,129,121,161,142,88,115,127,96,125];
y = [71,119,101,83,107,134,92];
[p,h,stats] = ranksum(x,y,'method','approximate')

%% examp9.4-8  分布拟合与检验
data = xlsread('体测成绩.xls');
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

%% examp9.4-9  列联表检验
[~,~,rawdata] = xlsread('汽车销售.xls');  % 读取原始数据
% 提取婚姻状况和车型数据，并将婚姻状况数据作为第一列，将车型数据作为第二列
data = rawdata(2:end,[7,3]);
myCrossTab(data);

%% examp9.5-1  方差分析
[x,y] = xlsread('高等数学成绩.xls');
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
VarNames = {'学院1','学院2','置信下限','组均值差','置信上限','p值'};
T1 = [gnames(c(:,1:2)),num2cell(c(:,3:end))];
T1 = cell2table(T1,'VariableNames',VarNames)

T2 = [gnames,num2cell(m)];
T2 = cell2table(T2,'VariableNames',{'学院','平均成绩','均值标准误'})

%% examp9.6-1  一元线性回归
data = xlsread('沪深股市收益率.xls');
x = data(:,5);
y = data(:,10);
figure;
plot(x, y, 'k.', 'Markersize', 15);
xlabel('深市收益率(x)');
ylabel('沪市收益率(y)');
%计算x和y的线性相关系数矩阵R
R = corrcoef(x, y)

mdl1 = fitlm(x,y)    % 模型求解

figure;
mdl1.plot;
xlabel('深市收益率(x)');
ylabel('沪市收益率(y)');
title('');
legend('原始散点','回归直线','置信区间');

xnew = [0.035,0.04]';
ynew = mdl1.predict(xnew)

Res = mdl1.Residuals;
Res_Stu = Res.Studentized;
id = find(abs(Res_Stu)>2);  
mdl2 = fitlm(x,y, 'Exclude',id)

figure;
mdl2.plot;
xlabel('深市收益率(x)');
ylabel('沪市收益率(y)');
title('');
legend('剔除异常数据后散点','回归直线','置信区间');

%% examp9.6-2 一元非线性回归
HeadData = xlsread('儿童颅脑发育情况指标.xls'); 
x = HeadData(:, 4); 
y = HeadData(:, 9);
figure;
plot(x, y, 'k.');
xlabel('年龄(x)');
ylabel('头围(y)');

HeadCirFun = @(beta, x)beta(1)*exp(beta(2)./(x+beta(3)));
beta0 = [53,-0.2604,0.6276];
nlm1 = fitnlm(x,y,HeadCirFun,beta0)

xnew = linspace(0,16,50)';
ynew = nlm1.predict(xnew);
figure;
plot(x, y, 'k.');
hold on;
plot(xnew, ynew, 'linewidth', 3);
xlabel('年龄(x)');
ylabel('头围(y)');
legend('原始数据散点','非线性回归曲线');

% 计算给定年龄处头围预测值、预测值的95%置信区间和观测值的95%预测区间
x0 = 10;
[yp,ypci1] = nlm1.predict(x0,'Prediction','curve')
[~,ypci2] = nlm1.predict(x0,'Prediction','observation')

cftool(x,y)

%% examp9.6-3 多元线性和广义线性回归
data = xlsread('人体耗氧能力测试.xls');
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
title('(a) 残差直方图');
xlabel('残差r');ylabel('f(r)');
subplot(1,2,2);
mmdl1.plotResiduals('probability');
title('(b) 残差正态概率图');
xlabel('残差');ylabel('概率');

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
legend('y的原始散点','5元线性回归拟合',...
    '3元线性回归拟合','完全二次回归拟合');
xlabel('y的观测序号'); 
ylabel('y'); 

mmdl4 = LinearModel.stepwise(X,y, 'poly22222')
yfitted = mmdl4.Fitted;
figure;
plot(y,'ko');
hold on
plot(yfitted,':','linewidth',2);
legend('y的原始散点','逐步回归拟合');
xlabel('y的观测序号');
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

%% examp9.6-4 多元非线性回归
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

%% examp9.7-2 Q型聚类分析
%***************************读取数据，并进行标准化***************************
[X,textdata] = xlsread('分地区居民人均消费支出.xls');
obslabel = textdata(4:end,1);
X = zscore(X);

%******************************* 分步聚类 **********************************
y = pdist(X);
Z = linkage(y,'ward');
H = dendrogram(Z,0,'orientation','right','labels',obslabel);
set(H,'LineWidth',2,'Color','k');
xlabel('标准化距离（Ward方法）')

%******************************* 聚类评价 **********************************
eva = evalclusters(X,'linkage','silhouette','KList',[2:6])

%*********************调用clusterdata函数进行一步聚类************************
id1 = clusterdata(X,'linkage','ward','maxclust',3);
obslabel(id1 == 1)
obslabel(id1 == 2)
obslabel(id1 == 3)

%******************************* K均值聚类 **********************************
startdata = X(1:3,:);
id2 = kmeans(X,3,'Start',startdata);
obslabel(id2 == 1)
obslabel(id2 == 2)
obslabel(id2 == 3)

%% examp9.7-3 R型聚类分析
%*************************读取数据，并转为距离向量***************************
[X,textdata] = xlsread('全国服装标准.xls');
y = 1 - X(tril(true(size(X)),-1))';

%***********************调用linkage函数创建系统聚类树************************
Z = linkage(y,'average');

%****************************** 绘制聚类树形图 *****************************
varlabel = textdata(2:end,1);
H = dendrogram(Z,0,'orientation','right','labels',varlabel);
set(H,'LineWidth',2,'Color','k');
xlabel('并类距离（类平均法）')


%% examp9.8-1 判别分析
% 1. 读取数据并创建数据表
T = readtable('蠓虫分类.xls','ReadRowNames',true);
T.Properties.VariableNames = {'x1','x2','y'};
T_train = T(1:15,:);

% 2. 训练分类器
ResponseVarName = 'y';
Mdl = fitcdiscr(T_train,ResponseVarName)

% 3. 用训练好的分类器进行预测
label = predict(Mdl,T)

% 4. 绘制分组散点图及分类线
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

%% examp9.9-1 主成分分析
% 1. 读取数据
[data,textdata] = xlsread('2016各地区月平均气温.xls');
cityname = textdata(2:end,1);

% 2. 主成分分析
[coeff,score,latent,~,explained] = pca(data)
cumsum(explained)

% 3. 绘制前两个主成分得分的散点图
figure;
plot(score(:,1),score(:,2),'o');
hold on
text(score(:,1),score(:,2),cityname)
xlabel('第一主成分得分（寒冷成分）')
ylabel('第二主成分得分（炎热成分）')