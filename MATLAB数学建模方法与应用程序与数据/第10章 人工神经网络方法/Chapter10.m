%--------------------------------------------------------------------------
%  ��10��  �����緽��
%--------------------------------------------------------------------------
% CopyRight��xiezhh

%% examp10.5-1 BP�����������
HeadData = xlsread('��ͯ­�Է������ָ��.xls'); 
x = HeadData(:, 4)'; 
y = HeadData(:, 9)';
rng(0)
net = fitnet(3);
trainedNet = train(net,x,y);
view(trainedNet)
xnew = linspace(0,18,50);
ynew = trainedNet(xnew);
figure;
plot(x,y,'.',xnew,ynew,'k')
xlabel('����(x)');
ylabel('ͷΧ(y)');
trainedNet.IW{1}
trainedNet.LW{2,1}
trainedNet.b

%% examp10.6-1 SOM�������
% 1. ��ȡ����
[data,TextData] = xlsread('2016��������ƽ������.xls','A2:M32');
ObsLabel = TextData(:,1);
data = data';
% 2. ����SOM������о���
net = selforgmap([3,1]);
trainedNet = train(net,data);
view(trainedNet)
figure;
plotsomtop(trainedNet)
y = trainedNet(data)
% 3. �鿴������
classid = vec2ind(y);
ObsLabel(classid == 1)  % �鿴��һ���а����ĳ���
ObsLabel(classid == 2)  % �鿴�ڶ����а����ĳ���
ObsLabel(classid == 3)  % �鿴�������а����ĳ���

%% examp10.7-1 BP����ģʽʶ��
[data1,textdata1] = xlsread('��Ԫ����ʶ��.xlsx','��¼A');
[data2,textdata2] = xlsread('��Ԫ����ʶ��.xlsx','��¼B');
[data3,textdata3] = xlsread('��Ԫ����ʶ��.xlsx','��¼C');
trainData = data1(:,3:end)';
n1 = size(trainData,2);
trainGroup = textdata1(2:end,2);
[Gid,Gname] = grp2idx(trainGroup);
Gid = full(ind2vec(Gid'));
net = patternnet(41);
net.divideParam.trainRatio = 85/100;
net.divideParam.valRatio = 15/100;
net.divideParam.testRatio = 0/100;
sampleData = data2(:,3:end)';
n2 = size(sampleData,2);
testData = data3(:,3:end)';
n3 = size(testData,2);

m = 20;
trainResult = zeros(7,n1,m);
sampleResult = zeros(7,n2,m);
testResult = zeros(7,n3,m);
for i = 1:m
    trainedNet = train(net,trainData,Gid);
    trainResult(:,:,i) = trainedNet(trainData);
    sampleResult(:,:,i) = trainedNet(sampleData);
    testResult(:,:,i) = trainedNet(testData);
end
trainResult = mean(trainResult,3);
sampleResult = mean(sampleResult,3);
testResult = mean(testResult,3);
figure;
plotconfusion(Gid,trainResult)
testGroup = Gname(vec2ind(testResult))
sampleGroup = Gname(vec2ind(sampleResult))

%% 10.8-1 ���������
% 1.��ξ���
data = xlsread('����.xlsx',1);
x = data(:,2:12);
d = pdist(x);
z = linkage(d,'ward');
dendrogram(z,0,'orientation','top');
set(gca,'XTickLabelRotation',-90);
result1 = clusterdata(x,'linkage','ward','maxclust',2);

% 2.���������
x2 = x';
net = selforgmap(2);
rng(1);
trainedNet = train(net,x2);
y = trainedNet(x2);
result2 = vec2ind(y);
id1 = find(result2(1:48) == 2)
id2 = find(result2(49:end) == 1)

% 3.ͳ���б�
T = readtable('����.xlsx','PreserveVariableNames',1);
T_train = T(:,[2:12,15]);
ResponseVarName = '��ʵ����';
Mdl = fitcdiscr(T_train,ResponseVarName);
result3 = Mdl.predict(T_train)'
x3 = [3,3,2,2,2,3,2,2,2,3,2; 1,1,2,1,1,2,2,2,1,1,2];
label1 = Mdl.predict(x3)

% 4.�������б�
data = xlsread('����.xlsx',1);
x = data(:,2:12)';
y = kron([1,0;0,1],ones(1,48));
xnew = [3,3,2,2,2,3,2,2,2,3,2; 1,1,2,1,1,2,2,2,1,1,2]';
net = patternnet(23);                % ����ģʽʶ�����磬1�����㣬����10���ڵ�
% ����ѵ�������и����֣�ѵ����������֤�����ԣ���ռ����
net.divideParam.trainRatio = 85/100;
net.divideParam.valRatio = 15/100;
net.divideParam.testRatio = 0/100;
rng(1);
trainedNet = train(net,x,y);
result4 = trainedNet(x);
figure;
plotconfusion(y,result4);
label2 = trainedNet(xnew)
