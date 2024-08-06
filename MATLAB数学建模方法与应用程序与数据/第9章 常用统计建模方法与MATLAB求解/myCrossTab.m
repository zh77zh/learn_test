function myCrossTab(data,varargin)
%   r*c��������麯�������øú���������Ļ����ʾ��������
%
%   myCrossTab(T)  ��������飬TΪ��ֵ����Ƶ���������е�Ԫ������ԭʼ���ݣ���
%         ��TΪԭʼ����ʱ����һ��Ϊ�б������ڶ���Ϊ�б�����
%
%   myCrossTab(T,Name1,Value1, ...)  
%         �óɶԳ��ֵĲ�����/����ֵ���Ƽ���ĸ������ԡ����õĲ�����/����ֵ���£�
%         'rowlabel' ---- �趨�б�����ǩ�������ֵΪ�ַ���Ԫ�����飬Ĭ��Ϊ{'r1','r2',...}
%         'collabel' ---- �趨�б�����ǩ�������ֵΪ�ַ���Ԫ�����飬Ĭ��Ϊ{'c1','c2',...}
%
%   Example:
%   T = [200,100;
%        300,100;
%        150,90];
%   rowlabel = {'A','B','C'};
%   collabel = {'��ʱ','����ʱ'};
%   myCrossTab(T,'rowlabel',rowlabel,'collabel',collabel)
%
%   CopyRight��xiezhh��л�л���,2016.11.27��д

if ~ismatrix(data) || numel(data) == 1
    warning('��һ���������ֻ������ֵ����Ƶ������Ԫ������ԭʼ���ݣ�');
    return;
elseif isnumeric(data)
    table = data;
elseif iscell(data)
    if size(data,2) ~= 2
        warning('ԭʼ����Ӧ�����е�Ԫ������');
        return;
    else
        [table] = crosstab(data(:,1),data(:,2));
        [~,rowlabel] = grp2idx(data(:,1));
        [~,collabel] = grp2idx(data(:,2));
    end
else
    warning('��һ���������ֻ������ֵ����Ƶ������Ԫ������ԭʼ���ݣ�');
    return;
end

if mod(numel(varargin),2) ~= 0
    warning('��������������ԣ����Ʋ���/����ֵӦ�ɶԳ���');
    return;
end

n = sum(table(:));
sz = size(table);
colsum = [sum(table),n];
expected = sum(table,2)*sum(table)/n;
chi2 = (table - expected).^ 2 ./ expected;
chi2s = sum(chi2(:));
df = (sz(1)-1)*(sz(2)-1);
p = 1-chi2cdf(chi2s,df);  

if isnumeric(data)
    rowlabel_dflts = strcat({'r'},num2str((1:sz(1))'));
    collabel_dflts = strcat({'c'},num2str((1:sz(2))'));
    names = {'rowlabel','collabel'};
    dflts = {rowlabel_dflts,collabel_dflts};
    [rowlabel,collabel] = ...
        internal.stats.parseArgsUser(names, dflts, 'userargs', varargin{:});
end
if (numel(rowlabel) ~= sz(1)) || (numel(collabel) ~= sz(2))
    warning('������ǩ����Ƶ����ƥ��');
    return;
end
fprintf('\n');
fprintf('-----------------------------------�����������-----------------------------------');
fprintf('\n');
fprintf('%4s ','');
fprintf(repmat('%10s ',[1,sz(2)]),collabel{:});
fprintf('%10s','�ϼ�');
fprintf('\n');
for j = 1:sz(1)
    fprintf('%4s ',rowlabel{j});
    fprintf([repmat('%10.2f ',[1,sz(2)]),'%10d'],[table(j,:),sum(table(j,:))]);
    fprintf('\n');
    fprintf('%5s ','');
    fprintf(repmat('%10.2f ',[1,sz(2)]),expected(j,:));
    fprintf('\n');
    fprintf('%5s ','');
    fprintf(repmat('%10.3f ',[1,sz(2)]),chi2(j,:));
    fprintf('\n');
    fprintf('\n');
end
fprintf('%s ','�ϼ�');
fprintf(repmat('%10d ',[1,sz(2)+1]),colsum);
fprintf('\n');
fprintf('������������ = %10.3f,  ���ɶ�DF = %d,  Pֵ = %5.4f',chi2s,df,p);
fprintf('\n');
end
