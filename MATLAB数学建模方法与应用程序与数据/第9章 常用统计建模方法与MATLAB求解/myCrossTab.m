function myCrossTab(data,varargin)
%   r*c列联表检验函数，调用该函数可在屏幕上显示检验结果。
%
%   myCrossTab(T)  列联表检验，T为数值矩阵（频数表）或两列的元胞矩阵（原始数据）。
%         当T为原始数据时，第一列为行变量，第二列为列变量。
%
%   myCrossTab(T,Name1,Value1, ...)  
%         用成对出现的参数名/参数值控制检验的各项属性。可用的参数名/参数值如下：
%         'rowlabel' ---- 设定行变量标签，其参数值为字符串元胞数组，默认为{'r1','r2',...}
%         'collabel' ---- 设定列变量标签，其参数值为字符串元胞数组，默认为{'c1','c2',...}
%
%   Example:
%   T = [200,100;
%        300,100;
%        150,90];
%   rowlabel = {'A','B','C'};
%   collabel = {'及时','不及时'};
%   myCrossTab(T,'rowlabel',rowlabel,'collabel',collabel)
%
%   CopyRight：xiezhh（谢中华）,2016.11.27编写

if ~ismatrix(data) || numel(data) == 1
    warning('第一个输入变量只能是数值矩阵（频数表）或元胞矩阵（原始数据）');
    return;
elseif isnumeric(data)
    table = data;
elseif iscell(data)
    if size(data,2) ~= 2
        warning('原始数据应是两列的元胞矩阵');
        return;
    else
        [table] = crosstab(data(:,1),data(:,2));
        [~,rowlabel] = grp2idx(data(:,1));
        [~,collabel] = grp2idx(data(:,2));
    end
else
    warning('第一个输入变量只能是数值矩阵（频数表）或元胞矩阵（原始数据）');
    return;
end

if mod(numel(varargin),2) ~= 0
    warning('输入参数个数不对，控制参数/参数值应成对出现');
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
    warning('变量标签数与频数表不匹配');
    return;
end
fprintf('\n');
fprintf('-----------------------------------列联表检验结果-----------------------------------');
fprintf('\n');
fprintf('%4s ','');
fprintf(repmat('%10s ',[1,sz(2)]),collabel{:});
fprintf('%10s','合计');
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
fprintf('%s ','合计');
fprintf(repmat('%10d ',[1,sz(2)+1]),colsum);
fprintf('\n');
fprintf('检验结果：卡方 = %10.3f,  自由度DF = %d,  P值 = %5.4f',chi2s,df,p);
fprintf('\n');
end
