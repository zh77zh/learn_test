function [Shortest_Route,Shortest_Length,R_best,L_best,L_ave] = ...
          acotsp(C,NC_max,m,Alpha,Beta,Rho,Q)
%=========================================================================
%  用蚁群算法求解旅行商问题
%  谢中华，天津科技大学
%  Email:xiezhh@tust.edu.cn
%  版权所有：谢中华（xiezhh）
%-------------------------------------------------------------------------
%  参数说明
%  C               n个城市的坐标，n×2的矩阵
%  NC_max          最大迭代次数
%  m               蚂蚁个数
%  Alpha           表征信息素重要程度的参数
%  Beta            表征启发式因子重要程度的参数
%  Rho             信息素蒸发系数
%  Q               信息素增加强度系数
%  Shortest_Route  最短路径
%  Shortest_Length 最短路径长度
%  R_best          各代最短路径
%  L_best          各代最短路径的长度
%  L_ave           各代路径的平均长度
%=========================================================================
% Example:
% m = 31; Alpha = 1; Beta = 5; Rho = 0.1; NC_max = 200; Q = 100;
% C = [116.38  39.92
%      ...];
% [Shortest_Route,Shortest_Length,R_best,L_best,L_ave]...
%    = acotsp(C,NC_max,m,Alpha,Beta,Rho,Q)
%=========================================================================

% 第一步：变量初始化
n = size(C,1); % n表示问题的规模（城市个数）
D = pdist2(C,C); % 各城市间距离矩阵
D = D + eps*eye(n); % 由于后面的启发因子要取倒数，把D的对角线元素用eps（浮点相对精度）表示

Eta = 1./D; % Eta为启发因子，这里设为距离的倒数
Tau = ones(n); % Tau为信息素矩阵
Tabu = zeros(m,n); % 存储并记录路径的生成
NC = 1; % 迭代计数器
R_best = zeros(NC_max,n); % 各代最佳路线
L_best = inf*ones(NC_max,1) ; % 各代最佳路线的长度
L_ave = zeros(NC_max,1) ; % 各代路线的平均长度

while NC <= NC_max  % 停止条件之一：达到最大迭代次数
    % 第二步：将m只蚂蚁随机放到n个城市上
    if m <= n
        Tabu(:,1) = randperm(n,m)';
    else
        Tabu(:,1) = [randperm(n)';randi(n,m-n,1)];
    end
    L = zeros(m,1);  % 各个蚂蚁的路径距离，开始距离为0，m*1的列向量
    
    % 第三步：m只蚂蚁按概率函数选择下一座城市，完成各自的周游
    for i = 1:m
        Visiting = Tabu(i,1);  % 第i只蚂蚁正在访问的城市
        UnVisited = 1:n;  % 第i只蚂蚁待访问的城市
        UnVisited(UnVisited == Visiting) = [];
        for j = 2:n  % 所在城市不计算
            % 下面计算访问城市的选择概率分布
            P = (Tau(Visiting,UnVisited).^Alpha).*(Eta(Visiting,UnVisited).^Beta);
            P = P/sum(P);
            % 按概率原则选取下一个城市
            if length(UnVisited) > 1
                NextVisit = randsample(UnVisited,1,true,P);
            else
                NextVisit = UnVisited;
            end
            Tabu(i,j) = NextVisit;
            % 计算第i只蚂蚁走过的距离，原距离加上当前城市到下一个城市的距离
            L(i) = L(i) + D(Visiting,NextVisit);
            Visiting = NextVisit;
            UnVisited(UnVisited == Visiting) = [];
        end
        % 第i只蚂蚁回到出发点，一轮下来后走过的距离
        L(i) = L(i) + D(Visiting,Tabu(i,1));
    end
       
    % 第四步：记录本次迭代最短路径
    [L_best(NC),pos] = min(L);  % 计算最短距离
    R_best(NC,:) = Tabu(pos,:);  % 此轮迭代后的最短路径
    L_ave(NC) = mean(L);  % 此轮迭代后的平均距离
       
    % 第五步：更新信息素
    Delta_Tau = zeros(n);  %开始时信息素增量为n*n的0矩阵
    for i = 1:m
        for j = 1:(n-1)
            % 第i只蚂蚁从第j个城市到第j+1个城市的信息素增量
            Delta_Tau(Tabu(i,j),Tabu(i,j+1)) = ...
                Delta_Tau(Tabu(i,j),Tabu(i,j+1))+Q/L(i);
        end
        % 第i只蚂蚁从第n个城市到第1个城市（出发地）的信息素增量
        Delta_Tau(Tabu(i,n),Tabu(i,1)) = Delta_Tau(Tabu(i,n),Tabu(i,1))+Q/L(i);
    end
    Tau = (1-Rho)*Tau + Delta_Tau;  % 考虑信息素挥发，更新后的信息素
   
    % 第六步：禁忌表清零
    Tabu = zeros(m,n);  % 直到最大迭代次数    
    NC = NC + 1;  % 迭代次数加1
end

% 第七步：输出结果
[Shortest_Length,pos] = min(L_best);  % 迭代结束后最短距离
Shortest_Route = R_best(pos,:); % 迭代结束后最短路径
figure;
subplot(1,2,1)  % 绘制第一个子图
BestRoute = Shortest_Route([1:end,1]);
plot(C(BestRoute,1),C(BestRoute,2),'-o')
hold on
text(C(BestRoute(1),1),C(BestRoute(1),2),'出发点')
title('TSP问题求解结果');
subplot(1,2,2)  % 绘制第二个子图
plot(L_best,'r')
hold on  % 图形保持
plot(L_ave,'b--')
title('各代平均距离和最短距离')   % 标题
legend('各代最短距离','各代平均距离') % 图例
end