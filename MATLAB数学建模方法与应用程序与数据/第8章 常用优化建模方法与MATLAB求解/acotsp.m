function [Shortest_Route,Shortest_Length,R_best,L_best,L_ave] = ...
          acotsp(C,NC_max,m,Alpha,Beta,Rho,Q)
%=========================================================================
%  ����Ⱥ�㷨�������������
%  л�л������Ƽ���ѧ
%  Email:xiezhh@tust.edu.cn
%  ��Ȩ���У�л�л���xiezhh��
%-------------------------------------------------------------------------
%  ����˵��
%  C               n�����е����꣬n��2�ľ���
%  NC_max          ����������
%  m               ���ϸ���
%  Alpha           ������Ϣ����Ҫ�̶ȵĲ���
%  Beta            ��������ʽ������Ҫ�̶ȵĲ���
%  Rho             ��Ϣ������ϵ��
%  Q               ��Ϣ������ǿ��ϵ��
%  Shortest_Route  ���·��
%  Shortest_Length ���·������
%  R_best          �������·��
%  L_best          �������·���ĳ���
%  L_ave           ����·����ƽ������
%=========================================================================
% Example:
% m = 31; Alpha = 1; Beta = 5; Rho = 0.1; NC_max = 200; Q = 100;
% C = [116.38  39.92
%      ...];
% [Shortest_Route,Shortest_Length,R_best,L_best,L_ave]...
%    = acotsp(C,NC_max,m,Alpha,Beta,Rho,Q)
%=========================================================================

% ��һ����������ʼ��
n = size(C,1); % n��ʾ����Ĺ�ģ�����и�����
D = pdist2(C,C); % �����м�������
D = D + eps*eye(n); % ���ں������������Ҫȡ��������D�ĶԽ���Ԫ����eps��������Ծ��ȣ���ʾ

Eta = 1./D; % EtaΪ�������ӣ�������Ϊ����ĵ���
Tau = ones(n); % TauΪ��Ϣ�ؾ���
Tabu = zeros(m,n); % �洢����¼·��������
NC = 1; % ����������
R_best = zeros(NC_max,n); % �������·��
L_best = inf*ones(NC_max,1) ; % �������·�ߵĳ���
L_ave = zeros(NC_max,1) ; % ����·�ߵ�ƽ������

while NC <= NC_max  % ֹͣ����֮һ���ﵽ����������
    % �ڶ�������mֻ��������ŵ�n��������
    if m <= n
        Tabu(:,1) = randperm(n,m)';
    else
        Tabu(:,1) = [randperm(n)';randi(n,m-n,1)];
    end
    L = zeros(m,1);  % �������ϵ�·�����룬��ʼ����Ϊ0��m*1��������
    
    % ��������mֻ���ϰ����ʺ���ѡ����һ�����У���ɸ��Ե�����
    for i = 1:m
        Visiting = Tabu(i,1);  % ��iֻ�������ڷ��ʵĳ���
        UnVisited = 1:n;  % ��iֻ���ϴ����ʵĳ���
        UnVisited(UnVisited == Visiting) = [];
        for j = 2:n  % ���ڳ��в�����
            % ���������ʳ��е�ѡ����ʷֲ�
            P = (Tau(Visiting,UnVisited).^Alpha).*(Eta(Visiting,UnVisited).^Beta);
            P = P/sum(P);
            % ������ԭ��ѡȡ��һ������
            if length(UnVisited) > 1
                NextVisit = randsample(UnVisited,1,true,P);
            else
                NextVisit = UnVisited;
            end
            Tabu(i,j) = NextVisit;
            % �����iֻ�����߹��ľ��룬ԭ������ϵ�ǰ���е���һ�����еľ���
            L(i) = L(i) + D(Visiting,NextVisit);
            Visiting = NextVisit;
            UnVisited(UnVisited == Visiting) = [];
        end
        % ��iֻ���ϻص������㣬һ���������߹��ľ���
        L(i) = L(i) + D(Visiting,Tabu(i,1));
    end
       
    % ���Ĳ�����¼���ε������·��
    [L_best(NC),pos] = min(L);  % ������̾���
    R_best(NC,:) = Tabu(pos,:);  % ���ֵ���������·��
    L_ave(NC) = mean(L);  % ���ֵ������ƽ������
       
    % ���岽��������Ϣ��
    Delta_Tau = zeros(n);  %��ʼʱ��Ϣ������Ϊn*n��0����
    for i = 1:m
        for j = 1:(n-1)
            % ��iֻ���ϴӵ�j�����е���j+1�����е���Ϣ������
            Delta_Tau(Tabu(i,j),Tabu(i,j+1)) = ...
                Delta_Tau(Tabu(i,j),Tabu(i,j+1))+Q/L(i);
        end
        % ��iֻ���ϴӵ�n�����е���1�����У������أ�����Ϣ������
        Delta_Tau(Tabu(i,n),Tabu(i,1)) = Delta_Tau(Tabu(i,n),Tabu(i,1))+Q/L(i);
    end
    Tau = (1-Rho)*Tau + Delta_Tau;  % ������Ϣ�ػӷ������º����Ϣ��
   
    % �����������ɱ�����
    Tabu = zeros(m,n);  % ֱ������������    
    NC = NC + 1;  % ����������1
end

% ���߲���������
[Shortest_Length,pos] = min(L_best);  % ������������̾���
Shortest_Route = R_best(pos,:); % �������������·��
figure;
subplot(1,2,1)  % ���Ƶ�һ����ͼ
BestRoute = Shortest_Route([1:end,1]);
plot(C(BestRoute,1),C(BestRoute,2),'-o')
hold on
text(C(BestRoute(1),1),C(BestRoute(1),2),'������')
title('TSP���������');
subplot(1,2,2)  % ���Ƶڶ�����ͼ
plot(L_best,'r')
hold on  % ͼ�α���
plot(L_ave,'b--')
title('����ƽ���������̾���')   % ����
legend('������̾���','����ƽ������') % ͼ��
end