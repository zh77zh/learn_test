function Problem = MakeObjFun(R,H,w,t,EdgeFun)
% ���㳤����ƽ����Ϻ�Բ���۵�����������Ƽӹ�����
% R�����Ӱ뾶
% L��ƽ��볤
% H�������۵���߶ȣ����ӵĺ�Ȳ��������ڣ�
% w��ľ�����
% t�����ӵĺ��
% EdgeFun�����ӱ�Ե����
% Problem���Ż�����
%    Ŀ�꺯����min ���۳���  &  max ת������
%    Լ��������
%          a + H < L < 5R
%          a + b < L < 5R
%          0 < Theta < pi/2
%          0 < Beta  < pi/2
%          (a+R)/2 < yA1 < R    % yA1�������ľ����˵� A �� y ����
%          zAi < zA1            % zAi�ǵ�i��ľ����˵� A �� z ���꣬

xi = -(R-w):2*w:R-w;    % ľ��������x����
yi = real(EdgeFun(xi)); % �����Եy����
a = yi(1);              % ���������ľ��һ�볤
ymaxi = max(yi);        % yi���ֵ
n = numel(xi);          % ľ����Ŀ
h = H;                  % ���ߣ����������Ӻ�ȣ�
Problem.objective = @ObjFcn;    % ����Ŀ�꺯��
Problem.nonlcon  = @NlinConFcn; % ���÷�����Լ������
Problem.x0 = [H+a,H+a-ymaxi];   % ���ó�ֵ
Problem.Aineq = [-1,0;1,0;-1,1];% ���Բ���ʽԼ��ϵ������
Problem.bineq = -[a+h;-5*R;ymaxi+t];  % ���Բ���ʽԼ����������
Problem.Aeq = [];               % ���Ե�ʽԼ��ϵ������
Problem.beq = [];               % ���Ե�ʽԼ����������
Problem.lb = [max([R,h+a]);0];  % �����½�
Problem.ub = [];                % �����Ͻ�
Problem.solver = 'fmincon';     % �������������
Problem.options = '';           % �����Ż�����

function fval = ObjFcn(Lb)
% Ŀ�꺯����min ���۳���/ת������
% Lb = [L,b]
% L��ƽ��볤
% b���ֽ�λ�ò���

L = Lb(1);
b = Lb(2);
J = 0;  % ת��������ʼֵ
yz0 = [a+sqrt((L-a)^2-h^2),h];  % ת���yz����
D = zeros(1,n);
for k = 1:n
    J = J + sum((yz0-[0,-t/2]).^2)*2*yi(k);  % �����ϸ�ľ����ת������
    % ƽչ״̬��y������һ��������������
    verti2 = [xi(k),L,-t/2;xi(k),yi(k),-t/2]; 
    % �����k������������н����Ҽ��������ֽ�ľ���
    [CosBetai,d] = SubFcn(yi(k),h,L,b);
    SinBetai = sqrt(1-CosBetai^2);  % �н�����
    Ai = repmat([xi(k),yi(k),0],[2,1]);   % ��ת�㣨������������
    Rot = [1,0,0;0,CosBetai,SinBetai;0,-SinBetai,CosBetai]; % ��ת����
    % �����k��������ת�����˵�����
    vert_Rot = (verti2-Ai)*Rot+Ai;
    xyzLeft = mean(vert_Rot);  % ��k��������������
    J = J + sum((yz0-xyzLeft(2:3)).^2)*(L-yi(k)); % y������һ�����ȵ�ת������
    D(k) = d;
    % ���¼���y�Ḻ��һ�����ȵ�ת������
    verti2 = [xi(k),-L,-t/2;xi(k),-yi(k),-t/2]; 
    Ai = repmat([xi(k),-yi(k),0],[2,1]);
    Rot = [1,0,0;0,CosBetai,-SinBetai;0,SinBetai,CosBetai];
    vert_Rot = (verti2-Ai)*Rot+Ai;
    xyzRight = mean(vert_Rot);
    J = J + sum((yz0-xyzRight(2:3)).^2)*(L-yi(k));
end

Li = abs(D-(L-b-yi));   % ��ľ���Ŀ��۳���
Ltotal = sum(Li);       % �����ܳ���
fval = Ltotal/J;        % Ŀ�꺯��ֵ�������ܳ���/ת������
end

%-------------------------------------------------------------
%  �Ӻ���
%-------------------------------------------------------------
function [CosBeta,d] = SubFcn(y,z,L,b)
    % ���ľ����ת�Ƕȵĸ�����ֵ;
    CosTheta = sqrt((L-a)^2-z^2)/(L-a);
    d = sqrt((y-a)^2+(L-a-b)^2-2*(y-a)*(L-a-b)*CosTheta);
    if y == a
        CosBeta = CosTheta;
    else
        CosBeta = -((y-a)^2+d^2-(L-a-b)^2)/(2*(y-a)*d);
    end
end

function [C,Ceq] = NlinConFcn(Lb)
    % ������Լ������
    L = Lb(1);  % ���Ӱ볤
    b = Lb(2);  % �ֽ�λ�ò���
    [CosBeta,di] = SubFcn(EdgeFun(0),h,L,b);  % ��ת�Ƕȵĸ�����ֵ
    Beta = acos(-CosBeta); % ��ת�Ƕ�
    SinTheta = h/(L-a);  % theta�ǵ�����ֵ
    Theta = asin(SinTheta); % theta��
    yA1 = a + (L-a)*cos(Theta);  % �����ľ����˵� A �� y ����(a+R)/2 < yA1 < R+t
    Ceq = di-(L-a-b);  % �����Ե�ʽԼ�������������νṹ��
    C = [-Beta;Beta-pi/2;-Theta;Theta-pi/2;yA1-R;(a+R)/2-yA1]; % �����Բ���ʽԼ��
    zLeftbottom = zeros(1,n);
    % ͨ��ѭ����������ŵ�z����
    for k = 1:n
        [CosBetai,di] = SubFcn(yi(k),h,L,b);
        SinBetai = sqrt(1-CosBetai^2);
        zLeftbottom(k) = (L-yi(k))*SinBetai;
    end
    % �����Բ���ʽԼ��
    C = [C;max(zLeftbottom(2:end-1))-zLeftbottom(1);di'+yi'-L];
end
end