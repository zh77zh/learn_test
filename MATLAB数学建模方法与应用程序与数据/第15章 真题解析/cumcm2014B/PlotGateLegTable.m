function Result = PlotGateLegTable(R,L,H,w,t,EdgeFun,b)
% �����۵����۵��Ķ���Ч��ͼ���������ľ�����۳���
% ����ԭ�������ӵ��±���
% R�����Ӱ뾶
% L��ƽ��볤
% H������չ����߶ȣ����ӵĺ�Ȳ��������ڣ�
% w��ľ�����
% t�����ӵĺ��
% EdgeFun�����ӱ�Ե����
% b���ֽ�λ�ò����������ľ����˵㵽�ֽ�ľ��룩
% Result��������㣬�����յ㣬���۳���
xi = -(R-w):2*w:R-w;    % ľ��������x����
yi = real(EdgeFun(xi)); % �����Եy����
a = yi(1);              % �����������ľ����һ�볤
if nargin < 7
    b = (L-a)/2;        % Ĭ�ϵ�b����
end
% ��û�������������۵����۵����̶�̬Ч��ͼ
if ~nargout
    figure('name','�۵����۵����̶�̬Ч��ͼ','numbertitle','off');
end

n = numel(xi);          % ľ������
fac = [1 2 3 4;2 6 7 3;4 3 7 8;1 5 8 4;1 2 6 5;5 6 7 8];% ����ľ��������Ķ�����
verti1 = cell(1,n);  % �����ϸ�ľ���Ķ�������
verti2 = verti1;     % ������ȵĶ�������
verti3 = verti1;     % �Ҳ����ȵĶ�������
h1 = zeros(1,n);
h2 = h1;
for i = 1:n
    vx = xi(i) + w*[-1 1 1 -1 -1 1 1 -1]'; % ����x����
    vz = t*[-1 -1 -1 -1 0 0 0 0]';         % ����z����
    vy1 = yi(i)*[1 1 -1 -1 1 1 -1 -1]';    % ����y����
    vy2 = [L,L,yi(i),yi(i),L,L,yi(i),yi(i)]'; % ����y����
    verti1{i} = [vx,vy1,vz];  % �����ϸ�ľ���Ķ�������
    patch('faces',fac,'vertices',verti1{i},'FaceColor',[1,1,0]); % ������ľ��
    verti2{i} = [vx,vy2,vz];  % ������ȵĶ�������
    % ���������
    h1(i) = patch('faces',fac,'vertices',verti2{i},'FaceColor',[0,0,1]);
    verti3{i} = [vx,-vy2,vz]; % �Ҳ����ȵĶ�������
    % ���Ҳ�����
    h2(i) = patch('faces',fac,'vertices',verti3{i},'FaceColor',[0,0,1]);
end
view([-60,24]);  % �����ӵ�λ��
xlabel('x');ylabel('y');zlabel('z');  % �������ǩ
axis equal; axis([-R,R,-L,L,-5,H+5]); % ������������ʾ����
set(gca,'ZDir','reverse','color','none'); % ����z����������
hold on;   % ͼ�α���
% �����ű�Ե��
hline = plot3(xi,L*ones(size(xi)),zeros(size(xi)),'r','linewidth',3);

hei = linspace(0,H,24);  % ������������
hei(hei>(L-a)) = [];     % ȥ������������
m = numel(hei);          % ���߸���
xTop = zeros(m,n);
yTop = xTop;
zTop = xTop;
D = xTop;
if ~nargout, f = getframe(gca); Frame = repmat(f.cdata,[1,1,1,m]);end
% ͨ��ѭ������ľ��λ��
for j = 1:m      % �����߽���ѭ��
    for k = 1:n  % ��ľ������ѭ��
        [CosBetai,di] = SubFcn(yi(k),hei(j));  % ����н����Ҽ�di
        SinBetai = real(sqrt(1-CosBetai^2));   % �н����� 
        Oi = repmat([xi(k),yi(k),0],[8,1]);    % ��������ϵ�����
        Rot = [1,0,0;0,CosBetai,SinBetai;0,-SinBetai,CosBetai]; % ��ת����
        vert_Rot = (verti2{k}-Oi)*Rot+Oi;  % �����������ת�󶥵�����
        set(h1(k),'vertices',vert_Rot);    % ����������ȵĶ�������
        xTop(j,k) = xi(k);  % ���ű�Ե��x����
        yTop(j,k) = yi(k)+(L-yi(k))*CosBetai; % ���ű�Ե��y����
        zTop(j,k) = (L-yi(k))*SinBetai;  % ���ű�Ե��z����
        D(j,k) = di;
        % �������ű�Ե�ߵ���ά����
        set(hline,'XData',xTop(j,:),'YData',yTop(j,:),'ZData',zTop(j,:))
        Oi = repmat([xi(k),-yi(k),0],[8,1]); % �Ҳ������ϵ�����
        Rot = [1,0,0;0,CosBetai,-SinBetai;0,SinBetai,CosBetai]; % ��ת����
        vert_Rot = (verti3{k}-Oi)*Rot+Oi;  % ���Ҳ�������ת�󶥵�����
        set(h2(k),'vertices',vert_Rot);  % �����Ҳ����ȵĶ�������
    end
    if j == 1
        % �������ű�Ե�߻���������
        hsurf = mesh(xTop([j;j],:),yTop([j;j],:),zTop([j;j],:),'FaceAlpha',0);
    else
        % �������ű�Ե�߻���������
        set(hsurf,'XData',xTop(1:j,:),'YData',yTop(1:j,:),'ZData',zTop(1:j,:));
    end
    if ~nargout, f = getframe(gca); Frame(:,:,:,j) = f.cdata;end
    pause(0.3);
    drawnow;
end

% ƽչ״̬�¸�ľ�����������յ��y����
yNotch = [(L-b)*ones(n,1),(D(end,:)+yi)'];
Result1 = [yNotch,abs(diff(yNotch,1,2))];  % ������㣬�����յ㣬���۳���
if nargout
    Result = Result1;  % ���������㣬�����յ㣬���۳��Ƚ��
else
    % ���Ƹ�ľ������λ��ʾ��ͼ
    figure('name','��ľ������λ��ʾ��ͼ','numbertitle','off');
    for i = 1:n
        vert = verti1{i}(5:8,:);
        % ��������ľ��ͶӰ
        patch('faces',1:4,'vertices',vert,'FaceColor',[0.5,0.5,0.5]);
        vert = verti2{i}(5:8,:);
        % �����������ͶӰ
        patch('faces',1:4,'vertices',vert,'FaceColor',[1,1,1]);
        vert = verti3{i}(5:8,:);
        % �����Ҳ�����ͶӰ
        patch('faces',1:4,'vertices',vert,'FaceColor',[1,1,1]);
    end
    view(2);
    hold on;
    plot(xi,yi,'r','linewidth',2);  % ���������Ե��
    plot(xi,-yi,'r','linewidth',2); % ���������Ե��
    plot(xi,[yNotch,-yNotch],'r','linewidth',2); % ���ƿ���λ����
    xlabel('x'); ylabel('y');  % �������ǩ
    axis equal;  axis([-R,R,-L,L]);  % ������������ʾ����
    
    % �����۵�����̬�仯���̵�ʾ��ͼ��ץȡ9����̬ͼƬ��
    figure('name','�۵�����̬�仯���̵�ʾ��ͼ');
    id = round(linspace(1,m,9));
    montage(Frame(:,:,:,id),'Size',[3,3]);
    % ��������д��excel�ļ�
    Result2 = [{'�������','�����յ�','���۳���'};num2cell(Result1)];
    xlswrite('��ľ�����������յ��y����.xlsx',Result2);
end

%-------------------------------------------------------------
%  �Ӻ���
%-------------------------------------------------------------
function [CosBeta,d] = SubFcn(y,z)
    % ���ľ����ת�Ƕȵ�����ֵ;
    CosTheta = sqrt((L-a)^2-z^2)/(L-a);
    d = sqrt((y-a)^2+(L-a-b)^2-2*(y-a)*(L-a-b)*CosTheta);
    if y == a
        CosBeta = CosTheta;
    else
        CosBeta = -((y-a)^2+d^2-(L-a-b)^2)/(2*(y-a)*d);
    end
end
end