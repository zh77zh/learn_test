function Example6_5_3_IcedWatermelon
% �����Ϸű�����೤ʱ����ܶ�͸��
% ���������Ϲ����е��¶���ɢ���⣨�ȴ������̣�
%       1   �6�8T     1  �6�8        �6�8T
%      ---.---- = ---.---( r^2.---)
%       a   �6�8t    r^2 �6�8r       �6�8r
%
%   ��ʼ������
%            t = 0, 0 �� r �� R, T = T0
%   ��ֵ������
%            t > 0, r = 0, �6�8T/�6�8r = 0
%            t > 0, r = R, h(T - Tinf) + k.�6�8T/�6�8r = 0
%
%   1. �������ϵ�������һ���뾶ΪR������������
%   2. ��������û�й�Ƥ��û�й��ѣ���ȿҲ����ȫ���ȵģ����ܶȡ�����������������
%      ȫ�Ǿ��ȵģ�һ�����Բ����������¶ȱ仯���仯
%   3. ������һ������ΪTinf ��Ļ���
%
%   August 4  2017, editted by ZhongHua-Xie, Tianjin University of Science and Technology.

% ��ز���
h = 5;        % �����뾲ֹ������Ķ�������ϵ�� [w/(m^2.K)]
k = 0.48;     % ���ϵĵ���ϵ�� [w/(m.K)]
p = 918;      % ���ϵ��ܶ� [kg/m^3]
Cp = 3990;    % ���� [J/(kg.K)]
a = k/(p*Cp); % ���ϵ�����ɢϵ�� [m^2/s]
T0 = 32;      % ��ʼ�¶� [��]
Tinf = 6;     % �ռ��¶� [��]��������������¶�

m = 2;                                     % �����е�m����
r = linspace(0,9,41)/100;                  % �����������r
t = linspace(0,10,41)*3600;                % ����ʱ������t
T = pdepe(m,@pdefun,@pdeic,@pdebc,r,t);    % �������

% ������ӻ�
figure;
surf(r*100,t/3600,T);                     % �����¶�T���ھ���r��ʱ��t����ά����
colorbar;
zlim([min(T(:))-2,max(T(:))+2]);
xlabel('���� r ��cm��');
ylabel('ʱ�� t ��h��')
zlabel('�¶� T ���棩'); 
view(116,33);
text(8,-0.5,34.5,'��ʼ���','Rotation',40)
text(0,5,26,'����')
text(7,2,10.5,'��Ƥ')

% �������ϵ�Ч�����¶���ɢ����
num = size(T,1);
theta = linspace(0,2*pi,61);
X = 100*r'*cos(theta);
Y = 100*r'*sin(theta);
Ti = T(1,:)'*ones(size(theta));
xmin = min(X(:))-1;
xmax = max(X(:))+1;
ymin = min(Y(:))-1;
ymax = max(Y(:))+1;
X_bac = linspace(xmin,xmax,50);
Y_bac = linspace(ymin,ymax,50);
[X_bac,Y_bac] = meshgrid(X_bac,Y_bac);
figure;
% ���Ʊ����¶���
surf(X_bac,Y_bac,Tinf*ones(size(X_bac)),...
    'FaceColor','interp','EdgeColor','none');
hold on;
% ������ʱ�����ϵ�Ч�¶�����
hs = surf(X,Y,Ti,'FaceColor','interp','EdgeColor','none');
axis equal;
axis([xmin,xmax,ymin,ymax])
xlabel('X');
ylabel('Y');
caxis([Tinf,max(T(:))]);
colorbar;
view(2);
ht = title(['��� ',num2str(t(1)/3600,'%2.1f'),' Сʱ']);

filename = '���������¶���ɢ����.gif';
f = getframe(gcf);
IM = frame2im(f);
[IM,map] = rgb2ind(IM,256);
imwrite(IM,map,filename,'gif', 'Loopcount',inf,'DelayTime',0.2);

% ͨ��ѭ�����Ƹ�ʱ�����ϵ�Ч�¶�����
for i = 2:num
    Ti = T(i,:)'*ones(size(theta));
    set(hs,'ZData',Ti);
    set(ht,'String',['��� ',num2str(t(i)/3600,'%2.1f'),' Сʱ']);
    pause(0.2);
    
    f = getframe(gcf);
    IM = frame2im(f);
    [IM,map] = rgb2ind(IM,256);
    imwrite(IM,map,filename,'gif','WriteMode','append','DelayTime',0.2);   
end

    function [c,f,s] = pdefun(r,t,T,dT)
        % ƫ΢�ַ��̺���
        c = 1/a;
        f = dT;
        s = 0;
    end

    function u0 = pdeic(r)
        % ��ֵ��������
        u0 = T0;
    end

    function [pa,qa,pb,qb] = pdebc(ra,Ta,rb,Tb,t)
        % ��ֵ��������
        pa = 0;
        qa = 1;
        pb = h*(Tb - Tinf);
        qb = k;
    end

end
