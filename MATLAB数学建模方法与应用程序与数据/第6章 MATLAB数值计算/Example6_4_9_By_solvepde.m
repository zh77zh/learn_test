%% ����parabolic�������ƫ΢�ַ����飨examp6.4-9��
% 1. ���������������̵�΢�ַ���ģ��
model6 = createpde(2);

% 2. �����������������
% ����3Ϊ���α�ţ�4�Ǳ������������򶥵�x����Ϊ[0,1,1,0],y����Ϊ[0,0,1,1]
R1 = [3,4,0,1,1,0,0,0,1,1]'; 
geom = R1;
ns = char('R1')';             % �����������
sf = 'R1';                    % �����������Ĺ�ʽ
g = decsg(geom,sf,ns);        % ����������
geometryFromEdges(model6,g);  % �����Զ����������
figure;
pdegplot(model6,'EdgeLabels','on');   % ����������򣬲���ʾ�߽��ǩ 
axis equal;
axis([-0.1,1.1,-0.1,1.1]);

% 3. ���ñ�ֵ������Dirichlet ��ֵ������
% u1���ұ߽��Ӧ�ı�ֵ����
applyBoundaryCondition(model6,'Edge',2,'u',1,'EquationIndex',1);
% u2����߽��Ӧ�ı�ֵ����
applyBoundaryCondition(model6,'Edge',4,'u',0,'EquationIndex',2);

% 4. ��������
Me = generateMesh(model6,'Hmax',0.1,'GeometricOrder','linear'); 

% 5. ���÷��̲���
m = 0;
d = 1;
c = [1/80;0;0;1/91;0;0];
a = 0;
specifyCoefficients(model6,'m',m,'d',d,'c',c,'a',a,'f',@ffun);

% 6. ���ó�ʼ����
u0 = [1;0];
setInitialConditions(model6,u0);

% 7. �������
tlist = linspace(0,2,31);  % ����ʱ������
u6 = solvepde(model6,tlist);  % �������

% 8. ������ӻ�
xi = Me.Nodes(1,:);  % �����x����
[xi,id] = sort(xi);
U = u6.NodalSolution;
U1 = squeeze(U(id,1,:))';
U2 = squeeze(U(id,2,:))';
figure;
surf(xi,tlist,U1);
shading interp;
title('u1(x,t)'); xlabel('Distance x'); ylabel('Time t');

figure;
surf(xi,tlist,U2);
shading interp;
title('u2(x,t)'); xlabel('Distance x'); ylabel('Time t');