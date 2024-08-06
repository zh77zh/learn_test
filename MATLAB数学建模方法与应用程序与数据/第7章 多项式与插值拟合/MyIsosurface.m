function MyIsosurface(X,Y,Z,V,value)
% ���ƺ���V = V(x,y,z)�ĵ�ֵ��ͼ
% MyIsosurface(X,Y,Z,V,value) ����X,Y,ZΪ��ά�����꣬VΪ��Ӧ�ĺ���ֵ��X,Y,Z,V
%                             ����ͬ��ģ�����飬value����ָ����ֵ���Ӧ����ֵ.
% CopyRight��xiezhh��л�л���  2012.2.15

cdata = smooth3(rand(size(V)),'box',5);  % ��ά����ƽ��
p = patch(isosurface(X,Y,Z,V,value));   % ���Ƶ�ֵ��
isonormals(X,Y,Z,V,p);  % �����ֵ�涥��ķ���
isocolors(X,Y,Z,cdata,p);  % �����ֵ����ɫ
% ��������ɫ��ʽΪ��ֵ��ɫ�����ñ��ߵ���ɫΪ��ɫ
set(p,'FaceColor','interp','EdgeColor','none');
view(3);   % ��ά�ӽ�
axis equal ;  % ������������ʾ������ͬ
% axis off;  % ����ʾ������
xlabel('x');ylabel('y');zlabel('z');
camlight; lighting phong;  % ���ù��պ͹���ģʽ 