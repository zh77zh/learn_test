function f = framp0(location, state)
% location  --- ���������ֶεĽṹ�����
%    location.x --- x����
%    location.y --- y����
%    location.z --- z����
% state  --- ���������ֶεĽṹ�����
%    state.u    --- ״̬����u
%    state.ux   --- u��x��һ��ƫ��
%    state.uy   --- u��y��һ��ƫ��
%    state.uz   --- u��z��һ��ƫ��
%    state.time --- ʱ�����

f = 5*sin(1.7*state.time).*(location.x.^2 + location.y.^2 <= 0.001);
end