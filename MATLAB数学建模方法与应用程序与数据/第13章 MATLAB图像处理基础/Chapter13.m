%--------------------------------------------------------------------------
%  ��13��  MATLABͼ�������
%--------------------------------------------------------------------------
% CopyRight��xiezhh

%% examp13.2-1
% ��ȡһ��jpeg��ʽͼ��
I = imread('MyHandwriting.jpg');
figure;
imshow(I);  % ��ʾͼ��

%% examp13.2-2
figure;
implay('rhinos.avi'); % ���Ŷ���

%% examp13.2-3
bingdundun;
frame = getframe;
I = frame2im(frame);
imwrite(I,'���ն�.jpg')

%% examp13.2-4
v = VideoReader('rhinos.avi');
I = readFrame(v);
figure
imshow(I)

%% examp13.3-1
% ͼ������
I = imread('football.jpg');
L1 = imresize(I,0.5,'bilinear');
L2 = imresize(I,[120,150],'bilinear');
figure
imshowpair(L1,L2,'montage')

%% examp13.3-2
% ͼ����ת
I = imread('football.jpg');
J = imrotate(I,30,'bilinear','crop');
figure
imshowpair(I,J,'montage')

%% examp13.3-3
I = imread('dd.png');
J = imadjust(I);
figure
imshowpair(I,J,'montage')

%% examp13.3-4
I = imread('jj.png');
J = histeq(I);
figure
imshowpair(I,J,'montage')

%% examp13.4-1
I1 = imread('moon.tif');
% Ԥ�����˲���
h = fspecial('unsharp');
I2 = imfilter(I1,h);
figure
imshowpair(I1,I2,'montage');

%% examp13.4-2
I = imread('Ӳ��.tif');
% ��ӽ�������
J = imnoise(I,'salt & pepper',0.02);
K = medfilt2(J);
figure
imshowpair(J,K,'montage');

%% examp13.4-3
I = imread('football.jpg'); 
J = imnoise(I,'salt & pepper',0.02);
figure
subplot(2,2,1)
imshow(J)
xlabel('��a������ͼ��')
F1 = fftn(J);
F2 = fftshift(F1);
S1 = uint8(20*log(abs(F2)));
subplot(2,2,2)
imshow(S1)
xlabel('��b�������ͼ��')
[m,n,k] = size(F2);
[x,y] = meshgrid(-1:2/(n-1):1,-1:2/(m-1):1);
r = 0.3;
id = sqrt(x.^2+y.^2) <= r;
id = repmat(id,[1,1,k]);
F2 = F2.*id;
S2 = uint8(20*log(abs(F2)));
subplot(2,2,3)
imshow(S2)
xlabel('��c��ȥ��������ͼ��')
F3 = ifftshift(F2);
K = uint8(real(ifftn(F3)));
subplot(2,2,4)
imshow(K)
xlabel('��d��ȥ���ͼ��')

%% examp13.5-1
I = imread('ct001.JPG');
BW = im2bw(I,0.3);
figure;
imshowpair(I,BW,'montage');

%% examp13.5-2
I = imread('printedtext.png');
figure;
imshow(I);
BW = imbinarize(I,'adaptive',...
    'ForegroundPolarity','dark',...
    'Sensitivity',0.4);
figure
imshow(BW);

%% examp13.5-3
I = imread('apple.bmp');
figure
imshow(I)
id = (I(:,:,1)-I(:,:,2))>10;
I(repmat(~id,[1,1,3])) = 0;
figure
imshow(I);

%% examp13.5-4
I = imread('0514.bmp');
I = rgb2gray(I);  % ���ͼ��ת�Ҷ�ͼ��
figure
subplot(2,2,1)
h_im = imshow(I);   % ��ʾͼ��
xlabel('��a��ԭʼͼ��')
h = imellipse;  % �ֶ�ѡ����Բ����Բ������
% h = imrect;  % ��������
% h = imfreehand;  % ��������  
wait(h);
id = ~createMask(h,h_im);
I(id) = 0;   % ͼ��Ԥ����ȥ������Ȥ����֮���ͼ�� 
subplot(2,2,2)
imshow(I)   % ��ʾͼ��
xlabel('��b������Ȥ����ͼ��')
BW = im2bw(I,0.7);  % ͼ���ֵ��
subplot(2,2,3)
imshow(BW);
xlabel('��c����ֵ��ͼ��')
BW = bwareaopen(BW,5000);  % �����ȥ������ͼ��
subplot(2,2,4)
imshow(BW);
xlabel('��d��ȥ�����ŵ��Ķ�ֵ��ͼ��')
S = regionprops(BW,'Area');
Area = S(1).Area

%% examp13.5-5
I = imread('rice.png');
BW = edge(I,'Sobel');
figure
imshowpair(I,BW,'montage');

%% examp13.5-6
I = imread('�������.png');
BW = imbinarize(I,'adaptive',...
    'ForegroundPolarity','dark',...
    'Sensitivity',0.05);
BW = ~BW;
BW = bwareaopen(BW,50);  % ȥ�����ŵ�
BW = imfill(BW,'holes');
figure
imshow(BW)
stats = regionprops('table',BW,'EquivDiameter',...
    'Area','centroid')
centers = stats.Centroid;
diameters = stats.EquivDiameter;
r = diameters/2;
figure
imshow(BW); 
hold on;
viscircles(centers,r);

%% examp13.6-1
%*******************************��ȡͼ������********************************
IM = imread('FittingLine.bmp');
whos  IM

%********************************ȥ�������*********************************
Red = IM(:,:,1);
Rrow = sum(Red,2);
[~,idrow] = mink(Rrow,2)

Rcol = sum(Red);
[~,idcol] = mink(Rcol,2)

% ��ȡ������ڲ���ͼ������
I = IM(min(idrow):max(idrow),min(idcol):max(idcol),:);
m = size(I, 1)
n = size(I, 2)
figure;
imshow(I)

%**************************��ȡͼ�������ϵ����������************************
BluePoints = (I(:,:,1) == 0 & I(:,:,2) == 0 & I(:,:,3) == 255);
[ypixel,xpixel] = find(BluePoints);
size(xpixel) 

%**************************����������ת��Ϊʵ������**************************
x_xishu = 1/(n-1);
y_xishu = 9/(m-1);
xreal = (xpixel-1)*x_xishu;
yreal = -11-(ypixel-1)*y_xishu;

%**********************************�������*********************************
fun = @(a,x)[a(1)+a(2)/2*(x-0.17).^2+a(3)/4*(x-0.17).^4];
a = nlinfit(xreal,yreal,fun,[0, 0, 0])

%*******************************�����ؽ�������******************************
yp = fun(a, xreal);
figure;
plot(xreal,yp);
xlabel('X');
ylabel('Y = f(X)');
text('Interpreter','latex',...
	'String',['$$-19.6749+\frac{22.2118}{2}(x-0.17)^2'...
    '+\frac{5.0905}{4}(x-0.17)^4$$'],'Position',[0.05, -12],...
	'FontSize',12);