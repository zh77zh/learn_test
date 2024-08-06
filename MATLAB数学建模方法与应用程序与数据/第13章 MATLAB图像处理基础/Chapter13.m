%--------------------------------------------------------------------------
%  第13章  MATLAB图像处理基础
%--------------------------------------------------------------------------
% CopyRight：xiezhh

%% examp13.2-1
% 读取一幅jpeg格式图像
I = imread('MyHandwriting.jpg');
figure;
imshow(I);  % 显示图像

%% examp13.2-2
figure;
implay('rhinos.avi'); % 播放动画

%% examp13.2-3
bingdundun;
frame = getframe;
I = frame2im(frame);
imwrite(I,'冰墩墩.jpg')

%% examp13.2-4
v = VideoReader('rhinos.avi');
I = readFrame(v);
figure
imshow(I)

%% examp13.3-1
% 图像缩放
I = imread('football.jpg');
L1 = imresize(I,0.5,'bilinear');
L2 = imresize(I,[120,150],'bilinear');
figure
imshowpair(L1,L2,'montage')

%% examp13.3-2
% 图像旋转
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
% 预定义滤波器
h = fspecial('unsharp');
I2 = imfilter(I1,h);
figure
imshowpair(I1,I2,'montage');

%% examp13.4-2
I = imread('硬币.tif');
% 添加椒盐噪声
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
xlabel('（a）加噪图像')
F1 = fftn(J);
F2 = fftshift(F1);
S1 = uint8(20*log(abs(F2)));
subplot(2,2,2)
imshow(S1)
xlabel('（b）振幅谱图像')
[m,n,k] = size(F2);
[x,y] = meshgrid(-1:2/(n-1):1,-1:2/(m-1):1);
r = 0.3;
id = sqrt(x.^2+y.^2) <= r;
id = repmat(id,[1,1,k]);
F2 = F2.*id;
S2 = uint8(20*log(abs(F2)));
subplot(2,2,3)
imshow(S2)
xlabel('（c）去噪后振幅谱图像')
F3 = ifftshift(F2);
K = uint8(real(ifftn(F3)));
subplot(2,2,4)
imshow(K)
xlabel('（d）去噪后图像')

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
I = rgb2gray(I);  % 真彩图形转灰度图像
figure
subplot(2,2,1)
h_im = imshow(I);   % 显示图像
xlabel('（a）原始图像')
h = imellipse;  % 手动选择椭圆（或圆）区域
% h = imrect;  % 矩形区域
% h = imfreehand;  % 自由区域  
wait(h);
id = ~createMask(h,h_im);
I(id) = 0;   % 图像预处理，去除感兴趣区域之外的图像 
subplot(2,2,2)
imshow(I)   % 显示图像
xlabel('（b）感兴趣区域图像')
BW = im2bw(I,0.7);  % 图像二值化
subplot(2,2,3)
imshow(BW);
xlabel('（c）二值化图像')
BW = bwareaopen(BW,5000);  % 按面积去除干扰图像
subplot(2,2,4)
imshow(BW);
xlabel('（d）去除干扰点后的二值化图像')
S = regionprops(BW,'Area');
Area = S(1).Area

%% examp13.5-5
I = imread('rice.png');
BW = edge(I,'Sobel');
figure
imshowpair(I,BW,'montage');

%% examp13.5-6
I = imread('区域分析.png');
BW = imbinarize(I,'adaptive',...
    'ForegroundPolarity','dark',...
    'Sensitivity',0.05);
BW = ~BW;
BW = bwareaopen(BW,50);  % 去除干扰点
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
%*******************************读取图像数据********************************
IM = imread('FittingLine.bmp');
whos  IM

%********************************去除坐标框*********************************
Red = IM(:,:,1);
Rrow = sum(Red,2);
[~,idrow] = mink(Rrow,2)

Rcol = sum(Red);
[~,idcol] = mink(Rcol,2)

% 提取坐标框内部的图像数据
I = IM(min(idrow):max(idrow),min(idcol):max(idcol),:);
m = size(I, 1)
n = size(I, 2)
figure;
imshow(I)

%**************************提取图中曲线上点的像素坐标************************
BluePoints = (I(:,:,1) == 0 & I(:,:,2) == 0 & I(:,:,3) == 255);
[ypixel,xpixel] = find(BluePoints);
size(xpixel) 

%**************************将像素坐标转换为实际坐标**************************
x_xishu = 1/(n-1);
y_xishu = 9/(m-1);
xreal = (xpixel-1)*x_xishu;
yreal = -11-(ypixel-1)*y_xishu;

%**********************************曲线拟合*********************************
fun = @(a,x)[a(1)+a(2)/2*(x-0.17).^2+a(3)/4*(x-0.17).^4];
a = nlinfit(xreal,yreal,fun,[0, 0, 0])

%*******************************绘制重建的曲线******************************
yp = fun(a, xreal);
figure;
plot(xreal,yp);
xlabel('X');
ylabel('Y = f(X)');
text('Interpreter','latex',...
	'String',['$$-19.6749+\frac{22.2118}{2}(x-0.17)^2'...
    '+\frac{5.0905}{4}(x-0.17)^4$$'],'Position',[0.05, -12],...
	'FontSize',12);