% 从地图上选取经纬度坐标
figure;
axesm('MapProjection','mercator');  % 建立墨卡托投影坐标系
Prince = shaperead('省界线.shp', 'UseGeoCoords', true);  % 读取省界线数据
plotm([Prince.Lat],[Prince.Lon],'b');  % 绘制省界线
[lat, lon] = inputm;
lat = [lat;lat(1)];
lon = [lon;lon(1)];
T = table(lat,lon,'VariableNames',{'纬度','经度'});
writetable(T,'河南省省界线经纬度坐标2.xlsx');