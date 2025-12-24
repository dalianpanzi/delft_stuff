% ncdisp('0811npw\map.nc');
% info = ncinfo("0811npw\map.nc");

node_x=ncread('0811npw\map.nc', 'mesh2d_node_x');
node_y=ncread('0811npw\map.nc', 'mesh2d_node_y');

face_x = ncread('0811npw\map.nc', 'mesh2d_face_x');  % 面中心经度
face_y = ncread('0811npw\map.nc', 'mesh2d_face_y');  % 面中心纬度
face_nodes=ncread("0811npw\map.nc",'mesh2d_face_nodes');

time = ncread('0811npw\map.nc', 'time');
start_time=datetime(2024,1,1,0,0,0);
datetime_vec=start_time+seconds(time);

fill_value = -999;
start_index = ncreadatt('0811npw\map.nc', 'mesh2d_face_nodes', 'start_index');

nFaces = size(face_nodes, 2);
nNodes = size(face_nodes, 1);  % 最大节点数，通常为4

wl = ncread('0811npw\map.nc', 'mesh2d_s1');
wl(wl==fill_value)=NaN;
[~,nwlt]=size(wl);

target_start=datetime(2024,1,16);
[~, start_idx]=min(abs(datetime_vec-target_start));
end_idx=nwlt-5;

sl_current=wl(:,start_idx:end_idx);
%绘制水位
% figure('Position',[100,100,600,1000])
% plot_faces=min(10000,nFaces);
% 
% for i =1:plot_faces
%     node_indices=face_nodes(:,i);
%     node_indices=node_indices(~isnan(node_indices));
%     
%     if ~isempty(node_indices)
%         x_coords=node_x(node_indices);
%         y_coords=node_y(node_indices);
%         patch(x_coords,y_coords,sl_current(i),...
%             'EdgeColor','none','FaceAlpha',0.8);
%         hold on
%     end
% end

%提取特定区域
lon_min=117; lon_max=128;
lat_min=24; lat_max=41;

in_region=(face_y>=lat_min&face_y<=lat_max&...
    face_x>=lon_min & face_x <= lon_max);

region_indices = find(in_region); 
region_face_x=face_x(in_region);
region_face_y=face_y(in_region);
region_sl=sl_current(in_region,:);
nRegionPoints = length(region_indices);

% figure
% scatter(region_face_x, region_face_y, 20, region_sl, 'filled')
% colorbar

[nFace,nTime]=size(region_sl);

time_analyze=datetime_vec(start_idx:end_idx);
startT=[2024,1,16,0,0,0];

valid_data_points=~any(isnan(region_sl), 2);
valid_indices = find(valid_data_points);
nValidPoints = length(valid_indices);

M2_amp_region = NaN(nRegionPoints, 1);  % 区域内所有点的振幅（初始为NaN）
M2_phase_region = NaN(nRegionPoints, 1); % 区域内所有点的相位（初始为NaN）
valid_count = 0;

for idx_counter = 1:nValidPoints
    i=valid_indices(idx_counter);
    el=region_sl(i,:)';
    lat=region_face_y(i);

    [amp,phase]=compute_tide_constants(el,lat,startT);
    M2_amp_region(i) = amp;
    M2_phase_region(i) = phase;
    if ~isnan(amp) && ~isnan(phase)
        valid_count = valid_count + 1;
    end

    if mod(i, 100) == 0
        fprintf('进度: %d/%d (%.1f%%)，有效点数: %d\n', ...
            i, nRegionPoints, i/nRegionPoints*100, valid_count);
    end
end

fprintf('\n计算完成！\n');
fprintf('总点数: %d，有效点数: %d (%.1f%%)\n', ...
    nRegionPoints, valid_count, valid_count/nRegionPoints*100);

%填色图
figure('Position',[100,100,600,1000]);
scatter(region_face_x,region_face_y,30,M2_amp_region,'filled');
colorbar;
caxis([min(M2_amp_region), max(M2_amp_region)]);
title('M2分潮振幅 - 散点图');
xlabel('经度 (°E)');
ylabel('纬度 (°N)');
xlim([lon_min, lon_max]);
ylim([lat_min, lat_max]);
grid on;

valid_points = ~isnan(region_face_x) & ~isnan(region_face_y) & ...
               ~isnan(M2_amp_region);

x_clean = region_face_x(valid_points);
y_clean = region_face_y(valid_points);
amp_clean = M2_amp_region(valid_points);

amp_clean_cm=amp_clean*100;  %原本得出的M2振幅是m

% 2. 创建插值函数
F_amp = scatteredInterpolant(x_clean, y_clean, amp_clean_cm, 'natural', 'none');

xi = linspace(lon_min, lon_max, 500);
yi = linspace(lat_min, lat_max, 600);
[XI, YI] = meshgrid(xi, yi);
ZI = F_amp(XI, YI);

try
    % 方法A: 使用alpha形状（更精确）
    shp = alphaShape(x_clean, y_clean, 0.1);
    in = inShape(shp, XI(:), YI(:));
catch
    % 方法B: 使用凸包（更简单）
    k = convhull(x_clean, y_clean);
    in = inpolygon(XI(:), YI(:), x_clean(k), y_clean(k));
end

mask = reshape(in, size(XI));
ZI(~mask) = NaN;  % 将区域外的值设为NaN

figure('Position',[100,100,800,1000]);

land = shaperead('landareas', 'UseGeoCoords', true);
geoshow(land, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'k', 'LineWidth', 0.5);
hold on;
% 绘制等值线
[C, h] = contour(XI, YI, ZI, 20, 'LineColor', 'k', 'LineWidth', 0.5);

% clabel(C,h,'FontSize', 8);
hold on;

%等值线变成整数，但效果不够好, 可跳过不使用
idx = 1; % 索引指针
while idx < size(C, 2)
    level = C(1, idx);          % 当前等值线的数值
    num_vertices = C(2, idx);   % 构成这条等值线的点的数量
   
    % 提取这条等值线上所有点的坐标
    x_data = C(1, idx+1:idx+num_vertices);
    y_data = C(2, idx+1:idx+num_vertices);
   
    % 选择在等值线的大约1/3长度处放置一个标签（避免在两端）
    if num_vertices >= 3
        label_pos = floor(num_vertices / 3);
        % 将数值四舍五入为整数，并创建文本
        % 如果你想直接舍去小数，可将round改为fix
        text(x_data(label_pos), y_data(label_pos), ...
             sprintf('%d', round(level)), ... % 关键：格式化为整数
             'FontSize', 8, ...
             'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'middle', ...
             'BackgroundColor', 'w', ...      % 白色背景让标签更清晰
             'Margin', 0.5, ...
             'EdgeColor', 'none');
    end
   
    % 移动索引指针到下一条等值线开始的位置
    idx = idx + num_vertices + 1;
end

% 添加颜色条
cbar = colorbar;
ylabel(cbar, '振幅 (cm)', 'FontSize', 12);
% hold on;
% land = shaperead('landareas', 'UseGeoCoords', true); % 读取全球陆地数据
% geoshow(land, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'k', 'LineWidth', 0.5); % 填充灰色并描边

% 设置图形属性
title('M2分潮振幅等值线图', ...
      'FontSize', 14, 'FontWeight', 'bold');
xlabel('经度 (°E)', 'FontSize', 12);
ylabel('纬度 (°N)', 'FontSize', 12);
grid on;
axis equal;
xlim([lon_min, lon_max]);
ylim([lat_min, lat_max]);


%方法单独创建一个文件
function [M2_amp, M2_phase]=compute_tide_constants(el,lat,start_time)
    if any(isnan(el))
        M2_amp=NaN;
        M2_phase=NaN;
        return;
    end
    interval_hours=1;
    [tc_name,tc_freq, tc_tidecon, ~]=t_tide(el, ...
        'interval',interval_hours,...
        'start time', start_time,...
        'rayleigh', ['M2';'S2';'N2';'K2';'K1';'O1';'P1';'Q1';],...
        'shallow',['MS4';'M4 ';'M6 ';],...
        'latitude',lat,...              
        'output','none');
    m2_idx=[];

    for i =1:size(tc_name,1)
        name=deblank(tc_name(i,:));
        if strcmp(strtrim(name),'M2')
            m2_idx=i;
            break;
        end
    end
    if ~isempty(m2_idx)
        M2_amp=tc_tidecon(m2_idx,1);
        M2_phase=tc_tidecon(m2_idx,3);
    else
        M2_phase=NaN;
        M2_amp=NaN;
    end
end

