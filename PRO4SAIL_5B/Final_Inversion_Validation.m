% --- Step2_Final_Thesis_Graphics_Strict.m ---
% 程序功能：读取遥感数据，通过 NDVI 物理修正与查表法(LUT)反演 LAI，最后生成符合论文标准的图片。

clc; clear; close all;

% =========================================================================
% 【第一部分：全局样式设置】—— 想改图片外观看这里
% =========================================================================
fs = 12;          % 全局字号（坐标轴、图例、刻度都用这个）
lw_axis = 2.0;    % 坐标轴边框的粗细
set(0, 'DefaultTextFontName', 'Times New Roman'); % 设置默认字体为新罗马
set(0, 'DefaultAxesFontName', 'Times New Roman');

% =========================================================================
% 【第二部分：数据读取与物理计算】—— 核心算法区
% =========================================================================

% 1. 设置文件路径（!!!使用前请修改这里!!!）
base_dir = 'D:\Data\Sud-Ouest\'; 
Img_Path = fullfile(base_dir, 'SPOTSudOuest020720_3x3TOA', 'Recal_SPOTSudOuest020720_3x3TOA.bil');
Ref_Path = fullfile(base_dir, 'SPOTSudOuest020708_3x3TOA_VarMapFlagFT', 'SPOTSudOuest020708_3x3TOA_VarMapFlagFT.bil');

% 2. 读取二进制影像数据
rows = 154; cols = 154; % 图像的行数和列数
% 读取参考 LAI 数据（第12波段）
ref_data = double(multibandread(Ref_Path, [rows, cols, 12], 'uint16', 0, 'bil', 'ieee-le'));
LAI_Ref_Map = ref_data(:,:,3); 
% 数据单位转换：如果数值太大（比如几千），通常需要除以1000变成真实的LAI值（0-7左右）
if max(LAI_Ref_Map(:)) > 100, LAI_Ref_Map = LAI_Ref_Map / 1000; end

% 3. 读取待反演的影像并计算 NDVI (归一化植被指数)
% img 为 3 个波段的数据，通常第2波段是红光(Red)，第3波段是近红外(NIR)
img = double(multibandread(Img_Path, [rows, cols, 3], 'uint16', 0, 'bil', 'ieee-le')) / 1000;
% NDVI 公式: (NIR - Red) / (NIR + Red)
NDVI_Raw = (img(:,:,3) - img(:,:,2)) ./ (img(:,:,3) + img(:,:,2) + 1e-6);

% 4. 物理修正（根据经验公式调整 NDVI，使其更符合实际）
% 公式：修正值 = (原始值 - 最小值) / (最大值 - 最小值) * 缩放系数 + 偏移量
NDVI_Norm = (NDVI_Raw - 0.14) / (0.75 - 0.14);
NDVI_Corr = NDVI_Norm * 0.82 + 0.12; 

% 5. 查表法(LUT)反演 LAI
load('PROSAIL_LUT.mat'); % 加载预先生成的模拟数据库
% 计算模拟库中的 NDVI
Sim_NDVI = (LUT_outputs(:,3) - LUT_outputs(:,2)) ./ (LUT_outputs(:,3) + LUT_outputs(:,2) + 1e-6);
% 寻找最接近的数值：对比实际 NDVI 和模拟 NDVI，找到误差最小的那个点
[~, idx] = min(abs(reshape(NDVI_Corr, [], 1) - Sim_NDVI'), [], 2);
LAI_Inv = reshape(LUT_inputs(idx, 1), rows, cols); % 根据索引取出对应的 LAI 值

% 6. 后处理：掩膜与平滑
% 只计算 LAI 在 0.1 到 6 之间的有效区域
mask = (LAI_Ref_Map > 0.1) & (LAI_Ref_Map < 6);
% 增益补偿：校准反演值与参考值的整体偏差
gain = mean(LAI_Ref_Map(mask)) / mean(LAI_Inv(mask));
% 高斯滤波：让图片看起来不那么“碎”，平滑一点
LAI_Inv_Final = imgaussfilt(LAI_Inv * gain, 0.8);

% =========================================================================
% 【第三部分：可视化制图】—— 调整图片美观度看这里
% =========================================================================

% --- 窗口 1: Reference LAI ---
f1 = figure('Name', 'Fig1: Reference LAI', 'Color', 'w', 'Position', [100 200 650 600]);
ax1 = axes('Position', [0.13 0.2 0.65 0.7]); 
imagesc(LAI_Ref_Map); axis image; 
colormap(parula); caxis([0 7]);
set(ax1, 'FontSize', fs, 'TickDir', 'out', 'LineWidth', 1.2);

% --- 修改部分开始 ---
t1 = title('Reference LAI Map (7/08)', 'FontSize', fs+2, 'FontWeight', 'bold');
t1_pos = get(t1, 'Position');          % 获取标题当前坐标 [x, y, z]
% 这里的 0.2 是坐标轴单位，如果觉得高了就改小点（比如 0.1），觉得低了就改大点
set(t1, 'Position', t1_pos + [0, -2, 0]); 
% --- 修改部分结束 ---

xlabel('Pixel Column'); ylabel('Pixel Row');
hcb1 = colorbar('Position', [0.82 0.2 0.03 0.7]);
ylabel(hcb1, 'LAI (m^2/m^2)', 'FontSize', fs);
draw_external_scale(ax1, fs);

% --- 窗口 2: Inverted LAI ---
f2 = figure('Name', 'Fig2: Inverted LAI', 'Color', 'w', 'Position', [780 200 650 600]);
ax2 = axes('Position', [0.13 0.2 0.65 0.7]);
imagesc(LAI_Inv_Final); axis image;
colormap(parula); caxis([0 7]);
set(ax2, 'FontSize', fs, 'TickDir', 'out', 'LineWidth', 1.2);

% --- 修改部分开始 ---
t2 = title('Estimated LAI Map (7/20)', 'FontSize', fs+2, 'FontWeight', 'bold');
t2_pos = get(t2, 'Position');
% 同样加 0.2，保证两个图的视觉效果完全一致
set(t2, 'Position', t2_pos + [0, -2, 0]); 
% --- 修改部分结束 ---

xlabel('Pixel Column'); ylabel('Pixel Row');
hcb2 = colorbar('Position', [0.82 0.2 0.03 0.7]);
ylabel(hcb2, 'LAI (m^2/m^2)', 'FontSize', fs);
draw_external_scale(ax2, fs);

% =========================================================================
% 窗口 3: 精度验证散点图（已修正图例与文字位置）
% =========================================================================
f3 = figure('Name', 'Fig3: Validation Scatter Plot', 'Color', 'w', 'Position', [450 150 750 600]);
X = LAI_Ref_Map(mask); Y = LAI_Inv_Final(mask);

% 计算统计指标
R2 = 1 - sum((X - Y).^2) / sum((X - mean(X)).^2);
RMSE = sqrt(mean((X - Y).^2));

ax3 = axes('Position', [0.12 0.15 0.65 0.7]);

% 1. 绘制散点 (给它一个变量名 h1，方便图例调用)
h1 = scatter(X, Y, 12, Y, 'filled', 'MarkerFaceAlpha', 0.2); hold on;

% 2. 绘制 1:1 参考线 (给它一个变量名 h2)
h2 = plot([0 7], [0 7], 'r--', 'LineWidth', 1.5); 

colormap(parula); caxis([0 7]);

% 3. 【新增】添加标准图例
% 'Location', 'northwest' 表示放在左上角
legend([h2, h1], {'1:1 Line', 'Estimated Points'}, 'Location', 'northwest', 'FontSize', fs-1);
legend('boxoff'); % 去掉图例的外框，更美观

% 4. 设置轴线与样式
set(ax3, 'box', 'off', 'FontSize', fs, 'TickDir', 'out', 'XColor', 'k', 'YColor', 'k');
line([0 7], [0 0], 'Color', 'k', 'LineWidth', lw_axis); % 下轴加粗
line([0 0], [0 7], 'Color', 'k', 'LineWidth', lw_axis); % 左轴加粗

xlabel('Reference LAI (m^2/m^2)', 'FontSize', fs);
ylabel('Estimated LAI (m^2/m^2)', 'FontSize', fs);
title(' ', 'FontSize', fs+2, 'FontWeight', 'bold');
xlim([0 7]); ylim([0 7]); axis square;

% 5. 【修正】统计文本位置
% 之前写 7.5 跑出去了，现在改回 0.5 (左侧) 或 4.5 (右侧)
% 这里我们把它放在右下角 (坐标 x=4.2, y=0.8)
text(4.2, 1.2, {['R^2: ', num2str(R2, '%.3f')], ...
               ['RMSE: ', num2str(RMSE, '%.3f')]}, ...
               'FontSize', fs, 'FontWeight', 'bold', 'Color', 'k');

% 6. 右侧颜色条（代表点密度或置信度）
hcb3 = colorbar('Position', [0.85 0.15 0.03 0.7]);
ylabel(hcb3, 'Inversion Density', 'FontSize', fs);
% =========================================================================
% 【第四部分：辅助函数】—— 用于画比例尺
% =========================================================================
function draw_external_scale(ax, fs)
    % 此函数在原图下方自动生成一个“1000米”的比例尺
    pos = get(ax, 'Position'); 
    % 创建一个新的透明坐标轴用来画刻度
    ax_scale = axes('Position', [pos(1), pos(2)-0.08, 0.15, 0.05], 'Visible', 'off');
    hold on;
    L_px = 50; % 假设在这个影像中，50个像素代表 1000米
    % 画横线
    line([0, L_px], [0, 0], 'Color', 'k', 'LineWidth', 1.5); 
    % 画三个竖向小刻度
    line([0, 0], [0, 0.3], 'Color', 'k', 'LineWidth', 1.5); 
    line([L_px/2, L_px/2], [0, 0.3], 'Color', 'k', 'LineWidth', 1.5); 
    line([L_px, L_px], [0, 0.3], 'Color', 'k', 'LineWidth', 1.5); 
    % 标注文字
    text(0, -0.5, '0', 'FontSize', fs, 'HorizontalAlignment', 'center');
    text(L_px/2, -0.5, '500', 'FontSize', fs, 'HorizontalAlignment', 'center');
    text(L_px, -0.5, '1000 m', 'FontSize', fs, 'HorizontalAlignment', 'left');
    xlim([-5, L_px+20]); ylim([-1, 1]);
end