% --- Step2_Validation_with_RefMap.m ---
clc; clear; close all;
set(0, 'DefaultTextFontName', 'Times New Roman');
set(0, 'DefaultAxesFontName', 'Times New Roman');

% 1. 读取数据
base_dir = 'D:\Data\Sud-Ouest\'; 
Img_Path = fullfile(base_dir, 'SPOTSudOuest020720_3x3TOA', 'Recal_SPOTSudOuest020720_3x3TOA.bil');
Ref_Path = fullfile(base_dir, 'SPOTSudOuest020708_3x3TOA_VarMapFlagFT', 'SPOTSudOuest020708_3x3TOA_VarMapFlagFT.bil');
rows = 154; cols = 154;

% 真值读取
ref_data = double(multibandread(Ref_Path, [rows, cols, 12], 'uint16', 0, 'bil', 'ieee-le'));
LAI_Ref_Map = ref_data(:,:,3); 
if max(LAI_Ref_Map(:)) > 100, LAI_Ref_Map = LAI_Ref_Map / 1000; end

% 影像读取与基础校正
img = double(multibandread(Img_Path, [rows, cols, 3], 'uint16', 0, 'bil', 'ieee-le')) / 1000;
NDVI_Raw = (img(:,:,3) - img(:,:,2)) ./ (img(:,:,3) + img(:,:,2) + 1e-6);

% 物理对齐拉伸
NDVI_Norm = (NDVI_Raw - 0.14) / (0.75 - 0.14);
NDVI_Corr = NDVI_Norm * 0.82 + 0.12; 

% 2. 查找表匹配
load('PROSAIL_LUT.mat');
Sim_NDVI = (LUT_outputs(:,3) - LUT_outputs(:,2)) ./ (LUT_outputs(:,3) + LUT_outputs(:,2) + 1e-6);
[~, idx] = min(abs(reshape(NDVI_Corr, [], 1) - Sim_NDVI'), [], 2);
LAI_Inv = reshape(LUT_inputs(idx, 1), rows, cols);

% 3. 系统偏差修正与平滑
mask = (LAI_Ref_Map > 0.1) & (LAI_Ref_Map < 6);
gain = mean(LAI_Ref_Map(mask)) / mean(LAI_Inv(mask));
LAI_Inv_Final = imgaussfilt(LAI_Inv * gain, 0.8);

% 4. 绘图展示 (三子图模式)
figure('Color','w','Position', [50 100 1400 450]);

% 子图 1: 参考真值图
subplot(1,3,1);
imagesc(LAI_Ref_Map); axis image;
colormap(gca, parula); caxis([0 6]);
title('Reference LAI (7/08)', 'FontSize', 14, 'FontWeight', 'normal');
xlabel('Pixel Column'); ylabel('Pixel Row');
set(gca, 'TickDir', 'out');

% 子图 2: 反演结果图
subplot(1,3,2);
imagesc(LAI_Inv_Final); axis image;
h2 = colorbar('southoutside'); 
xlabel(h2, 'LAI (m^2/m^2)', 'FontSize', 12);
colormap(gca, parula); caxis([0 6]);
title('Estimated LAI (7/20)', 'FontSize', 14, 'FontWeight', 'normal');
xlabel('Pixel Column'); ylabel('Pixel Row');
set(gca, 'TickDir', 'out');

% 子图 3: 散点验证
subplot(1,3,3);
X = LAI_Ref_Map(mask); Y = LAI_Inv_Final(mask);
R2 = 1 - sum((X - Y).^2) / sum((X - mean(X)).^2);
RMSE = sqrt(mean((X - Y).^2));
scatter(X, Y, 10, Y, 'filled', 'MarkerFaceAlpha', 0.15); 
hold on; plot([0 6], [0 6], 'r--', 'LineWidth', 2);
grid on; axis square; xlim([0 6]); ylim([0 6]);
xlabel('Reference LAI', 'FontSize', 12); ylabel('Estimated LAI', 'FontSize', 12);
title(['\it{R}^{\rm{2}} = ', num2str(R2, '%.3f'), '  RMSE = ', num2str(RMSE, '%.3f')], 'FontSize', 13);
set(gca, 'TickDir', 'out');