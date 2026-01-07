%% 基于 PRO4SAIL 的参数反演 - LAI 
clear; clc; close all;

% 屏蔽优化器警告
warning('off', 'optim:lsqncommon:NotEnoughEquations');

% =========================================================
% 1. 数据准备与地理提取
% =========================================================
[img, R] = readgeoraster('SPOT_TOA.tif'); 
img = double(img);
[rows, cols, bands] = size(img);

% 动态提取影像中心经纬度 (用于表格显示)
[x_w, y_w] = intrinsicToWorld(R, cols/2, rows/2);
try
    [lat_center, lon_center] = projinv(projcrs(32631), x_w, y_w);
catch
    lat_center = 25.0; lon_center = 45.0; 
end

% =========================================================
% 2. 参数设置 (反演 LAI，固定 Cab, Cw, Cdm)
% =========================================================
% 反演变量仅为 LAI: x = [LAI]
x0 = 2.0;    % 初始猜测值
lb = 0.01;   % 下界
ub = 8.0;    % 上界

% 固定值 (根据你之前的提取结果设置)
fixed.Cab   = 40;     % ug/cm2
fixed.Cw    = 0.01;   % g/cm2
fixed.Cdm   = 0.009;  % g/cm2
fixed.B     = 0.5;
fixed.SMp   = 30;
fixed.Cs    = 0.1;
fixed.Cca   = 10;
fixed.Cant  = 0;
fixed.N     = 1.5;
fixed.LIDFa = -0.35;
fixed.LIDFb = -0.15;
fixed.TypeLidf = 1;
fixed.Car   = 8;
fixed.Ant   = 0;
fixed.Cbrown = 0;
fixed.hspot = 0.01;
fixed.tts   = 30;
fixed.tto   = 0;
fixed.psi   = 0;
fixed.psoil = 1;

% 加载光谱基准
data = dataSpec_PDB; 
rsoil0 = fixed.psoil*data(:,11) + (1-fixed.psoil)*data(:,12);
skyl = 0.847 - 1.61*sin((90-fixed.tts)*pi/180) + 1.04*sin((90-fixed.tts)*pi/180)^2;

% 初始化结果矩阵
map_LAI  = zeros(rows, cols);
map_RMSE = zeros(rows, cols);

% =========================================================
% 3. 执行反演 (逐像元计算 LAI)
% =========================================================
fprintf('开始反演 LAI... 正在处理 %d 个像元...\n', rows*cols);
tic;
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'none');
band_idx = [146, 246, 441]; % Green, Red, NIR

% 若要提速，可将此处的 for 改为 parfor
for i = 1:rows
    for j = 1:cols
        obs_pixel = squeeze(img(i,j,:));
        if any(isnan(obs_pixel)) || any(obs_pixel <= 0), continue; end
        
        % 目标函数定义 (只优化 LAI)
        obj_fun = @(x) cost_func_lai_only(x, obs_pixel, band_idx, rsoil0, data(:,9), data(:,10), skyl, fixed);
        [x_fit, resnorm] = lsqnonlin(obj_fun, x0, lb, ub, options);
        
        map_LAI(i,j)  = x_fit(1);
        map_RMSE(i,j) = sqrt(resnorm / bands);
    end
end
toc;

% =========================================================
% 4. 绘制 LAI 分布图 (灰度图)
% =========================================================
figure('Color', 'w', 'Name', 'LAI Inversion Result');
imagesc(map_LAI);
colormap(gray); 
h = colorbar;
ylabel(h, 'Leaf Area Index (m^2/m^2)');
title('PROSAIL 反演生成的 LAI 空间分布图');
axis image; 
% 自动优化对比度
if max(map_LAI(:)) > 0
    caxis([quantile(map_LAI(map_LAI>0), 0.05), quantile(map_LAI(map_LAI>0), 0.95)]);
end

% =========================================================
% 5. 构造多列垂直对比表格
% =========================================================
[r_find, c_find] = find(map_LAI > 0);
num_to_show = min(15, length(r_find)); 

ParamLabels = {
    'fitted parameters(will appear AFTER running the model)';
    'B'; 'lat'; 'lon'; 'SMp'; 
    'Cab (ug cm-2)'; 'Cw (g cm-2)'; 'Cdm (g cm-2)'; 
    'Cs (a.u)'; 'Cca (ug cm-2)'; 'Cant'; 
    'N (dimensionless)'; 'LAI'; 'LIDFa'; 'LIDFb'; 
    'RMSE (mod-meas spectra)'
};

final_cell_data = ParamLabels;

for k = 1:num_to_show
    r = r_find(k); c = c_find(k);
    current_col = {
        sprintf('Pixel_%d', k);
        fixed.B; lat_center; lon_center; fixed.SMp;
        fixed.Cab; fixed.Cw; fixed.Cdm;
        fixed.Cs; fixed.Cca; fixed.Cant;
        fixed.N; map_LAI(r,c); fixed.LIDFa; fixed.LIDFb;
        map_RMSE(r,c)
    };
    final_cell_data = [final_cell_data, current_col];
end

ResultTable = cell2table(final_cell_data);
disp('>>> 拟合参数结果提取 (LAI 反演版) <<<');
disp(ResultTable);
writetable(ResultTable, 'SAIL_LAI.xlsx', 'WriteVariableNames', false);

% =========================================================
% 辅助函数: 仅针对 LAI 进行优化
% =========================================================
function err = cost_func_lai_only(x, obs, idx, rsoil, Es, Ed, skyl, f)
    % x(1) 此时代表 LAI
    [rdot, rsot, ~, ~] = PRO4SAIL(f.N, f.Cab, f.Car, f.Ant, f.Cbrown, ...
                                 f.Cw, f.Cdm, f.LIDFa, f.LIDFb, f.TypeLidf, ...
                                 x(1), f.hspot, f.tts, f.tto, f.psi, rsoil);
    
    % 结合直射与散射光的反射率
    PARdiro = (1-skyl) * Es;
    PARdifo = (skyl) * Ed;
    resv = (rdot .* PARdifo + rsot .* PARdiro) ./ (PARdiro + PARdifo);
    
    err = resv(idx) - obs(:); 
end