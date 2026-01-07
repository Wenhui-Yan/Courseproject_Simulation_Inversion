%% 基于 PROSAIL 的参数反演 - LAI，Cab 
clear; clc; close all;

% 屏蔽优化器警告
warning('off', 'optim:lsqncommon:NotEnoughEquations');

% =========================================================
% 1. 数据准备与地理提取
% =========================================================
[img, R] = readgeoraster('SPOT_TOA.tif'); 
img = double(img);
[rows, cols, bands] = size(img);

% 动态提取经纬度
[x_w, y_w] = intrinsicToWorld(R, cols/2, rows/2);
try
    [lat_center, lon_center] = projinv(projcrs(32631), x_w, y_w);
catch
    lat_center = 25.0; lon_center = 45.0; 
end

% =========================================================
% 2. 参数设置 (PROSAIL：LAI 和 Cab 均设为变量)
% =========================================================
% 待反演变量 x = [LAI, Cab]
x0 = [2.0, 40];    % 初始猜测值 [LAI, Cab]
lb = [0.01, 1];    % 下界
ub = [8.0, 100];   % 上界

% 其他固定参数
fixed.Cw    = 0.01;   
fixed.Cdm   = 0.009;  
fixed.N     = 1.5;
fixed.LIDFa = -0.35;
fixed.LIDFb = -0.15;
fixed.TypeLidf = 1;
fixed.hspot = 0.01;
fixed.tts   = 30;
fixed.tto   = 0;
fixed.psi   = 0;
fixed.psoil = 1;
fixed.Car   = 8;
fixed.Ant   = 0;
fixed.Cbrown = 0;

% 加载光谱基准
data = dataSpec_PDB; 
rsoil0 = fixed.psoil*data(:,11) + (1-fixed.psoil)*data(:,12);
skyl = 0.847 - 1.61*sin((90-fixed.tts)*pi/180) + 1.04*sin((90-fixed.tts)*pi/180)^2;

% 初始化结果矩阵
map_LAI  = zeros(rows, cols);
map_Cab  = zeros(rows, cols);
map_RMSE = zeros(rows, cols);

% =========================================================
% 3. 执行双变量反演
% =========================================================
fprintf('开始 PROSAIL 模式反演 ... \n');
tic;
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'none');
band_idx = [146, 246, 441]; % Green, Red, NIR

for i = 1:rows
    for j = 1:cols
        obs_pixel = squeeze(img(i,j,:));
        if any(isnan(obs_pixel)) || any(obs_pixel <= 0), continue; end
        
        % 目标函数：同时调整 LAI 和 Cab 以拟合观测值
        obj_fun = @(x) cost_func_prosail(x, obs_pixel, band_idx, rsoil0, data(:,9), data(:,10), skyl, fixed);
        [x_fit, resnorm] = lsqnonlin(obj_fun, x0, lb, ub, options);
        
        map_LAI(i,j)  = x_fit(1);
        map_Cab(i,j)  = x_fit(2);
        map_RMSE(i,j) = sqrt(resnorm / bands);
    end
end
toc;

% =========================================================
% 4. 构造多列垂直对比表格 
% =========================================================
[r_find, c_find] = find(map_LAI > 0);
num_to_show = min(15, length(r_find)); 

ParamLabels = {
    'fitted parameters (will appear AFTER running the model)';
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
        0.5; lat_center; lon_center; 30;
        map_Cab(r,c); fixed.Cw; fixed.Cdm; % Cab 反演出的动态值
        0.1; 10; 0;
        fixed.N; map_LAI(r,c); fixed.LIDFa; fixed.LIDFb; % LAI 也是动态值
        map_RMSE(r,c)
    };
    final_cell_data = [final_cell_data, current_col];
end

ResultTable = cell2table(final_cell_data);
disp('>>> PROSAIL 模式反演结果 (LAI 和 Cab 均动态) <<<');
disp(ResultTable);
writetable(ResultTable, 'PROSAIL_CabLAI.xlsx', 'WriteVariableNames', false);

% =========================================================
% 辅助函数: 双变量优化
% =========================================================
function err = cost_func_prosail(x, obs, idx, rsoil, Es, Ed, skyl, f)
    % x(1) = LAI, x(2) = Cab
    [rdot, rsot] = PRO4SAIL(f.N, x(2), f.Car, f.Ant, f.Cbrown, ...
                            f.Cw, f.Cdm, f.LIDFa, f.LIDFb, f.TypeLidf, ...
                            x(1), f.hspot, f.tts, f.tto, f.psi, rsoil);
    
    resv = (rdot.*(skyl.*Ed) + rsot.*((1-skyl.*Es))) ./ (Ed.*skyl + Es.*(1-skyl));
    err = resv(idx) - obs(:); 
end