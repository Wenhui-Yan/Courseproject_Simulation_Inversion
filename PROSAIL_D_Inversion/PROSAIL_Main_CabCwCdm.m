%% 基于 PRO4SAIL 的参数反演 - 多列垂直对比表输出
clear; clc; close all;

% 屏蔽警告
warning('off', 'optim:lsqncommon:NotEnoughEquations');

% =========================================================
% 1. 数据准备与地理信息动态提取
% =========================================================
[img, R] = readgeoraster('SPOT_TOA.tif'); 
img = double(img);
[rows, cols, bands] = size(img);

% 动态提取影像中心经纬度 (用于填充表格中的 lat/lon 行)
[x_w, y_w] = intrinsicToWorld(R, cols/2, rows/2);
try
    % 假设投影为 UTM 31N (EPSG:32631)
    [lat_center, lon_center] = projinv(projcrs(32631), x_w, y_w);
catch
    lat_center = 25.0; lon_center = 45.0; % 容错处理
end

% =========================================================
% 2. 参数与变量初始化
% =========================================================
% 反演变量 [Cab, Cw, Cdm]
x0 = [40, 0.01, 0.009]; 
lb = [1, 0.0001, 0.0001]; 
ub = [100, 0.1, 0.1];

% 固定参数定义
fixed.B = 0.5; fixed.SMp = 30; fixed.Cs = 0.1; fixed.Cca = 10;
fixed.Cant = 0; fixed.N = 1.5; fixed.LAI = 2.0;
fixed.LIDFa = -0.35; fixed.LIDFb = -0.15; fixed.TypeLidf = 1;
fixed.Car = 8; fixed.Ant = 0; fixed.Cbrown = 0;
fixed.hspot = 0.01; fixed.tts = 30; fixed.tto = 0; fixed.psi = 0; fixed.psoil = 1;

% 加载光谱数据
data = dataSpec_PDB;
rsoil0 = fixed.psoil*data(:,11) + (1-fixed.psoil)*data(:,12);
skyl = 0.847 - 1.61*sin((90-fixed.tts)*pi/180) + 1.04*sin((90-fixed.tts)*pi/180)^2;

% 用于存储结果的临时矩阵
map_Cab = zeros(rows, cols);
map_Cw = zeros(rows, cols);
map_Cdm = zeros(rows, cols);
map_RMSE = zeros(rows, cols);

% =========================================================
% 3. 执行反演 (建议先跑小块区域测试，如 i=1:5)
% =========================================================
fprintf('开始反演计算... \n');
tic;
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'none');
band_idx = [146, 246, 441];

% 如果你的电脑支持并行，可以将 for 改为 parfor 提速
for i = 1:rows
    for j = 1:cols
        obs_pixel = squeeze(img(i,j,:));
        if any(isnan(obs_pixel)) || any(obs_pixel <= 0), continue; end
        
        obj_fun = @(x) cost_func_final(x, obs_pixel, band_idx, rsoil0, data(:,9), data(:,10), skyl, fixed);
        [x_fit, resnorm] = lsqnonlin(obj_fun, x0, lb, ub, options);
        
        map_Cab(i,j) = x_fit(1);
        map_Cw(i,j) = x_fit(2);
        map_Cdm(i,j) = x_fit(3);
        map_RMSE(i,j) = sqrt(resnorm / bands);
    end
end
toc;

% =========================================================
% 4. 构造"多列垂直对比表"
% =========================================================
% 寻找前 15 个有效的反演结果
[r_find, c_find] = find(map_Cab > 0);
num_to_show = min(15, length(r_find)); % 确定展示多少列

% 初始化表格数据 (第一列是参数名)
ParamLabels = {
    'fitted parameters(will appear AFTER running the model)';
    'B'; 'lat'; 'lon'; 'SMp'; 
    'Cab (ug cm-2)'; 'Cw (g cm-2)'; 'Cdm (g cm-2)'; 
    'Cs (a.u)'; 'Cca (ug cm-2)'; 'Cant'; 
    'N (dimensionless)'; 'LAI'; 'LIDFa'; 'LIDFb'; 
    'RMSE (mod-meas spectra)'
};

% 创建一个 cell 数组用于合并所有列
final_cell_data = ParamLabels;

for k = 1:num_to_show
    r = r_find(k); c = c_find(k);
    
    % 当前像元的列数据
    current_col = {
        sprintf('Pixel_%d', k); % 第一行标题下列名
        fixed.B; lat_center; lon_center; fixed.SMp;
        map_Cab(r,c); map_Cw(r,c); map_Cdm(r,c);
        fixed.Cs; fixed.Cca; fixed.Cant;
        fixed.N; fixed.LAI; fixed.LIDFa; fixed.LIDFb;
        map_RMSE(r,c)
    };
    
    % 横向合并到总表
    final_cell_data = [final_cell_data, current_col];
end

% 转换为 Table 并导出
ResultTable = cell2table(final_cell_data);
disp('>>> 垂直多列对比表 (展示前15个有效像元) <<<');
disp(ResultTable);

% 保存到 Excel，不带变量名以保持我们的自定义标题
writetable(ResultTable, 'PROSAIL_CabCwCdm.xlsx', 'WriteVariableNames', false);

% =========================================================
% 辅助函数
% =========================================================
function err = cost_func_final(x, obs, idx, rsoil, Es, Ed, skyl, f)
    [rdot, rsot, ~, ~] = PRO4SAIL(f.N, x(1), f.Car, f.Ant, f.Cbrown, ...
                                 x(2), x(3), f.LIDFa, f.LIDFb, f.TypeLidf, ...
                                 f.LAI, f.hspot, f.tts, f.tto, f.psi, rsoil);
    resv = (rdot.*(skyl.*Ed) + rsot.*((1-skyl).*Es)) ./ (Ed.*skyl + Es.*(1-skyl));
    err = resv(idx) - obs(:); 
end