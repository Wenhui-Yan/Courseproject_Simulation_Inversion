% --- Step1_LUT_Precision_Alignment.m ---
clc; clear;
data = dataSpec_P5B; 
num_samples = 15000; 

% --- 核心改进：等步长采样确保覆盖范围 ---
LAI_vec = linspace(0.01, 8, num_samples)'; 
Cab = linspace(20, 80, num_samples)';  % 提高叶绿素上限
psoil = linspace(0, 1, num_samples)';  % 土壤亮度全覆盖

% 其他参数保持相对稳健
LIDFa = -0.35; LIDFb = -0.15; % 锁定球型分布，减少干扰
Cw = 0.015; Cm = 0.012; N = 1.5; 
Car = 8; Cbrown = 0; TypeLidf = 1; hspot = 0.1; 
tts = 35; tto = 0; psi = 0; 
skyl = 0.25; % 针对 7/20 较重的大气影响，调高散射比

LUT_outputs = zeros(num_samples, 3); 
LUT_inputs = LAI_vec; 

disp('正在生成高密度物理 LUT...');
for i = 1:num_samples
    rsoil0 = psoil(i)*data(:,10) + (1-psoil(i))*data(:,11);
    [rdot, rsot, ~, ~] = PRO4SAIL(N, Cab(i), Car, Cbrown, Cw, Cm, LIDFa, LIDFb, TypeLidf, LAI_vec(i), hspot, tts, tto, psi, rsoil0);
    resv = rdot * skyl + rsot * (1 - skyl);
    LUT_outputs(i, 2) = mean(resv(610-400+1 : 680-400+1)); % Red
    LUT_outputs(i, 3) = mean(resv(790-400+1 : 890-400+1)); % NIR
end
save('PROSAIL_LUT.mat', 'LUT_inputs', 'LUT_outputs');
disp('高密度 LUT 生成成功！');