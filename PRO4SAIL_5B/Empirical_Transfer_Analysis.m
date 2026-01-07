% --- Empirical_Transfer_Analysis.m ---
% 在运行完 Step 2 后接着运行此段代码

% 1. 提取有效的数据点对
X_data = NDVI_Corrected(mask); % 经过拉伸修正后的 NDVI
Y_ref = LAI_Ref(mask);         % 参考真值 LAI

% 2. 执行经验拟合 (典型的 LAI = a * exp(b * NDVI) 模型)
% 为了方便拟合，我们取对数处理：ln(LAI) = ln(a) + b * NDVI
valid_idx = (Y_ref > 0.1); 
fit_coeffs = polyfit(X_data(valid_idx), log(Y_ref(valid_idx)), 1);

b = fit_coeffs(1);
a = exp(fit_coeffs(2));

% 3. 生成经验拟合曲线
ndvi_range = 0.1:0.01:0.9;
lai_fit = a * exp(b * ndvi_range);

% 4. 可视化经验传递函数
figure('Color', 'w');
scatter(X_data, Y_ref, 5, 'filled', 'MarkerFaceAlpha', 0.1); hold on;
plot(ndvi_range, lai_fit, 'r', 'LineWidth', 3);
grid on;
xlabel('Corrected NDVI'); ylabel('Reference LAI');
title(['经验传递函数拟合: LAI = ', num2str(a, '%.3f'), ' * e^{(', num2str(b, '%.3f'), ' * NDVI)}']);
legend('像元点分布', '经验传递曲线');