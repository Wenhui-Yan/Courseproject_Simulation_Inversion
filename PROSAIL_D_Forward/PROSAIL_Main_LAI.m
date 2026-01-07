%% ========================================================================
% PRO4SAIL：不同 LAI 观测角反射率仿真 (800 nm)
% ========================================================================
clear; clc; close all;

%% ---------- 1. 基本参数设置 ----------
LAI_all = 0.5:0.5:4;     % LAI 范围
nLAI = numel(LAI_all);
vza_all = -80:5:80;
nVZA = numel(vza_all);
nColor = 256;
cmap = parula(nColor);

% 读取数据（确保 dataSpec_PDB.m 在路径下）
data = dataSpec_PDB;
lambda = data(:,1);
Es = data(:,9);
Ed = data(:,10);
Rsoil1 = data(:,11);
Rsoil2 = data(:,12);

psoil = 1; % 干土
rsoil0 = psoil*Rsoil1 + (1-psoil)*Rsoil2;

% 其他固定参数 (Cab 固定为 40, 观察 LAI 的影响)
N = 1.5; Cab = 40; Car = 8; Ant = 0; Cbrown = 0; Cw = 0.01; Cm = 0.009;
LIDFa = 30; LIDFb = 0; TypeLidf = 2;
tts = 30; hspot = 0.01;

% 选定波长：800 nm (近红外对 LAI 最敏感)
[~, idx800] = min(abs(lambda - 800));
rd = pi/180;
skyl = 0.847 - 1.61*sin((90-tts)*rd) + 1.04*sin((90-tts)*rd).^2;
PARdir = (1-skyl).*Es;
PARdif = skyl.*Ed;

%% ---------- 2. 核心计算循环 ----------
refl_pri = zeros(nLAI, nVZA);
refl_perp = zeros(nLAI, nVZA);

fprintf('正在计算 LAI 曲线，请稍候...\n');
for i = 1:nLAI
    LAI = LAI_all(i);
    for k = 1:nVZA
        vza = vza_all(k);
        tto = abs(vza);
        
        % 主平面计算 (psi = 0/180)
        psi_pri = (vza < 0)*0 + (vza >= 0)*180;
        [rdot1, rsot1] = PRO4SAIL(N, Cab, Car, Ant, Cbrown, Cw, Cm, ...
            LIDFa, LIDFb, TypeLidf, LAI, hspot, tts, tto, psi_pri, rsoil0);
        resv_pri = (rdot1.*PARdif + rsot1.*PARdir) ./ (PARdir + PARdif);
        refl_pri(i,k) = resv_pri(idx800);
        
        % 垂直主平面计算 (psi = 90)
        [rdot2, rsot2] = PRO4SAIL(N, Cab, Car, Ant, Cbrown, Cw, Cm, ...
            LIDFa, LIDFb, TypeLidf, LAI, hspot, tts, tto, 90, rsoil0);
        resv_perp = (rdot2.*PARdif + rsot2.*PARdir) ./ (PARdir + PARdif);
        refl_perp(i,k) = resv_perp(idx800);
    end
end
fprintf('计算完成，开始作图。\n');

%% ---------- 3. 绘图：主平面 ----------
figure('Color','w', 'Name', 'LAI Principal Plane'); 
hold on;
for i = 1:nLAI
    idxCol = round((LAI_all(i)-min(LAI_all))/(max(LAI_all)-min(LAI_all))*(nColor-1)) + 1;
    plot(vza_all, refl_pri(i,:), 'Color', cmap(idxCol,:), 'LineWidth', 1.2);
end

% 样式美化
xlabel('View zenith angle ($^\circ$)', 'Interpreter', 'latex');
ylabel('Reflectance', 'Interpreter', 'latex');
title('$\mathrm{Reflectance\ with\ varying\ LAI\ at\ 800\,nm\ (Principal\ Plane)}$', ...
      'Interpreter','latex','FontSize',10);

xlim([-80 80]); 
ylim([0.35 0.68]); 
set(gca, 'FontName','Times New Roman', 'FontSize', 9, 'Box','on');

% 布局与 Colorbar
drawnow; ax = gca;
ax.Position = [0.12 0.15 0.68 0.75];
colormap(cmap);
c = colorbar;
c.Position = [ax.Position(1)+ax.Position(3)+0.05, ax.Position(2), 0.02, ax.Position(4)];

% 色带设置
c.Limits = [0 1]; 
c.Ticks = linspace(0, 1, nLAI);
c.TickLabels = compose('%.1f', LAI_all); 
c.Label.String = '$\mathrm{LAI}$';       
c.Label.Interpreter = 'latex';
set(c, 'FontName', 'Times New Roman', 'FontSize', 9);

% 导出
exportgraphics(gcf, 'LAI_PROSAIL_Principal.png', 'Resolution', 600);


%% ---------- 4. 绘图：垂直主平面 ----------
figure('Color','w', 'Name', 'LAI Perpendicular Plane'); 
hold on;
for i = 1:nLAI
    idxCol = round((LAI_all(i)-min(LAI_all))/(max(LAI_all)-min(LAI_all))*(nColor-1)) + 1;
    plot(vza_all, refl_perp(i,:), 'Color', cmap(idxCol,:), 'LineWidth', 1.2);
end

xlabel('View zenith angle ($^\circ$)', 'Interpreter', 'latex');
ylabel('Reflectance', 'Interpreter', 'latex');
title('$\mathrm{Reflectance\ with\ varying\ LAI\ at\ 800\,nm\ (Perpendicular\ Plane)}$', ...
      'Interpreter','latex','FontSize',10);

xlim([-80 80]); 
ylim([0.38 0.55]); 
set(gca, 'FontName','Times New Roman', 'FontSize', 9, 'Box','on');

drawnow; ax = gca;
ax.Position = [0.12 0.15 0.68 0.75];
colormap(cmap);
c = colorbar;
c.Position = [ax.Position(1)+ax.Position(3)+0.05, ax.Position(2), 0.02, ax.Position(4)];

% 色带
c.Limits = [0 1];
c.Ticks = linspace(0, 1, nLAI);
c.TickLabels = compose('%.1f', LAI_all);
c.Label.String = '$\mathrm{LAI}$';
c.Label.Interpreter = 'latex';
set(c, 'FontName', 'Times New Roman', 'FontSize', 9);

% 导出 
exportgraphics(gcf, 'LAI_PROSAIL_Perpendicular.png', 'Resolution', 600);