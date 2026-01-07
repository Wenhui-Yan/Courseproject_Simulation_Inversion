% ********************************************************************************
%% Plot of Reflectance and Transmittance with varying Cab (PROSPECT-D)
%% Author: WenhuiYan
%% Affiliation: CAU
%% Email: wenhuiyan233@163.com
% ********************************************************************************

clear; clc; close all;

%% ---- Parameter range ----
Cab_all = 0:10:100;     % Chlorophyll content (ug/cm^2)
nC = length(Cab_all);

%% ---- Continuous colormap ----
nColor = 256;                      
cmap = parula(nColor);             

%% ---- Figure ----
figure('Color', 'w'); 
hold on;

%% ---- Loop over Cab ----
for i = 1:nC
    Cab = Cab_all(i);

    % ---- Fixed leaf parameters ----
    N      = 1.2;
    Car    = 10;
    Anth   = 1;
    Cbrown = 0;
    Cw     = 0.015;
    Cm     = 0.009;

    % ---- PROSPECT-D ----
    LRT = prospect_DB(N, Cab, Car, Anth, Cbrown, Cw, Cm);

    wl = LRT(:,1);
    R  = LRT(:,2);
    T  = LRT(:,3);

    % ---- Color mapping ----
    idx = round( (Cab - min(Cab_all)) / (max(Cab_all) - min(Cab_all)) * (nColor-1) ) + 1;
    col = cmap(idx,:);

    % ---- Plot reflectance ----
    plot(wl, R, 'Color', col, 'LineWidth', 0.8);

    % ---- Plot mirrored transmittance ----
    plot(wl, 1 - T, 'Color', col, 'LineWidth', 0.8);
end

%% ---- Axes settings ----
x_margin = 0.02 * (2500 - 400);
xlim([400 - x_margin, 2500 + x_margin]);
xticks([400 700 1000 1300 1600 1900 2200 2500]);

% Left y-axis
yyaxis left
ylim([0 1]);
yticks(0:0.2:1);
ylabel('Reflectance');
set(gca, 'YColor', 'k');

% Right y-axis
yyaxis right
ylim([0 1]);
set(gca, 'YDir','reverse');
yticks(0:0.2:1);
ylabel('Transmittance');
set(gca, 'YColor', 'k');

xlabel('Wavelength (nm)');

%% ---- Title ----
title('$\mathrm{Reflectance\ and\ Transmittance\ with\ varying\ } C_{ab}\ \mathrm{(PROSPECT\!\!-\!D)}$', ...
      'Interpreter','latex','FontSize',9);

set(gca, 'FontName','Times New Roman', 'FontSize',9, 'Box','on');

%% ---- Layout Adjustment for Colorbar ---- 
drawnow;
ax = gca;
ax.Position = [0.10 0.12 0.70 0.80];

colormap(cmap);
c = colorbar;

% Colorbar position
cb_width = 0.018;
cb_gap   = 0.06;
c.Position = [ ...
    ax.Position(1) + ax.Position(3) + cb_gap, ...
    ax.Position(2), ...
    cb_width, ...
    ax.Position(4) ...
];

% ---- Colorbar settings ----
c.Limits = [0 1];
c.Ticks = linspace(0,1,nC);
c.TickLabels = compose('%d', Cab_all);
c.Label.String = '$C_{ab}\ (\mathrm{\mu g/cm^2})$';
c.Label.Interpreter = 'latex';
set(c, 'FontName', 'Times New Roman', 'FontSize', 9);

%% ---- Save ----
filename = 'Cab_Prospect_D.png';
exportgraphics(gcf, filename, 'Resolution', 600);
