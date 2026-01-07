% ********************************************************************************
%% Plot of Reflectance and Transmittance with varying Cw (SIP-Leaf Model)
%% Author: WenhuiYan
%% Affiliation: CAU
%% Email: wenhuiyan233@163.com
% ********************************************************************************

clear; clc; close all;

%% ---- Parameter range ----
Cw_all = 0.004:0.004:0.032;     % Leaf water content (Cw) range in cm
nC = length(Cw_all);

%% ---- Continuous colormap ----
nColor = 256;                      
cmap = parula(nColor);             

%% ---- Figure ----
figure('Color', 'w'); 
hold on;

%% ---- Fixed leaf parameters ----
Cab    = 60;
Car    = 20;
Anth   = 0;
Cbrown = 0;
Cm     = 0.009;

%% ---- Loop over Cw ----
for i = 1:nC
    Cw = Cw_all(i);

    % ---- SIP-Leaf model ----
    % LRT columns:
    % 1 = Wavelength (nm)
    % 2 = Single scattering albedo
    % 3 = Reflectance
    % 4 = Transmittance
    LRT = SIP_Model(Cab, Car, Anth, Cbrown, Cw, Cm);

    wl = LRT(:,1);
    R  = LRT(:,3);
    T  = LRT(:,4);

    % ---- Color mapping ----
    idx = round( (Cw - min(Cw_all)) / (max(Cw_all) - min(Cw_all)) * (nColor-1) ) + 1;
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
ylabel('Leaf reflectance');
set(gca, 'YColor', 'k');

% Right y-axis
yyaxis right
ylim([0 1]);
set(gca, 'YDir','reverse');
yticks(0:0.2:1);
ylabel('Leaf transmittance');
set(gca, 'YColor', 'k');

xlabel('Wavelength (nm)');

%% ---- Title ----
title('$\mathrm{Reflectance\ and\ Transmittance\ with\ varying\ } C_{w}\ \mathrm{(SIP\!-\!Leaf)}$', ...
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

%% ---- Colorbar settings ----
c.Limits = [0 1];
c.Ticks = linspace(0,1,nC);
c.TickLabels = compose('%.3f', Cw_all);
c.Label.String = '$C_{w}\ (\mathrm{cm})$';
c.Label.Interpreter = 'latex';
set(c, 'FontName', 'Times New Roman', 'FontSize', 9);

%% ---- Save ----
filename = 'Cw_SIP_Leaf.png';
exportgraphics(gcf, filename, 'Resolution', 600);
