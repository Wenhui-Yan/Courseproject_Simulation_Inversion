% ********************************************************************************
%% Plot of Reflectance and Transmittance with varying Cw (PROSPECT-D)
%% Description:
%   Generates leaf reflectance and transmittance spectra using the PROSPECT-D
%   model for a range of equivalent water thickness (Cw).
% ********************************************************************************

clear; clc; close all;

%% ---- Parameter range ----
Cw_all = 0.004:0.004:0.032;     % Leaf water content range in cm
nC = length(Cw_all);

%% ---- Continuous colormap ----
nColor = 256;                      
cmap = parula(nColor);             

%% ---- Figure ----
figure('Color', 'w'); 
hold on;

%% ---- Plot each curve ----
for i = 1:nC
    Cw = Cw_all(i);

    % ---- Leaf biochemical parameters (fixed) ----
    N      = 1.2;
    Cab    = 60;
    Car    = 20;
    Anth   = 0;
    Cbrown = 0;
    Cm     = 0.009;

    % ---- Run PROSPECT-D ----
    LRT = prospect_DB(N, Cab, Car, Anth, Cbrown, Cw, Cm);

    wl = LRT(:,1);     % wavelength (nm)
    R  = LRT(:,2);     % reflectance
    T  = LRT(:,3);     % transmittance

    % ---- Map Cw to continuous colormap ----
    idx = round( (Cw - min(Cw_all)) / (max(Cw_all) - min(Cw_all)) * (nColor-1) ) + 1;
    col = cmap(idx,:);

    % ---- Plot reflectance (bottom) ----
    plot(wl, R, 'Color', col, 'LineWidth', 0.8);

    % ---- Plot transmittance (top, mirrored) ----
    plot(wl, 1 - T, 'Color', col, 'LineWidth', 0.8);
end

%% ---- Axes settings ----
x_margin = 0.02 * (2500 - 400);   % Leave 2% margin at both ends of x-axis
xlim([400 - x_margin, 2500 + x_margin]);
xticks([400 700 1000 1300 1600 1900 2200 2500]);

% Left y-axis (Reflectance)
yyaxis left
ylim([0 1]);
yticks(0:0.2:1);
ylabel('Reflectance');
set(gca, 'YColor', 'k');

% Right y-axis (Transmittance)
yyaxis right
ylim([0 1]);
yticks(0:0.2:1);
set(gca, 'YDir', 'reverse');    % Mirror effect
ylabel('Transmittance');
set(gca, 'YColor', 'k');

xlabel('Wavelength (nm)');

%% ---- Title ----
title('$\mathrm{Reflectance\ and\ Transmittance\ with\ varying\ } C_w\ \mathrm{(PROSPECT\!\!-\!D)}$', ...
      'Interpreter','latex','FontSize',9);     

set(gca, 'FontName','Times New Roman', 'FontSize',9, 'Box','on');

%% ---- Layout Adjustment for Colorbar ---- 
drawnow;
ax = gca;

% Adjust main axis to leave space for colorbar
ax.Position = [0.10 0.12 0.70 0.80];

% Add colorbar
colormap(cmap);
c = colorbar;

cb_width = 0.018;
cb_gap   = 0.06;

c.Position = [ ...
    ax.Position(1) + ax.Position(3) + cb_gap, ...
    ax.Position(2), ...
    cb_width, ...
    ax.Position(4) ...
];

% ---- Colorbar styling ----
c.Limits       = [0 1];
c.Ticks        = linspace(0,1,nC);
c.TickLabels   = compose('%.3f', Cw_all);
c.Label.String = '$C_w\ (\mathrm{cm})$';
c.Label.Interpreter = 'latex';
set(c, 'FontName', 'Times New Roman', 'FontSize', 9);

%% ---- Save figure ---- 
filename = 'Cw_Prospect_D.png'; 
exportgraphics(gcf, filename, 'Resolution', 600);
