%% VCRS with Ejector Pin
clearvars; clc;

% Input Parameters
refrigerant = 'R744';
T_evap_array = [275.15, 277.15, 279.15, 281.15, 283.15]; 
T_cond = 308.15;       
Q_target = 3517;      
eta_isen = 0.85;   
eta_n = 0.85;
PLR = 1.194;
dT_superheat = 0;
dT_subcooling = 0;      
omega = 0.68; 
P_gc = 110e5; 
COP_paper = [3.4, 3.55, 3.75, 3.95, 4.15];

n = length(T_evap_array);
COP_code = zeros(1, n);

for i = 1:n
    T_evap = T_evap_array(i);
    P_evap = py.CoolProp.CoolProp.PropsSI('P', 'T', T_evap, 'Q', 1, refrigerant);
    h4 = py.CoolProp.CoolProp.PropsSI('H', 'T', T_evap, 'Q', 1, refrigerant);
    P_sep = P_evap * PLR;
    h6 = py.CoolProp.CoolProp.PropsSI('H', 'P', P_sep, 'Q', 0, refrigerant);
    h7 = h6;
    m_evap = Q_target / (h4 - h7);
    m_comp = m_evap / omega; 
    P1 = P_sep;
    h1 = py.CoolProp.CoolProp.PropsSI('H', 'P', P1, 'Q', 1, refrigerant);
    s1 = py.CoolProp.CoolProp.PropsSI('S', 'P', P1, 'Q', 1, refrigerant);
    h2s = py.CoolProp.CoolProp.PropsSI('H', 'P', P_gc, 'S', s1, refrigerant);
    h2 = h1 + (h2s - h1) / eta_isen;
    W_total = m_comp * (h2 - h1);
    COP_code(i) = Q_target / W_total;
end
for i=1:n
    fprintf("COP:%.3f\n", COP_code(i))
end


%% Plotting

x_label_text = 'Evaporator Temperature (Kelvin)'; 
y_label_text = 'COP';                         

legend_label_1 = 'Research Paper Values';            
legend_label_2 = 'Code Generated Values';              

condenser_temp_text = 'Condenser Temperature = 308.15 K'; 
output_filename = 'COP_Validation.png';   

fig = figure('Color', 'w', 'Units', 'pixels', 'Position', [100, 100, 200, 200]);
t = tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ax1 = nexttile;
hold(ax1, 'on');

% Plot lines
p1 = plot(ax1, T_evap_array, COP_paper, '-ro','LineWidth', 1.5, 'MarkerSize', 7,'MarkerFaceColor', 'w');
    
p2 = plot(ax1, T_evap_array, COP_code, '-bo','LineWidth', 1.5, 'MarkerSize', 7,'MarkerFaceColor', 'w');

% Axis styling (Inward ticks, clean grid)
set(ax1, 'Box', 'on', 'Color', 'w', 'LineWidth', 1,'XColor', 'k', 'YColor', 'k','FontSize', 10, 'FontName', 'Arial','TickDir', 'in','GridColor', [0.8 0.8 0.8], 'GridAlpha', 1); 
grid(ax1, 'on');

% Axis labels
xlabel(ax1, x_label_text, 'FontSize', 11, 'FontWeight', 'bold', 'Color', 'k');
ylabel(ax1, y_label_text, 'FontSize', 11, 'FontWeight', 'bold', 'Color', 'k');

xlim(ax1, [min(T_evap_array) - 2, max(T_evap_array) + 2]);
ylim(ax1, 'padded'); 

lgd = legend([p1, p2], legend_label_1, legend_label_2, 'Location', 'northeast');
title(lgd, condenser_temp_text, 'FontSize', 9, 'FontWeight', 'normal');
set(lgd, 'FontSize', 9, 'TextColor', 'k', 'EdgeColor', [0.6 0.6 0.6], 'Color', 'w');
hold(ax1, 'off');

% Export at 600 DPI for high-quality journal submission
exportgraphics(fig, output_filename, 'Resolution', 600, 'BackgroundColor', 'white');