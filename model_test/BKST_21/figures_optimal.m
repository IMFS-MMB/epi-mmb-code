


% C_fig includes lambda_p (age specific and universal), and time outside
% variables (age specific, universal, partial equilibrium, and benchmark)
% C_fig{1, 1}: lambda_p age specific
% C_fig{1, 2}: lambda_p universal
% C_fig{2, 1}: time outside (age specific)
% C_fig{2, 2}: time outside (universal)
% C_fig{2, 3}: time outside (partial equilibrium)
% C_fig{2, 4}: time outside (benchmark, or some other scenario)

%% figure preliminaries

close all

% figure preliminaries
clear pos_vec

fig_pos = [10, 50, 920, 375];

figure
set(gcf, 'position', fig_pos)

for i = 1:2
    handle = subplot(1, 2, i);
    pos_vec(i, :) = get(handle, 'position');    
end

close

% adjusts some properties of subplot
pos_vec(:, 2) = pos_vec(:, 2) * 1.75;
pos_vec(:, 4) = pos_vec(:, 4) * 0.85; % height

color_y = [0, 0.4470, 0.7410]; % blue
color_o = [255 178 102]./255; % orange

%% figure 1

if all(size(C_fig) == [2, 4])
    figure
    set(gcf, 'position', fig_pos)

    % subplot 1: lambdas
    subplot('position', pos_vec(1, :))
    grid on
    hold on
    plot(C_fig{1, 1}(1, :), 'linewidth', 2, 'color', color_y)
    plot(C_fig{1, 1}(2, :), 'linewidth', 2, 'color', color_o)
    title('$ \lambda_p $', 'interpreter', 'latex', 'fontsize', 14)

    % subplot 2: time outside
    subplot('position', pos_vec(2, :))
    grid on
    hold on
    plot(NaN, 'linewidth', 10, 'color', color_y)
    plot(NaN, 'linewidth', 10, 'color', color_o)
    plot(C_fig{2, 1}(1, :), 'linewidth', 2, 'color', color_y)
    plot(C_fig{2, 1}(2, :), 'linewidth', 2, 'color', color_o)
    plot(C_fig{2, 4}(1, :), '--', 'linewidth', 2, 'color', color_y)
    plot(C_fig{2, 4}(2, :), '--', 'linewidth', 2, 'color', color_o)
    plot(C_fig{2, 3}(1, :), ':', 'linewidth', 2, 'color', color_y, 'MarkerSize', 5)
    plot(C_fig{2, 3}(2, :), ':', 'linewidth', 2, 'color', color_o, 'MarkerSize', 5)
    myplot = findobj('Type', 'line');
    title('Time outside (weekly hours)', 'interpreter', 'latex', 'fontsize', 14)

    leg = {'Young', 'Old', 'Optimal', 'Benchmark', 'Partial equil.'};

    legend([myplot(8) myplot(7) myplot(6) myplot(4) myplot(2)], leg, ...
        'Orientation', 'horizontal', ...
        'Position', [0.5 0.05 0 0], ... % [x y width height]
        'FontSize', 12, 'interpreter', 'latex')

    % saves
    if exist('filename_figures_optimal_suffix', 'var')==1 && filename_figures_optimal_suffix~=""
        filename = strcat('figures/optimal1_', filename_figures_optimal_suffix);
        saveas(gcf, filename, 'epsc')
    end

end

%% figure 1b (without partial equilibrium)

if all(size(C_fig) == [2, 4])
    figure
    set(gcf, 'position', fig_pos)

    % subplot 1: lambdas
    subplot('position', pos_vec(1, :))
    grid on
    hold on
    plot(C_fig{1, 1}(1, :), 'linewidth', 2, 'color', color_y)
    plot(C_fig{1, 1}(2, :), 'linewidth', 2, 'color', color_o)
    title('$ \lambda_p $', 'interpreter', 'latex', 'fontsize', 14)

    % subplot 2: time outside
    subplot('position', pos_vec(2, :))
    grid on
    hold on
    plot(NaN, 'linewidth', 10, 'color', color_y)
    plot(NaN, 'linewidth', 10, 'color', color_o)
    plot(C_fig{2, 1}(1, :), 'linewidth', 2, 'color', color_y)
    plot(C_fig{2, 1}(2, :), 'linewidth', 2, 'color', color_o)
    plot(C_fig{2, 4}(1, :), '--', 'linewidth', 2, 'color', color_y)
    plot(C_fig{2, 4}(2, :), '--', 'linewidth', 2, 'color', color_o)
%     plot(C_fig{2, 3}(1, :), ':', 'linewidth', 2, 'color', color_y, 'MarkerSize', 5)
%     plot(C_fig{2, 3}(2, :), ':', 'linewidth', 2, 'color', color_o, 'MarkerSize', 5)
    myplot = findobj('Type', 'line');
    title('Time outside (weekly hours)', 'interpreter', 'latex', 'fontsize', 14)

    leg = {'Young', 'Old', 'Optimal', 'Benchmark'};

    legend([myplot(6) myplot(5) myplot(4) myplot(2)], leg, ...
        'Orientation', 'horizontal', ...
        'Position', [0.5 0.05 0 0], ... % [x y width height]
        'FontSize', 12, 'interpreter', 'latex')

    % saves
    if exist('filename_figures_optimal_suffix', 'var')==1 && filename_figures_optimal_suffix~=""
        filename = strcat('figures/optimal1b_', filename_figures_optimal_suffix);
        saveas(gcf, filename, 'epsc')
    end

end

%% figure 2

if all(size(C_fig) == [2, 4])
    figure
    set(gcf, 'position', fig_pos)

    % subplot 1: lambdas
    subplot('position', pos_vec(1, :))
    grid on
    hold on
    plot(C_fig{1, 1}(1, :), 'linewidth', 2, 'color', color_y)
    plot(C_fig{1, 1}(2, :), 'linewidth', 2, 'color', color_o)
    plot(C_fig{1, 2}(1, :), 'linewidth', 2, 'color', [0.5 0.5 0.5])
    title('$ \lambda_p $', 'interpreter', 'latex', 'fontsize', 14)

    % subplot 2: time outside
    subplot('position', pos_vec(2, :))
    grid on
    hold on
    plot(NaN, 'linewidth', 10, 'color', color_y)
    plot(NaN, 'linewidth', 10, 'color', color_o)
    plot(C_fig{2, 1}(1, :), 'linewidth', 2, 'color', color_y)
    plot(C_fig{2, 1}(2, :), 'linewidth', 2, 'color', color_o)
    plot(C_fig{2, 2}(1, :), '--', 'linewidth', 2, 'color', color_y)
    plot(C_fig{2, 2}(2, :), '--', 'linewidth', 2, 'color', color_o)
    myplot = findobj('Type', 'line');
    title('Time outside (weekly hours)', 'interpreter', 'latex', 'fontsize', 14)

    leg = {'Young', 'Old', 'Age-specific optimal', 'Universal optimal'};

    legend([myplot(6) myplot(5) myplot(4) myplot(2)], leg, ...
        'Orientation', 'horizontal', ...
        'Position', [0.5 0.05 0 0], ... % [x y width height]
        'FontSize', 12, 'interpreter', 'latex')

    % saves
    if exist('filename_figures_optimal_suffix', 'var')==1 && filename_figures_optimal_suffix~=""
        filename = strcat('figures/optimal2_', filename_figures_optimal_suffix);
        saveas(gcf, filename, 'epsc')
    end
end

%% figure 3

if all(size(C_fig) == [2, 2])

    figure
    set(gcf, 'position', fig_pos)

    % subplot 1: lambdas
    subplot('position', pos_vec(1, :))
    grid on
    hold on
    plot(C_fig{1, 1}(1, :), 'linewidth', 2, 'color', color_y)
    plot(C_fig{1, 1}(2, :), 'linewidth', 2, 'color', color_o)
    plot(C_fig{1, 2}(1, :), '--', 'linewidth', 2, 'color', color_y)
    plot(C_fig{1, 2}(2, :), '--', 'linewidth', 2, 'color', color_o)
    title('$ \lambda_p $', 'interpreter', 'latex', 'fontsize', 14)

    % subplot 2: time outside
    subplot('position', pos_vec(2, :))
    grid on
    hold on
    plot(NaN, 'linewidth', 10, 'color', color_y)
    plot(NaN, 'linewidth', 10, 'color', color_o)
    plot(C_fig{2, 1}(1, :), 'linewidth', 2, 'color', color_y)
    plot(C_fig{2, 1}(2, :), 'linewidth', 2, 'color', color_o)
    plot(C_fig{2, 2}(1, :), '--', 'linewidth', 2, 'color', color_y)
    plot(C_fig{2, 2}(2, :), '--', 'linewidth', 2, 'color', color_o)
    myplot = findobj('Type', 'line');
    title('Time outside (weekly hours)', 'interpreter', 'latex', 'fontsize', 14)

    leg = {'Young', 'Old', 'Optimal (with quarantine)', 'Optimal (without quarantine)'};

    legend([myplot(6) myplot(5) myplot(4) myplot(2)], leg, ...
        'Orientation', 'horizontal', ...
        'Position', [0.5 0.05 0 0], ... % [x y width height]
        'FontSize', 12, 'interpreter', 'latex')

    % saves
    if exist('filename_figures_optimal_suffix', 'var')==1 && filename_figures_optimal_suffix~=""
        filename = strcat('figures/optimal3_', filename_figures_optimal_suffix);
        saveas(gcf, filename, 'epsc')
    end

end















