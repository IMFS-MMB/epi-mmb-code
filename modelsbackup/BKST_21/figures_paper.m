
% if filename_figures_paper_suffix is defined, saves one figure to
% folder 'figures'
% figure name is 'aggregates_' (or 'choices_h_') and a suffix given by
% variable filename_figures_paper_suffix

close all

color_y = [0, 0.4470, 0.7410]; % blue
% color_o = [0.8500, 0.3250, 0.0980]; % red
color_o = [255 178 102]./255; % orange


% color_y = [35/250, 55/250, 59/250]; % metropolis beamer theme
% color_o = [235/250, 129/250, 27/250]; % metropolis beamer theme

if exist('t0_fig', 'var')==1 && exist('t1_fig', 'var')==1
    t0 = t0_fig;
    t1 = t1_fig;
else
    t0 = 1; % first t to plot
    t1 = T; % last t to plot
end

%-------------------------------------------------------------------------%
%                               Aggregates                                %
%-------------------------------------------------------------------------%

if exist('S_fig', 'var')==1
    X = cell(6, 2);
    
    if flag_epidemiological==1
        i_vec = [2, 1];
        leg = {'Young', 'Old', 'Young (epidem.)', 'Old (epidem.)'};
    else
        i_vec = [1, 2];
        
        if flag_Pi==1
            leg = {'Young', 'Old', 'Young (optimal policy)', 'Old (optimal policy)'};
        else
            leg = {'Young', 'Old', 'Young (benchmark)', 'Old (benchmark)'};
        end
    end
else
    X = cell(6, 1);
    leg = {'Young', 'Old'};
    i_vec = 1;
end

for i = 1:length(i_vec)
    if length(i_vec)==1 || length(i_vec)>1 && i==1
        X{1, i_vec(i)} = M_fi;
        X{2, i_vec(i)} = M_i;
        X{3, i_vec(i)} = M_s;
        X{4, i_vec(i)} = M_i_all;
        X{5, i_vec(i)} = M_d;
        X{6, i_vec(i)} = gdp;
    else
        X{1, i_vec(i)} = S_fig.M_fi;
        X{2, i_vec(i)} = S_fig.M_i;
        X{3, i_vec(i)} = S_fig.M_s;
        X{4, i_vec(i)} = S_fig.M_i_all;
        X{5, i_vec(i)} = S_fig.M_d;
        X{6, i_vec(i)} = S_fig.gdp;
    end
end

fig_pos = [10, 50, 1280, 600];
% titles = {'Measure, fever-infected ($M_{fi}$)', 'Measure, infected ($M_i$)', ...
%     'Measure, serious symptoms ($M_s$)', 'Measure, all infected ($M_{fi} + M_i + M_s$)', ...
%     'Measure, deceased ($M_d$)', 'GDP'};
titles = {{'Measure, fever-infected', '$M_t(f_i,a)$'}, {'Measure, infected', '$M_t(i,a)$'}, ...
    {'Measure, hospitalized', '$M_t(h,a)$'}, {'Measure, all infected', '$M_t(f_i,a) + M_t(i,a) + M_t(h,a)$'}, ...
    {'Measure, deceased ', '$M_t(d,a)$'}, 'GDP'};
linestyles = {'-', '--'};
colors = {color_y, color_o};

if flag_epidemiological==1
    log_vec = [1, 1, 1, 1, 1, 0];
else
    log_vec = zeros(1, 6);
end

% make_fig(fig_pos, X, t0, t1, linestyles, colors, titles, leg, log_vec)

% preliminaries
clear pos_vec

figure
set(gcf, 'position', fig_pos)

for i = 1:6
    handle = subplot(2, 3, i);
    pos_vec(i, :) = get(handle, 'position');
end

close

% adjusts some properties of subplot
pos_vec(1:3, 2) = pos_vec(1:3, 2) * 1;
pos_vec(4:6, 2) = pos_vec(4:6, 2) * 1.075;
pos_vec(:, 4) = pos_vec(:, 4) * 0.9; % height

% starts figure
figure
set(gcf, 'position', fig_pos);

% loops over subplots
for i = 1:size(X, 1)
%     subplot(n_rows, n_cols, i)
    subplot('position', pos_vec(i, :))
    hold on

    % some lines of code if we want to plot variable in logs
    if log_vec(i)==1
        if contains(titles{i}, 'deceased')==1
            ylim_min_vec = [];
            for j = 1:size(X, 2)
                ylim_min_vec = [ylim_min_vec; X{i, j}(:)];
            end
            ylim_min = min(ylim_min_vec(ylim_min_vec > 0));
        else
            ylim_min_vec = [];
            for j = 1:size(X, 2)
                ylim_min_vec = [ylim_min_vec; X{i, j}(:, 2)];
            end
            ylim_min = min(ylim_min_vec);

            if ylim_min==0
                log_vec(i) = 0;
            end
        end
    end

    % loops over benchmark/counterfactual/epidemiological
    for j = 1:size(X, 2)
        Y = X{i, j};

        % loops over ages
        for g = 1:size(Y, 1)
            Y0 = Y(g, :);
            t1m = min(t1, length(Y0));
            plot(Y0(1, t0:t1m), linestyles{j}, 'linewidth', 2, ...
                'color', colors{g})
        end
    end

    % changes y limits if variable is in logs
    if log_vec(i)==1
        set(gca, 'YScale', 'log')
        ylim([ylim_min, inf])
        newy = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1];
        set(gca,'YTick', newy);
    end

    % write "log scale" in the y label
    if log_vec(i)==1
        ylabel('log scale', 'interpreter', 'latex')
    end

    xlabel('Weeks', 'interpreter', 'latex')

    % title and grid
    title(titles{i}, 'interpreter', 'latex', 'FontSize', 14)
    grid on

    % adds legend
    if i==1
        legend(leg, ...
            'Orientation', 'horizontal', ...
            'Position', [0.5 0.025 0 0], ... % [x y width height]
            'FontSize', 12, 'interpreter', 'latex')
    end
end

if exist('filename_figures_paper_suffix', 'var')==1 && filename_figures_paper_suffix~=""
    filename = strcat('figures/aggregates_', filename_figures_paper_suffix);
    saveas(gcf, filename, 'epsc')
end

%-------------------------------------------------------------------------%
%                            Choices (healthy)                            %
%-------------------------------------------------------------------------%

if exist('S_fig', 'var')==1
    X = cell(3, 2);
    
    if flag_epidemiological==1
        i_vec = [2, 1];
        leg = {'Young', 'Old', 'Young (epidem.)', 'Old (epidem.)'};
    else
        i_vec = [1, 2];
        
        if flag_Pi==1
            leg = {'Young', 'Old', 'Young (optimal policy)', 'Old (optimal policy)'};
        else
            leg = {'Young', 'Old', 'Young (benchmark)', 'Old (benchmark)'};
        end
    end
else
    X = cell(6, 1);
    leg = {'Young', 'Old'};
    i_vec = 1;
end

for i = 1:length(i_vec)
    if length(i_vec)==1 || length(i_vec)>1 && i==1
        X{1, i_vec(i)} = 112 * (d_h + v_h);
        X{2, i_vec(i)} = 112 * n_h;
        X{3, i_vec(i)} = 112 * l_h;
    else
        X{1, i_vec(i)} = 112 * (S_fig.d_h + S_fig.v_h);
        X{2, i_vec(i)} = 112 * (S_fig.n_h);
        X{3, i_vec(i)} = 112 * (S_fig.l_h);
    end
end

% preliminaries
clear pos_vec

fig_pos = [10, 50, 1280, 325];

figure
set(gcf, 'position', fig_pos)

for i = 1:3
    handle = subplot(1, 3, i);
    pos_vec(i, :) = get(handle, 'position');    
end

close

% adjusts some properties of subplot
pos_vec(:, 2) = pos_vec(:, 2) * 2.15;
pos_vec(:, 4) = pos_vec(:, 4) * 0.8; % height

% makes figures
linestyles = {'-', '--'};
colors = {color_y, color_o};
titles = {'Time at home, $d_t(s,a) + v_t(s,a)$', 'Time at work, $n_t(s,a)$', 'Leisure outside, $\ell_t(s,a)$'};

figure
set(gcf, 'Position', fig_pos);

for i = 1:3 % 3 choices
    subplot('position', pos_vec(i, :))
    hold on
    
    for k = 1:size(X, 2) % scenarios
        for j = 1:2 % age
            t1m = min(t1, length(X{i, k}(j, :)));
            plot(X{i, k}(j, t0:t1m), linestyles{k}, 'color', colors{j}, 'LineWidth', 2)
        end
    end
    
    title(titles{i}, 'interpreter', 'latex', 'FontSize', 14)
    
    if i==1
        ylabel('Weekly hours', 'interpreter', 'latex')
    end
    
    grid on
    
    xlabel('Weeks', 'interpreter', 'latex')
    
    if i==1
        legend(leg, ...
            'Orientation', 'horizontal', ...
            'Position', [0.5 0.06 0 0], ... % [x y width height]
            'FontSize', 12, 'interpreter', 'latex')
    end
end

% saves
if exist('filename_figures_paper_suffix', 'var')==1 && filename_figures_paper_suffix~=""
    filename = strcat('figures/choices_h_', filename_figures_paper_suffix);
    saveas(gcf, filename, 'epsc')
end

%-------------------------------------------------------------------------%
%                             M_i_all and GDP                             %
%-------------------------------------------------------------------------%

if exist('S_fig', 'var')==1
    X = cell(2, 2);
    
    if flag_epidemiological==1
        i_vec = [2, 1];
        leg = {'Young', 'Old', 'Young (epidem.)', 'Old (epidem.)'};
    else
        i_vec = [1, 2];
        
        if flag_Pi==1
            leg = {'Young', 'Old', 'Young (optimal policy)', 'Old (optimal policy)'};
        else
            leg = {'Young', 'Old', 'Young (benchmark)', 'Old (benchmark)'};
        end
    end
else
    X = cell(2, 1);
    leg = {'Young', 'Old'};
    i_vec = 1;
end

for i = 1:length(i_vec)
    if length(i_vec)==1 || length(i_vec)>1 && i==1
        X{1, i_vec(i)} = M_i_all;
        X{2, i_vec(i)} = gdp;
    else
        X{1, i_vec(i)} = S_fig.M_i_all;
        X{2, i_vec(i)} = S_fig.gdp;
    end
end

% preliminaries
clear pos_vec

fig_pos = [10, 50, 800, 350];

figure
set(gcf, 'position', fig_pos)

for i = 1:2
    handle = subplot(1, 2, i);
    pos_vec(i, :) = get(handle, 'position');    
end

close

% adjusts some properties of subplot
pos_vec(:, 2) = pos_vec(:, 2) * 2;
pos_vec(:, 4) = pos_vec(:, 4) * 0.8; % height

% makes figures
linestyles = {'-', '--'};
colors = {color_y, color_o};
% titles = {'Measure, all infected ($M_{fi} + M_i + M_s$)', 'GDP'};
titles = {{'Measure, all infected ', '$M_t(f_i,a) + M_t(i,a) + M_t(h,a)$'}, 'GDP'};

figure
set(gcf, 'Position', fig_pos);

for i = 1:2 % 2 plots
    subplot('position', pos_vec(i, :))
    hold on
    
    for k = 1:size(X, 2) % scenarios
        N = size(X{i}, 1);
        for j = 1:N % age
            t1m = min(t1, length(X{i, k}(j, :)));
            plot(X{i, k}(j, t0:t1m), linestyles{k}, 'color', colors{j}, 'LineWidth', 2)
        end
    end
    
    title(titles{i}, 'interpreter', 'latex', 'FontSize', 14)
    
    xlabel('Weeks', 'interpreter', 'latex')
    
    grid on
    
    if i==1
        legend(leg, ...
            'Orientation', 'horizontal', ...
            'Position', [0.5 0.06 0 0], ... % [x y width height]
            'FontSize', 12, 'interpreter', 'latex')
    end
end

% saves
if exist('filename_figures_paper_suffix', 'var')==1 && filename_figures_paper_suffix~=""
    filename = strcat('figures/M_i_all_gdp_', filename_figures_paper_suffix);
    saveas(gcf, filename, 'epsc')
end


