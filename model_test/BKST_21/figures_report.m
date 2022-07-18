
% makes figures for the report (pdf file)
%
% if variable filename_figures_report_suffix is defined, this code saves
% some figures to report/figures

close all

% color_y = [0, 0.4470, 0.7410]; % blue
% % color_o = [0.8500, 0.3250, 0.0980]; % red
% color_o = [255 178 102]./255; % orange

color_y = [35/250, 55/250, 59/250]; % metropolis beamer theme
color_o = [235/250, 129/250, 27/250]; % metropolis beamer theme

fig_pos = [0, 0, 1000, 871.33]; % bigger figure

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
    X = cell(9, 2);
    
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
    X = cell(9, 1);
    leg = {'Young', 'Old'};
    i_vec = 1;
end

for i = 1:length(i_vec)
    if length(i_vec)==1 || length(i_vec)>1 && i==1
        X{1, i_vec(i)} = M_h;
        X{2, i_vec(i)} = M_i;
        X{3, i_vec(i)} = M_fi;
        X{4, i_vec(i)} = M_fh;
        X{5, i_vec(i)} = M_s;
        X{6, i_vec(i)} = M_r;
        X{7, i_vec(i)} = M_d;
        X{8, i_vec(i)} = gdp;
        X{9, i_vec(i)} = Pi;
    else
        X{1, i_vec(i)} = S_fig.M_h;
        X{2, i_vec(i)} = S_fig.M_i;
        X{3, i_vec(i)} = S_fig.M_fi;
        X{4, i_vec(i)} = S_fig.M_fh;
        X{5, i_vec(i)} = S_fig.M_s;
        X{6, i_vec(i)} = S_fig.M_r;
        X{7, i_vec(i)} = S_fig.M_d;
        X{8, i_vec(i)} = S_fig.gdp;
        X{9, i_vec(i)} = S_fig.Pi;
    end
end

titles = {'$M_h$', '$M_i$', '$M_{fi}$', '$M_{fh}$', '$M_s$', '$M_r$', '$M_d$', 'GDP', '$\Pi$'};
linestyles = {'-', '--'};
colors = {color_y, color_o};

make_fig(fig_pos, X, t0, t1, linestyles, colors, titles, leg)

if exist('filename_figures_report_suffix', 'var')==1
    filename = strcat('report/figures/aggregates_', filename_figures_report_suffix);
    saveas(gcf, filename, 'epsc')
end

%-------------------------------------------------------------------------%
%                            Choices (healthy)                            %
%-------------------------------------------------------------------------%

if exist('S_fig', 'var')==1
    X = cell(7, 2);
    
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
    X = cell(7, 1);
    leg = {'Young', 'Old'};
    i_vec = 1;
end

for i = 1:length(i_vec)
    if length(i_vec)==1 || length(i_vec)>1 && i==1
        X{1, i_vec(i)} = c_h;
        X{2, i_vec(i)} = n_h;
        X{3, i_vec(i)} = v_h;
        X{4, i_vec(i)} = l_h;
        X{5, i_vec(i)} = x_h;
        X{6, i_vec(i)} = d_h;
        X{7, i_vec(i)} = v_h ./ (n_h + v_h);
    else
        X{1, i_vec(i)} = S_fig.c_h;
        X{2, i_vec(i)} = S_fig.n_h;
        X{3, i_vec(i)} = S_fig.v_h;
        X{4, i_vec(i)} = S_fig.l_h;
        X{5, i_vec(i)} = S_fig.x_h;
        X{6, i_vec(i)} = S_fig.d_h;
        X{7, i_vec(i)} = S_fig.v_h ./ (S_fig.n_h + S_fig.v_h);
    end
end

titles = {'$c_h$', '$n_h$', '$\nu_h$', '$\ell_h$', '$x_h$', '$d_h$', '$\nu_h/(n_h+\nu_h)$'};
linestyles = {'-', '--'};
colors = {color_y, color_o};

make_fig(fig_pos, X, t0, t1, linestyles, colors, titles, leg)

if exist('filename_figures_report_suffix', 'var')==1
    filename = strcat('report/figures/choices_h_', filename_figures_report_suffix);
    saveas(gcf, filename, 'epsc')
end

%-------------------------------------------------------------------------%
% time outside in hours
%-------------------------------------------------------------------------%

if exist('S_fig', 'var')==1
    A = 112 * (1 - (v_h + d_h));
    B = 112 * (1 - (S_fig.v_h + S_fig.d_h));

    figure
    hold on
    plot(A(1, :), '-', 'color', color_y, 'linewidth', 2)
    plot(A(2, :), '-', 'color', color_o, 'linewidth', 2)
    plot(B(1, :), '--', 'color', color_y, 'linewidth', 2)
    plot(B(2, :), '--', 'color', color_o, 'linewidth', 2)
    grid on
    title('Time outside (weekly hours)', 'interpreter', 'latex', 'fontsize', 14)
    
    leg = {'Young', 'Old', 'Young (benchmark)', 'Old (benchmark)'};
    legend(leg, ...
        'Orientation', 'horizontal', ...
        'Position', [0.5 0.0275 0 0], ... % [x y width height]
        'FontSize', 10)

    if exist('filename_figures_report_suffix', 'var')==1
        filename = strcat('report/figures/time_outside_', filename_figures_report_suffix);
        saveas(gcf, filename, 'epsc')
    end
end

%-------------------------------------------------------------------------%
%
%-------------------------------------------------------------------------%

figure
hold on
plot(lambda_p_h(1, :), 'color', color_y, 'linewidth', 2)
plot(lambda_p_h(2, :), 'color', color_o, 'linewidth', 2)
grid on
legend('Young', 'Old')
title('$\lambda_p(h)$', 'interpreter', 'latex')

if exist('filename_figures_report_suffix', 'var')==1
    filename = strcat('report/figures/lambda_p_h_', filename_figures_report_suffix);
    saveas(gcf, filename, 'epsc')
end

%-------------------------------------------------------------------------%
%                                delta_vec                                %
%-------------------------------------------------------------------------%

figure
hold on
plot(delta_vec(1, :), 'color', color_y, 'linewidth', 2)
plot(delta_vec(2, :), 'color', color_o, 'linewidth', 2)
grid on
legend('Young', 'Old')
title('$\delta$', 'interpreter', 'latex')

if exist('filename_figures_report_suffix', 'var')==1
    filename = strcat('report/figures/delta_', filename_figures_report_suffix);
    saveas(gcf, filename, 'epsc')
end

%-------------------------------------------------------------------------%
%                           Internal functions                            %
%-------------------------------------------------------------------------%

function make_fig(fig_pos, X, t0, t1, linestyles, colors, titles, leg)
    
    % starts figure
    figure
    set(gcf, 'position', fig_pos);
    
    % loops over subplots
    for i = 1:size(X, 1)
        subplot(3, 3, i)
        hold on
        
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

        % title and grid
        title(titles{i}, 'interpreter', 'latex', 'FontSize', 14)
        grid on
        
        % adds legend
        if i==1
            legend(leg, ...
                'Orientation', 'horizontal', ...
                'Position', [0.5 0.9875 0 0], ... % [x y width height]
                'FontSize', 10)
        end
    end
    
    
    
end

