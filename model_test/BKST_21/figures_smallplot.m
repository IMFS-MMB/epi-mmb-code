
% makes plots for slides
% 
% if filename_figures_smallplot_suffix is defined, saves one figure to
% folder 'figures'
% figure name is 'smallplot_' and a suffix given by variable
% filename_figures_smallplot_suffix

close all

% color_y = [0, 0.4470, 0.7410]; % blue
% color_o = [0.8500, 0.3250, 0.0980]; % red
% color_o = [255 178 102]./255; % orange

color_y = [35/250, 55/250, 59/250]; % metropolis beamer theme
color_o = [235/250, 129/250, 27/250]; % metropolis beamer theme

fig_pos = [10, 50, 800, 500];

if exist('t0_fig', 'var')==1 && exist('t1_fig', 'var')==1
    t0 = t0_fig;
    t1 = t1_fig;
else
    t0 = 1; % first t to plot
    t1 = T; % last t to plot
end

%-------------------------------------------------------------------------%
%                               Small plot                                %
%-------------------------------------------------------------------------%

if exist('S_fig', 'var')==1
    X = cell(4, 2);
    
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
    X = cell(4, 1);
    leg = {'Young', 'Old'};
    i_vec = 1;
end

for i = 1:length(i_vec)
    if length(i_vec)==1 || length(i_vec)>1 && i==1
        X{1, i_vec(i)} = M_i_all;
        X{2, i_vec(i)} = M_d;
        X{3, i_vec(i)} = gdp;
        X{4, i_vec(i)} = d_h + v_h;
    else
        X{1, i_vec(i)} = S_fig.M_i_all;
        X{2, i_vec(i)} = S_fig.M_d;
        X{3, i_vec(i)} = S_fig.gdp;
        X{4, i_vec(i)} = S_fig.d_h + S_fig.v_h;
    end
end

titles = {'Measure, all infected ($M_{fi} + M_{i} + M_{s}$)', 'Measure, deceased ($M_d$)', 'GDP', 'Time at home ($d_h + \nu_h$)'};
linestyles = {'-', '--'};
colors = {color_y, color_o};

if flag_epidemiological==1
    log_vec = [1, 1, 0, 0];
else
    log_vec = zeros(1, 4);
end

make_fig(fig_pos, X, t0, t1, linestyles, colors, titles, leg, log_vec)

if exist('filename_figures_smallplot_suffix', 'var')==1 && filename_figures_smallplot_suffix~=""
    filename = strcat('figures/smallplot_', filename_figures_smallplot_suffix);
    saveas(gcf, filename, 'epsc')
end

%-------------------------------------------------------------------------%
%                           Internal functions                            %
%-------------------------------------------------------------------------%

function make_fig(fig_pos, X, t0, t1, linestyles, colors, titles, leg, log_vec)
    
    % starts figure
    figure
    set(gcf, 'position', fig_pos);
    
    % loops over subplots
    for i = 1:size(X, 1)
        subplot(2, 2, i)
        hold on
        
        % some lines of code if we want to plot variable in logs
        if log_vec(i)==1
            ylim_min_vec = [];
            for j = 1:size(X, 2)
                ylim_min_vec = [ylim_min_vec; X{i, j}(:, 2)];
            end
            ylim_min = min(ylim_min_vec);
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
            newy = [0.00001, 0.0001, 0.001, 0.01, 0.1];
            set(gca,'YTick', newy);
        end
        
        % write "log scale" in the y label
        if log_vec(i)==1
            ylabel('log scale', 'interpreter', 'latex')
        end

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
    
    
    
end

