% this code uses the variable filename_results_table
% uses the following scripts to make figures:
% figures_report, figures_smallplot, figures_paper

% clears some variables
%clear S_table S_Pi S_fig

lambda_p_i = zeros(n_age, T);
lambda_p_h = zeros(n_age, T);
lambda_p_r = zeros(n_age, T);
lambda_p_f = zeros(n_age, T);

xi_p = zeros(n_age, 1);

descs = [];

flag_append = 0;

%-------------------------------------------------------------------------%
%                                Benchmark                                %
%-------------------------------------------------------------------------%

disp('=========')
disp('Benchmark')
disp('=========')

% runs and saves results to csv
equilibrium
results_table

% Susceptibles = frac_young.*(M_h(i_young,:)+M_fh(i_young,:))+frac_old.*(M_h(i_old,:)+M_fh(i_old,:)); 
% Susceptibles=Susceptibles';
% Infected = frac_young.*(M_i(i_young,:)+M_fi(i_young,:)+M_s(i_young,:))+frac_old.*(M_i(i_old,:)+M_fi(i_old,:)++M_s(i_old,:));
% Infected=Infected';
% %Infected = frac_young.*N_c_s(i_young,:)+frac_old.*N_c_s(i_old,:); %infected stock
% Recovered = frac_young.*(M_r(i_young,:))+frac_old.*(M_r(i_old,:)); 
% Recovered=Recovered';
% Deaths = frac_young.*(M_d(i_young,:))+frac_old.*(M_d(i_old,:)); %note: Covid deaths only
% Deaths=Deaths';
% %Deaths = frac_young.*(M_d(i_young,:)+M_dn(i_young,:))+frac_old.*(M_d(i_old,:)+M_dn(i_old,:)); %note: Covid deaths and natural deaths
% 
% Consumption = frac_young.*(c_r(i_young,:).*M_r(i_young,:)+c_h(i_young,:).*M_h(i_young,:)+c_i(i_young,:).*(M_i(i_young,:)+M_s(i_young,:))+c_f(i_young).*(M_fh(i_young,:)+M_fi(i_young,:)))...
%     +frac_old.*(c_r(i_old,:).*M_r(i_old,:)+c_h(i_old,:).*M_h(i_old,:)+c_i(i_old,:).*(M_i(i_old,:)+M_s(i_old,:))+c_f(i_old).*(M_fh(i_old,:)+M_fi(i_old,:)));
% Consumption=Consumption';
% Labour = frac_young.*((n_r(i_young,:)+v_r(i_young,:)).*M_r(i_young,:)+(n_h(i_young,:)+v_h(i_young,:)).*M_h(i_young,:)+(n_i(i_young,:)+v_i(i_young,:)).*(M_i(i_young,:)+M_s(i_young,:))+(n_f(i_young,:)+v_f(i_young,:)).*(M_fh(i_young,:)+M_fi(i_young,:)));
% Labour=Labour';
% Output = gdp';
% 
% Consumption_ss = frac_young.*c_r(i_young,1)+frac_old.*c_r(i_old,1);
% Labour_ss = frac_young.*(n_r(i_young,1)+v_r(i_young,1));
% Output_ss = w*(n_r(i_young,1)+tau_r(i_young,1));
% Susceptibles_ss=0;
% Infected_ss=0;
% Recovered_ss=0;
% Deaths_ss=0;
% 
% save('simulated_results.mat','Consumption','Labour','Output','Deaths','Susceptibles','Infected','Recovered');
% save('simulated_results_ss.mat','Consumption_ss','Labour_ss','Output_ss','Deaths_ss','Susceptibles_ss','Infected_ss','Recovered_ss');
% 



%{
% saves some variables of benchmark
S_table = struct('M_d', M_d, 'gdp', gdp, 'gdp_pc', gdp_pc, 't_peak_I', t_peak_I);
S_Pi = struct('Pi', Pi, 'I', I, 'M_s', M_s);
S_fig = struct('M_h', M_h, 'M_i', M_i, 'M_fi', M_fi, 'M_fh', M_fh, 'M_s', M_s, 'M_r', M_r, 'M_d', M_d, 'M_dn', M_dn, 'Pi', Pi, 'c_h', c_h, 'n_h', n_h, 'v_h', v_h, 'l_h', l_h, 'x_h', x_h, 'd_h', d_h, 'c_f', c_f, 'n_f', n_f, 'v_f', v_f, 'l_f', l_f, 'x_f', x_f, 'd_f', d_f, 'c_i', c_i, 'n_i', n_i, 'v_i', v_i, 'l_i', l_i, 'x_i', x_i, 'd_i', d_i, 'gdp', gdp, 'M_i_all', M_i_all);
S_c_eq = struct('c_h_term', c_h_term, 'x_h_term', x_h_term, 'n_h_term', n_h_term, 'v_h_term', v_h_term, 'l_h_term', l_h_term, 'd_h_term', d_h_term, 'flag_fake_young', flag_fake_young, 'c_i_term', c_i_term, 'x_i_term', x_i_term, 'n_i_term', n_i_term, 'v_i_term', v_i_term, 'l_i_term', l_i_term, 'd_i_term', d_i_term, 'Pi', Pi, 'c_r', c_r, 'x_r', x_r, 'n_r', n_r, 'v_r', v_r, 'l_r', l_r, 'd_r', d_r, 'c_i', c_i, 'x_i', x_i, 'n_i', n_i, 'v_i', v_i, 'l_i', l_i, 'd_i', d_i, 'c_h', c_h, 'x_h', x_h, 'n_h', n_h, 'v_h', v_h, 'l_h', l_h, 'd_h', d_h, 'c_f', c_f, 'x_f', x_f, 'n_f', n_f, 'v_f', v_f, 'l_f', l_f, 'd_f', d_f, 'delta_vec', delta_vec, 'xi_p', xi_p);
descs = [descs, "Benchmark"];
%-------------------------------------------------------------------------%
%                             Epidemiological                             %
%-------------------------------------------------------------------------%
disp('=======')
disp('Epidem.')
disp('=======')
flag_epidemiological = 1;
% runs and saves results to csv
equilibrium
results_table
% figures
suffix = strcat(filename_results_table, '_benchmark');
filename_figures_paper_suffix = suffix;
filename_figures_report_suffix = suffix;
filename_figures_smallplot_suffix = suffix;
figures_paper
figures_report
figures_smallplot
% turns off flag epidemiological
flag_epidemiological = 0;
descs = [descs, "Epidem."];
%-------------------------------------------------------------------------%
%                             Age ext. partial                            %
%-------------------------------------------------------------------------%
disp('================')
disp('Age ext. partial')
disp('================')
flag_Pi = 1;
flag_fake_young = 1;
% runs and saves results to csv
equilibrium
results_table
flag_Pi = 0;
flag_fake_young = 0;
descs = [descs, "Age ext. partial"];
%-------------------------------------------------------------------------%
%                             Age ext. general                            %
%-------------------------------------------------------------------------%
disp('================')
disp('Age ext. general')
disp('================')
flag_fake_young = 1;
% runs and saves results to csv
equilibrium
results_table
flag_fake_young = 0;
descs = [descs, "Age ext. general"];
%-------------------------------------------------------------------------%
%                             Selective mixing                            %
%-------------------------------------------------------------------------%
disp('================')
disp('Selective mixing')
disp('================')
% sets selective mixing
zeta = zeta_data;
% runs and saves results to csv
equilibrium
results_table
% figures
suffix = strcat(filename_results_table, '_selmix');
filename_figures_paper_suffix = suffix;
filename_figures_report_suffix = suffix;
filename_figures_smallplot_suffix = suffix;
figures_paper
figures_report
figures_smallplot
% back without selective mixing
zeta = 0;
descs = [descs, "Sel. mix."];
%-------------------------------------------------------------------------%
%                                 Testing                                 %
%-------------------------------------------------------------------------%
disp('=======')
disp('Testing')
disp('=======')
testing_intensities = [0.25, 0.5, 0.75, 1];
for i_testing_intensity = 1:length(testing_intensities)
    testing_intensity = testing_intensities(i_testing_intensity);
    
    for i_who = 1:3
        if i_who==1 % all
            xi_p(:) = testing_intensity;
            who_string = "a";
        elseif i_who==2 % young
            xi_p(i_young) = testing_intensity;
            xi_p(i_old) = 0;
            who_string = "y";
        else % old
            xi_p(i_young) = 0;
            xi_p(i_old) = testing_intensity;
            who_string = "o";
        end
        desc = strcat("T-", who_string, "-", num2str(fix(100*testing_intensity)));
        
%         % uncomment this if I want to skip some simulations
%         I_want = ["T-a-100"];
%         if ismember(desc, I_want)==0
%             continue
%         end
        
        disp(desc)
        
        equilibrium
        results_table
        
        % figures
        if desc=="T-a-100"
            suffix = strcat(filename_results_table, '_', desc);
            filename_figures_paper_suffix = suffix;
            filename_figures_report_suffix = suffix;
            filename_figures_smallplot_suffix = suffix;
            figures_paper
            figures_report
            figures_smallplot
        end
        
        descs = [descs, desc];
    end
end
%-------------------------------------------------------------------------%
%                                Quarantine                               %
%-------------------------------------------------------------------------%
disp('==========')
disp('Quarantine')
disp('==========')
if alpha1==0
    new_lambdas = [8.13814414450348, 19.8648352905615]; % 75%, 90% increase in d (calibration without teleworking)
else
    new_lambdas = [2.45096990297207, 32.3598330005704]; % 75%, 90% increase in d+v (calibration with teleworking)
end
testing_intensities = [0.5, 1]; % 50%, 100% testing
descs_new_lambda = ["75", "90"];
for i_new_lambda = 1:length(new_lambdas)
    new_lambda = new_lambdas(i_new_lambda);
    desc_new_lambda = descs_new_lambda(i_new_lambda);
    
    lambda_p_i(:,:) = max(new_lambda - lambda_i, 0);
    lambda_p_h(:,:) = 0;
    lambda_p_r(:,:) = 0;
    lambda_p_f(:,:) = 0;
    
    for i_testing_intensity = 1:length(testing_intensities)
        testing_intensity = testing_intensities(i_testing_intensity);
        
        for i_who = 1:3
            if i_who==1 % all
                xi_p(:) = testing_intensity;
                desc_who = "a";
            elseif i_who==2 % young
                xi_p(i_young) = testing_intensity;
                xi_p(i_old) = 0;
                desc_who = "y";
            else % old
                xi_p(i_young) = 0;
                xi_p(i_old) = testing_intensity;
                desc_who = "o";
            end
            
            desc = strcat("Q", desc_new_lambda, "-", desc_who, "-", num2str(floor(100*testing_intensity)), "t");
            
%             % uncomment this if I want to skip some simulations
%             I_want = ["Q90-a-50t"];
%             if ismember(desc, I_want)==0
%                 continue
%             end
            
            disp(desc)
        
            % runs and saves results to csv
            equilibrium
            results_table
            
            % figures
            if desc=="Q90-a-50t"
                suffix = strcat(filename_results_table, '_', desc);
                filename_figures_paper_suffix = suffix;
                filename_figures_report_suffix = suffix;
                filename_figures_smallplot_suffix = suffix;
                figures_paper
                figures_report
                figures_smallplot
            end
            
            descs = [descs, desc];
        end
    end
end
xi_p(:) = 0; % no more testing
%-------------------------------------------------------------------------%
%                             Shelter-at-home                             %
%-------------------------------------------------------------------------%
disp('===============')
disp('Shelter-at-home')
disp('===============')
if alpha1==0
    new_lambdas = [1.00667809548158, 2.93688139007933, 8.13814414450348, 19.8648352905615]; % 25, 50, 75%, 90% increase in d (calibration without teleworking)
else
    new_lambdas = [0.383675452740587, 1.06851029113554, 2.45096990297207, 32.3598330005704]; % 25%, 50%, 75%, 90% increase in d+v (calibration with teleworking)
end
durations = [4, 8, 12, 26, 35, 78]; % weeks
descs_new_lambda = ["25", "50", "75", "90"];
for i_new_lambda = 1:length(new_lambdas)
    new_lambda = new_lambdas(i_new_lambda);
    desc_new_lambda = descs_new_lambda(i_new_lambda);
    for i_duration = 1:length(durations)
        duration = durations(i_duration);
        
        for i_who = 1:3 % all, y, o
            lambda_p_i(:,:) = 0;
            lambda_p_h(:,:) = 0;
            lambda_p_r(:,:) = 0;
            lambda_p_f(:,:) = 0;
            
            if i_who==1
                lambda_p_i(:, 1:duration) = max(new_lambda - lambda_i, 0);
                lambda_p_h(:, 1:duration) = new_lambda;
                lambda_p_r(:, 1:duration) = new_lambda;
                lambda_p_f(:, 1:duration) = new_lambda;
                desc_who = "a";
            elseif i_who==2
                lambda_p_i(i_young, 1:duration) = max(new_lambda - lambda_i, 0);
                lambda_p_h(i_young, 1:duration) = new_lambda;
                lambda_p_r(i_young, 1:duration) = new_lambda;
                lambda_p_f(i_young, 1:duration) = new_lambda;
                desc_who = "y";
            elseif i_who==3
                lambda_p_i(i_old, 1:duration) = max(new_lambda - lambda_i, 0);
                lambda_p_h(i_old, 1:duration) = new_lambda;
                lambda_p_r(i_old, 1:duration) = new_lambda;
                lambda_p_f(i_old, 1:duration) = new_lambda;
                desc_who = "o";
            end
            
            desc = strcat("SH", desc_new_lambda, "-", desc_who, "-", num2str(floor(duration)));
            
%             % uncomment this if I want to skip some simulations
%             I_want = ["SH90-a-26"];
%             if ismember(desc, I_want)==0
%                 continue
%             end
            
            disp(desc)
            % runs and saves results
            equilibrium
            results_table
            
            % figures
            I_want = ["SH90-a-26", "SH90-o-26"];
            if ismember(desc, I_want)==1
                suffix = strcat(filename_results_table, '_', desc);
                filename_figures_paper_suffix = suffix;
                filename_figures_report_suffix = suffix;
                filename_figures_smallplot_suffix = suffix;
                figures_paper
                figures_report
                figures_smallplot
            end
            
            descs = [descs, desc];
        end
    end
end
% clears lambda_p variables
lambda_p_i(:, :) = 0;
lambda_p_h(:, :) = 0;
lambda_p_r(:, :) = 0;
lambda_p_f(:, :) = 0;
% exports descriptions to a csv
filename = strcat('tables/', filename_results_table, '_descs.csv');
descs = cellstr(descs);
writetable(cell2table(descs), filename, 'writevariablenames', 0)
%}
