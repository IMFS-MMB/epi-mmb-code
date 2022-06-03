function main_epimmb(modellist, macrovariablelist, horizon, re_simulate, inf_ini)
% function starts
%====================================================================================================
% navigate to the root folder (if exists)
try
    cd(p.path.root)
end
p.path.root = convertCharsToStrings(pwd);
p.path.working = fullfile(p.path.root, 'working') ;
p.path.share.dyModules =  fullfile(p.path.root, "shared", "dyModules");
num_vars=7;
num_vars_with_ss=3;
try
    rmdir(p.path.working,'s')
end
mkdir working
maxhorizon = 700;

result_mat = [];
num_model = length(modellist);
num_macrovar = length(macrovariablelist);

for ind_model = 1: num_model
    modelname = modellist(ind_model);
    disp(modelname);
    p.path.models = fullfile(p.path.root, 'models');
    t.path.model = fullfile(p.path.models, modelname);
    t.path.modeldymodules = fullfile(t.path.model, "dyModules");
    cd(t.path.model);
        
    if re_simulate == 0
        load simulated_results_base
    else
        % copy everything to working folder
        copyfile( t.path.model ,  p.path.working)
        try
            copyfile( t.path.modeldymodules ,  p.path.working)
        end
            %copyfile( p.path.share.dyModules, p.path.working)

        cd(p.path.working)
        
        %load json file to see how to make adjustments
        fname = strcat(modelname,".json"); 
        fid = fopen(fname); 
        raw = fread(fid,inf); 
        str = char(raw'); 
        fclose(fid); 
        jcode = jsondecode(str);
        
        if inf_ini(1)==1
            helper=inf_ini(2);
            save inf_ini.mat helper
        else
            helper=jcode.initial_value;
            save inf_ini.mat helper
        end
        if isfile('modelmasterscript.m')
            modelmasterscript
            ress_load=load('simulated_results.mat');
            ress=nan(num_vars+num_vars_with_ss,size(ress_load.Consumption,1));
            ress(1,:)=ress_load.Consumption;
            ress(2,:)=ress_load.Consumption; %for SS name
            ress(3,:)=ress_load.Deaths;
            ress(4,:)=ress_load.Infected;           
            ress(5,:)=ress_load.Labour;
            ress(6,:)=ress_load.Labour; %for SS name
            ress(7,:)=ress_load.Output;
            ress(8,:)=ress_load.Output;  %for SS name
            ress(9,:)=ress_load.Recovered;
            ress(10,:)=ress_load.Susceptibles;
            ress_ss=nan(num_vars+num_vars_with_ss,1);
            ress_ss(1)=ress_load.Consumption_ss;
            ress_ss(5)=ress_load.Labour_ss;
            ress_ss(7)=ress_load.Output_ss;
            ressnames=fieldnames(ress_load); %sorted list including SS
        else
            eval(strcat("dynare ", modelname))
            load simulated_results
        end
        
    end
    disp("results")
    if jcode.Code_type=="Dynare"
        for ind_macrovar = 1: num_macrovar
            macrovar = macrovariablelist(ind_macrovar);
            series_pos = find(strcmp(results.M_.endo_names,macrovar));
            if length(series_pos) > 0
                series_level = results.oo_.endo_simul(series_pos,:);
                series_ss = results.oo_.steady_state(series_pos,:);
                
                if modelname=="LFA_21"
                  if series_ss == 0 | macrovar=="Susceptibles"| macrovar=="Infected"| macrovar=="Recovered"| macrovar=="Deaths";
                    series = series_level;
                    if length(series) < maxhorizon
                        series(end+1:maxhorizon)=nan;
                    end
                    result_mat(ind_model,ind_macrovar,:) = series(1:maxhorizon);

                else
                    series = (series_level - series_ss)/series_ss;
                    if length(series) < maxhorizon
                        series(end+1:maxhorizon)=nan;
                    end                    
                    result_mat(ind_model,ind_macrovar,:) = series(2:maxhorizon);
                end
            else
                if series_ss == 0 | macrovar=="Susceptibles"| macrovar=="Infected"| macrovar=="Recovered"| macrovar=="Deaths";
                    series = series_level;
                    if length(series) < maxhorizon
                        series(end+1:maxhorizon)=nan;
                    end
                    result_mat(ind_model,ind_macrovar,:) = series(1:maxhorizon);

                else
                    series = (series_level - series_ss)/series_ss;
                    if length(series) < maxhorizon
                        series(end+1:maxhorizon)=nan;
                    end                    
                    result_mat(ind_model,ind_macrovar,:) = series(1:maxhorizon);
                end
                end
           
            else
                 result_mat(ind_model,ind_macrovar,:) = nan(1,maxhorizon);
            end

        end
    elseif jcode.Code_type=="Matlab"
        for ind_macrovar = 1: num_macrovar
            macrovar = macrovariablelist(ind_macrovar);
            series_pos = find(strcmp(ressnames,macrovar));
            
            result_mat(ind_model,ind_macrovar,:) = nan(1,maxhorizon);
            
            if length(series_pos) > 0
                %S.(F{idx})
                series_level = ress(series_pos,:);
                try
                    series_ss = ress_ss(series_pos);
                end
                
                if series_ss == 0 | macrovar=="Susceptibles"| macrovar=="Infected"| macrovar=="Recovered"| macrovar=="Deaths";
                    series = series_level;
                    if length(series) < maxhorizon
                        series(end+1:maxhorizon)=nan;
                    end
                    result_mat(ind_model,ind_macrovar,:) = series(1:maxhorizon);

                else
                    series = (series_level - series_ss)/series_ss;
                    if length(series) < maxhorizon
                        series(end+1:maxhorizon)=nan;
                    end                    
                    result_mat(ind_model,ind_macrovar,:) = series(1:maxhorizon);
                end
            end

        end    
    else
        display("Invalid Code Type! Choose Matlab or Dynare or code it yourself!")
        pause(4)
    end
    
    if convertCharsToStrings(pwd) == p.path.working
        cd(p.path.root)
        try
            rmdir(p.path.working,'s')
            mkdir working
        end
    end
    
end

%% plot block
%========================================
figure('units','normalized','outerposition',[0 0 1 1])
if num_macrovar>3
   subplot_rownum =  ceil(num_macrovar/3);
   subplot_colnum = 3;
else
   subplot_rownum = 1;
   subplot_colnum = num_macrovar;
end

for ind_var_plot = 1:num_macrovar
    subplot(subplot_rownum,subplot_colnum,ind_var_plot)
    for ind_model_plot = 1:num_model
        tmp_legendname = strrep(modellist(ind_model_plot),'_','\_');
        tmp_series = squeeze(result_mat(ind_model_plot,ind_var_plot,:));
        if modellist(ind_model_plot) =="LFA_21" | modellist(ind_model_plot) =="JPV_21"   % model simulations start from the second period
            if macrovariablelist(ind_var_plot) == "Infected" % tailored for specific variables
                plot(1:horizon, tmp_series(1:(horizon)), 'Linewidth', 2 , 'DisplayName',tmp_legendname), hold on
            else
            plot(1:horizon, tmp_series(2:(horizon+1)), 'Linewidth', 2 , 'DisplayName',tmp_legendname), hold on
            legend
            end
        else
            if ~all(isnan(tmp_series))
                plot(1:horizon, tmp_series(1:horizon), 'Linewidth', 2 , 'DisplayName',tmp_legendname), hold on
                legend
            else
                plot(1:horizon, tmp_series(1:horizon), 'Linewidth', 2 , 'DisplayName','unavailable' ), hold on
            end
        end
    end
   
    hold off
    title(macrovariablelist(ind_var_plot))
end

%function ends
cd(p.path.root)


end