function result_mat=ocp_epimmb(modelname, macrovariablelist, inf_ini,dy_root)
% function starts
%====================================================================================================
% navigate to the root folder (if exists)
try
    cd(ppp.path.root)
end
ppp.path.root = convertCharsToStrings(pwd);
ppp.path.working = fullfile(ppp.path.root, "working");
ppp.path.share.dyModules =  fullfile(ppp.path.root, "shared", "dyModules") ;
%num_vars=10;
num_vars = length(macrovariablelist);
try
    rmdir(ppp.path.working,'s')
end
mkdir working
maxhorizon = 700;

num_macrovar = length(macrovariablelist);
result_mat = NaN(num_macrovar,maxhorizon);

disp(modelname);
ppp.path.models = fullfile(ppp.path.root, "models");
ttt.path.model = fullfile(ppp.path.models, modelname);
ttt.path.modeldymodules = fullfile(ttt.path.model, "dyModules") ;
cd(ttt.path.model);


% copy everything to working folder
copyfile( ttt.path.model ,  ppp.path.working)
try
    copyfile( ttt.path.modeldymodules ,  ppp.path.working)
end
    %copyfile( ppp.path.share.dyModules, ppp.path.working)

cd(ppp.path.working)

%load json file to see how to make adjustments
fname = strcat(modelname,".json"); 
fid = fopen(fname); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 
jcode = jsondecode(str);

%Remove because every model must be run with Dy 5.1
%{
if (jcode.Dynare_version ~= "0") && (jcode.Dynare_version ~= "??")
    dy_path=dy_root+jcode.Dynare_version+"/matlab";
    addpath(dy_path);
elseif jcode.Dynare_version == "??"
    disp('Dynare path ??  Please check which is the correct version');
end
%}

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
    ress_load_ss=load('simulated_results_ss.mat');
    ress=nan(num_vars,size(ress_load.Susceptibles,1));
    ress_ss=nan(num_vars,1);
    ressnames=fieldnames(ress_load); %sorted list including SS
    ress_ssName = fieldnames(ress_load_ss);
    
      
    % retriving results from .mat
    for iii = 1:length(fieldnames(ress_load))
        ressName_tmp = char(ressnames(iii));
        eval(strcat("ress(iii,:) = ress_load." , ressName_tmp, ";"))
        
        ress_ssName_tmp = char(ress_ssName(iii));
        eval(strcat("ress_ss(iii) = ress_load_ss." , ress_ssName_tmp, ";"))
    end
%     ress(1,:)=ress_load.Consumption;
%     ress(2,:)=ress_load.Deaths;
%     ress(3,:)=ress_load.Infected;           
%     ress(4,:)=ress_load.Labour;
%     ress(5,:)=ress_load.Output;
%     ress(6,:)=ress_load.Recovered;
%     ress(7,:)=ress_load.Susceptibles;
%     ress(8,:)=ress_load.Interest;
%     ress(9,:)=ress_load.Inflation;
%     ress(10,:)=ress_load.Investment;
% 
%     ress_ss(1)=ress_load_ss.Consumption_ss;
%     ress_ss(2)=ress_load_ss.Deaths_ss;
%     ress_ss(3)=ress_load_ss.Infected_ss;
%     ress_ss(4)=ress_load_ss.Labour_ss;
%     ress_ss(5)=ress_load_ss.Output_ss;
%     ress_ss(6)=ress_load_ss.Recovered_ss;
%     ress_ss(7)=ress_load_ss.Susceptibles_ss;
%     ress_ss(8)=ress_load_ss.Interest_ss;
%     ress_ss(9)=ress_load_ss.Inflation_ss;
%     ress_ss(10)=ress_load_ss.Investment_ss;


else
    eval(strcat("dynare ", modelname))
    load simulated_results
end

disp("results")
if jcode.Code_type=="Dynare"
    for ind_macrovar = 1: num_macrovar
        macrovar = macrovariablelist(ind_macrovar);
        series_pos = find(strcmp(results.M_.endo_names,macrovar));
        if length(series_pos) > 0
            series_level = results.oo_.endo_simul(series_pos,:);
            series_ss = results.oo_.steady_state(series_pos,:);

            if modelname=="LFA_22"
                  if series_ss == 0                  
                        series = 100*series_level;
                   
                    if series==zeros(1,length(series))
                        series=nan(1,length(series));
                    end
                    if length(series) < maxhorizon
                        series(end+1:maxhorizon)=nan;
                    end
                    result_mat(ind_macrovar,1:(maxhorizon-1)) = series(2:maxhorizon);

                  else
                    if macrovar=="Susceptibles"| macrovar=="Infected"| macrovar=="Recovered"| macrovar=="Deaths";
                        series = 100*series_level;
                    else
                        if  macrovar=="Interest"| macrovar=="Inflation";
                        series = 100*(series_level- series_ss);
                        else
                        series = 100*(series_level - series_ss)/series_ss;
                        end
                    end
                    if series==zeros(1,length(series))
                        series=nan(1,length(series));
                    end
                    if length(series) < maxhorizon
                        series(end+1:maxhorizon)=nan;
                    end                    
                    result_mat(ind_macrovar,1:(maxhorizon-1)) = series(2:maxhorizon);
                end
            else
                if series_ss == 0 
                        series = 100*series_level;
                    
                    if series==zeros(1,length(series))
                        series=nan(1,length(series));
                    end
                    if length(series) < maxhorizon
                        series(end+1:maxhorizon)=nan;
                    end
                    result_mat(ind_macrovar,:) = series(1:maxhorizon);

                else
                    if macrovar=="Susceptibles"| macrovar=="Infected"| macrovar=="Recovered"| macrovar=="Deaths";
                        series = 100*series_level;
                    else
                        if  macrovar=="Interest"| macrovar=="Inflation";
                        series = 100*(series_level- series_ss);
                        else
                        series = 100*(series_level - series_ss)/series_ss;
                        end
                    end
                    if series==zeros(1,length(series))
                        series=nan(1,length(series));
                    end
                    if length(series) < maxhorizon
                        series(end+1:maxhorizon)=nan;
                    end                    
                    result_mat(ind_macrovar,:) = series(1:maxhorizon);
                end
                end
        else
             result_mat(ind_macrovar,:) = nan(1,maxhorizon);
        end

    end
elseif jcode.Code_type=="Matlab"
    for ind_macrovar = 1: num_macrovar
        macrovar = macrovariablelist(ind_macrovar);
        series_pos = find(strcmp(ressnames,macrovar));
        ss_pos = find(strcmp(ress_ssName,strcat(macrovar,"_ss")));
        result_mat(ind_macrovar,:) = nan(1,maxhorizon);

        if length(series_pos) > 0
            %S.(F{idx})
            series_level = ress(series_pos,:);
            series_ss = ress_ss(ss_pos);
            
            if exist('series_ss','var') == 1
                if series_ss == 0;
                   
                       series = 100* series_level;
                   
                    if series==zeros(1,length(series))
                        series=nan(1,length(series));
                    end
                    if length(series) < maxhorizon
                        series(end+1:maxhorizon)=nan;
                    end
                    result_mat(ind_macrovar,:) = series(1:maxhorizon);

                else
                    if  macrovar=="Interest"| macrovar=="Inflation";
                        series = 100*(series_level- series_ss);
                    else
                        series = 100*(series_level - series_ss)/series_ss;
                    end   
                    if series==zeros(1,length(series))
                        series=nan(1,length(series));
                    end
                    if length(series) < maxhorizon
                        series(end+1:maxhorizon)=nan;
                    end                    
                    result_mat(ind_macrovar,:) = series(1:maxhorizon);
                    end

                end
           
        else
            result_mat(ind_macrovar,:) = nan(1,maxhorizon);
        end

    end    
else
    display("Invalid Code Type! Choose Matlab or Dynare or code it yourself!")
    pause(4)
end

try
    rmpath(dy_path);
end
if convertCharsToStrings(pwd) == ppp.path.working
    cd(ppp.path.root)
    try
        rmdir(ppp.path.working,'s')
        mkdir working
    end
end

end

