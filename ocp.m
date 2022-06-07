clear all; clc; close all;
dy_path='c:/dynare/5.1/matlab';
dy_root="c:/dynare/";
%dy_path='c:/dynare/4.5.7/matlab';
addpath(dy_path); 
savepath=pwd;
addpath(pwd);
addpath(fullfile(pwd,'frequency'));

Time=100;

if(size(getenv('MODEL_NAME'), 1))
    list_models = { getenv('MODEL_NAME') };
else
    cd models
    modelspath=pwd;
    all = dir;
    all = all(3:end);

    list_models = {};
    for i =1:size(all,1)
        if all(i).isdir
            list_models{end+1}=all(i).name;
        end
    end
    cd ..
end

%list_models=list_models(6);
% main_epimmb(  [modellist (string)], ...
%               [macrovariablelist (string)],...
%               horizon (number), re_simulate (0=no, 1=yes), ...
%               [indicator_shock_size(0=model specific, 1=common) common_shock_size])
%main_epimmb(["HHK_20","AJRT_20","ACS_21","ERT_21_NK","KUX_20"],...
shocks=zeros(4,2); % number shocks x shock size
shocks(1,:)=[0 0]; %model specific
shocks(2,:)=[1 0.0005]; %low
shocks(3,:)=[1 0.001 ]; %medium
shocks(4,:)=[1 0.0025]; %large
shocklist={'Model-specific Initial Infections','Low Initial Infections','Medium Initial Infections','High Initial Infections'};
for index_m=1:size(list_models,2)
    for index_s=1:size(shocks,1)
        modelname=list_models{index_m};
        shockname=shocklist{index_s};
        str_modelname=string(modelname);
        str_shockname=string(shockname);
        json_filename=str_modelname+"-"+str_shockname+".output.json";        
        %cd(list_models{index_m})
        results_mat=ocp_epimmb(list_models(index_m),["Consumption","Labour","Output","Susceptibles","Infected","Recovered","Deaths"],[shocks(index_s,1) shocks(index_s,2)],dy_root);

        if(~exist('results', 'file'))
          mkdir results;
        end

        cd results
        jso.model=modelname;
        jso.shock=shockname;
        jso.data.Consumption=results_mat(1,1:Time+1); 
        jso.data.Labour=results_mat(2,1:Time+1); 
        jso.data.Output=results_mat(3,1:Time+1); 
        jso.data.Susceptibles=results_mat(4,1:Time+1); 
        jso.data.Infected=results_mat(5,1:Time+1); 
        jso.data.Recovered=results_mat(6,1:Time+1); 
        jso.data.Deaths=results_mat(7,1:Time+1); 
        json_print=jsonencode(jso);
        new_string = strrep(json_print, ',', ',\n');
        % add a return character after curly brackets:
        new_string = strrep(new_string, '{', '{\n');
        fid=fopen(json_filename,'w'); 
        fprintf(fid, new_string); 
        fclose('all'); 
        cd ..
        
    end
end



    
