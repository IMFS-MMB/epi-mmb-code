function unionvariablelist = unionvariable(model_list,rootdir)

 for iii =1:length(model_list)    
    modelName = model_list{iii};
    cd models
    eval(strcat("cd ",modelName,";" ))
    fname = strcat(modelName,".json"); 
    fid = fopen(fname); 
    raw = fread(fid,inf); 
    str = char(raw'); 
    fclose(fid); 
    jcode = jsondecode(str);
    allvarnames_tmp = convertCharsToStrings(jcode.variables);

% %% To check if there are some typos in the json file in case spoted
%     check1=ismember("Susceptible",allvarnames_tmp);
% 
%     if check1 ==1
%         modelName
%         display("Susceptible is in the list")
%     end      
        
    if iii == 1
        allvarnames = allvarnames_tmp;
    end
    allvarnames = union(allvarnames,allvarnames_tmp);

    cd(rootdir)
 end
 unionvariablelist = allvarnames;