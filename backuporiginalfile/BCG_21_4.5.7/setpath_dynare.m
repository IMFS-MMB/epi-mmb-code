

restoredefaultpath

if ~ismac 
    dir1 = 'C:\dynare\4.5.7\matlab';
    dir2 = '..\toolkit_files';
    
else 
    dir1 = '/Applications/Dynare/4.5.7/MATLAB';
    dir2 = '../occbin_20140630/toolkit_files';
    
end

path(path,dir1);
path(path,dir2);

dynare_config