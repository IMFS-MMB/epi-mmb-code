function [ys,params,check] = two_sector_steadystate(ys,exo,M_,options_)


%global M_  



check = 0;

% real economy parameters
paramfile_one_sector_vcu


% export parameters
nparams = size(M_.param_names,1);
for icount = 1:nparams
    eval(['M_.params(icount) = ' M_.param_names{icount} ';'])
end
params = M_.params;
% transfer parameters to Dynare.
% send chip to the parameter list
%M_.params(strmatch('chip',M_.param_names,'exact')) = chip;
%M_.params(strmatch('pbss',M_.param_names,'exact')) = pbss;
%M_.params(strmatch('puss',M_.param_names,'exact')) = puss;
%M_.params(strmatch('pcss',M_.param_names,'exact')) = pcss;




nvars = M_.endo_nbr;

ys = zeros(nvars,1);



 for i_indx = 1:nvars
     eval(['ys(i_indx)=',M_.endo_names{i_indx,:},';']) 
 end


