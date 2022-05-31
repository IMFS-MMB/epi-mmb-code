%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get IRFS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dset,history,history1] = makeirfsecondorder(M_,oo_,nperiods,shock,order,startingpoint1,startingpoint2)

lgy_ = M_.endo_names;
dr_ = oo_.dr;
ys_ = oo_.dr.ys;


    nfwrd = M_.nfwrd;
    nstatic =M_.nstatic; 

% for older version of Dynare:    
%    nfwrd = dr_.nfwrd;
%    nstatic = dr_.nstatic;

if ~exist('startingpoint1','var')
[history, history1, h0] = mymkirf2(dr_,nstatic,nfwrd,ys_,nperiods,shock,0,order,0);  % returns matrix containing IRFs for all variables
else
    if ~exist('startingpoint2','var')
        error('Also need the starting point for the second-order pruned chain')
    end
    [history, history1, h0] = mymkirf2(dr_,nstatic,nfwrd,ys_,nperiods,shock,0,order,0,...
                                startingpoint1,startingpoint2); 
   
end
reordered = lgy_(dr_.order_var,:);          % the vector reordered contains names of variables
                                            % corresponding to rows of
                                            % history

% this loop assigns each row of history to its variable
% it also assigns steady state values
% IRFs are stored under varname_irf  where varname is the variable name
% SS values are stored under varname
% NB: SS values are scalar, IRFs are column vectors.

% history = history(dr_.order_var,:);
% history = history';

for indxi = 1:M_.endo_nbr
    eval(['dset.',deblank(reordered(indxi,:)),'_irf=transpose(history(indxi,:));']);
    eval(['dset.',deblank(lgy_(indxi,:)),'_ss= ys_(indxi);']);
    eval(['dset.',deblank(reordered(indxi,:)),'_h0=transpose(h0(indxi,:));']);
end

history = history(dr_.inv_order_var,:);
history1 = history1(dr_.inv_order_var,:);



