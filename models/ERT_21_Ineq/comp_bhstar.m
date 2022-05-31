%this script computes the path of assets such that the assets FOC holds
disp(' ');
disp('Computing bhstar. This may take a while...');
disp(' ');

U_idx=strmatch('U',char(M_.endo_names),'exact');
C_idx=[]; 
lambh_idx=strmatch('lambh',char(M_.endo_names),'exact');
residEuler_idx=strmatch('residEuler',char(M_.endo_names),'exact');

oo_.exo_simul(:,1)=0; %set diffbhstar=0 for all t

options_.verbosity=0;
[oo_.endo_simul, oo_.deterministic_simulation] = sim1(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, M_, options_);

opts_fmincon=optimoptions('fmincon','Display','iter','ConstraintTolerance',1e-7,'StepTolerance',1e-7,'TolFun',1e-7,...
    'MaxFunctionEvaluations',2000,'UseParallel',true,'FiniteDifferenceStepSize',1e-3);

if load_ini_guess_diffbhstar==1     
        load ini_guess_diffbhstar;
        guess=ini_guess_diffbhstar;    
else
    guess=zeros(horzz,1);
end


LB=[];UB=[];%lower and upper bounds

sol = fmincon(@get_welf,guess,[],[],[],[],LB,UB,@nonlinconstraint,opts_fmincon,oo_.endo_simul, oo_.exo_simul, oo_.steady_state, M_, options_,U_idx,lambh_idx,rstar,betta,horzz,C_idx,residEuler_idx);
 
ini_guess_diffbhstar=sol;
save ini_guess_diffbhstar ini_guess_diffbhstar;

options_.verbosity=1;
[welf,cc,ceq,oo_.endo_simul, oo_.deterministic_simulation,oo_.exo_simul]=get_welf(sol,oo_.endo_simul, oo_.exo_simul, oo_.steady_state, M_, options_,U_idx,lambh_idx,rstar,betta,horzz,C_idx,residEuler_idx);
dyn2vec(M_, oo_, options_);


function [welf,cc,ceq,oo_endo_simul, oo_deterministic_simulation,oo_exo_simul]=get_welf(guess,oo_endo_simul, oo_exo_simul, oo_steady_state, M_, options_,U_idx,lambh_idx,rstar,betta,horzz,C_idx,residEuler_idx)

oo_exo_simul(2:2+horzz-1,1)=guess(1:horzz); %set bhstareps

[oo_endo_simul, oo_deterministic_simulation] = sim1(oo_endo_simul, oo_exo_simul, oo_steady_state, M_, options_); %simulate model

welf=-oo_endo_simul(U_idx,2)'; %get welfare (flip sign since we are using fmincon). 

%ceq=0, equality constraint
ceq=oo_endo_simul(residEuler_idx,2:horzz+1); %impose Euler equation for assets, i.e. drive Euler error to zero

%cc<0, inequality constraint.
cc=[];

end


function [cc,ceq] = nonlinconstraint(guess,oo_endo_simul, oo_exo_simul, oo_steady_state, M_, options_,U_idx,lambh_idx,rstar,betta,horzz,C_idx,residEuler_idx)
%get value of nonlinear constraint (i.e. Euler equation) must hold
%with equality (ceq)
[~,cc,ceq]=get_welf(guess,oo_endo_simul, oo_exo_simul, oo_steady_state, M_, options_,U_idx,lambh_idx,rstar,betta,horzz,C_idx,residEuler_idx);
end

