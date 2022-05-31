% Simple Planning Problem - finite differences
% Fernando Alvarez, David Argente, Francesco Lippi --  August 2020


%close all
%tic
lockdown_case=0;
% setting for saving data and figures
savedata=0; 

% time discount parameters
r=  0.05; % annual interest rate
dur_virus   = 1; % number of days until cured (cured = vaccine + treatment)
nu = 1/dur_virus ; % daily discount with vaccine + treatment

% propagation parameters
beta    = parameters(R,1) ;  % propagation rate Sdot = beta S * X / N
ggamma  = parameters(R,2) ;  % avg 20 days to get out of infection state (either recovered or dead)
theta   = parameters(R,3) ;  % effectiveness of transmission    
alpha   = parameters(R,4) ;  % parameter of the Testing-Tracing cost c(T ,S,X) = alpha/2 (T (S+X)/X)^2;
zeta    = parameters(R,5) ;  % effectiveness of tracing, zeta=0, random sampling, zeta=1 perfect tracing

% parameter for death rate function  Daily Death (t) =  (ggammad +  kappa*I(i)^alphad) * I(t) ;  
varphi  = parameters(R,6) ;  %   ; % baseline case mortality is 0.68%  adjusted by duration of desease
kappa   = parameters(R,7) ; %   ; % benchmark 0.034      % non-linear death rate (units adjusted by duration of desease), kappa=0 to allow for TTQ protocol increases by 1pp for every 10pp increase in Infected

% Lockdown, testing, tracing, quarantine
tau = parameters(R,8) ; % tau =1 (test) tau=0 (no test) 
l_max = parameters(R,9) ;  % bench = 0.75;
T_max = parameters(R,10) ;   % maximun testing-tracing rate, T_max = 1 to allow for TTQ protocol

% Welfare cost parameters
w   =  parameters(R,11) ; % normalized wage
% value statistical life
vsl = parameters(R,12) ;

%Start 300 grid points
N_S = parameters(R,13) ;
k = parameters(R,14) ; % ratio of Delta_S = k Delta_X
N_X = (N_S-1)*k + 1 ;

Delta_S = 1/(N_S-1);
Delta_X = 1/(N_X-1);

S        = linspace(0,1,N_S); 
X        = linspace(0,1,N_X); 

max_dt_mat = ones(N_S,N_X)/(r+nu); 
for i=2:N_S-1
    for j=2:N_X - (i-1)*k
        max_dt_mat(i,j) = 1/ ( r+nu+ beta* S(i)*X(j)*(1/Delta_S+1/Delta_X) + (ggamma+T_max)* X(j) / Delta_X );
        if S(i)+X(j)>1
            disp('shit')
        end
    end
end

max_dt = min(min(max_dt_mat));

dt = 0.95 * max_dt;

if dt > max_dt
    disp(['dt =  ',num2str(dt)]) ;
    disp('dt is too large');
    disp(['dt must be smaller than ',num2str(max_dt)]);   
end

V = zeros(N_S,N_X);
V_initial = V;
V_next = V ;
L = zeros(N_S,N_X);



for j=1:N_X
     V_initial(1,j) = vsl*X(j) * ( varphi*ggamma/(r+nu+ggamma) + kappa*ggamma*X(j)/(r+nu+2*ggamma) ) ;
end
for i=2:N_S-1
    for j=2:N_X+1-i
        V_initial(i,j) = w*l_max*( tau*(S(i)+X(j)) +1-tau )/(r+nu) +  vsl*X(j) * ( varphi*ggamma/(r+nu+ggamma) + kappa*ggamma*X(j)/(r+nu+2*ggamma) ) ;    
    end
end

V = V_initial; 
V_next = V_initial; 


diffmax=10; 
diffmax_policy= diffmax;   
iteration=0;
L_optimal = zeros(N_S,N_X);
L_old     = L_optimal ;
T_optimal = zeros(N_S,N_X);
T_old     = L_optimal ;
check_foc = zeros(N_S,N_X);
Diff_Der_I_m_S = zeros(N_S,N_X);
Diff_Der_X     = zeros(N_S,N_X);
Convex_HBJ   = zeros(N_S,N_X); 


%%% Set to NaN to values outside feasible set. 
for i = 1: N_S 
    for j = 1:N_X
        if S(i)+X(j) > 1 %(i+j-2)*Delta>1
            V(i,j)         = NaN;
            V_initial(i,j) = NaN;
            V_next(i,j)   =  NaN;
            Convex_HBJ(i,j) = NaN;
            Diff_Der_I_m_S(i,j) = NaN;
        end
    end
end

N_non_boundary = N_S*N_X- sum(sum(isnan(V_initial))) - N_X - N_S + 1 ;

diff_v = 0; 
s=0;

Ind_Cost_TTQ = (w+vsl*varphi*ggamma)/(r+nu+ggamma);


while  diffmax>0.000001*r || (diffmax_policy > 0.000001) || iteration < 2/dt 
    iteration=iteration+1   ;
    for i= 2: N_S-1
        for j=2 : N_X-1-(i-1)*k
            Diff_Der_I_m_S(i,j)  = 1/Delta_X * (V(i,j+1) - V(i,j)) - 1/Delta_S * (V(i,j) - V(i-1,j)) ;
            
            Diff_Der_X(i,j) = 1/Delta_X * (V(i,j) - V(i,j-1)) ;
            T_foc = ( Diff_Der_X(i,j) - Ind_Cost_TTQ) * ( X(j)/ (S(i)+X(j)) )^(2*(1-zeta)) * 1/alpha;
            T_foc = min( T_max, max( T_foc , 0 ) );  
            T_optimal(i,j) = T_foc;
            T = T_optimal(i,j) ;
            
            if  Diff_Der_I_m_S(i,j) > 0  %%% right hand side of HJB is convex
                Convex_HBJ(i,j)= 1 ;
                Lfoc = - w*( tau*(S(i)+X(j)) + 1-tau ) / ( 2*theta^2 * beta * S(i)*X(j)* Diff_Der_I_m_S(i,j) ) + 1/theta ;  %%% use foc
                check_foc(i,j) =  w *( tau*(S(i)+X(j)) + 1-tau ) - 2*theta*beta*S(i)*X(j)*(1-theta*Lfoc) * Diff_Der_I_m_S(i,j) ;
                L_optimal(i,j) = min(l_max, max(Lfoc, 0));
                L = L_optimal(i,j) ;
                V_next(i,j)  =  L *  w * (tau *( S(i)+X(j)) + 1-tau ) * dt  +  X(j)*(varphi+kappa*X(j))*ggamma * vsl * dt ...
                    +  Ind_Cost_TTQ * T * dt + alpha/2* ( T* ( (S(i)+X(j))/X(j) )^(1-zeta) )^2 * dt  ...
                    + (1-dt*(r+nu)) * ( 1- ( beta*S(i)*X(j)*(1-theta*L)^2*(dt/Delta_S+dt/Delta_X) + (ggamma*X(j)+T)*(dt/Delta_X) ) /(1-dt*(r+nu)) ) * V(i,j) ...
                    + (1-dt*(r+nu)) * ( beta*S(i)*X(j)*(1-theta*L)^2 * (dt/Delta_X)/(1-dt*(r+nu)) ) * V(i,j+1)   ...
                    + (1-dt*(r+nu)) * ( beta*S(i)*X(j)*(1-theta*L)^2 * (dt/Delta_S)/(1-dt*(r+nu)) ) * V(i-1,j)  ...
                    + (1-dt*(r+nu)) * ( (ggamma*X(j)+T)* (dt/Delta_X)/(1-dt*(r+nu)) ) * V(i,j-1) ;
                diff_v = diff_v + abs(V_next(i,j)-V(i,j)) ;
                s=s+1;
            else %%% RHS of HJB is concave
                Convex_HBJ(i,j)= 0 ;
                L = 0;
                V_0  =  L *  w * (tau *( S(i)+X(j)) + 1-tau ) * dt  +  X(j)*(varphi+kappa*X(j))*ggamma * vsl * dt ...
                    +  Ind_Cost_TTQ * T * dt + alpha/2* ( T* ( (S(i)+X(j))/X(j) )^(1-zeta) )^2 * dt...        
                    + (1-dt*(r+nu)) * ( 1- ( beta*S(i)*X(j)*(1-theta*L)^2*(dt/Delta_S+dt/Delta_X) + (ggamma*X(j)+T)*(dt/Delta_X) ) /(1-dt*(r+nu)) ) * V(i,j) ...
                    + (1-dt*(r+nu)) * ( beta*S(i)*X(j)*(1-theta*L)^2 * (dt/Delta_X)/(1-dt*(r+nu)) ) * V(i,j+1)   ...
                    + (1-dt*(r+nu)) * ( beta*S(i)*X(j)*(1-theta*L)^2 * (dt/Delta_S)/(1-dt*(r+nu)) ) * V(i-1,j)  ...
                    + (1-dt*(r+nu)) * ( (ggamma*X(j)+T)* (dt/Delta_X)/(1-dt*(r+nu)) ) * V(i,j-1) ;
                L = l_max;
                V_max  =  L *  w * (tau *( S(i)+X(j)) + 1-tau ) * dt  +  X(j)*(varphi+kappa*X(j))*ggamma * vsl * dt ...
                    +  Ind_Cost_TTQ * T * dt + alpha/2* ( T* ( (S(i)+X(j))/X(j) )^(1-zeta) )^2 * dt  ...        
                    + (1-dt*(r+nu)) * ( 1- ( beta*S(i)*X(j)*(1-theta*L)^2*(dt/Delta_S+dt/Delta_X) + (ggamma*X(j)+T)*(dt/Delta_X) ) /(1-dt*(r+nu)) ) * V(i,j) ...
                    + (1-dt*(r+nu)) * ( beta*S(i)*X(j)*(1-theta*L)^2 * (dt/Delta_X)/(1-dt*(r+nu)) ) * V(i,j+1)   ...
                    + (1-dt*(r+nu)) * ( beta*S(i)*X(j)*(1-theta*L)^2 * (dt/Delta_S)/(1-dt*(r+nu)) ) * V(i-1,j)  ...
                    + (1-dt*(r+nu)) * ( (ggamma*X(j)+T)* (dt/Delta_X)/(1-dt*(r+nu)) ) * V(i,j-1) ;
                if V_max < V_0
                    L_optimal(i,j) = l_max;
                    V_next(i,j) = V_max ;
                    diff_v = diff_v + abs(V_next(i,j)-V(i,j)) ;
                    s=s+1;
                else
                    L_optimal(i,j) = 0 ;
                    V_next(i,j) = V_0;
                    diff_v = diff_v + abs(V_next(i,j)-V(i,j)) ;
                    s=s+1;
                end
            end
        end
        
        j = N_X - (i-1)*k ;
        Diff_Der_I_m_S(i,j)  = 1/Delta_S * ( V(i-1,j+k) -V(i,j) ) ;
        
         Diff_Der_X(i,j) = 1/Delta_X * (V(i,j) - V(i,j-1)) ;
         T_foc = ( Diff_Der_X(i,j) - Ind_Cost_TTQ) * ( X(j)/ (S(i)+X(j)) )^(2*(1-zeta)) * 1/alpha;
         T_foc = min( T_max, max( T_foc , 0 ) );  
         T_optimal(i,j) = T_foc;
         T = T_optimal(i,j) ;
        
        %%% optimal L
        if  Diff_Der_I_m_S(i,j) > 0  %%% right hand side of HJB is convex
            Convex_HBJ(i,j) =1;
            Lfoc = - w*( tau*(S(i)+X(j)) + 1-tau ) / ( 2*theta^2 * beta * S(i)*X(j)* Diff_Der_I_m_S(i,j) ) + 1/theta ;  %%% use foc
            check_foc(i,j) =  w *( tau*(S(i)+X(j)) + 1-tau ) - 2*theta*beta*S(i)*X(j)*(1-theta*Lfoc) * Diff_Der_I_m_S(i,j) ;
            L_optimal(i,j) = min(l_max, max(Lfoc, 0));
            L = L_optimal(i,j) ;
            V_next(i,j)  =  L *  w * (tau *( S(i)+X(j)) + 1-tau ) * dt  +  X(j)*(varphi+kappa*X(j))*ggamma * vsl * dt ...       
                +  Ind_Cost_TTQ * T * dt + alpha/2* ( T * ( (S(i)+X(j))/X(j) )^(1-zeta) )^2 * dt ...
                + (1-dt*(r+nu)) * ( 1- ( beta*S(i)*X(j)*(1-theta*L)^2*(dt/Delta_S) + (ggamma*X(j)+T)*dt/Delta_X) /(1-dt*(r+nu)) ) * V(i,j) ...
                + (1-dt*(r+nu)) * ( beta*S(i)*X(j)*(1-theta*L)^2 * (dt/Delta_S)/(1-dt*(r+nu)) ) * V(i-1,j+k)   ...
                + (1-dt*(r+nu)) * (  (ggamma*X(j)+T) * (dt/Delta_X)/(1-dt*(r+nu)) )  * V(i,j-1) ;
            diff_v = diff_v + abs(V_next(i,j)-V(i,j)) ;
            s=s+1;
        else %%% RHS of HJB is concave
            Convex_HBJ(i,j)= 0 ;
            L = 0;
            V_0  =  L *  w * (tau *( S(i)+X(j)) + 1-tau ) * dt  +  X(j)*(varphi+kappa*X(j))*ggamma * vsl * dt ...
                +  Ind_Cost_TTQ * T * dt + alpha/2* ( T* ( (S(i)+X(j))/X(j) )^(1-zeta) )^2 * dt ...        
                + (1-dt*(r+nu)) * ( 1- ( beta*S(i)*X(j)*(1-theta*L)^2*(dt/Delta_S) + (ggamma*X(j)+T)*dt/Delta_X) /(1-dt*(r+nu)) ) * V(i,j) ...
                + (1-dt*(r+nu)) * ( beta*S(i)*X(j)*(1-theta*L)^2 * (dt/Delta_S)/(1-dt*(r+nu)) ) * V(i-1,j+k)   ...
                + (1-dt*(r+nu)) * (  (ggamma*X(j)+T) * (dt/Delta_X)/(1-dt*(r+nu)) )  * V(i,j-1) ;
            L = l_max;
            V_max  =  L *  w * (tau *( S(i)+X(j)) + 1-tau ) * dt  +  X(j)*(varphi+kappa*X(j))*ggamma * vsl * dt ...
                +  Ind_Cost_TTQ * T * dt + alpha/2* ( T* ( (S(i)+X(j))/X(j) )^(1-zeta) )^2 * dt...        
                + (1-dt*(r+nu)) * ( 1- ( beta*S(i)*X(j)*(1-theta*L)^2*(dt/Delta_S) + (ggamma*X(j)+T)*dt/Delta_X) /(1-dt*(r+nu)) ) * V(i,j) ...
                + (1-dt*(r+nu)) * ( beta*S(i)*X(j)*(1-theta*L)^2 * (dt/Delta_S)/(1-dt*(r+nu)) ) * V(i-1,j+k)   ...
                + (1-dt*(r+nu)) * ( (ggamma*X(j)+T) * (dt/Delta_X)/(1-dt*(r+nu)) )  * V(i,j-1) ;
            if V_max < V_0
                L_optimal(i,j) = l_max;
                V_next(i,j) = V_max ;
                diff_v = diff_v + abs(V_next(i,j)-V(i,j)) ;
                s=s+1;
            else
                L_optimal(i,j) = 0 ;
                V_next(i,j) = V_0;
                diff_v = diff_v + abs(V_next(i,j)-V(i,j)) ;
                s=s+1;
            end
            %end
        end
    end
    
    diffmax_policy = 1/2* norm(L_optimal-L_old)/(1+norm(L_old)) + 1/2*norm(T_optimal-T_old)/(1+norm(T_old));
    diffmax = r*diff_v / N_non_boundary ;
    diff_v = 0;
    s=0;
    V = V_next ;
    L_old = L_optimal ;
    T_old = T_optimal ;  
    
end

%toc;

V_adj = V ; 
V_ini_adj = V_initial;

policy_L = L_optimal;
policy_T = T_optimal;

Diff_Der_X_adj = Diff_Der_X;


for i = 1: N_S 
    for j = 1:N_X
        if S(i)+X(j) > 1  
            policy_L(i,j)=NaN;
            policy_T(i,j)=NaN;
            V_adj(i,j)=NaN;
            V_ini_adj(i,j)=NaN;
            Diff_Der_X_adj(i,j) = NaN;
        end
    end
end


Dt = 1/(10*365); 
 
%Simulation with policy function
Nyears = 2 ; % for simulations and figures%
time = round(iteration*dt/Dt) ; 
period=[1:1:time];

% setup vars
Stime=[]; Xtime=[]; Rtime=[]; Ntime=[]; Dtime=[];  
Stime_Lck=[]; Xtime_Lck=[]; Rtime_Lck=[]; Ntime_Lck=[];  Dtime_Lck=[]; 
Qtime=[];


% initial conditions
%Xtime(1)=0.01;
helper=load('inf_ini.mat');
Xtime(1)=helper.helper;

Stime(1)=1-Xtime(1)-0.02;  

%%% put it in a grid
[X0,X0_indx] = min(abs(X-Xtime(1)));
[S0,S0_indx] = min(abs(S-Stime(1)));

Xtime(1) = X(X0_indx);
Stime(1) = S(S0_indx);


Rtime(1)=1-Xtime(1)-Stime(1);
Dtime(1)=0;
Ntime(1)=Xtime(1)+Rtime(1)+Stime(1);

Xtime_Lck(1)=Xtime(1); 
Stime_Lck(1)=Stime(1); 
Rtime_Lck(1)=Rtime(1);
Dtime_Lck(1)=0;
Ntime_Lck(1) = Xtime_Lck(1)+ Rtime_Lck(1) + Stime_Lck(1);
Qtime(1) = 0; 

[X0,X0_indx] = min(abs(X-Xtime(1)));
[S0,S0_indx] = min(abs(S-Stime(1)));

Lck(1) = policy_L(S0_indx,X0_indx) ;  
TTQ(1) = policy_T(S0_indx,X0_indx) ;  


for t=2:time 
    % benchmark uncontrolled process
    Stime(t) =  Stime(t-1)*(1 - Dt*beta*Xtime(t-1)*(1-0) );
    Xtime(t) =  Xtime(t-1)*(1 - Dt*ggamma + Dt*beta*Stime(t-1)*(1-0))   ;

    fatality =  varphi*ggamma + kappa*ggamma*Xtime(t-1);
    Dtime(t) =  Dtime(t-1) + Xtime(t-1)*fatality*Dt ;
    Rtime(t) =  Rtime(t-1) + Xtime(t-1)*Dt*(ggamma-fatality) ;
    Ntime(t) =  Stime(t) + Xtime(t) + Rtime(t); 
    
    % Controlled process
    Stime_Lck(t) =  Stime_Lck(t-1)*(1 - Dt*beta*Xtime_Lck(t-1)*(1-theta*Lck(t-1))^2) ;
    Xtime_Lck(t) =  Xtime_Lck(t-1)*(1 - Dt*ggamma + Dt*beta*Stime_Lck(t-1)*(1-theta*Lck(t-1))^2) - Dt*TTQ(t-1) ;
    % Stock in quarantine 
    Qtime(t)  = Qtime(t-1) + TTQ(t-1)*Dt;
    
    fatality = varphi*ggamma + kappa*ggamma*Xtime_Lck(t-1);
    Dtime_Lck(t) =  Dtime_Lck(t-1) + Xtime_Lck(t-1)*fatality*Dt ;
    Rtime_Lck(t) =  Rtime_Lck(t-1) + Xtime_Lck(t-1)*Dt*(ggamma+TTQ(t-1)-fatality)  ;
    Ntime_Lck(t) =  Stime_Lck(t)   + Xtime_Lck(t) + Rtime_Lck(t);
    
    % Update for policy
    [Xt,Xt_indx] = min(abs(X-Xtime_Lck(t) ));
    [St,St_indx] = min(abs(S-Stime_Lck(t) ));
    
    % Update for no policyf
    [Xt_nopol,Xt_nopol_indx] = min(abs(X-Xtime(t) ));
    [St_nopol,St_nopol_indx] = min(abs(S-Stime(t) ));     
    
    %%% linearly interpolate L and T
    
                %%% compute the derivatives as finite differences, two
                %%% sided when possible
                if St_indx > 1 & St_indx < N_S
                    Der_L_S = (policy_L(St_indx+1,Xt_indx) - policy_L(St_indx-1,Xt_indx))/(2*Delta_S);
                    Der_T_S = (policy_T(St_indx+1,Xt_indx) - policy_T(St_indx-1,Xt_indx))/(2*Delta_S);
                elseif St_indx == 1
                    Der_L_S = (policy_L(St_indx+1,Xt_indx) - policy_L(St_indx,Xt_indx))/(Delta_S);
                    Der_T_S = (policy_T(St_indx+1,Xt_indx) - policy_T(St_indx,Xt_indx))/(Delta_S);
                else
                    Der_L_S = (policy_L(St_indx,Xt_indx) - policy_L(St_indx-1,Xt_indx))/(Delta_S);
                    Der_T_S = (policy_T(St_indx,Xt_indx) - policy_T(St_indx-1,Xt_indx))/(Delta_S);
                end
                if Xt_indx > 1 & Xt_indx < N_X
                    Der_L_X = (policy_L(St_indx,Xt_indx+1) - policy_L(St_indx,Xt_indx-1))/(2*Delta_X);
                    Der_T_X = (policy_T(St_indx,Xt_indx+1) - policy_T(St_indx,Xt_indx-1))/(2*Delta_X);
                elseif Xt_indx == 1
                    Der_L_X = (policy_L(St_indx,Xt_indx+1) - policy_L(St_indx,Xt_indx))/(Delta_X);
                    Der_T_X = (policy_T(St_indx,Xt_indx+1) - policy_T(St_indx,Xt_indx))/(Delta_X);
                else
                    Der_L_X = (policy_L(St_indx,Xt_indx) - policy_L(St_indx,Xt_indx-1))/(Delta_X);
                    Der_T_X = (policy_T(St_indx,Xt_indx) - policy_T(St_indx,Xt_indx-1))/(Delta_X);
                end
                                
    Lck(t) = policy_L(St_indx,Xt_indx) + Der_L_S*(Stime_Lck(t)-S(St_indx)) + Der_L_X*(Xtime_Lck(t)-X(Xt_indx)) ;
    TTQ(t) = policy_T(St_indx,Xt_indx) + Der_T_S*(Stime_Lck(t)-S(St_indx)) + Der_T_X*(Xtime_Lck(t)-X(Xt_indx)) ;
    
    Lck(t) = min( l_max, max( Lck(t), 0)) ;
    TTQ(t) = min( T_max, max( TTQ(t), 0)) ;

    
    clear Der_L_S Der_L_X Der_T_S Der_T_X 
    
end

Frac_Pop_Lck= ( tau*(Stime_Lck +  Xtime_Lck) + 1-tau ).*Lck;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmax = min( round(365/2) , time) ;
tmax = min(Nyears* round(365/2) , time) ;

rho = exp(-(r+nu)*Dt);

for t=1:time
    rho_vec(1,t)=rho^(t-1);
end
 
Dtime_flow(1)=Dtime_Lck(1,1); Dtime_flow2(1)=Dtime_Lck(1,1);
for t=2:time
    Dtime_flow(1,t)             = Dtime_Lck(t)- Dtime_Lck(t-1);
    Dtime_flow2(1,t)            = Dt*(varphi*ggamma + kappa*ggamma*Xtime_Lck(t-1))*Xtime_Lck(t-1);
    Dtime_flow_nopolicy(1,t)    = Dtime(t)-Dtime(t-1);
    Dtime_flow2_nopolicy(1,t)   = Dt*(varphi*ggamma + kappa*ggamma*Xtime(t-1))*Xtime(t-1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Xtime0 = Xtime(1);

V_ini =  V(S0_indx,X0_indx);
V_last     =  V_initial(St_indx,Xt_indx);


    if R<15
PV_Lck_loss   =        Dt * sum(( tau*(Stime_Lck+Xtime_Lck)+1-tau)*w .* Lck .* rho_vec)  ...
                     + Dt * sum( (w+vsl*varphi*ggamma)/(r+nu+ggamma) * TTQ .* rho_vec) + ....
                     + Dt * sum( alpha/2 * ( TTQ .* ((Xtime_Lck+Stime_Lck)./Xtime_Lck ).^(1-zeta) ).^2  .*  rho_vec );
PV_Death_loss =      Dt * sum( (varphi+ kappa*Xtime_Lck) * ggamma .* Xtime_Lck * vsl .*rho_vec);
PV_loss       =      PV_Lck_loss + PV_Death_loss + rho^time*V_last;
    else   
PV_Lck_loss   =        Dt * w * const * sum( ( ( tau*(Stime_Lck+Xtime_Lck)+1-tau) .* Lck ).^2 *1/2 .* rho_vec)  ...
                     + Dt * sum( (w+vsl*varphi*ggamma)/(r+nu+ggamma) * TTQ .* rho_vec) + ....
                     + Dt * sum( alpha/2 * ( TTQ .* ((Xtime_Lck+Stime_Lck)./Xtime_Lck ).^(1-zeta) ).^2  .*  rho_vec );

PV_Death_loss =      Dt * sum( (varphi+ kappa*Xtime_Lck) * ggamma .* Xtime_Lck * vsl .*rho_vec);
PV_loss       =      PV_Lck_loss + PV_Death_loss + rho^time*V_last;
    end
    

flow_value_fct_opt = r * V_ini ;
flow_pv_optimal    = (r * PV_loss)*100;
output_flowcost = (r * PV_Lck_loss )*100 ;
PV_nopolicy_loss  = 0 ;
PV_nopolicy_Death_loss= sum(Dtime_flow_nopolicy * vsl .*rho_vec);
PV_loss_nopolicy= PV_nopolicy_loss + PV_nopolicy_Death_loss ;

V_nopolicy_last      =  V_initial(St_nopol_indx,Xt_nopol_indx);
flow_pv_nopolicy     = (r * (PV_loss_nopolicy + rho^time*V_nopolicy_last))*100;


TABLE1(R,:) =[flow_pv_optimal,output_flowcost,flow_pv_nopolicy];


%save(['SPP_row' num2str(R) '.mat']);
%no control

%need to check this one!!
if lockdown_case==0
    Consumption = w.*(Stime+Xtime+Rtime);
    Labour = (Stime+Xtime+Rtime);
    Output = w.*(Stime+Xtime+Rtime);
    Susceptibles = Stime; 
    Infected = Xtime;
    Recovered = Rtime;
    Deaths = Dtime;
elseif lockdown_case==1
    Consumption = w.*((1-Lck).*(Stime_lck+Xtime_lck)+Rtime_lck);
    Labour = ((1-Lck).*(Stime_lck+Xtime_lck)+Rtime_lck);
    Output = w.*((1-Lck).*(Stime_lck+Xtime_lck)+Rtime_lck);
    Susceptibles = Stime_lck; 
    Infected = Xtime_lck;
    Recovered = Rtime_lck;
    Deaths = Dtime_lck;
end
Consumption_ss= w;
Labour_ss= 1;
Output_ss= w;
%save results!
save('simulated_results.mat','Consumption','Labour','Output','Susceptibles','Infected','Recovered','Deaths',...
    'Consumption_ss','Labour_ss','Output_ss');
