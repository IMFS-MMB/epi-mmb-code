% A Simple Planning Problem for COVID-19 Lockdown
% Fernando Alvarez, David Argente, Francesco Lippi - May 2020

%clear all
%close all

%tic
lockdown_case=0;
% Setting for saving data and figures
savedata   = 0 ; 
initialize = 1 ; % == 0 if already have VF

% RATES - ALLL MEASURED AS PER-YEAR RATES 
r  = 0.05;    % annual interest rate
nu = 1 / 1.5; % annual vaccine arrival rate  

% Propagation RATES (ANNUAL)
beta    = 0.20*365;    %  doubling every 0.7/(beta/365)  days  
ggamma  = 365/18  ;    %  avg 18 DAYS to get out of infection state (either recovered or dead)

% Parameter for death rate function  Daily Death (t) =  (varphi +  kappa*I(i)) * I(t) ;  
varphi  = 0.01 ; % baseline mortality is 1%  adjusted by duration of desease
kappa   = 0.05 ; % benchmark 0.05 ;  % non-linear death rate (units adjusted by duration of desease) increases by 1pp for every 10pp increase in Infected

% Welfare cost parameters
w     = 1   ;  % ANNUAL GDP per capita 
vsl   = 20  ;  % value of  stistical life as multiple of ANNUAL earnings (since w=1 ANNUAL OUTPUT)

% discrete time period 
dt    = 1/365; % discrete time period in daily model; used to map annual rates into daily rates

% policy vector
theta      =  0.50;      % effectiveness of lockdown    
L_max      =  0.70;  % bench = 0.70;
L_gridsize =  50; 
Lvec       =  linspace(0,L_max,L_gridsize);

% Initial conditions for graphs
%I0 = 0.01; % S(1) = 1- I(1)
helper=load('inf_ini.mat');
I0=helper.helper;

% S infinity if uncontrolled
s_inf=fzero(@(s) log(s)-beta/ggamma*(s-1), 0.3 );
%disp(['long run stationary level S is s_infinity =',num2str(s_inf,3)])

%disp(['Params :  Lmax=',num2str(L_max,3), ';  gamma   =',num2str(ggamma,3),';  beta  =',num2str(beta,3),';  theta=',num2str(theta,3)...
%    ,';  vsl =',num2str(vsl,3) ,';  varphi  =',num2str(varphi,3),';  kappa =',num2str(kappa,3)])

%  grid points
gridsize_S = 200;% bench 200
gridsize_I = gridsize_S ;% bench 200
S          = linspace(0,1,gridsize_S); 
I          = linspace(0,1,gridsize_I); 

% generate S, I domain 
[I_mat,S_mat] = ndgrid(I',S') ; % 

if initialize==1
    
    V =  zeros(gridsize_S,gridsize_I)  ;
    
    for s=1:gridsize_S
        
        I_feasible_indx = find(I<=1-S(s));  
        
        for i=1:max(I_feasible_indx)
            if i>1
                I_feasible_indx = find(I<=1-S(s));
                for i=1:max(I_feasible_indx)
                    V(s,i) =   vsl  *(I(i) + S(s))*   ggamma*( varphi / (r + nu + ggamma)  +  kappa*I(i) / (r + nu + 2*ggamma)  ) ; % PV in daily output units w=1
                end
            else % i ==1
                V(s,i) = 0;
            end
        end
    end
    
    V_guess = V; 
    
    Vnew = V_guess;
    %{ 
    %value function
    figure(30)
    h2=surf(S,I, V_guess'*r ); hold on;
    set(h2,'Linestyle','none')
    th=title('Value function  guess (annual flow units)');set(th,'fontsize',18)
    xlabel('Susceptible','FontSize',14);
    ylabel('Infected','FontSize',14);
    zlabel('value','FontSize',14);
    view(-55,25)
   %} 
    policy=zeros(gridsize_S,gridsize_I);  policy_indx=policy; policy_L=policy_indx; 
   iteration=0;
end

diffmax=10; diffmax_policy= diffmax; 
%Value function iteration

while diffmax>0.00001 || (diffmax_policy>0.00001) 
    iteration=iteration+1;
    for s=1:gridsize_S % S-sgrid
        
        I_feasible_indx = find(I<=1-S(s)); %find(I<=1-S(s));
    
        for i=1:max(I_feasible_indx)  % I-sgrid
            
            % mortality function
            fct_phi = ggamma *(varphi + kappa*I(i) ) ;  
 
            if i>1  && s>1 
            
            Sprime_vec =  S(s)*(1 - beta*I(i)*(1-theta*Lvec).^2*dt ) ;
            Iprime_vec =  I(i)*(1  + ( beta*S(s)*(1-theta*Lvec).^2 - ggamma)*dt ) ;

            v_candidate = (w*Lvec*(S(s)+I(i)) + fct_phi*I(i)*vsl )*dt +  exp(-(nu+r)*dt) * interp2(I_mat',S_mat',V,Iprime_vec,Sprime_vec) ;    % PV as multiple of daily values
           
            [Vnew(s,i) policy_indx(s,i)]=min(v_candidate); 
                        
            elseif i == 1 && s>1  
                
             Vnew(s,i) = 0; 
             policy_indx(s,i)=1 ;
             
            else   % S=0 ,  any I 
                
            Vnew(s,i) = vsl *  I(i) * ggamma* (   varphi / (r +nu +ggamma)   +  kappa * I(i) / (r + nu +2*ggamma)       ); 
            
            policy_indx(s,i)=1;       
            
            end
             
            policy_L(s,i) = Lvec(policy_indx(s,i)) ; 

            v_candidate=[];  
       
        end%i
        
        I_feasible_indx=[];
    end%s

    diffmax_policy= norm(policy_indx-policy)/norm(policy);
    policy=policy_indx; 
    diffmax = norm(( Vnew -V  ))/ norm (V)   ;
    V = Vnew;
    Vnew=[];
    
end      

%toc    
 
% Some figures
policy_L_adj=policy_L; 
V_adj = V ; 

for s = 1: gridsize_S
    I_unfeasible_indx = find(I>1-S(s));
    policy_L_adj(s,I_unfeasible_indx )=NaN;
    V_adj(s,I_unfeasible_indx )=NaN;
    I_unfeasible_indx =[];
end
%{
figure(1)
h2=surf(S,I,policy_L_adj');
set(h2,'Linestyle','none')
th=title('Policy function');set(th,'fontsize',18)
xlabel('Susceptible','FontSize',14); 
ylabel('Infected','FontSize',14); 
zlabel('Policy','FontSize',14); view(0,90)


% value function
figure(3)
h2=surf(S,I, V_adj'*r ); hold on;
set(h2,'Linestyle','none')
th=title('Value function (annual flow units)');set(th,'fontsize',18)
xlabel('Susceptible','FontSize',14); 
ylabel('Infected','FontSize',14); 
zlabel('value','FontSize',14); 
view(-55,25)

% figures font size
aaf  = 12 ;   ttf  = 14;  llf  = 14; 


figure(4)
subplot(1,2,1)
h2=surf(S,I,policy_L_adj'); hold on;
set(h2,'Linestyle','none')
th=title('Policy function (heat map)');set(th,'fontsize',ttf)
xlabel('Suceptible','FontSize',aaf); 
ylabel('Infected','FontSize',aaf); 
view(0,90)

subplot(1,2,2)
h2=surf(S,I, V_adj'*r ); hold on;
set(h2,'Linestyle','none')
th=title('Value function (annual flow units)');set(th,'fontsize',ttf)
xlabel('Suceptible','FontSize',aaf); 
ylabel('Infected','FontSize',aaf); 
zlabel('value','FontSize',aaf); 
view(-55,25)
%}
%Simulation with policy function

Nyears = 2 ; % for simulations and figures

time= 365*Nyears ; 
period=[1:1:time];

Itime=zeros(time,1);
%Xtime(1)=helper.helper;
Itime(1)=I0;

Stime=zeros(time,1); 
%Itime(2:time)=Stime(2:time); % Itime(1)  set as parameter on top
Rtime=zeros(time,1); 
Ntime=zeros(time,1);
Dtime=zeros(time,1); 

% initial conditions 
Stime(1)= 1-Itime(1)-0.02; %??? why the -0.02?
Lck=[];

Stime_Lck=zeros(time,1); 
Itime_Lck=zeros(time,1); 
Rtime_Lck=zeros(time,1); 
Ntime_Lck=zeros(time,1);  
Dtime_Lck=Dtime;

% Same initial conditions
Itime_Lck(1)=Itime(1); 
Stime_Lck(1)=Stime(1); 

Ntime_Lck(1) = Itime_Lck(1)+ Rtime_Lck(1) + Stime_Lck(1);
Lck(1) = interp2(I_mat',S_mat',policy_L,Itime(1),Stime(1)) ;  

for t=2:time 
    % benchmark uncontrolled process
    Stime(t) =  Stime(t-1)*(1 - beta*Itime(t-1)*(1-0)*dt );
    Itime(t) =  Itime(t-1)*(1 - ggamma*dt + beta*Stime(t-1)*(1-0)*dt) ;

    % fatality unctrolled
    fct_phi = ggamma*(varphi + kappa*Itime(t-1));

    Dtime(t) =  Dtime(t-1) + Itime(t-1)*fct_phi*dt ;
    Rtime(t) =  Rtime(t-1) + Itime(t-1)*(ggamma-fct_phi)*dt ;
    Ntime(t) =  Stime(t)   + Itime(t) + Rtime(t); 
    
    % Controlled process
    % fatality ctrolled
    fct_phi_Lck = ggamma*(varphi + kappa*Itime_Lck(t-1) );

    Stime_Lck(t) =  Stime_Lck(t-1)*(1 - beta*Itime_Lck(t-1)*( 1-theta*Lck(t-1))^2*dt) ;
    Itime_Lck(t) =  Itime_Lck(t-1)*(1 - ggamma*dt + beta*Stime_Lck(t-1)*(1-theta*Lck(t-1))^2*dt) ;
    
    Dtime_Lck(t) =  Dtime_Lck(t-1) + Itime_Lck(t-1)*fct_phi_Lck*dt  ;
    Rtime_Lck(t) =  Rtime_Lck(t-1) + Itime_Lck(t-1)*(ggamma-fct_phi_Lck )*dt ;
    
    Ntime_Lck(t) =  Stime_Lck(t)   + Itime_Lck(t) + Rtime_Lck(t);
    
    Lck(t) = interp2(I_mat',S_mat',policy_L,Itime_Lck(t),Stime_Lck(t)) ;  
end

Frac_Pop_Lckdwn= (Stime_Lck + Itime_Lck).*Lck';
%{
figure(7)
plot(period,Stime, 'b' , 'LineWidth', 3); hold on;
plot(period,Stime_Lck, 'r' , 'LineWidth', 3);  hold off;
lh=legend('Susceptible (no  control)','Suceptible (control)');set(lh,'fontsize',18)
xh=xlabel('Time','FontSize',14); 
yh=ylabel('Share of the Population','FontSize',14); 

figure(8)
plot(period,Itime, 'b' , 'LineWidth', 3); hold on;
plot(period,Itime_Lck, 'r' , 'LineWidth', 3);  hold off;
lh=legend('Infected (no  control)','Infected (control)');set(lh,'fontsize',18)
xh=xlabel('Time','FontSize',16); 
yh=ylabel('Share of the Population','FontSize',16); 

figure(9)
plot(period,Rtime, 'b' , 'LineWidth', 3); hold on;
plot(period,Rtime_Lck, 'r' , 'LineWidth', 3);  hold off;
lh=legend('Recovered (no  control)','Recovered (control)');set(lh,'fontsize',18)
xh=xlabel('Time','FontSize',14); 
yh=ylabel('Share of the Population','FontSize',14); 

figure(10)
plot(period,Dtime, 'b' , 'LineWidth', 3); hold on;
plot(period,Dtime_Lck, 'r' , 'LineWidth', 3);  hold off;
lh=legend('Dead (no  control)','Dead (control)');set(lh,'fontsize',18)
xh=xlabel('Time','FontSize',14); 
yh=ylabel('Share of the Population','FontSize',14); 

figure(11)
plot(period,Lck, 'k' , 'LineWidth', 3); hold on;
lh=legend('% in lockdown');set(lh,'fontsize',18)
xh=xlabel('Time (days)','FontSize',16); 
yh=ylabel('Share of Population in Lockdown','FontSize',16); 
vline([15 30 45 60 90 120 150 180 365]);
%}
%%%%%%%%%%%%%%%%%%%%%%%%

tmax =  round(365*Nyears ) ;

aaf  = 12 ;  
ttf  = 14; 
llf  = 14; 
%{
figure(12)

subplot(3,1,1)
plot(period(1:tmax),100*Lck(1:tmax), 'k' , 'LineWidth', 3); hold on;
plot(period(1:tmax),100*Frac_Pop_Lckdwn(1:tmax), 'b' , 'LineWidth', 3);  hold off;
th=title('Lockdown policy');set(th,'fontsize',ttf)
lh=legend('Lockdown rate','% of population'); set(lh,'fontsize',llf)
xh=xlabel('Time (days)','FontSize',aaf); 
yh=ylabel('Population (%)','FontSize',aaf); 
vline([15 30 45 60 90 120 150 180 ]);

subplot(3,1,2)
plot(period(1:tmax),100*Itime(1:tmax), 'b' , 'LineWidth', 3); hold on;
plot(period(1:tmax),100*Itime_Lck(1:tmax), 'r' , 'LineWidth', 3);  hold off;
th=title('Infected (%)');set(th,'fontsize',ttf)
lh=legend('Infected (no  control)','Infected (control)');set(lh,'fontsize',llf)
xh=xlabel('Time (days)','FontSize',aaf); 
yh=ylabel('Population (%)','FontSize',aaf); 

subplot(3,1,3)
plot(period(1:tmax),100*Dtime(1:tmax), 'b' , 'LineWidth', 3); hold on;
plot(period(1:tmax),100*Dtime_Lck(1:tmax), 'r' , 'LineWidth', 3);  hold off;
th=title('Dead (%)');set(th,'fontsize',ttf)
lh=legend('Dead (no  control)','Dead (control)');set(lh,'fontsize',llf)
xh=xlabel('Time','FontSize',aaf); 
yh=ylabel('Population (%)','FontSize',aaf); 
%}
% Checking value function

daily_tot_disc = exp(-(nu+r)*dt) ; 
time = iteration;

for t=1:time
    daily_tot_disc_vec(t,1)=daily_tot_disc^(t-1);
end
 
Dtime_flow(1)=Dtime_Lck(1,1); 
for t=2:time
    Dtime_flow(t,1) = Dtime_Lck(t)- Dtime_Lck(t-1);
    Dtime_flow_nopolicy(t,1) = Dtime(t)-Dtime(t-1);
end

V_initial = interp2(I_mat',S_mat',V,Itime_Lck(1),Stime_Lck(1));
V_last    = interp2(I_mat',S_mat',V_guess,Itime_Lck(iteration),Stime_Lck(iteration)) ;

% present value
PV_Death_loss  = sum(Dtime_flow((1:time))*vsl.*daily_tot_disc_vec);
PV_Lck_loss    = sum((Stime_Lck(1:time)+Itime_Lck(1:time))*w*dt.*Lck(1:time)'.*daily_tot_disc_vec);
PV_loss        = PV_Lck_loss + PV_Death_loss + daily_tot_disc^iteration*V_last  ;

output_flow_loss=r*PV_Lck_loss      ;
flow_valuefunctio_optimal=r*V_initial;
flow_pv_optimal          =r*PV_loss;

% No policy
PV_nopolicy_GDPloss       = 0;
PV_nopolicy_Death_loss =sum(Dtime_flow_nopolicy(1:time)*vsl.*daily_tot_disc_vec);
PV_loss_nopolicy = PV_nopolicy_GDPloss + PV_nopolicy_Death_loss ;
flow_pv_nopolicy = r*(PV_loss_nopolicy );
%{
% Phase diagram pic 
fake_policy=gridsize_I*ones(time,1);
figure(110)
h2=surf(S,I,policy_L_adj'); 
view(2), shading interp; hold on;
plot3(Stime_Lck(1:time),Itime_Lck(1:time),fake_policy(1:time), '-.c' , 'LineWidth', 4); hold on;
plot3(Stime(1:time),Itime(1:time),fake_policy(1:time), '<-r' , 'LineWidth', 4); hold off;
set(h2,'Linestyle','none')
xlabel('Susceptible','FontSize',14); 
ylabel('Infected','FontSize',14);   
zlabel('Policy','FontSize',14); 
view(0,90)
%}
if lockdown_case==0
    Consumption_day = w.*(Stime+Itime+Rtime);
    Labour_day = (Stime+Itime+Rtime);
    Output_day = w.*(Stime+Itime+Rtime);
    Susceptibles_day = Stime; 
    Infected_day = Itime;
    Recovered_day = Rtime;
    Deaths_day = Dtime;
elseif lockdown_case==1
    Consumption = w.*((1-Lck).*(Stime_lck+Itime_lck)+Rtime_lck);
    Labour = ((1-Lck).*(Stime_lck+Itime_lck)+Rtime_lck);
    Output = w.*((1-Lck).*(Stime_lck+Itime_lck)+Rtime_lck);
    Susceptibles = Stime_lck; 
    Infected = Itime_lck;
    Recovered = Rtime_lck;
    Deaths = Dtime_lck;
end
Consumption=day2week(Consumption_day);
Labour=day2week(Labour_day);
Output=day2week(Output_day);
Susceptibles=day2week(Susceptibles_day);
Infected=day2week(Infected_day);
Recovered=day2week(Recovered_day);
Deaths=day2week(Deaths_day);


Consumption_ss= w;
Labour_ss= 1;
Output_ss= w;
Susceptibles_ss=0;
Infected_ss=0;
Recovered_ss=0;
Deaths_ss=0;
%save results!
save('simulated_results.mat','Consumption','Labour','Output','Susceptibles','Infected','Recovered','Deaths');
save('simulated_results_ss.mat','Consumption_ss','Labour_ss','Output_ss','Deaths_ss','Susceptibles_ss','Infected_ss','Recovered_ss');
