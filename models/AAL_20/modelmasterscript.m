% Simple Planning Problem - finite differences
% Fernando Alvarez, David Argente, Francesco Lippi --  August 2020

%clc
%clear all
%
%close all
easy=1;
% TABLE 1
beta    = 0.13 * 365 ;  % propagation rate Sdot = beta S * X / N
ggamma  = 1/18 * 365   ;  % avg 20 days to get out of infection state (either recovered or dead)
theta   = 0.5   ;      % effectiveness of transmission    
alpha   = 0.02    ;       % parameter of the Testing-Tracing cost c(T ,S,X) = alpha/2 (T (S+X)/X)^2;
zeta    = 1;             % effectiveness of tracing, zeta=0, random sampling, zeta=1 perfect tracing

% parameter for death rate function  Daily Death (t) =  (ggammad +  kappa*I(i)^alphad) * I(t) ;  
varphi  = 0.0068  ; %   ; % baseline case mortality is 0.68%  adjusted by duration of desease
kappa   = 0.034 ; %   ; % benchmark 0.034      % non-linear death rate (units adjusted by duration of desease), kappa=0 to allow for TTQ protocol increases by 1pp for every 10pp increase in Infected

% Lockdown, testing, tracing, quarantine
tau = 0    ; % tau =1 (test) tau=0 (no test), tau=0 to allow for TTQ protocol
l_max = 0.7;  % bench = 0.75;
T_max = 0   ;   % maximun testing-tracing rate, T_max = 1 to allow for TTQ protocol

% Welfare cost parameters
w   =      1    ; % normalized wage
% value statistical life
vsl = 40 * w    ;

%Start 300 grid points
N_S = 300 ;
k = 5; % ratio of Delta_S = k Delta_X
const = 1; % quadratic constant

% tau = 0
row1 = [beta, ggamma, 0.3, alpha, zeta, varphi, kappa, tau, l_max, T_max, w, vsl, N_S, k,const];
row2 = [beta, ggamma, theta, alpha, zeta, varphi, kappa, tau, l_max, T_max, w, vsl, N_S, k,const]; % benchmark
% alternative vls
row3 = [beta, ggamma, theta, alpha, zeta, varphi, kappa, tau, l_max, T_max, w, 50 * w , 400, k,const];
row4 = [beta, ggamma, theta, alpha, zeta, varphi, kappa, tau, l_max, T_max, w, 70 * w , 650, k,const];
% constant fatality rate kappa=0
row5 = [beta, ggamma, 0.3, alpha, zeta, varphi, 0, tau, l_max, T_max, w, vsl, N_S, k,const];
row6 = [beta, ggamma, 0.5, alpha, zeta, varphi, 0, tau, l_max, T_max, w, vsl, N_S, k,const];
% antibody test tau=1
row7 = [beta, ggamma, theta, alpha, zeta, varphi, kappa, 1, l_max, T_max, w, vsl, N_S, k,const];
row8 = [beta, ggamma, theta, alpha, zeta, varphi, kappa, 1, l_max, T_max, w, 50 * w, 400, k,const];
row9 = [beta, ggamma, theta, alpha, zeta, varphi, kappa, 1, l_max, T_max, w, 70 * w, 650, k,const];
row10 = [beta, ggamma, theta, alpha, zeta, varphi, 0, 1, l_max, 1, w, vsl, N_S, k,const];
row11 = [beta, ggamma, theta, alpha, zeta, varphi, 0, 1, l_max, 1, w, 50 * w, 400, k,const];
row12 = [beta, ggamma, theta, alpha, zeta, varphi, 0, 1, l_max, 1, w, 70 * w, 650, k,const];
% less pessimistic parameter values
row13 = [0.1 * 365, ggamma, theta, alpha, zeta, varphi, kappa, tau, l_max, T_max, w, vsl, N_S, k,const];
row14 = [beta, ggamma, theta, alpha, zeta, 0.005, kappa, tau, l_max, T_max, w, vsl, N_S, k,const];
%quadratic losses
row15 = [beta, ggamma, theta, alpha, zeta, varphi, kappa, tau, l_max, T_max, w, vsl, N_S, k, 5.5347];
row16 = [beta, ggamma, theta, alpha, zeta, varphi, kappa, 1, l_max, T_max, w, vsl, N_S, k, 6.9245];
    
    

parameters = [row1;row2;row3;row4;row5;row6;row7;row8;row9;row10;row11;row12;row13;row14;row15;row16];
if easy==1
    COVID19_interp
else
    R=2;%benchmark
    SolvePlanningProblem
end
%{
selector = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16];

for i = 1:length(selector)
    R=selector(i)
    tic
    if R<15
      SolvePlanningProblem
    else
      SolvePlanningProblem_quadratic 
    end
    toc
end

clearvars -except TABLE1

TABLE1
%}