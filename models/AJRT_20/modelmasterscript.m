job_name='Baseline_pi4med_C1L1B4_d0510_05';

pi1_shr_target=1/6;
pi2_shr_target=1/6;
pi3_shr_target=4/6;
pi4_level=0.000001;

helper=load('inf_ini.mat');
I0in=helper.helper;
%I0in=0.001;

alphain = 0.53; % Home bias
phiin = 0.8; % infected labor loss
delta_muin=0.5;
delta_nuin=1.0;

alphain = 0.53; % Home bias
phiin = 0.8; % infected labor loss

pid_accel = 0.05;

Run_Benchmarks;
%Run_Planner;
%Run_Nash;
