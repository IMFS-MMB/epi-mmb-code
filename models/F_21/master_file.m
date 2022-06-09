%% Master file for Faria-e-Castro (2021) "Fiscal Policy during a Pandemic"
% Journal of Economic Dynamics & Control
% Federal Reserve Bank of St. Louis

% Baseline version of the model
cd main
calibrate_model
multipliers_normal_times
crisis_experiments

% APPENDIX B - Keynesian Supply Shocks as in Guerrieri-Lorenzoni-Straub-Werning 2020
cd ../keynesian_ss
calibrate_model
multipliers_normal_times
crisis_experiments

% APPENDIX C - Version of the model where borrowers also consume
cd ../borr_consume
calibrate_model
multipliers_normal_times
crisis_experiments