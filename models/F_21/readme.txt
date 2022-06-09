# fiscal_policy_pandemic
Replication code for Faria-e-Castro (2021) "Fiscal Policy during a Pandemic", forthcoming in the *Journal of Economic Dynamics & Control*

A stable version of the working paper version can be found [here](https://research.stlouisfed.org/wp/more/2020-006). 
If not, please check my [website](www.fariaecastro.net).

Code written on MATLAB R2020a, requires Dynare (works with v4.6.3). 

# Instructions
Download this folder and run `master_file.m`. This should generate all the figures and most tables in the draft.

There are three folders, all with the same structure:
1. Files in the `main` folder correspond to the main analysis inthe body of the paper.
2. Files in the `keynesian_ss` folder correspond to the Keynesian supply shock analysis of Appendix B.
3. Files in the `borr_consume` folder correspond to the model extension presented in Appendix C. 

The structure in each folder is as follows:
1. `calibrate_model.m` sets up the model's calibration, loads the fiscal policy data, and estimates the sequence of pandemic shocks to match the path of the unemployment rate. It generates Figures 1-2. 
2. `multipliers_normal_times.m` computes fiscal multipliers for "normal times", in the absence of shocks. It generates Table 4.
3. `crisis_experiments.m` computes IRFs for the pandemic shock as well as for fiscal policy, as well as estimating fiscal multipliers under the pandemic and for the CARES Act. It generates Figures 3-12, and Tables 3 and 5. 
