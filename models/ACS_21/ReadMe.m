%%%%%%%%%%%% ACS (2021). How to replicate paper Figures 3 and 4

%%%%%%%%%%%% 1) Run SS_banchmark to store steady state.

%%%%%%%%%%%% 2) Run nonfarm_YoverH (using hpfilter and us_YoverH_nnfarm) to
%%%%%%%%%%%% generate and store empirical path of labor productivity.

%%%%%%%%%%%% 3) Run acs_benchmark on dynare 4.4.3 (use plotter in
%%%%%%%%%%%% plot_agg_result, go_calibrate and calibrate_pi) to store
%%%%%%%%%%%% results and simulate benchmark response.

%%%%%%%%%%%% 4) Run plot_paper for benchmark figure 3 in the paper.

%%%%%%%%%%%% 5) Run plot_prod_paper for decomposition productivity (Fig.4).