This packet offers replication codes for "Social Distancing and Supply Disruptions in a Pandemic"
by Martin Bodenstein, Giancarlo Corsetti and Luca Guerrieri


Link paper and codes: 
http://www.lguerrieri.com/published-research.html


The programs require a local installation of Dynare --- the original code used version 4.5.7
The programs also require the OccBin toolbox, which can be downloaded from http://www.lguerrieri.com/occbin_20140630.zip

Before running the programs modify setpath_dynare.m to point to the local installation of Dynare and OccBin.

This packet replicates the figures the model-based figures, starting with Figure 3. Figures 1 and 2 can be replicated using a companion packet 
for the regression analysis.
--------------------------------------------------------------------------------------------------------------------------
EPI-MMB: For the replication to work with Dynare 4.6 and 5.1, we set in = 1e-12 in one_sector_vcu_binding.mod and two_sector_vcu_binding.mod , the original code set in=0
--------------------------------------------------------------------------------------------------------------------------
Figure 3, 4, and 5 --- call_two_sector_path.m

Figure 6  --- call_two_sector_social_distancing2.m, set case_switch=1 on line 28 

Figure 7  --- call_two_sector_social_distancing2.m, set case_switch=5 on line 28

Figure A.1 --- call_two_sector_path_unemp.m

Figure A.2 --- call_two_sector_beta_sensitivity.M

Figure A.3 --- call_two_sector_social_distancing2.m, set case_switch=7 on line 28

Figure A.4 --- call_two_sector_social_distancing2.m, set case_switch=8 on line 28
