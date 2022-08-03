% Contacts at work: taken from "Projecting contact matrices in 177 
% geographical regions: an update and comparison with empirical data
% for the COVID-19 era" (2020) by Kiesha Prem, Kevin van Zandvoort,
% Petra Klepac, Rosalind M Eggo, Nicholas G Davies, Alex R Cook, Mark Jit
%
% US work contacts are the average across the Polymod countries adjusted
% by the labor force participation in each age group.

%  0-4  5-9  10-14 15-19 20-24  25-29    30-34   35-39  40-44    45-49  50-54  55-59   60-64   65-69   70-74   75+
Contact_work_USA = [...
    0    0    0 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000
    0    0    0 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000
    0    0    0 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000
    0    0    0 0.12040 0.20985 0.13723 0.12792 0.11888 0.13380 0.10898 0.06976 0.03879 0.01791 0.00510 0.00009 0.00000
    0    0    0 0.14144 0.59548 0.64892 0.54157 0.59653 0.45503 0.37650 0.31335 0.16870 0.10728 0.01393 0.00094 0.00026
    0    0    0 0.14136 0.64437 1.23383 0.88285 0.86171 0.80428 0.58255 0.52196 0.28258 0.16271 0.02639 0.00083 0.00025
    0    0    0 0.07639 0.44647 0.83702 1.09531 0.96030 0.86946 0.73335 0.48217 0.32966 0.15235 0.02704 0.00034 0.00029
    0    0    0 0.15232 0.37025 0.78595 0.84157 1.16257 1.14267 0.81693 0.64836 0.30231 0.10794 0.02032 0.00076 0.00050
    0    0    0 0.09611 0.44029 0.78533 0.92624 0.99095 1.22657 1.01817 0.79766 0.32428 0.15990 0.02379 0.00085 0.00107
    0    0    0 0.12411 0.31001 0.60990 0.78194 0.85944 0.89830 0.90574 0.63593 0.37919 0.13158 0.02669 0.00089 0.00050
    0    0    0 0.09930 0.27601 0.65725 0.76609 0.79753 1.09108 1.08727 0.90036 0.51185 0.16820 0.01868 0.00094 0.00080
    0    0    0 0.08054 0.22447 0.46947 0.65271 0.61526 0.79466 0.64762 0.65446 0.47987 0.19151 0.01958 0.00048 0.00048
    0    0    0 0.02226 0.18184 0.36853 0.40594 0.47961 0.52393 0.52853 0.44017 0.39815 0.13389 0.02344 0.00048 0.00086
    0    0    0 0.01392 0.04447 0.12966 0.11957 0.04819 0.10963 0.09766 0.12191 0.07699 0.05380 0.00565 0.00050 0.00024
    0    0    0 0.00080 0.00119 0.00175 0.00295 0.00333 0.00296 0.00254 0.00169 0.00327 0.00151 0.00034 0.00002 0.00003
    0    0    0 0.00000 0.00329 0.00386 0.00388 0.00696 0.00230 0.00085 0.00085 0.00000 0.00000 0.00021 0.00000 0.00000 ...
    ];

%  0-4  5-9  10-14 15-19 20-24  25-29    30-34   35-39  40-44    45-49  50-54  55-59   60-64   65-69   70-74   75+
Lfp_USA = [0	0	0	0.30611	0.70694	0.82237	0.82653	0.82695	0.82919	0.82072	0.79186	0.72625	0.5779	0.19975	0.01 0.01];


% if the model does not feature labor force participation explicitly, the
% aggregate contacts are:
% aggregating to 3 age groups: young (0-19), middle (20-64), old (65+)
Contact_work_agg_3 = ...
[sum(sum(Contact_work_USA(1:4,1:4),2).*popsharebyage(1:4)')/sum(popsharebyage(1:4)),...
 sum(sum(Contact_work_USA(1:4,5:13),2).*popsharebyage(1:4)')/sum(popsharebyage(1:4)),...
 sum(sum(Contact_work_USA(1:4,14:16),2).*popsharebyage(1:4)')/sum(popsharebyage(1:4));...
 sum(sum(Contact_work_USA(5:13,1:4),2).*popsharebyage(5:13)')/sum(popsharebyage(5:13)),...
 sum(sum(Contact_work_USA(5:13,5:13),2).*popsharebyage(5:13)')/sum(popsharebyage(5:13)),...
 sum(sum(Contact_work_USA(5:13,14:16),2).*popsharebyage(5:13)')/sum(popsharebyage(5:13));...
 sum(sum(Contact_work_USA(14:16,1:4),2).*popsharebyage(14:16)')/sum(popsharebyage(14:16)),...
 sum(sum(Contact_work_USA(14:16,5:13),2).*popsharebyage(14:16)')/sum(popsharebyage(14:16)),...
 sum(sum(Contact_work_USA(14:16,14:16),2).*popsharebyage(14:16)')/sum(popsharebyage(14:16))];

% aggregate
Contact_work_agg_1 = ...
[sum(sum(Contact_work_USA(1:16,1:16),2).*popsharebyage(1:16)')/sum(popsharebyage(1:16))];


