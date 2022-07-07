Replication Package
***************************************************************************


Title: The Euro Area's pandemic recession: A DSGE-based interpretation
Authors: Roberta Cardani, Olga Croitorov, Massimo Giovannini, Philipp Pfeiffer, Marco Ratto, Lukas Vogel

Corresponding author: Marco Ratto

Link to paper:
https://ec.europa.eu/jrc/en/publication/eur-scientific-and-technical-research-reports/euro-area-s-pandemic-recession-dsge-interpretation

***************************************************************************

The package has been tested with:
MATLAB 2017b
Dynare 4.6.4 (https://www.dynare.org/release/) 
  
***************************************************************************


To replicate the paper results, please proceed as follows:

1)  Add your Dynare path to the Matlab search path:

        addpath [YOUR DYNARE PATH] \dynare_4.6.4\matlab\

2)  gemc.dyn runs and replicates* the irfs:

    - gemc\Output (IRFs using deterministic simulations, see Figures 2-4)

*NOTE: in the paper we use occbin toolkit which will be available only in the upcoming dynare version (4.7), hence we provide IRFs using deterministic simulations, rather than stochastic, which numerically slightly differ.

- Run the model by typing in the command window:
    dynare gemc console nointeractive



