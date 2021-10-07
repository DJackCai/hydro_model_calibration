# Markov-Chain-Monte-Carlo-in-hydrology
This repository contains documents relevant to the use of Markov Chain Monte Carlo (MCMC) algorithm in the calibration of hydrological model to both 
1) obtain sets of optimal parameter values to attain improve model performance and 
2) quantify the associated parameter uncertainty in parameter estimates.

The R codes with "MISDc" covers the R version of the MISDc model designed by Brocca et al. (2011) and Brocca et al. (2012).
This is a parsimonious hydrological model originally to predict flood events. Its layered structure facilitates the assimilation of satellite data 
or satellite-based soil moisture product to the model, as well as joint calibration using satellite soil moisture data.

A modified version of model that is written as a hydromad (HMD) object in hydromad package in R to help conduct efficient sensitivity analysis using 
the Sensitivity package without the need to explicitly write up the sampling algorithm. It will also be convenient to utilise the parameter optimisation, 
model evaluation and visualisation of outputs in hydromad package.


