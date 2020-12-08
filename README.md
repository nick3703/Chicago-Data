# Chicago-Data

Data to Reproduce Analysis of "A Class of Spatially Correlated Self-Exciting Statistical Models", by Clark and Dixon

BayesINGARCHChi.R contains Rstan code for conducting Bayesian inference for an INGARCH model for Chicago Burglaries, equation 26 of the manuscript.  


BayesSPINGARCHChi.R contains Rstan code for conducting Bayesian inference for a SPINGARCH model for Chicago Burglaries, equation 24 of the manuscript.

Both of these files are self-contained with both Stan code as well as R script to pull in data. In order to use the code for other locations, you must provide the spatial neighborhood matrix, a matrix of crimes where each row is a location and each column is a day, and any covariates used for large-scale spatial structure.

