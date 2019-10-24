## Code for the analysis of CAV dataset

1. cav-functions.R: Contains all the functions written to analyze CAV (both using the proposed method as well as parametric and non-parametric alternatives.
2. cav-analysis.R: Runs the uniformization-based MCMC analysis of the CAV data. Generates RDS files containing MCMC runs for both parameters and the latent multi-state survival process
3. posterior-plots.R: Generates traceplots, posterior parameter distribution plots, and posterior survival curves.
4. cav-altmodels-analysis.R: Runs the parametric and non-parametric analysis of CAV.  Generates the comparison plots of the estimated survival curves.
