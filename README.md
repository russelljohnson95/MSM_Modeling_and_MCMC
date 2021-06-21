# MSM_Modeling_and_MCMC
Workflow for performing Markov Chain Monte Carlo analyses on OpenSim musculoskeletal models, and includes an implementation of the MCMC algorithm with a mass-spring-damper model. 

# DRAM_MCMC_Matlab: 
This is the Matlab implementation of the DRAM MCMC algorithm used with this project. It was obtained from this github: https://github.com/mjlaine/mcmcstat and documented further here: https://mjlaine.github.io/mcmcstat/#orgcdeadeb. Presented here with light edits, mostly for formatting purposes. Users will need to make sure these files are in their Matlab path to run MCMC for either the . 

# Mass_Spring_Damper_Model: 
This folder runs an MCMC analysis to recover the parameters of a mass-spring-damper system, with a variable-stiffness spring. 

**deriv:** A function that finds the derivative of a time series trajectory, used to find velocity from position data of the mass. 

**effective_sample_sizeCalc:** A script that calculates the R-hat and effective sample size of the MCMC results 

**Figure_Manuscript_fromChain:** A script that plots seen in the manuscript from the MCMC results 

**rank_plot_SMD:** A function that calculates and plots the rank of the results, called by the script Spring_mass_Damper_mcmc_parallel.m

**rank_plot_SMD_param:** A function that calculates and plots the rank of the results for a single parameter at a time, called by the script Figure_Manuscript_fromChain.m

**Spring_mass_Damper:** A script that will generate data for the mass-spring-damper system. In this file, you can set up the "true" parameters of the system, and simulate the oscillation of the mass over time. 

**Spring_mass_Damper_mcmc_parallel:** This is the main script that sets up and then calls the MCMC algorithm to run. It will also plot some of the results afterwards. 
