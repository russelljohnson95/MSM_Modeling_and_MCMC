# MSM_Modeling_and_MCMC
This repository presents a workflow for performing Markov Chain Monte Carlo (MCMC) analyses with OpenSim musculoskeletal models. It also includes an implementation of the MCMC algorithm with a mass-spring-damper model. 

## Prerequisites
This project _requires_ the appropriate set up of: 
1. DRAM MCMC MATLAB Algorithm (see below)  
2. [Parallel Computing Toolbox](https://www.mathworks.com/products/parallel-computing.html) in MATLAB
3. For the Arm 16 Model,  OpenSim 3.3 with MATLAB Scripting API

**_IMPORTANT:_** You must have [OpenSim 3.3](https://simtk.org/frs/download_confirm.php/file/4580/OpenSim-3.3.0-win64VC13P.exe?group_id=91) installed  on your computer to run the model associated with this code, with the MATLAB Scripting API set up. [Click here to learn how to set up the API](https://simtk-confluence.stanford.edu/display/OpenSim33/Scripting+with+Matlab). 

This set up does *NOT* work with OpenSim 4.0 or later due to scripting changes. 

{TODO: add reference to paper}

## DRAM_MCMC_Matlab: 
The MATLAB implementation of the DRAM MCMC Algorithm used with this project. Users will need to make sure the files in this folder are in their MATLAB path to run MCMC for either the Mass Spring Damper Model or Arm 16 Model.

These files were obtained from Marko Laine's [`mcmcstat` Repository](https://github.com/mjlaine/mcmcstat) on Sept. 3rd, 2019, and are documented further [here](https://mjlaine.github.io/mcmcstat/#orgcdeadeb).

They are presented here with light edits, mostly for formatting purposes. 

## Mass_Spring_Damper_Model: 
This folder runs an MCMC analysis to recover the parameters of a mass-spring-damper system, with a variable-stiffness spring. See Supplementary Material paper for full description of model and algorithm. 

### Project Files Description

* **Mass_Spring_Damper_mcmc_parallel:** This is the _main_ script that sets up and then calls the MCMC algorithm to run. It will also plot some of the results afterwards. 

* **Mass_Spring_Damper_Simulator:** If you want to generate new reference data, this script that will generate kinematics for the mass-spring-damper system. In this script, you can set up the "true" parameters of the system, and simulate the oscillation of the mass over a set period of time. 

* **deriv:** A function that finds the derivative of a time series trajectory, used to find velocity from position data of the mass. 

* **effective_sample_sizeCalc:** A script that calculates the R-hat and effective sample size of the MCMC results.

* **Figure_Manuscript_fromChain:** A script that plots seen in the manuscript from the MCMC results.

* **rank_plot_SMD:** A function that calculates and plots the rank of the results, called by the script `Spring_mass_Damper_mcmc_parallel.m`.

* **rank_plot_SMD_param:** A function that calculates and plots the rank of the results for a single parameter at a time, called by the script `Figure_Manuscript_fromChain.m`.

#### **Also Included:**
* Position data (with and without noise) used to track mass oscillations in MCMC and used for plotting.
* Sample MCMC results `chain_20210526T111546.mat` and `results_20210526T111546.mat`.

## Arm_16_Model 
This folder runs an MCMC analysis to find the plausible muscle forces for a elbow flexion model (`arm16_millard_rigidtendon.osim`) in OpenSim. This uses compact radial basis functions (CRBF's) to represent muscle excitations, that then get forward integrated using OpenSim to simulate the motion of the model during each iteration. 

### Project Files Description

* **Arm16_CRBF_6musc_mcmc_parallel:** This is the _main_ script that sets up and then calls the MCMC algorithm to run. It will also plot some of the results afterwards. 

* **Arm16_Figure_FromChains:** A script that plots seen in the manuscript from the MCMC results for the elbow model.

* **Arm16_SimManager_controls_CRBF_6musc:** This function uses the muscle excitation signals to simulate the motion of the Elbow Model via OpenSim. It returns the positions and velocities of the elbow joint angle to the script `Arm16_CRBF_6musc_parallel.m`.

* **Arm16_SimManager_controls_CRBF_6musc_wForce:** This function uses the muscle excitation signals to simulate the motion of the Elbow Model via OpenSim. It returns the positions and velocities of the elbow joint angle AND the muscle forces from each of the muscles in the model to the script `Arm16_Figure_FromChains.m`.

* **CRBF_excit:** This function takes the CRBF amplitudes from the MCMC run and calculates the muscle excitation signal for a muscle. It uses compact radial basis functions then transforms it via a logit transform to put them onto a scale from 0 to 1. It gets called by both `Arm16_SimManager_controls_CRBF_6musc.m` and `Arm16_SimManager_controls_CRBF_6musc_wForce.m`.

* **Effective_Sample_Size_Elbow:**  A script that calculates the R-hat and effective sample size of the MCMC results for the elbow model.

* **rank_plot_Arm16_10CRBVs_fxn:** A function that helps calculate the rank plots of each muscle, called by `Arm16_Figure_FromChains.m` and `Arm16_CRBF_6musc_parallel.m`.

* **rank_plot_muscle_10node:** A function that calculates the rank plots of each muscle, called by `rank_plot_Arm16_10CRBVs_fxn.m`.

* **arm16_millard_rigidtendon.osim:** Arm model used for this project with 6 muscles and 1 mechanical degree-of-freedom at the elbow. 

* **arm16_pert4_:** Reference data for the MCMC analysis. 

TODO: Not included (right now) are the sample data from the MCMC run for this project - the file size is too big for GitHub... 
