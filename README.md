# Codes for the *A model-based approach for historical borrowing, with an application to neovascular age-related macular degeneration* paper

This repository contains the codes used to run and reproduce the analyses done in the *A model-based approach for historical borrowing, with an application to neovascular age-related macular degeneration* paper.

1. Construct a Monolix dataset for the analysis, containing data from six historical nAMD trials (ANCHOR, AVENUE, HARBOR, MARINA, PIER, STAIRWAY). Data cannot be shared externally. 
- `R/PlotSummaryData.R` provides code for exploratory plots summarizing the data.
2. Fit two non-linear mixed effect models (referred to as the PK/BCVA and PK/BCVA-TRIAL models) to the data. Code:
- `Monolix/StructModels/PKBCVAtrans2.txt` contains the Monolix structural model, which is the same for both the PK/BCVA and PK/BCVA-TRIAL models.
- `Monolix/ModelsRes/pkbcva_trans2.mlxtran` and `Monolix/ModelsRes/pkbcva_trans2_trial.mlxtran` are the Monolix [mlxtran](https://mlxtran.lixoft.com/) files for the PK/BCVA and PK/BCVA-TRIAL models. mlxtran files fully specify a model in Monolix.
- `Monolix/ModelsRes/pkbcva_trans2/` and `Monolix/ModelsRes/pkbcva_trans2_trial/` are the Monolix results folders for the PK/BCVA and PK/BCVA-TRIAL models
3. Validate the models, producing goodness of fit plots based on both the entire data-set and a cross-validation strategy at the trial level (of note, PIER is never excluded in the cross-validation procedure). Codes used are:
-  Model parameter for the PK/BCVA and PK/BCVA-TRIAL models are fitted excluding one parameters. The mlxtran files are in `Monolix/ModelsRes/pkbcva_trans2_cv_*.mlxtran` and the results folders are `Monolix/ModelsRes/pkbcva_trans2_cv_*/`. `*` denotes the excluded trial, *i.e.* one of ANC, AVE, HAR, MAR, and STA.
- `R/PlotModels.R` provides standard goodness of fit plots based on the full model and data.
- `R/OOSSim.R` runs the out of sample predictions for all individuals in a out-of-sample trial, via *simulx*, based on the *.mlxtran* out-of-sample files.
- `R/OOSSimID.R` runs the out of sample predictions for a single individual, via *simulx*, based on the *.mlxtran* out-of-sample files.
- `R/PlotOOS.R` plots goodness of fit plots for the out-of-sample predictions.
4. Simulate the baseline characteristics and dosing regimen information of two hypothetical future nAMD trials, one investigating a treatment-naive and another a pre-treated population. Codes:
- `R/MasterTrialSim.R` defines all scenarios for trial simulation (different sample and effect sizes, useful in analysis step number 6), as well as baseline characteristics and dosing information.
5. Obtain model-informed priors (based on the PK/BCVA and PK/BCVA-TRIAL models), and compare them with other priors typically used for historical borrowing. Codes:
- `R/PriorSim.R` based on the baseline characteristics + dosing info (defined in `R/MasterTrialSim.R`), simulates longitudinal profiles of individuals in treatment-naive and pre-treated trials.
- `R/PlotsPrior.R` based on longitudinal simulations, derive model-informed prior. Also compares model-informed prior with standard meta-analytic and power priors. The power prior is implemented in the `Stan/pow_pri.stan` file.
6. Evaluate the operating characteristics of the hypothetical pre-treated trial, comparing the performance of different priors and sample sizes.
- `R/SimData.R` creates simulated data for all the scenario of interest.
- `R/FunDataSim.R` provides helper function for simulating the future trial.
- `R/RunTrialSim.R` runs trial simulations, *i.e.* estimate trial parameters based on different priors and sample sizes. The trial simulations are based upon:  1) `Stan/trial_anal.stan` which carries out Bayesian inference for current trial data (assuming the normal cross-sectional analysis model), given specification of a prior (either weakly informative, model informed, model informed power prior and model informed mixture prior). 2) `R/FunTrialAnalysis.R` a collection of function to analyse the results of a Bayesian trial, and report the most important outputs.
- `R/PlotsTrialSim.R` provides plots summarizing the trial simulations.

Of note, a Monolix licence is needed in order to run these codes.
