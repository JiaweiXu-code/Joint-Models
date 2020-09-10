# Software for paper "Bayesian Design of Clinical Trials Using Joint Models for Longitudinal and Time-to-event Data" by Xu et al.

-------------------------------------------------------------------------------------------------------------------------------------

All programs are setup to be executed on a Linux computing cluster using R (3.6.0) or SAS (9.4). The paths referenced in all programs will need to be updated for the code to work. Once all paths are updated, one can use the SLURM scheduler shell scripts to submit jobs on a SLURM-based computing cluster. 


--------------------------------------------------- RUN ORDER ------------------------------------------------

[1] ssd.sh -- Determine the desired sample size such that a specified number of events is obtained in a specified interval of time on average. This process is performed based on massive (e.g., 1,000) simulations and calls the R program "ssd.R" which performs a single simulation. The "ssd.R" code requires the following inputs:

   (1)  integer      - seed     - seed for a single simulation (e.g., 1-1,000 for sample size determination and 1-4,000 for data generation)
   
   (2)  integer      - v        - number of events 
   
   (3)  double       - k        - ratio between number of enrolled patients and event total 
   
   (4)  integer      - nlinear  - number of intervals for peicewise linear trajectory 
   
   (5)  integer      - measures - number of scheduled longitudinal measurements 
   
   (6)  integer      - npieces  - number of intervals for piecewise constant baseline hazard function 
   
   (7)  double       - eta      - period of time (years) for patient enrollment 
   
   (8)  double       - maxt     - maximum follow-up time (years) under ideal cases (i.e., no dropout or censoring) 
   
   (9)  double       - censor   - period of time (years) for dropout 
   
   (10) double       - censorp  - dropout probability 
   
   (11) double       - interval - length of interval for discretization to simulate survival time 
   
   (12) vector (dbl) - pknots   - knots placement for piecewise linear trajectory 
   
   (13) vector (dbl) - knots    - knots placement for piecewise constant baseline hazard function 
   
   (14) double       - sigma    - standard deviation of measurement error
   
   (15) vector (dbl) - p        - allocation probabilities for treatment and baseline covariates 
   
   (16) double       - mSigma   - standard deviation or covariance matrix of random effects
   
   (17) vector (dbl) - gamma    - regression coefficients of intercept, treatment indicator and baseline covariates in longitudinal process 
   
   (18) vector (dbl) - plinear  - coefficients for time covariates in longitudinal process
   
   (19) vector (dbl) - pinter   - coefficients for interactions between treatment and time covariates in longitudinal process
   
   (20) vector (dbl) - llambda  - piecewise constant baseline hazards in log scale
   
   (21) double       - beta     - association parameter
   
   (22) vector (dbl) - alpha    - direct effects of treatment and baseline convariate on time-to-event endpoint
   
   (23) vector (dbl) - tpoints  - scheduled time points (years) at which longitudinal outcomes are measured
   
   (24) double       - mk       - ratio between standard deviations of random slope and random intercept when used in random effect misspecification

OUTPUTS: Seed, trial duration, sample size, event total



[2] mean.R -- Compute the average trial duration, sample size and event total based on the 1,000 simulations.

[3] data_array.R -- Simulate the longitudinal and time-to-event data by calling the R program "dataset.R" which performs a single simulation and fit the joint models to the corresponding dataset by calling SAS programs starting with "JM" and "SJM". The "dataset.R" code requires the same inputs as for "ssd.R". The example code for "dataset.R" simulates 50 datasets in a single run and event total automatically increases once the number of datasets hits the requirement (e.g., 4,000).

OUTPUTS (from "dataset.R" of simulated data):

   (1)  integer - ID        - patient ID

   (2)  double - measure    - longitudinal measures

   (3)  double - time       - measurement time (years) since enrollment

   (4)  binary - trt        - treatment indicator (1 = treated, 0 = control)

   (5)  double - toltime    - measurement time (years) since trial started

   (6)  binary - nodegrp    - lymph node group, binary baseline covariate (1 = more severe group, 0 = less severe group)

   (7)  integer - sim       - ID for simulated data

   (8)  double - s          - follow-up time (years) since enrollment

   (9)  double - sr         - follow-up time (years) since trial started

   (10) double - r0         - enrollment time (years)

   (11) binary - censorship - indicator for censoring (1 = event, 0 = censored)
   
   
   
[4] power.sas -- Estimate the Bayesian type I error rate or power by event total using the estimated avarage hazard ratios based on joint models.

[5] logrank.sas -- Estimate the type I error rate or power by event total based on log-rank test.

[6] PH.R -- Estimate the type I error rate or power by event total based on piecewise proportional hazard model. The average parameter estimates are also computed.



----------------------------------------- Folder MainResults-Figures2&3 --------------------------------------

Description: This folder contains data generation and model fitting programs for results in Section 3.1. 

dataset.R -- Simulate data where patient-level heterogeneity is based on a random intercept. A 4-component piecewise linear trajectory and a 5-component piecewise constant baseline hazard function is assumed. 

JM.sas -- Fit the proposed joint model with a random intercept, a 4-component piecewise linear trajectory and a 5-component piecewise constant baseline hazard function. 

SJM.sas -- Fit the simplified joint model that omits the random effects with a 4-component piecewise linear trajectory and a 5-component piecewise constant baseline hazard function. 

The ratios between sample size and event total used in the paper for data generation can be found below:

Table1: Ratios (k) used for data generation in Figure 2.

| α(first element in alpha) = 0.0 | γ (pinter) = {0.0,0.0,0.0,0.0} | γ (pinter) = {0.2,0.2,0.2,0.2} | γ (pinter) = {0.4,0.3,0.2,0.1} |
| --- | --- | --- | --- |
| β(beta) = -0.15 | 2.75 | 2.80 | 2.80 |
| β(beta) = -0.30 | 2.55 | 2.70 | 2.70 |
| β(beta) = -0.45 | 2.40 | 2.60 | 2.60 |

Table2: Ratios (k) used for data generation in Figure 3.

| α(first element in alpha) = -0.2 | γ (pinter) = {0.0,0.0,0.0,0.0} | γ (pinter) = {0.2,0.2,0.2,0.2} | γ (pinter) = {0.4,0.3,0.2,0.1} |
| --- | --- | --- | --- |
| β(beta) = -0.15 | 2.95 | 3.05 | 3.00 |
| β(beta) = -0.30 | 2.75 | 2.90 | 2.90 |
| β(beta) = -0.45 | 2.70 | 2.80 | 2.75 |


----------------------------------------- Folder Trajectory-Misspecification --------------------------------------

Description: This folder contains data generation and model fitting programs for results in Section 3.2. 

dataset-leveloff.R  -- Simulate data where the trajectory for treated group increases initially but levels off over time. A 6-component piecewise linear trajectory is assumed.

dataset-concave.R   -- Simulate data where the trajectory for treated group increases initially and subsequently decreases. A 6-component piecewise linear trajectory is assumed.

JM3.sas             -- Fit a joint model with random intercept. A 3-component piecewise linear trajectory is assumed.

JM6.sas             -- Fit a joint model with random intercept. A 6-component piecewise linear trajectory is assumed.

SJM3.sas            -- Fit a simplified joint model without random effects. A 3-component piecewise linear trajectory is assumed.

SJM6.sas            -- Fit a simplified joint model without random effects. A 6-component piecewise linear trajectory is assumed.

Trajectory-curves.R -- Compute the correct and average fitted trajectories for treated and control groups.

Note: Patient-level heterogeneity is based on a random intercept. A 5-component piecewise constant baseline hazard function is used for both data generation and model fitting.


--------------------------------------------- Folder RE-Misspecification ------------------------------------------

Description: This folder contains data generation and model fitting programs for results in Section 3.3. 

dataset.R  -- Simulate data where patient-level heterogeneity is based on a random intercept and a random slope. The two random effects are assumed to be independent. Ratio between standard deviations of random slope and random intercept "mk" varies in {0.25,0.50}.

JM.sas     -- Fit a joint model with a random intercept.

SJM.sas    -- Fit a simplified joint model without random effects.

TrueJM.sas -- Fit a joint model with a random intercept and a random slope.

Note: A 4-component piecewise linear trajectory and a 5-component piecewise constant baseline hazard function is used for both data generation and model fitting.


--------------------------------------------- Folder DropoutProb-AppendixB ----------------------------------------

Description: This folder contains data generation and model fitting programs for results in Appendix B of the Supplementary Materials. 

dataset.R  -- Simulate data where patient-level heterogeneity is based on a random intercept. The dropout probability "censorp" varies in {0.05,0.10,0.20}.

JM.sas     -- Fit a joint model with a random intercept.

SJM.sas    -- Fit a simplified joint model without random effects.

Note: A 4-component piecewise linear trajectory and a 5-component piecewise constant baseline hazard function is used for both data generation and model fitting.


------------------------------------------- Folder SurvivalCurves-AppendixD ---------------------------------------

Description: This folder contains programs to plot Figure2 in Appendix D of the Supplementary Materials. 

survival_curves.R -- Compute the estimated survival curves based on parameter estimates from proposed joint model.

Surv_plot.sas     -- Plot survival curves for treated and control groups using the estimated survival probabilities from "survival_curves.R".


------------------------------------------- Folder Heterogeneity-AppendixE ----------------------------------------

Description: This folder contains data generation and model fitting programs for results in Appendix E of the Supplementary Materials. 

dataset.R     -- Simulate data where patient-level heterogeneity is based on a random intercept. The standard deviation of random intercept varies in {0.000,0.356,0.712}.

JM.sas        -- Fit a joint model with a random intercept.

SJM.sas       -- Fit a simplified joint model without random effects.

estimates.sas -- Compute the average parameter estimates based on the joint models. 

Note: A linear time trajectory and constant baseline hazard function is used for both data generation and model fitting.


--------------------------------------------- Folder Overfitting-AppendixF ----------------------------------------

Description: This folder contains data generation and model fitting programs for results in Appendix F of the Supplementary Materials. 

dataset.R          -- Simulate data where patient-level heterogeneity is based on a random intercept. A linear time trajectory and constant baseline hazard function is assumed.

JM-none.sas        -- Fit a joint model with a random intercept. A linear trajectory and a constant baseline hazard function is used. 

JM-BL.sas          -- Fit a joint model with a random intercept. A linear trajectory and a piecewise constant baseline hazard function is used. 

JM-trajectory.sas  -- Fit a joint model with a random intercept. A piecewise linear trajectory and a constant baseline hazard function is used. 

JM-both.sas        -- Fit a joint model with a random intercept. A piecewise linear trajectory and a piecewise constant baseline hazard function is used. 

SJM-none.sas       -- Fit a simplified joint model without random effects. A linear trajectory and a constant baseline hazard function is used. 

SJM-BL.sas         -- Fit a simplified joint model without random effects. A linear trajectory and a piecewise constant baseline hazard function is used. 

SJM-trajectory.sas -- Fit a simplified joint model without random effects. A piecewise linear trajectory and a constant baseline hazard function is used. 

SJM-both.sas       -- Fit a simplified joint model without random effects. A piecewise linear trajectory and a piecewise constant baseline hazard function is used. 

