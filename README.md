# Software for paper "Bayesian Design of Clinical Trials Using Joint Models for Longitudinal and Time-to-event Data" by Xu et al.

-------------------------------------------------------------------------------------------------------------------------------------

All programs are setup to be executed on a Linux computing cluster using R (3.6.0) or SAS (9.4). The paths referenced in all programs will need to be updated for the code to work. Once all paths are updated, one can use the SLURM scheduler shell scripts to submit jobs on a SLURM-based computing cluster. 


-------------------------- RUN ORDER ------------------------------

[1] ssd.sh -- Determine the desired sample size such that a specified number of events is obtained in a specified interval of time on average. This process is performed based on 1,000 simulations and calls the R program "ssd.R" which performs a single simulation. The "ssd.R" code requires the following inputs:

   (1)  integer - seed - seed for a single simulation
   (2)  integer - v    - number of events
   (3)  double  - k    - ratio between number of enrolled patients and event total
   (4)  integer - nlinear - number of intervals for peicewise linear trajectory
   (5)  integer - measures
   (6)  integer - npieces 
   (7)  double  - eta
   (8)  double  - maxt 
   (9)  double  - censor 
   (10) double  - censorp = 0.05, 
   (11) double  - interval = 0.001, 
   (12) vector (dbl) - pknots = c(0.25,0.75,1.25), 
   (13) vector (dbl) - knots = c(1.91, 2.43, 3.00, 3.80), 
   (14) double      - sigma = 0.65642581, 
   (15) vector (dbl) - p = c(0.5,0.5), 
   (16) double - mSigma = 0.71202084, 
   (17) vector (dbl) - gamma = c(0.26634693,0,-0.03393003),
   (18) vector (dbl) - plinear = c(-0.31869067,
   (19) vector (dbl) - pinter = c(0.0,
   (20) vector (dbl) - llambda = c(-3.61388568,
   (21) double - beta = -0.15,
   (22) vector (dbl) - alpha = c(-0.2,
   (23) vector (dbl) - tpoints 




-------------- Folder MainResults-Figures2&3 -----------------

Description: Contains all programs and scripts to generate the results in Section 3.1.


		 
dataset.R     -- R code for dataset generation. Data structure or paramter settings can be different under different scenarios. 
                 Thus, the correspoding dataset.R code will be found under different folders.

JM.sas        -- Any SAS program starts with "JM" fits a joint models with a random intercept, with default of a 4-component piecewise linear time trajectory and a 5-component                      piecewise constant baseline hazard function. Interpretation of different settings can be found in the README file under the correspoding folder.

SJM.sas       -- Any SAS program starts with "SJM" fits a simplified joint model without random effects, with default of a 4-component piecewise linear time trajectory and a 5-                    component piecewise constant baseline hazard function. Interpretation of different settings can be found in the README file under the correspoding folder.

PH.R          -- R code that fits a piecewise proportional hazard models. This code is for all cases where a piecewise PH model is fitted

logrank.sas   -- SAS program that runs a log-rank tests. This program is for all cases where a log-rank test is applied.

power.sas     -- SAS program that computes the estimated Bayesian type I error rate and power curves after obtaining avarage hazard ratios based on all simulated datasets.
                 This program is for all cases where Bayesian type I error rate or power should be computed based on joint models. 

estimates.sas -- SAS program that computes the average parameter estimates based on the joint models. This program is for all cases where joint models with/without are fitted.


-------------------------------------------------------------------------------------------------
Inputs in "dataset.R" for data simulation --
 
v        -- number of events

k        -- ratio between number of patients and event total

nlinear  -- number of intervals for peicewise linear trajectory

measures -- number of scheduled longitudinal measurements

npieces  -- number of intervals for piecewise constant baseline hazard function

eta      -- period of time (years) for patient enrollment

maxt     -- maximum follow-up time (years) under ideal cases (i.e., no dropout or censoring)

censor   -- period of time (years) for dropout

censorp  -- dropout probability

interval -- length of interval for discretization to simulate survival time

pknots   -- knots placement for piecewise linear trajectory

knots    -- knots placement for piecewise constant baseline hazard function

sigma    -- standard deviation of measurement error

p        -- allocation probabilities for treatment and baseline covariates 

mSigma   -- standard deviation or covariance matrix of random effects

gamma    -- regression coefficients of intercept, treatment indicator and baseline covariates in longitudinal process 

plinear  -- coefficients for time covariates in longitudinal process

pinter   -- coefficients for interactions between treatment and time covariates in longitudinal process

llambda  -- piecewise constant baseline hazards in log scale

beta     -- association parameter

alpha    -- direct effects of treatment and baseline convariate on time-to-event endpoint

tpoints  -- scheduled time points (years) at which longitudinal outcomes are measured

mk       -- ratio between standard deviations of random slope and random intercept when used in random effect misspecification



--------------------------------------------------------------------------------------------------
Outputs from "dataset.R" of simulated data --

ID         -- patient ID

measure    -- longitudinal measures

time       -- measurement time (years) since enrollment

trt        -- treatment indicator (1 = treated, 0 = control)

toltime    -- measurement time (years) since trial started

nodegrp    -- lymph node group, binary baseline covariate (1 = more severe group, 0 = less severe group)

sim        -- ID for simulated data

s          -- follow-up time (years) since enrollment

sr         -- follow-up time (years) since trial started

r0         -- enrollment time (years)

censorship -- indicator for censoring (1 = event, 0 = censored)
