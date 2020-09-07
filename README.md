# Programs used in "Bayesian Design of Clinical Trials Using Joint Models for Longitudinal and Time-to-event Data"

ssd.R         -- R code for sample size determination by varying the ratio between sample size and event total to obtain specified number of events in a specified interval of time                  on average. This code should be run before each data generation process to determine the desired sample size.
		 
dataset.R     -- R code for dataset generation. Data structure or paramter settings can be different under different scenarios. 
                 Thus, the correspoding dataset.R code will be found under different folders.

JM.sas        -- Any SAS program starts with "JM" fits a joint models with a random intercept, with default of a 4-component piecewise linear time trajectory and a 5-component                      piecewise constant baseline hazard function. Interpretation of different settings can be found in the README file under the correspoding folder.

SJM.sas       -- Any SAS program starts with "SJM" fits a simplified joint model without random effects, with default of a 4-component piecewise linear time trajectory and a 5-                    component piecewise constant baseline hazard function. Interpretation of different settings can be found in the README file under the correspoding folder.

PH.R          -- R code that fits a piecewise proportional hazard models. This code is for all cases where a piecewise PH model is fitted

logrank.sas   -- SAS program that runs a log-rank tests. This program is for all cases where a log-rank test is applied.

power.sas     -- SAS program that computes the estimated Bayesian type I error rate and power curves after obtaining avarage hazard ratios based on all simulated datasets.
                 This program is for all cases where Bayesian type I error rate or power should be computed based on joint models. 

estimates.sas -- SAS program that computes the average parameter estimates based on the joint models. This program is for all cases where joint models with/without are fitted.



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



Outputs from "dataset.R" of simulated data --

ID         -- patient ID

measure    -- longitudinal measures

time       -- measurement time (years) since enrollment

trt        -- treatment indicator (1 = treated, 0 = control)

toltime    -- measurement time (years) since trial started

nodegrp    -- lymph node group, binary baseline covariate (1 = more severe group, 0 = less severe group)

sim        -- ID for simulated data

t          -- follow-up time (years) since enrollment

sr         -- follow-up time (years) since trial started

r0         -- enrollment time (years)

censorship -- indicator for censoring (1 = event, 0 = censored)
