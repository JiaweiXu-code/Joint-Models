# Programs used in "Bayesian Design of Clinical Trials Using Joint Models for Longitudinal and Time-to-event Data"

ssd.R       -- R code for sample size determination by varying the ratio between sample size and event total to obtain specified number of events in a specified interval of time                  on average. This code should be run before each data generation process to determine the desired sample size.
		 
dataset.R   -- R code for dataset generation. Data structure or paramter settings can be different under different scenarios. 
               Thus, the correspoding dataset.R code will be found under different folders.

JM.sas      -- Any SAS program starts with "JM" fits a joint models with a random intercept, with default of a 4-component piecewise linear time trajectory and a 5-component                      piecewise constant baseline hazard function. Interpretation of different settings can be found in the README file under the correspoding folder.

SJM.sas     -- Any SAS program starts with "SJM" fits a simplified joint model without random effects, with default of a 4-component piecewise linear time trajectory and a 5-                    component piecewise constant baseline hazard function. Interpretation of different settings can be found in the README file under the correspoding folder.

PH.R        -- R code that fits a piecewise proportional hazard models. This code is for all cases where a piecewise PH model is fitted

logrank.sas -- SAS program that runs a log-rank tests. This program is for all cases where a log-rank test is applied.

power.sas   -- SAS program that computes the estimated Bayesian type I error rate and power curves after obtaining avarage hazard ratios based on all simulated datasets.
               This program is for all cases where Bayesian type I error rate or power should be computed based on joint models. 


