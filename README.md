# Programs used in "Bayesian Design of Clinical Trials Using Joint Models for Longitudinal and Time-to-event Data"

PH.R         -- R code for fitting piecewise proportional hazard models. This code is for all cases where a piecewise PH model is fitted

logrank.sas  -- SAS program for fitting log-rank tests. This program is for all cases where a log-rank test is applied.

power.sas    -- SAS program for computing the estimated Bayesian type I error rate and power curves after obtaining avarage hazard ratios based on all simulated datasets.
                This program is for all cases where Bayesian type I error rate or power should be computed based on joint models. 

ssd.R        -- R code for sample size determination by varying the ratio between sample size and event total to obtain specified number of events in a specified interval of time                 on average. This code should be run before each data generation process to determine the desired sample size.

