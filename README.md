# Algorithm

Multi-level network meta-regression with a survival outcome

# Maintainer

Ruofan Jia 1155165249@link.cuhk.edu.hk

# Introduction

R functions and codes for the manuscript Multi-level network meta-regression with a survival outcome and an application to NSCLC.

In biopharmaceutical studies, it is often of interest to compare the effects of two treatments (say, A and B), but data that contains a direct comparison is not available. However, there may exist studies that compare them against another competitor (say, C). In this case, an ``indirect'' comparison of treatments A versus B through C needs to be conducted. In this article, we consider the scenario where individual-level data is available for the comparison of A versus C, but only aggregated data, for example from a publication, is available for the comparison of B versus C. For such analysis, multi-level network meta regression (ML-NMR) is advantageously needed since it can combine evidence from multiple trials with either survival IPD or AgD and can compare the treatments of interest in any target population. Most of the existing ML-NMR studies have been focused on binary and continuous outcomes, while, relatively, research on censored survival outcomes remains limited with only one study modeling the marginal likelihood for aggregated data. Here, we aim to extend the ML-NMR for time-to-event outcomes. We consider multiple popular parametric survival models and develop Bayesian estimation approaches built on mean and median survival. 

# Usage


surv_mlnmr: function to get the estimates from Weibull and Gompertz survival model

surv_mlnmr_trans: function to get the estimates from log normal and log logistic survival model

# References

1. Bennett I, Gregory J, Smith S, Birnie R. Roche/MAIC. url: 10.5281/zenodo.6624152. 2022-06-08;version: 0.3.0.

2. Gsteiger S, Howlett N, Ashlee B. gemtcPlus: Provides a suite of extension functions for NMA using the ‘gemtc‘ package. 2022. R package version
1.0.0.

3. Wiecek W and Pafitis S. certara/survivalnma. 2019.	R package version 0.1.1.

4. Remiro-Azócar A, Heath A, Baio G. Methods for population adjustment with limited access to individual patient data: a review and simulation
study. Research synthesis methods. 2021;12(6):750–775.
