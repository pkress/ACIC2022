# ACIC2022

This repository contains the code and notes for my Mathematica ACIC 2022 submission (see competition details [here](https://acic2022.mathematica.org/)). 

## Method: TMLE + BART

I decided to use TMLE + BART estimation based on its excellent performance in a previous causal inference contest (See [Dorie et al, 2017](https://arxiv.org/pdf/1707.02641.pdf)). 

## TMLE Overview

My use of the TMLE is based on the textbook: van der Laan MJ, Rose S. _Targeted learning: causal inference for observational and experimental data._ Springer; New York: 2011. Using this method for the ATT is summarized in [this](https://stats.stackexchange.com/questions/520472/can-targeted-maximum-likelihood-estimation-find-the-average-treatment-effect-on) stack-exchange post. Additional helpful TMLE resources can be found [here](https://www.khstats.com/blog/tmle/tutorial/) (simple explanation for ATE) and [here](https://ehsanx.github.io/TMLEworkshop/tmle.html#step-5-estimate-epsilon) (a little more in depth for ATE).

TMLE is Targeted Maximum Likelihood Estimation, a doubly robust estimation method that allows for asymptotic inference on "black-box" machine learning outcome estimates. It is doubly robust because if either the treatment surface OR the outcome surface are correctly specified, than the TMLE estimator will produce valid inference for effect estimates. 

The estimator effectively works in two steps: 
- Estimation of the outcome and treatment probability
- Updating those estimates to maximize the bias-variance tradeoff for the target estimand (in our case, SATT)

We can then summarize the difference in the updated treated and untreated estimates for samples of interest with valid inference. 

### Estimation: BART

I use Bayesian Additive Regression Trees (BART) to estimate each individual-year's outcome (healthcare spending) based on available covariates. I also use BART to estimate each medical practice's probability of being treated. 

BART is a Bayesian MCMC procedure to estimate an outcome (discrete, binary, or continuous) from covariates. It is similar to random forest/boosting in that it uses ensambles of small, simple models to create a sophisticated aggreage model. However, as a Bayesian MCMC procedure, it estimateds the posterior distribution of the outcome for individual units based on independent draws from the MCMC (once stationary). The mean of these draws is the expected outcome for a given unit, conditional on that unit's covariates. 

### Update the estimates

We now update this initial fit towards the optimal bias–variance tradeoff for the parameter of interest (the ATT) instead of the overall density. This is where the TMLE comes into play. This is a process that iteratively updates the initial estimates based on the actual outcome/treatment estimates and clever covariates (similar to the propensity score for the outcome estimate, and similar to the difference between estimated outcomes and ATT for the treatment estimate). 

Once we've updated the estimates to target the ATT, we then estimate the effects for each subsample using the targeted ATT estimate, and use the SD of the efficient influence function to estimate the upper and lower bounds. 





