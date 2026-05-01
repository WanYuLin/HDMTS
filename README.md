# HDMTS
High‐dimensional mediation testing with stepwise regression

Author: Wan-Yu Lin, Institute of Health Data Analytics and Statistics, National Taiwan University

If you use this program, please cite: A Systematic Evaluation of High-Dimensional Mediation Methods in Epigenome-Wide Studies: Implications for Method Selection. Thank you.

Suppose you put the three files (HDMTS.R, toydata.csv, toydata2.csv) under "H:/". The following is R code:

setwd("H:/")

aa <- read.csv("toydata.csv")

M <- aa[,5:ncol(aa)]    # mediators

X <- aa[,2]    # exposure

Y <- aa[,1]    # outcome

COV <- aa[,3:4]    # covariates

FDRcut <- 0.05    # the desired false discovery rate of mediator identification

VIFcut <- 5.0    # the largest acceptable variance inflation factor (VIF) of the mediator-outcome model

nboot <- 1000    # the number of bootstrap samples for calculating bootstrap confidence intervals

bootCI <- 0.95    # the confidence level of the bootstrap confidence intervals 

source("HDMTS.R")

HDMTS(X, M, Y, COV = COV, FDRcut = 0.05, VIFcut = 5.0, nboot = 1000, bootCI = 0.95)

aa <- read.csv("toydata2.csv")

M <- aa[,5:ncol(aa)]    # mediators

X <- aa[,2]    # exposure

Y <- aa[,1]    # outcome

COV <- NULL    # No covariates

HDMTS(X, M, Y, COV = NULL, FDRcut = 0.05, VIFcut = 5.0, nboot = 1000, bootCI = 0.95)

# Outputs:

alpha: the association between the exposure and the mediator

p_alpha: the two-sided p-value for testing H0: alpha=0 vs. H1: alpha!=0

beta: the association between the mediator and the outcome

p_beta: the two-sided p-value for testing H0: beta=0 vs. H1: beta!=0

Pmax: the maximum of p_alpha and p_beta

Mediation.Effect: alpha*beta

Mediation.Prop: mediation proportion = mediation effects/total effects = mediation effects/(direct effects + the sum of mediation effects from all mediators)

p_Bonferroni: Bonferroni correction on Pmax

p_FDR: Benjamini-Hochberg correction on Pmax

CI_lower: the lower bound of the bootstrap confidence intervals 

CI_upper: the upper bound of the bootstrap confidence intervals 

p_Boot: Bootstrap p-value

p_Boot_FDR: Benjamini-Hochberg correction on p_Boot






