# HDMTS
High‐dimensional mediation testing with stepwise regression

Author: Wan-Yu Lin, Institute of Health Data Analytics and Statistics, National Taiwan University

If you use this program, please cite: High‐dimensional mediation testing with stepwise regression: Exploring epigenome-wide mediation effects. Thank you.

Suppose you put the three files (HDMTS.R, toydata.csv, toydata2.csv) under "H:/". The following is R code:

setwd("H:/")

aa <- read.csv("toydata.csv")

M <- aa[,5:ncol(aa)]    # mediators

X <- aa[,2]    # exposure

Y <- aa[,1]    # outcome

COV <- aa[,3:4]    # covariates

source("HDMTS.R")

HDMTS(X, M, Y, COV = COV, FDRcut = 0.05, VIFcut = 5.0, nboot = 1000, bootCI = 0.95)

aa <- read.csv("toydata2.csv")

M <- aa[,5:ncol(aa)]    # mediators

X <- aa[,2]    # exposure

Y <- aa[,1]    # outcome

COV <- NULL    # No covariates

HDMTS(X, M, Y, COV = NULL, FDRcut = 0.05, VIFcut = 5.0, nboot = 1000, bootCI = 0.95)

