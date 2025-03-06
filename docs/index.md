---
layout: full
homepage: true
disable_anchors: true
description: Multiple Outcome Mendelian Randomization
---
## METEOR Overview
![METEOR\_pipeline](METEOR.png)
METEOR is an R package for efficient statistical inference of multi-outcomes mendelian randomization analysis, METEOR utilizes a set of correlated SNPs, self-adaptively accounts for the sample structure of both exposure and outcomes, the uncertainty that these correlated SNPs may exhibit multiple pleiotropic effects. The term ‘self-adaptive’ represents that METEOR is able to automatically infer the sample structure and the probability that a SNP has specific pleiotropic effect from the data at hand. METEOR places the inference of the causal effects into a likelihood-framework and relies on a scalable sampling-based algorithm to obtain calibrated $p$-values, freely available at https://github.com/Liye222/METEOR. 

Installation
------------
You can install the released version of MAPLE from Github with the following code. This package is supported for Windows 10/11, and Linux. The package has been tested on the following systems:
* Windows 10, 11
* Linux: Ubuntu (22.04.4)

### Dependencies 
* R version >= 3.6.0
* R packages: R packages: Rcpp, RcppArmadillo, RcppDist, dplyr, magrittr, readr, parallel


### 1. Install `devtools` if necessary
```r
install.packages('devtools')
```

### 2. Install `METEOR`
```r
devtools::install_github('Liye222/METEOR')
```
### 3. Load package
```r
library(METEOR)
```

### Issues
All feedback, bug reports and suggestions are warmly welcomed! Please make sure to raise issues with a detailed and reproducible example and also please provide the output of your sessionInfo() in R! 


How to use `METEOR`
-------------------

The METEOR User Manual: [here](https://github.com/Liye222/METEOR/blob/main/docs/METEOR_user_manual.pdf).

