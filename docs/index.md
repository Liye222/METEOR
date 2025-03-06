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
You can install the released version of GIFT from Github with the following code, for more installation details or solutions that might solve related issues (specifically MacOS system) see the [link](https://yuanzhongshang.github.io/GIFT/documentation/02_installation.html).

### Dependencies 
* R version >= 4.0.0.
* R packages: Rcpp, RcppArmadillo, parallel


### 1. Install `devtools` if necessary
```r
install.packages('devtools')
```

### 2. Install `GIFT`
```r
devtools::install_github('yuanzhongshang/GIFT')
```
### 3. Load package
```r
library(GIFT)
```

This package is supported for Windows 10, MAC and Linux. The package has been tested on the following systems:
- Windows 10
- MAC: OSX (10.14.1)
- Linux: Ubuntu (16.04.6)

### Issues
All feedback, bug reports and suggestions are warmly welcomed! Please make sure to raise issues with a detailed and reproducible example and also please provide the output of your sessionInfo() in R! 

How to cite `GIFT`
-------------------
Liu, L., Yan, R., Guo, P. et al. Conditional transcriptome-wide association study for fine-mapping candidate causal genes. Nat Genet 56, 348–356 (2024).
[https://doi.org/10.1038/s41588-023-01645-y](https://doi.org/10.1038/s41588-023-01645-y)

How to use `GIFT`
-------------------
Example Analysis with GIFT: [here](https://yuanzhongshang.github.io/GIFT/documentation/04_GIFT_Example.html).

The GIFT Manual: [here](https://github.com/yuanzhongshang/GIFT/blob/main/docs/GIFT%20manual.pdf).

The genome-wide eQTL summary statistics from GEUVADIS data
-------------------
This data is availible in the dropbox: [here](https://www.dropbox.com/scl/fo/4nqcmkblerspfmva5stwf/ANHZU_kX2AlveEEbx9DKbZU?rlkey=qjcxprlk83t7pw8ka2ne2v4w9&dl=0).

The correlation matrix among gene expressions for each chromosome from GEUVADIS data is also availible in the dropbox: [here](https://www.dropbox.com/home/GEUVADIS/correlation%20matrix).
