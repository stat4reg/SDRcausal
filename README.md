# SDRcausal
The package implements the methods described in Ghosh et al. (2020): Sufficient Dimension Reduction for Feasible and Robust Estimation of Average Causal Effect. It uses a semiparametric locally efficient dimension reduction approach to assess both the treatment assignment mechanism and the average responses in both treated and nontreated groups. It then integrates all results through imputation, inverse probability weighting and doubly robust augmentation estimators.

# References
Ghosh, T., Ma, Y. and de Luna, X. (2020). Sufficient Dimension Reduction for Feasible and Robust Estimation of Average Causal Effect. Statistica Sinica, Vol. 31, (2) : 821-842. DOI: 10.5705/ss.202018.0416. ArXiv version with Supplementary material: https://arxiv.org/abs/1811.01992  

Ghasempour, M. and de Luna, X. (2021). SDRcausal: an R package for causal inference based on sufficient dimension reduction. ArXiv: https://arxiv.org/abs/2105.02499

# Installing
To install and load this package in R from GitHub, run the following commands:
```
install.packages("devtools")
library(devtools) 
install_github("stat4reg/SDRcausal")
library(SDRcausal)
```
## Mac OS
Mac user may need to run the following command in a terminal window before installing the package:
```
xcode-select --install
```
OpenMP is not available on Mac OS and only one thread can be used. 
