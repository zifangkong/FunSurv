
## **FunSurv**

### Funcational data analysis for recurrent event data

***FunSurv*** provides a novel framework for the joint analysis of
recurrent event and longitudinal data using a functional regression
approach. This package applies functional data analysis techniques to
model the time-varying effects of longitudinal predictors on recurrent
event processes, addressing the complex dependence structure between the
two data types.

### Installation

You can install and load **FunSurv** from CRAN using

``` r
install.packages("FunSurv")
library(FunSurv)
```

You can install **FunSurv** from Github with:

``` r
## install.packages("devtools")
devtools::install_github("zifangkong/FunSurv", ref = "master")
```

### Citation

Cite ***FunSurv*** with `citation("FunSurv")`.

``` r
citation("FunSurv")
```

### Online documentation

Add a vignette

### References:

Yau, K. K. W., & McGilchrist, C. A. (1998). ML and REML estimation in
survival analysis with time dependent correlated frailty. *Statistics in
Medicine*, **17**(11), 1201-1213.

Yao, F., Müller, H. G., & Wang, J. L. (2005). Functional linear
regression analysis for longitudinal data. *The Annals of Statistics*,
**33**(6), 2873–2903.

Yao, Fang. Functional principal component analysis for longitudinal and
survival data. *Statistica Sinica*, (2007): 965-983.

Kong, D., Ibrahim, J. G., Lee, E., & Zhu, H. (2018). FLCRM: Functional
linear cox regression model. *Biometrics*, **74**(1), 109-117.
