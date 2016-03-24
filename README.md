# pedRSP
Pedigree Relatedness Summary Parameters (pedRSP) R package.

This package is for a paper accepted in Statistics in Biosciences:

```{ref}
Raffa J.D. and Thompson E.A. (2016) Power and Effective Study Size in Heritability Studies.  Statistics in Biosciences.
```
[Available Here:](http://link.springer.com/article/10.1007/s12561-016-9143-2)(http://link.springer.com/article/10.1007/s12561-016-9143-2)

An earlier version is available as a Tech Report:

```{ref}
University of Washington, Department of Statistics, Technical Report no. 630 (389 KB)
    Power and Effective Study Size Based on Approximations to the Expected Likelihood Ratio Test in Heritability Studies 
    Jesse Raffa and Elizabeth Thompson 
    October 2014
```

[PDF:](http://www.stat.washington.edu/research/reports/2014/tr630.pdf) [(http://www.stat.washington.edu/research/reports/2014/tr630.pdf)](http://www.stat.washington.edu/research/reports/2014/tr630.pdf)

The development of this software was supported in part by NIH grant R37-GM046255.

# To Install

```{r}
install.packages("devtools");
library(devtools);
install_github("jraffa/pedRSP");
library(pedRSP);
```

# To Get Started

```{r}
?computeRSPs
```
