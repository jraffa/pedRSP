# pedRSP
Pedigree Relatedness Summary Parameters (pedRSP) R package.
For more details, see the tech report:

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
