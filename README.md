# ModelEvolvR

## Overview
This package is a basic implementation of feature selection through an abstraction of the process of evolution. It can be applied to any model that accepts a formula as input. The user must however define their own fitness function that defines how well a model performs relative to other models. See the vignette [vignette](https://htmlpreview.github.io/?https://github.com/HunterGleason/ModelEvolvR/blob/master/vignettes/introduction.html) for a quick tutorial showing how to implement this package in R. 

## Installation
**Using dev tools:**
```
remotes::install_github("HunterGleason/ModelEvolvR",build_vignettes=TRUE")
library(ModelEvolvR)
```




