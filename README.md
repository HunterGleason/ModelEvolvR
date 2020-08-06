# ModelEvolvR

## Overview
This package is a basic implementation of feature selection through an abstraction of the process of evolution. It can be applied to any model that accepts a formula as input. The user must however define their own fitness function that defines how well a model performs relative to other models. Below is a quick tutorial showing how to implement this package in R. 

## Installation
**Using dev tools:**
```
devtools::install_github("HunterGleason/ModelEvolvR)
library(ModelEvolvR)
```

## An Example

**Create a data frame**
```
x1<-c(1:100)+sample(-10:10,100,replace = T)
x2<-c(101:200)+sample(-10:10,length(x1),replace = T)
x3<-c(301:400)+sample(-10:10,length(x1),replace = T)
x4 <- sample(1:400,length(x1))
x5 <- sample(1:400,length(x1))
x6 <- sample(1:400,length(x1))
x7 <- sample(1:400,length(x1))
x8 <- sample(1:400,length(x1))
x9 <- sample(1:400,length(x1))
x10 <- sample(1:400,length(x1))

x_preds<-cbind(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10)

#Define a known function to test feature selection
y<-0.5*x1+.25*x2+x3

DF<-as.data.frame(cbind(y,x_preds))

```
**Run GA feature selection**
```
x_names<-colnames(DF)[-1]
y_name<-colnames(DF)[1]

#Example of fitness function
fit_fun<-function(FORMULA,DF)
{

  mod<-lm(formula = FORMULA, data = DF)

  bic<-BIC(mod)
  
  return(bic)

}

max_terms<-ncol(DF)-1

pop_size<-100

percentile<-0.25

max_mute_rate<-.3

generations<-2000

maximize<-FALSE

evolved<-ModelEvolvR::evolve(y_name,x_names,DF,fit_fun,max_terms,pop_size,percentile,max_mute_rate,generations,maximize)
```
**Explore Results**
```
plot(evolved[[3]], type='l', xlab="Generation",ylab = "Pop. Mean Fitness (BIC)")

best_model<-evolved[[1]][[which.min(evolved[[2]])]]

summary(lm(best_model,DF))
```


