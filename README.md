# pollination-networks-implicit-feedback
Latent factor models to predict plant-pollinator interactions using positive and negative implicit feedback 

**Introduction**

This is the R implementation of the paper [Predicting Links in Plant-Pollinator Interaction Networks Using Latent Factor Models with Implicit Feedback](https://www.aaai.org/ocs/index.php/AAAI/AAAI18/paper/view/17131/15762).


**Requirements**

R version: 3.4.2

R Packages
- reshape (for building a matrix)
- proxy (for similarity measurements)

## How to run
```
source("main.R")
```

For all years (2011-2015) to be evaluated
```
EvaluateAllYears() 
```

For a year to be evaluated
```
EvaluateOneYear(11)
```

The final results are saved in the "Results" folder.
