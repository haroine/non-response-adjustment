# non-response-adjustment
Functions developed at INSEE to produce non-response adjusted estimators using the Homogeneous Response Groups method (Haziza and Beaumont, 2007).

## Short example

A short example of non-response adjustment using table _dataPop_ from package Icarus.

First, we load the package Icarus and the functions for non-response adjustment:

```
install.packages("icarus")
library(icarus)
source("fonctions_grh.R")
```

Then we fit a simple logistic model to predict non-response:

```
samplePop <- dataPop[dataPop$weight > 0,] ## Restriction to units selected in the sample
modelTest <- glm(responding ~ Y1 + Y11 + Y12 + X , 
                 data = samplePop, family = "binomial")
summary(modelTest)
```

Then we can compute the Homogeneous Response Groups (HRG) using this model and obtain the non-response adjusted weights:

```
NRAweights1 <- ajouterPoidsGRH(samplePop, modelTest,
                              colPoids="weight", colRepondant = "responding")

```
In the dataframe _NRAweights1_ we will thus find the column of NRA weights (_POIDS\_CNR_ and other useful information, such as the ratio between adjusted and initial weights, etc.).

Another possibility is to directly use a column of computed response probabilities (for example coming from a different learning model than the logistic regression). Here, we have the real non-response probabilities that were used to generate the sample, so we might as well use it:

```
NRAweights2 <- ajouterPoidsGRH(samplePop, modelGRH=NULL, pHat=1-(samplePop$simul_nr),
                               colPoids="weight", colRepondant = "responding")
```

## Notes

Please note that some function names are still in French and that unit tests are not implemented yet. These functions come with no guarantee, use at your own risk.
