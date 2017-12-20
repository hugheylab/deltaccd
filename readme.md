# deltaccd
`deltaccd` is a package for inferring progression of the circadian clock using clock gene co-expression. `deltaccd` works even if the samples are not labeled with time of day and do not cover the entire circadian cycle.

For details about the method and to see how we used it to analyze circadian gene expression in human cancer, check out [Shilts et al. (bioRxiv)](https://dx.doi.org/10.1101/130765) and the [accompanying results](https://figshare.com/s/2eaf11e88642418f7e81).

## Install using drat
```R
install.packages('drat')
drat::addRepo('hugheylab')
install.packages('deltaccd', type='source')
```
You can update the package by calling `drat::addRepo('hugheylab')`, then `update.packages`.

## Install using devtools
```R
install.packages('devtools')
devtools::install_github('hugheylab/deltaccd')
```
You can then update the package using these same two lines.

## Install using docker
You can use a pre-built [docker image](https://hub.docker.com/r/hugheylab/hugheyverse), which has all dependencies already installed:
```
docker pull hugheylab/hugheyverse
```

## Getting started
See the example in the documentation:
```R
library('deltaccd')
?calcDeltaCCD
```
