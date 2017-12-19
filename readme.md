# deltaccd

`deltaccd` is a package for inferring progression of the circadian clock using clock gene co-expression. `deltaccd` works even if the samples are not labeled with time of day and do not cover the entire circadian cycle.

For details about the method and to see how we used it to analyze circadian gene expression in human cancer, check out [Shilts et al. (bioRxiv)](https://dx.doi.org/10.1101/130765) and the [accompanying results](https://figshare.com/s/2eaf11e88642418f7e81).

## Installation using drat
This is the preferred method.
```R
install.packages('drat')
drat::addRepo('hugheylab')
install.packages('deltaccd')
```

You can then update the package using `update.packages`.

## Installation using devtools
For the initial installation:
```R
install.packages('devtools')
devtools::install_github('hugheylab/deltaccd')
```

To update the package:
```R
devtools::install_github('hugheylab/deltaccd')
```

## Installation using docker
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
