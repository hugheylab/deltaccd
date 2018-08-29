# deltaccd
`deltaccd` is a package for inferring progression of the circadian clock using clock gene co-expression. `deltaccd` works even if the samples are not labeled with time of day and do not cover the entire circadian cycle.

For details about the method and to see how we used it to analyze circadian gene expression in human cancer, check out [Shilts et al. (PeerJ)](https://doi.org/10.7717/peerj.4327) and the [accompanying results](https://doi.org/10.6084/m9.figshare.4906745).

## Installation
First install drat.
```R
install.packages('drat')
```

Then add the following line to your `.Rprofile` file (located at "~/.Rprofile"), which gets run every time R starts. See [here](https://csgillespie.github.io/efficientR/3-3-r-startup.html#r-startup) for details.
```R
drat::addRepo('hugheylab')
```

Now you can install the package.
```R
install.packages('deltaccd', type = 'source')
```
You can update the package using `update.packages()`.

## Docker
You can also use a pre-built [docker image](https://hub.docker.com/r/hugheylab/hugheyverse), which has all dependencies installed.
```bash
docker pull hugheylab/hugheyverse
```

## Getting started
See the example in the documentation:
```R
library('deltaccd')
?calcDeltaCCD
```
