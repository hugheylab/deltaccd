# deltaccd

[![CircleCI](https://circleci.com/gh/hugheylab/deltaccd.svg?style=shield)](https://circleci.com/gh/hugheylab/deltaccd)

`deltaccd` infers the progression of the circadian clock using clock gene co-expression. `deltaccd` works even if the samples are not labeled with time of day and do not cover the entire circadian cycle.

For details about the method and to see how we used it to analyze circadian gene expression in human cancer, check out [Shilts et al. (PeerJ)](https://doi.org/10.7717/peerj.4327) and the [accompanying results](https://doi.org/10.6084/m9.figshare.4906745).

## Installation

If you use RStudio, go to Tools -> Global Options... -> Packages -> Add... (under Secondary repositories), then enter:

- Name: hugheylab
- Url: https://hugheylab.github.io/drat/

You only have to do this once. Then you can install or update the package by entering:

```R
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('deltaccd')
```

Alternatively, you can install or update the package by entering:

```R
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('deltaccd', site_repository = 'https://hugheylab.github.io/drat/')
```

There's also a [docker image](https://hub.docker.com/r/hugheylab/hugheyverse), which has all dependencies installed.

```bash
docker pull hugheylab/hugheyverse
```

## Usage

See the example in the documentation for `calcDeltaCCD()`, as well as the full [reference documentation](reference/index.html).
