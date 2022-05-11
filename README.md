# deltaccd

[![check-deploy](https://github.com/hugheylab/deltaccd/workflows/check-deploy/badge.svg)](https://github.com/hugheylab/deltaccd/actions)
[![codecov](https://codecov.io/gh/hugheylab/deltaccd/branch/master/graph/badge.svg)](https://codecov.io/gh/hugheylab/deltaccd)
[![Netlify Status](https://api.netlify.com/api/v1/badges/ddd35b6b-7210-442f-83e7-7115b23d9585/deploy-status)](https://app.netlify.com/sites/jovial-lovelace-7e335a/deploys)
[![CRAN Status](https://www.r-pkg.org/badges/version/deltaccd)](https://cran.r-project.org/package=deltaccd)
[![drat version](https://raw.githubusercontent.com/hugheylab/drat/gh-pages/badges/deltaccd_drat_badge.svg)](https://github.com/hugheylab/drat/tree/gh-pages/src/contrib)

`deltaccd` infers the progression of circadian rhythms using gene co-expression. `deltaccd` works even if the samples are not labeled with time of day and do not cover the entire circadian cycle.

For details about the method and to see how we used it to analyze circadian gene expression in human cancer, check out [Shilts et al. (PeerJ)](https://doi.org/10.7717/peerj.4327) and the [accompanying results](https://doi.org/10.6084/m9.figshare.4906745).

## Installation

### Option 1: CRAN

```r
install.packages('deltaccd')
```

### Option 2: Hughey Lab Drat Repository

1. Install [`BiocManager`](https://cran.r-project.org/package=BiocManager).

    ```r
    if (!requireNamespace('BiocManager', quietly = TRUE))
      install.packages('BiocManager')
    ```

1. If you use RStudio, go to Tools → Global Options... → Packages → Add... (under Secondary repositories), then enter:

    - Name: hugheylab
    - Url: https://hugheylab.github.io/drat/

    You only have to do this once. Then you can install or update the package by entering:

    ```r
    BiocManager::install('deltaccd')
    ```

    Alternatively, you can install or update the package by entering:

    ```r
    BiocManager::install('deltaccd', site_repository = 'https://hugheylab.github.io/drat/')
    ```

## Usage

See the example in the documentation for `calcDeltaCCD()`, as well as the full [reference documentation](https://deltaccd.hugheylab.org/reference/index.html).
