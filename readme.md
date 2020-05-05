# deltaccd
`deltaccd` infers the progression of the circadian clock using clock gene co-expression. `deltaccd` works even if the samples are not labeled with time of day and do not cover the entire circadian cycle.

For details about the method and to see how we used it to analyze circadian gene expression in human cancer, check out [Shilts et al. (PeerJ)](https://doi.org/10.7717/peerj.4327) and the [accompanying results](https://doi.org/10.6084/m9.figshare.4906745).

## Installation
First add the hugheylab repository to your repos. There are multiple ways to do this.

If you use RStudio, go to Tools -> Global Options... -> Packages -> Add... (under Secondary repositories), then enter the following values.

- Name: hugheylab
- Url: https://hugheylab.github.io/drat/

You only have to do this once.

Alternatively, enter the following command each time you want to install or update the package.
```R
options(repos = c(getOption('repos'), 'https://hugheylab.github.io/drat/'))
```

Now you can install the package.
```R
install.packages('deltaccd', type = 'source')
```
You can update the package using `update.packages()`.

There's also a pre-built [docker image](https://hub.docker.com/r/hugheylab/hugheyverse), which has all dependencies installed.
```bash
docker pull hugheylab/hugheyverse
```

## Usage
See the example in the documentation:
```R
library('deltaccd')
?calcDeltaCCD
```
