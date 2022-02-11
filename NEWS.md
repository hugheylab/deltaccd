# deltaccd 1.0.2
* Updated examples.

# deltaccd 1.0.1
* Downsampled genes in GSE19188 to reduce package size.
* Fixed recent bug in deltaCCD calculation.

# deltaccd 0.2.7
* `calcCCD` and `calcDeltaCCD` now produce .csv lists of genes with zero group-wise variance or which are missing from the reference.

# deltaccd 0.2.6
* Fixed bug in calcCCD.

# deltaccd 0.2.5
* Replaced `globalVariables` call to get around `R CMD check` warning with null variable assignment and using the rlang .data variable for ggplot2.

# deltaccd 0.2.4
* Added tests for all functions.

# deltaccd 0.2.3
* Modified `calcCCD` and `calcDeltaCCD` to throw errors for genes with zero group-wise variance.
* `calcCCD` and `calcDeltaCCD` now require all genes from the reference to be present in the expression matrix passed to `emat`.

# deltaccd 0.2.2
* Updated gene symbol.

# deltaccd 0.2.1
* Fixed bugs in `plotRefHeatmap()`.
* Ordered genes in human blood reference to make heatmap more interpretable.

# deltaccd 0.2.0
* Modified `calcCorr()`, `calcCCD()`, `calcDeltaCCD()` to use data.table.
* Added optional `scale` argument to `calcCCD()` and `calcDeltaCCD()` to scale CCD by number of gene pairs.
* Added reference correlations for human blood.

# deltaccd 0.1.2
* Fixed obvious bug in `plotHeatmap()`.

# deltaccd 0.1
* Added `pkgdown` site.
* Updated documentation.
