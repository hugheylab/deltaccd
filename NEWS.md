# deltaccd 0.2.1
* Modified `calcCCD` and `calcDeltaCCD` to throw errors for genes with zero group-wise variance.
* `calcCCD` and `calcDeltaCCD` now require all genes from the reference to be present in the expression matrix passed to `emat`.

# deltaccd 0.2.0
* Modified `calcCorr()`, `calcCCD()`, `calcDeltaCCD()` to use data.table.
* Added optional `scale` argument to `calcCCD()` and `calcDeltaCCD()` to scale CCD by number of gene pairs.
* Added reference correlations for human blood.

# deltaccd 0.1.2
* Fixed obvious bug in `plotHeatmap()`.

# deltaccd 0.1
* Added `pkgdown` site.
* Updated documentation.
