calcCorrSimple = function(df) {
	geneNames = setdiff(colnames(df), 'group')
	df1 = data.frame(cor(as.matrix(df[,geneNames]), method=method), gene1 = geneNames,
					 stringsAsFactors=FALSE, check.names=FALSE) %>%
		tidyr::gather(-gene1, key=gene2, value=rho) %>%
		dplyr::filter(gene1!=gene2)}


calcCorr = function(ematNow, groupVec) {
	df = data.frame(t(ematNow), group = groupVec, stringsAsFactors=FALSE, check.names=FALSE) %>%
		dplyr::group_by(group) %>%
		dplyr::do(calcCorrSimple(.)) %>%
		dplyr::ungroup() %>%
		dplyr::mutate(gene1 = factor(gene1, rownames(ematNow)),
						  gene2 = factor(gene2, rev(rownames(ematNow))))}


calcColorLimits = function(vals, vLow=-1, vMid=0, vHigh=1, cLow='#e66101', cMid='#f7f7f7', cHigh='#5e3c99') {
	valRange = seq(0, 1, length.out=201)
	colorScale = scales::div_gradient_pal(low=cLow, mid=cMid, high=cHigh)(valRange)

	minVal = (min(vals, na.rm=TRUE) - vLow) / (vHigh - vLow)
	idxLow = which.min(abs(minVal - valRange))

	maxVal = (max(vals, na.rm=TRUE) - vLow) / (vHigh - vLow)
	idxHigh = which.min(abs(maxVal - valRange))
	return(colorScale[c(idxLow, idxHigh)])}


plotHeatmapSimple = function(ggObj, cLims) {
	p = ggObj +
		ggplot2::geom_tile(ggplot2::aes(x=gene1, y=gene2, fill=rho)) +
		ggplot2::labs(x='Gene', y='Gene') +
		ggplot2::scale_fill_gradient2(low=cLims[1], mid='#f7f7f7', high=cLims[2], na.value='grey80') +
		ggplot2::theme_light() +
		ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90, vjust=0.5, hjust=1),
							strip.text = ggplot2::element_text(color='black'),
							axis.text = ggplot2::element_text(color='black'),
							legend.margin = ggplot2::margin(t=0, r=0, b=0, l=0, unit='cm'))}


#' Visualize gene co-expression.
#'
#' `plotHeatmap` creates heatmaps of the co-expression (Spearman correlation)
#' between pairs of selected genes in a dataset.
#'
#' @param geneNames Vector indicating the subset of genes in the rownames of `emat` for
#' which to calculate the correlations in expression.
#' @param emat Matrix of expression values, where each row corresponds to a gene and
#' each column corresponds to a sample. The elements of `geneNames` should be present
#' in the rownames of `emat`.
#' @param groupVec Optional vector indicating the group to which group each sample belongs.
#' If not provided, the function assumes all samples belong to the same group.
#'
#' @return A `ggplot` object, which can be saved using `\link[ggplot2]{ggsave}()`.
#' Heatmap colors will be directly comparable to any heatmaps created by this function
#' or by `\link{plotRefHeatmap}()`.
#'
#' @examples
#' \dontrun{
#' library('deltaccd')
#' library('doParallel')
#' library('doRNG')
#'
#' registerDoParallel(cores=2)
#' set.seed(35813)
#'
#' refCor = getRefCor()
#' ccdResult = calcCCD(refCor, GSE19188$emat, GSE19188$groupVec, dopar=TRUE)
#' deltaCcdResult = calcDeltaCCD(refCor, GSE19188$emat, GSE19188$groupVec, 'non-tumor', dopar=TRUE)
#'
#' pRef = plotRefHeatmap(refCor)
#' pTest = plotHeatmap(rownames(refCor), GSE19188$emat, GSE19188$groupVec)
#' }
#'
#' @seealso `\link{calcCCD}`, `\link{calcDeltaCCD}`, `\link{plotRefHeatmap}`
#'
#' @export
plotHeatmap = function(geneNames, emat, groupVec=NULL) {
	method = 'spearman'

	ematNow = emat[geneNames,]
	if (nrow(ematNow)<2) {
		stop('Fewer than two genes in the supplied vector are in the expression matrix.')
	} else if (nrow(ematNow) < length(geneNames)) {
		warning(sprintf('%d gene(s) in the supplied vector is/are not in the expression matrix.',
							 length(geneNames)-nrow(ematNow)))}

	if (is.null(groupVec)) {
		groupVec = rep('all', ncol(emat))
	} else if (length(groupVec)!=ncol(emat)) {
		stop('Length of groupVec does not match the number of columns in emat.')
	} else if (min(table(groupVec)) < 3) {
		stop('Each unique group in groupVec must have at least three samples.')}

	df = calcCorr(ematNow, groupVec)
	cLims = calcColorLimits(df$rho)
	p = plotHeatmapSimple(ggplot2::ggplot(df) + ggplot2::facet_wrap(~ group), cLims)}


#' Visualize the reference pattern of gene co-expression.
#'
#' `plotRefHeatmap()` creates a heatmap of the reference correlation matrix
#' for gene co-expression.
#'
#' @param refCor Correlation matrix, such as comes from `\link{getRefCor}()`.
#'
#' @return A `ggplot` object, which can be saved using `\link[ggplot2]{ggsave}()`.
#' Heatmap colors will be directly comparable to any heatmaps created by this function
#' or by `\link{plotHeatmap}()`.
#'
#' @examples
#' \dontrun{
#' library('deltaccd')
#' library('doParallel')
#' library('doRNG')
#'
#' registerDoParallel(cores=2)
#' set.seed(35813)
#'
#' refCor = getRefCor()
#' ccdResult = calcCCD(refCor, GSE19188$emat, GSE19188$groupVec, dopar=TRUE)
#' deltaCcdResult = calcDeltaCCD(refCor, GSE19188$emat, GSE19188$groupVec, 'non-tumor', dopar=TRUE)
#'
#' pRef = plotRefHeatmap(refCor)
#' pTest = plotHeatmap(rownames(refCor), GSE19188$emat, GSE19188$groupVec)
#' }
#'
#' @seealso `\link{getRefCor}`, `\link{plotHeatmap}`
#'
#' @export
plotRefHeatmap = function(refCor) {
	if (any(rownames(refCor)!=colnames(refCor)) || !isSymmetric(refCor)) {
		stop('refCor must be a correlation matrix, with identical rownames and colnames.')}

	df = data.frame(refCor, gene1 = geneNames, stringsAsFactors=FALSE, check.names=FALSE) %>%
		tidyr::gather(-gene1, key=gene2, value=rho) %>%
		dplyr::filter(gene1!=gene2) %>%
		dplyr::mutate(gene1 = factor(gene1, rownames(ematNow)),
						  gene2 = factor(gene2, rev(rownames(ematNow))))

	cLims = calcColorLimits(df$rho)
	p = plotHeatmapSimple(ggplot2::ggplot(df), cLims)}
