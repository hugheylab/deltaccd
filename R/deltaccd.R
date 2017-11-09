#' @importFrom foreach foreach
#' @importFrom foreach "%do%"
#' @importFrom foreach "%dopar%"
#' @importFrom doRNG "%dorng%"

globalVariables('ii')


#' Retrieve the reference correlation matrix for clock gene co-expression.
#'
#' The reference matrix is based on a fixed-effects meta-analysis of eight circadian
#' transcriptome datasets from mice, as described in \url{https://dx.doi.org/10.1101/130765}.
#'
#' @param species Currently either 'human' or 'mouse'. Only affects the row and column names
#' of the correlation matrix, not the actual values.
#' @param useEntrezGeneId If FALSE, row and column names of correlation matrix
#' will correspond to gene symbols (e.g., PER2).
#'
#' @return A matrix of Spearman correlation values.
#'
#' @examples
#' \dontrun{
#' library('deltaccd')
#' library('doParallel')
#' library('doRNG')
#'
#' registerDoParallel(cores=2)
#' set.seed(35811)
#'
#' refCor = getRefCor()
#' ccdResult = calcCCD(refCor, GSE19188$emat, GSE19188$groupVec, dopar=TRUE)
#' deltaCcdResult = calcDeltaCCD(refCor, GSE19188$emat, GSE19188$groupVec, 'non-tumor', dopar=TRUE)
#' }
#'
#' @seealso \code{\link{GSE19188}}, \code{\link{calcCCD}}, \code{\link{calcDeltaCCD}}
#'
#' @export
getRefCor = function(species='human', useEntrezGeneId=TRUE) {
	if (species=='human') {
		if (useEntrezGeneId) {
			colname = 'entrez_hs'
		} else {
			colname = 'symbol_hs'}
	} else if (species=='mouse') {
		if (useEntrezGeneId) {
			colname = 'entrez_mm'
		} else {
			colname = 'symbol_mm'}}
	ref = refCorMouseEntrez
	rownames(ref) = clockGenes[[colname]]
	colnames(ref) = rownames(ref)
	return(ref)}


calcDist = function(r1, r2) sqrt(sum((r1-r2)^2, na.rm=TRUE))


calcCCDSimple = function(ref, emat, method='spearman') {
	corVecRef = ref[upper.tri(ref)]
	corMatTest = cor(t(emat), method=method)
	corVecTest = corMatTest[upper.tri(corMatTest)]
	return(calcDist(corVecRef, corVecTest))}


#' Calculate clock correlation distance (CCD).
#'
#' \code{calcCCD} quantifies the similarity of gene co-expression between a reference
#' and a test dataset. Statistical significance is calculated using permutation of the genes.
#'
#' @param refCor Correlation matrix to be used as the reference, such as comes
#' from \code{\link{getRefCor}}. Should contain Spearman correlation values.
#' @param emat Matrix of expression values, where each row
#' corresponds to a gene and each column corresponds to a sample. The rownames and colnames
#' of \code{refCor} should be present in the rownames of \code{emat}. For the p-value
#' calculation, it is important that \code{emat} include all measured genes, not just those
#' in \code{refCor}.
#' @param groupVec Optional vector indicating the group to which group each sample belongs.
#' If not provided, the function assumes all samples belong to the same group.
#' @param refEmat Optional expression matrix for calculating co-expression for the reference,
#' with the same organization as \code{emat}. Only used if \code{refCor} is not provided.
#' @param nPerm Number of permutations for assessing statistical significance.
#' @param geneNames Optional vector indicating a subset of genes in \code{refCor}, \code{emat},
#' and/or \code{refEmat} to use for calculating the CCD.
#' @param dopar Logical indicating whether to process features in parallel. Prior to calling
#' \code{calcCCD}, use \code{\link[doParallel]{registerDoParallel}} to register the parallel
#' backend, followed by \code{\link{set.seed}} to make the p-values reproducible.
#'
#' @return A data frame with columns for group name, CCD, and p-value.
#'
#' @examples
#' \dontrun{
#' library('deltaccd')
#' library('doParallel')
#' library('doRNG')
#'
#' registerDoParallel(cores=2)
#' set.seed(35811)
#'
#' refCor = getRefCor()
#' ccdResult = calcCCD(refCor, GSE19188$emat, GSE19188$groupVec, dopar=TRUE)
#' deltaCcdResult = calcDeltaCCD(refCor, GSE19188$emat, GSE19188$groupVec, 'non-tumor', dopar=TRUE)
#' }
#'
#' @seealso \code{\link{getRefCor}}, \code{\link{calcDeltaCCD}}
#'
#' @export
calcCCD = function(refCor, emat, groupVec=NULL, refEmat=NULL, nPerm=1000, geneNames=NULL, dopar=FALSE) {
	method = 'spearman'
	doOp = ifelse(dopar, `%dorng%`, `%do%`)

	if (missing(refCor)) {
		if (is.null(refEmat)) {
			stop('Either refCor or refEmat must be supplied.')}
		refCor = cor(t(refEmat[geneNames,]), method=method)
	} else if (any(rownames(refCor)!=colnames(refCor)) || !isSymmetric(refCor)) {
		stop('refCor must be a correlation matrix, with identical rownames and colnames.')}

	geneNames = rownames(refCor)[rownames(refCor) %in% rownames(emat)]
	refNow = refCor[geneNames, geneNames]
	if (length(geneNames)<2) {
		stop('Fewer than two genes in the reference are in the expression matrix.')
	} else if (length(geneNames) < nrow(refCor)) {
		warning(sprintf('%d gene(s) in reference is/are not in the expression matrix.',
							 nrow(refCor)-length(geneNames)))}

	if (is.null(groupVec)) {
		groupVec = rep('all', ncol(emat))
	} else if (length(groupVec)!=ncol(emat)) {
		stop('Length of groupVec does not match the number of columns in emat.')
	} else if (min(table(groupVec)) < 3) {
		stop('Each unique group in groupVec must have at least three samples.')}

	nComb = choose(nrow(emat), length(geneNames))

	if (nPerm>1) {
		result = foreach(groupNow=sort(unique(groupVec)), .combine=rbind) %do% {
			ccdObs = calcCCDSimple(refNow, emat[geneNames, groupVec==groupNow], method=method)

			ccdRand = doOp(foreach(ii=1:nPerm, .combine=c), {
				genesNow = rownames(emat)
				genesNowRand = genesNow[sample.int(length(genesNow))]
				idxRand = genesNowRand %in% geneNames
				calcCCDSimple(refNow, emat[idxRand, groupVec==groupNow], method=method)})

			pvalue = statmod::permp(sum(ccdRand <= ccdObs), nperm=nPerm, total.nperm=nComb,
											twosided=FALSE, method='approximate')
			data.frame(group = groupNow, CCD = ccdObs, Pvalue = pvalue, stringsAsFactors=FALSE)}
	} else {
		result = foreach(groupNow=sort(unique(groupVec)), .combine=rbind) %do% {
			ccdObs = calcCCDSimple(refNow, emat[geneNames, groupVec==groupNow], method=method)
			data.frame(group = groupNow, CCD = ccdObs, Pvalue = NA, stringsAsFactors=FALSE)}}

	return(result)}


calcDeltaCCDSimple = function(ref, emat, idx, method='spearman') {
	corVecRef = ref[upper.tri(ref)]
	corMat0 = cor(t(emat[,!idx]), method=method)
	corVec0 = corMat0[upper.tri(corMat0)]
	corMat1 = cor(t(emat[,idx]), method=method)
	corVec1 = corMat1[upper.tri(corMat1)]
	d = calcDist(corVecRef, corVec1) - calcDist(corVecRef, corVec0)
	return(d)}


makePerms = function(idx, nPerm=1000, dopar=FALSE) {
	doOp = ifelse(dopar, `%dorng%`, `%do%`)
	doOp(foreach(ii=1:nPerm, .combine=rbind), {
		sample(idx, length(idx))})}


#' Calculate delta clock correlation distance.
#'
#' \code{calcDeltaCCD} calculates the difference between the clock correlation distances (CCDs),
#' relative to a reference, for two groups of samples. Statistical significance is calculated
#' using permutation of the samples that belong to either of those two groups.
#'
#' @param refCor Correlation matrix to be used as the reference, such as comes
#' from \code{\link{getRefCor}}. Should contain Spearman correlation values.
#' @param emat Matrix of expression values, where each row
#' corresponds to a gene and each column corresponds to a sample. The rownames and colnames
#' of \code{refCor} should be present in the rownames of \code{emat}. For the p-value
#' calculation, it is important that \code{emat} include all measured genes, not just those
#' in \code{refCor}.
#' @param groupVec Vector indicating the group to which group each sample belongs. It is ok
#' for groupVec to have more than two groups.
#' @param groupNormal Value indicating the group in groupVec that corresponds to normal or
#' healthy. Other groups will be compared to this group.
#' @param refEmat Optional expression matrix for calculating co-expression for the reference,
#' with the same organization as \code{emat}. Only used if \code{refCor} is not provided.
#' @param nPerm Number of permutations for assessing statistical significance.
#' @param geneNames Optional vector indicating a subset of genes in \code{refCor}, \code{emat},
#' and/or \code{refEmat} to use for calculating the CCD.
#' @param dopar Logical indicating whether to process features in parallel. Prior to calling
#' \code{calcCCD}, use \code{\link[doParallel]{registerDoParallel}} to register the parallel
#' backend, followed by \code{\link{set.seed}} to make the p-values reproducible.
#'
#' @return A data frame with columns for group 1, group 2, deltaCCD, and p-value. In each row,
#' the deltaCCD is the CCD of group 2 minus the CCD of group 1, so group 1 corresponds
#' to \code{groupNormal}.
#'
#' @examples
#' \dontrun{
#' library('deltaccd')
#' library('doParallel')
#' library('doRNG')
#'
#' registerDoParallel(cores=2)
#' set.seed(35811)
#'
#' refCor = getRefCor()
#' ccdResult = calcCCD(refCor, GSE19188$emat, GSE19188$groupVec, dopar=TRUE)
#' deltaCcdResult = calcDeltaCCD(refCor, GSE19188$emat, GSE19188$groupVec, 'non-tumor', dopar=TRUE)
#' }
#'
#' @seealso \code{\link{getRefCor}}, \code{\link{calcCCD}}
#'
#' @export
calcDeltaCCD = function(refCor, emat, groupVec, groupNormal, refEmat=NULL,
								nPerm=1000, geneNames=NULL, dopar=FALSE) {
	method = 'spearman'
	doOp = ifelse(dopar, `%dorng%`, `%do%`)

	if (missing(refCor)) {
		if (is.null(refEmat)) {
			stop('Either refCor or refEmat must be supplied.')}
		refCor = cor(t(refEmat[geneNames,]), method=method)
	} else if (any(rownames(refCor)!=colnames(refCor)) || !isSymmetric(refCor)) {
		stop('refCor must be a correlation matrix, with identical rownames and colnames.')}

	geneNames = rownames(refCor)[rownames(refCor) %in% rownames(emat)]
	refNow = refCor[geneNames, geneNames]
	if (length(geneNames)<2) {
		stop('Fewer than two genes in the reference are in the expression matrix.')
	} else if (length(geneNames) < nrow(refCor)) {
		warning(sprintf('%d gene(s) in reference is/are not in the expression matrix.',
							 nrow(refCor)-length(geneNames)))}

	if (length(groupVec)!=ncol(emat)) {
		stop('Length of groupVec does not match the number of columns in emat.')
	} else if (!(groupNormal %in% groupVec)) {
		stop('The supplied value for groupNormal is not present in groupVec.')
	} else {
		tt = table(groupVec)
		if (length(tt)<2) {
			stop('groupVec contains only one unique group.')
		} else if (min(tt)<3) {
			stop('Each unique group in groupVec must have at least three samples.')}}

	result = data.frame(group1 = groupNormal, group2 = setdiff(sort(unique(groupVec)), groupNormal),
							  stringsAsFactors=FALSE)

	if (nPerm>1) {
		resultTmp = foreach(group2Now=result$group2, .combine=rbind) %do% {
			idx1 = groupVec %in% c(groupNormal, group2Now)
			idx2 = groupVec[idx1]==group2Now
			ematNow = emat[geneNames, idx1]
			deltaCcdObs = calcDeltaCCDSimple(refNow, ematNow, idx2, method=method)

			idxPerm = makePerms(idx2, nPerm=nPerm, dopar=dopar)
			deltaCcdRand = doOp(foreach(ii=1:nrow(idxPerm), .combine=c), {
				calcDeltaCCDSimple(refNow, ematNow, idxPerm[ii,], method=method)})

			nComb = choose(length(idx2), sum(idx2))
			pvalue = statmod::permp(sum(deltaCcdRand >= deltaCcdObs), nperm=nPerm, total.nperm=nComb,
											twosided=FALSE, method='approximate')
			data.frame(DeltaCCD = deltaCcdObs, Pvalue = pvalue, stringsAsFactors=FALSE)}
	} else {
		resultTmp = foreach(group2Now=result$group2, .combine=rbind) %do% {
			idx1 = groupVec %in% c(groupNormal, group2Now)
			idx2 = groupVec[idx1]==group2Now
			ematNow = emat[geneNames, idx1]
			deltaCcdObs = calcDeltaCCDSimple(refNow, ematNow, idx2, method=method)
			data.frame(DeltaCCD = deltaCcdObs, Pvalue = NA, stringsAsFactors=FALSE)}}

	result = cbind(result, resultTmp)
	return(result)}


#' Gene expression data for GSE19188.
#'
#' Data of gene expression measured by microarray for tumor and non-tumor samples
#' from human non-small cell lung cancer. The data is used in examples for the \code{deltaccd} package.
#'
#' @format A list with two objects:
#' \describe{
#'    \item{emat}{matrix of normalized expression values, where each row corresponds to a gene
#'    (rownames correspond to Entrez gene IDs) and each column corresponds to a sample}
#'    \item{groupVec}{vector of condition (tumor or non-tumor) for each sample}
#' }
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19188}
#'
#' @seealso \code{\link{getRefCor}}, \code{\link{calcCCD}}, \code{\link{calcDeltaCCD}}
"GSE19188"
