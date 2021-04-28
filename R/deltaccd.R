#' @importFrom foreach foreach %do% %dopar%
#' @importFrom doRNG %dorng%
#' @importFrom data.table data.table
NULL


globalVariables(c('ii', 'groupNow', 'group', '.', 'gene1', 'gene2',
                  'rho', 'group2Now', 'geneNames', 'ematNow'))


#' Retrieve the reference correlation matrix for clock gene co-expression.
#'
#' The reference matrix is based on a fixed-effects meta-analysis of eight
#' circadian transcriptome datasets from mice, as described in
#' [Shilts et al. 2018](https://doi.org/10.7717/peerj.4327).
#'
#' @param species Currently either 'human' or 'mouse'. Only affects the row and
#'   column names of the correlation matrix, not the actual values.
#' @param useEntrezGeneId If `FALSE`, row and column names of correlation matrix
#'   will correspond to gene symbols (e.g., PER2).
#'
#' @return A matrix of Spearman correlation values.
#'
#' @examples
#' \dontrun{
#' library('deltaccd')
#' library('doParallel')
#' library('doRNG')
#'
#' registerDoParallel(cores = 2)
#' set.seed(35813)
#'
#' refCor = getRefCor()
#' ccdResult = calcCCD(refCor, GSE19188$emat, GSE19188$groupVec, dopar = TRUE)
#' deltaCcdResult = calcDeltaCCD(refCor, GSE19188$emat, GSE19188$groupVec,
#'                               'non-tumor', dopar = TRUE)
#'
#' pRef = plotRefHeatmap(refCor)
#' pTest = plotHeatmap(rownames(refCor), GSE19188$emat, GSE19188$groupVec)
#' }
#'
#' @seealso [GSE19188], [plotRefHeatmap()], [calcCCD()], [calcDeltaCCD()]
#'
#' @export
getRefCor = function(species = 'human', useEntrezGeneId = TRUE) {
  if (species == 'human') {
    if (useEntrezGeneId) {
      colname = 'entrez_hs'
    } else {
      colname = 'symbol_hs'}
  } else if (species == 'mouse') {
    if (useEntrezGeneId) {
      colname = 'entrez_mm'
    } else {
      colname = 'symbol_mm'}}
  ref = refCorMouseEntrez
  rownames(ref) = clockGenes[[colname]]
  colnames(ref) = rownames(ref)
  return(ref)}


calcDist = function(r1, r2) sqrt(sum((r1 - r2)^2, na.rm = TRUE))


calcCCDSimple = function(ref, emat, method = 'spearman', scale = FALSE) {
  
  nPairs = choose(ncol(ref), 2)
  
  corVecRef = ref[upper.tri(ref)]
  corMatTest = stats::cor(t(emat), method = method)
  corVecTest = corMatTest[upper.tri(corMatTest)]
  
  ccd = calcDist(corVecRef, corVecTest)
  
  if(isTRUE(scale)) {
    
    ccd = ccd/nPairs}
  
  return(ccd)}


#' Calculate clock correlation distance (CCD).
#'
#' Quantify the similarity of gene co-expression between a reference and a test
#' dataset. Statistical significance is calculated using permutation of the
#' genes.
#'
#' @param refCor Correlation matrix to be used as the reference, such as comes
#'   from [getRefCor()]. Should contain Spearman correlation values.
#' @param emat Matrix of expression values, where each row corresponds to a
#'   gene and each column corresponds to a sample. The rownames and colnames of
#'   `refCor` should be present in the rownames of `emat`. For the p-value
#'   calculation, it is important that `emat` include all measured genes, not
#'   just those in `refCor`.
#' @param groupVec Optional vector indicating the group to which group each
#'   sample belongs. If not provided, the function assumes all samples belong to
#'   the same group.
#' @param refEmat Optional expression matrix for calculating co-expression for
#'   the reference, with the same organization as `emat`. Only used if `refCor`
#'   is not provided.
#' @param nPerm Number of permutations for assessing statistical significance.
#' @param geneNames Optional vector indicating a subset of genes in `refCor`,
#'   `emat`, and/or `refEmat` to use for calculating the CCD.
#' @param dopar Logical indicating whether to process features in parallel. Make
#'   sure to register a parallel backend first.
#'
#' @return A data.table with columns for group name, CCD, and p-value.
#'
#' @examples
#' \dontrun{
#' library('deltaccd')
#' library('doParallel')
#' library('doRNG')
#'
#' registerDoParallel(cores = 2)
#' set.seed(35813)
#'
#' refCor = getRefCor()
#' ccdResult = calcCCD(refCor, GSE19188$emat, GSE19188$groupVec, dopar = TRUE)
#' deltaCcdResult = calcDeltaCCD(refCor, GSE19188$emat, GSE19188$groupVec,
#'                               'non-tumor', dopar = TRUE)
#'
#' pRef = plotRefHeatmap(refCor)
#' pTest = plotHeatmap(rownames(refCor), GSE19188$emat, GSE19188$groupVec)
#' }
#'
#' @seealso [getRefCor()], [calcDeltaCCD()], [plotHeatmap()]
#'
#' @export
calcCCD = function(refCor, emat, groupVec = NULL, refEmat = NULL, nPerm = 1000,
                   geneNames = NULL, dopar = FALSE, scale = FALSE) {
  method = 'spearman'
  doOp = ifelse(dopar, `%dorng%`, `%do%`)

  if (missing(refCor)) {
    if (is.null(refEmat)) {
      stop('Either refCor or refEmat must be supplied.')}
    refCor = stats::cor(t(refEmat[geneNames,]), method = method)
  } else if (any(rownames(refCor) != colnames(refCor)) || !isSymmetric(refCor)) {
    stop('refCor must be a correlation matrix, with identical rownames and colnames.')}

  geneNames = rownames(refCor)[rownames(refCor) %in% rownames(emat)]
  refNow = refCor[geneNames, geneNames]
  if (length(geneNames) < 2) {
    stop('Fewer than two genes in the reference are in the expression matrix.')
  } else if (length(geneNames) < nrow(refCor)) {
    warning(sprintf('%d gene(s) in reference is/are not in the expression matrix.',
                    nrow(refCor) - length(geneNames)))}

  if (is.null(groupVec)) {
    groupVec = rep('all', ncol(emat))
  } else if (length(groupVec) != ncol(emat)) {
    stop('Length of groupVec does not match the number of columns in emat.')
  } else if (min(table(groupVec)) < 3) {
    stop('Each unique group in groupVec must have at least three samples.')}

  nComb = choose(nrow(emat), length(geneNames))

  if (nPerm>1) {
    
    result = foreach(groupNow = sort(unique(groupVec)), .combine = rbind) %do% {
      
      ccdObs = calcCCDSimple(refNow, emat[geneNames, groupVec == groupNow],
        method = method, scale = scale)

      ccdRand = doOp(foreach(ii = 1:nPerm, .combine = c), {
        genesNow = rownames(emat)
        genesNowRand = genesNow[sample.int(length(genesNow))]
        idxRand = genesNowRand %in% geneNames
        calcCCDSimple(refNow, emat[idxRand, groupVec == groupNow],
                      method = method, scale = scale)})

      pvalue = statmod::permp(sum(ccdRand <= ccdObs), nperm = nPerm,
        total.nperm = nComb, twosided = FALSE, method = 'approximate')
      
      data.table(group = groupNow, CCD = ccdObs, Pvalue = pvalue)}
    
  } else {
    
    result = foreach(groupNow = sort(unique(groupVec)), .combine = rbind) %do% {
      
      ccdObs = calcCCDSimple(refNow, emat[geneNames, groupVec == groupNow],
        method = method, scale = scale)
      
      data.table(group = groupNow, CCD = ccdObs, Pvalue = NA)}}

  return(result)}


calcDeltaCCDSimple = function(ref, emat, idx, method = 'spearman'
  , scale = FALSE) {
  
  d0 = calcCCDSimple(ref, emat[!idx, ], method = method, scale = scale)
  d1 = calcCCDSimple(ref, emat[idx, ], method = method, scale = scale)
  
  d = d1 - d0
  
  return(d)}


makePerms = function(idx, nPerm = 1000, dopar = FALSE) {
  doOp = ifelse(dopar, `%dorng%`, `%do%`)
  doOp(foreach(ii = 1:nPerm, .combine = rbind), {
    sample(idx, length(idx))})}


#' Calculate delta clock correlation distance.
#'
#' Calculate the difference between the clock correlation distances (CCDs),
#' relative to a reference, for two groups of samples. Statistical significance
#' is calculated using permutation of the samples that belong to either of those
#' two groups.
#'
#' @param refCor Correlation matrix to be used as the reference, such as comes
#'   from [getRefCor()]. Should contain Spearman correlation values.
#' @param emat Matrix of expression values, where each row corresponds to a gene
#'   and each column corresponds to a sample. The rownames and colnames of
#'   `refCor` should be present in the rownames of `emat`. For the p-value
#'   calculation, it is important that `emat` include all measured genes, not
#'   just those in `refCor`.
#' @param groupVec Vector indicating the group to which group each sample
#'   belongs. It's ok for groupVec to have more than two groups.
#' @param groupNormal Value indicating the group in groupVec that corresponds to
#'   normal or healthy. Other groups will be compared to this group.
#' @param refEmat Optional expression matrix for calculating co-expression for
#'   the reference, with the same organization as `emat`. Only used if `refCor`
#'   is not provided.
#' @param nPerm Number of permutations for assessing statistical significance.
#' @param geneNames Optional vector indicating a subset of genes in `refCor`,
#'   `emat`, and/or `refEmat` to use for calculating the CCD.
#' @param dopar Logical indicating whether to process features in parallel. Make
#'   sure to register a parallel backend first.
#'
#' @return A data.table with columns for group 1, group 2, deltaCCD, and
#'   p-value. In each row, the deltaCCD is the CCD of group 2 minus the CCD of
#'   group 1, so group 1 corresponds to `groupNormal`.
#'
#' @examples
#' \dontrun{
#' library('deltaccd')
#' library('doParallel')
#' library('doRNG')
#'
#' registerDoParallel(cores = 2)
#' set.seed(35813)
#'
#' refCor = getRefCor()
#' ccdResult = calcCCD(refCor, GSE19188$emat, GSE19188$groupVec, dopar = TRUE)
#' deltaCcdResult = calcDeltaCCD(refCor, GSE19188$emat, GSE19188$groupVec,
#'                               'non-tumor', dopar = TRUE)
#'
#' pRef = plotRefHeatmap(refCor)
#' pTest = plotHeatmap(rownames(refCor), GSE19188$emat, GSE19188$groupVec)
#' }
#'
#' @seealso [getRefCor()], [calcCCD()], [plotHeatmap()]
#'
#' @export
calcDeltaCCD = function(refCor, emat, groupVec, groupNormal, refEmat = NULL,
  nPerm = 1000, geneNames = NULL, dopar = FALSE, scale = FALSE) {
  
  method = 'spearman'
  doOp = ifelse(dopar, `%dorng%`, `%do%`)

  if (missing(refCor)) {
    if (is.null(refEmat)) {
      stop('Either refCor or refEmat must be supplied.')}
    refCor = stats::cor(t(refEmat[geneNames,]), method = method)
  } else if (any(rownames(refCor) != colnames(refCor)) || !isSymmetric(refCor)) {
    stop('refCor must be a correlation matrix, with identical rownames and colnames.')}

  geneNames = rownames(refCor)[rownames(refCor) %in% rownames(emat)]
  refNow = refCor[geneNames, geneNames]
  if (length(geneNames) < 2) {
    stop('Fewer than two genes in the reference are in the expression matrix.')
  } else if (length(geneNames) < nrow(refCor)) {
    warning(sprintf('%d gene(s) in reference is/are not in the expression matrix.',
                    nrow(refCor) - length(geneNames)))}

  if (length(groupVec) != ncol(emat)) {
    stop('Length of groupVec does not match the number of columns in emat.')
  } else if (!(groupNormal %in% groupVec)) {
    stop('The supplied value for groupNormal is not present in groupVec.')
  } else {
    tt = table(groupVec)
    if (length(tt) < 2) {
      stop('groupVec contains only one unique group.')
    } else if (min(tt) < 3) {
      stop('Each unique group in groupVec must have at least three samples.')}}

  result = data.table(group1 = groupNormal,
                      group2 = setdiff(sort(unique(groupVec)), groupNormal))

  if (nPerm>1) {
    resultTmp = foreach(group2Now = result$group2, .combine = rbind) %do% {
      idx1 = groupVec %in% c(groupNormal, group2Now)
      idx2 = groupVec[idx1] == group2Now
      ematNow = emat[geneNames, idx1]
      deltaCcdObs = calcDeltaCCDSimple(refNow, ematNow, idx2,
                                       method = method, scale = scale)

      idxPerm = makePerms(idx2, nPerm = nPerm, dopar = dopar)
      deltaCcdRand = doOp(foreach(ii = 1:nrow(idxPerm), .combine = c), {
        calcDeltaCCDSimple(refNow, ematNow, idxPerm[ii,], method = method
                           , scale = scale)})

      nComb = choose(length(idx2), sum(idx2))
      pvalue = statmod::permp(sum(deltaCcdRand >= deltaCcdObs),
                              nperm = nPerm, total.nperm = nComb,
                              twosided = FALSE, method = 'approximate')
      data.table(DeltaCCD = deltaCcdObs, Pvalue = pvalue)}
    
  } else {
    
    resultTmp = foreach(group2Now = result$group2, .combine = rbind) %do% {
      idx1 = groupVec %in% c(groupNormal, group2Now)
      idx2 = groupVec[idx1] == group2Now
      ematNow = emat[geneNames, idx1]
      deltaCcdObs = calcDeltaCCDSimple(refNow, ematNow, idx2,
                                       method = method, scale = scale)
      data.table(DeltaCCD = deltaCcdObs, Pvalue = NA)}}

  result = cbind(result, resultTmp)
  return(result)}


#' Gene expression data for GSE19188.
#'
#' Data of gene expression measured by microarray for tumor and non-tumor
#' samples from human non-small cell lung cancer.
#'
#' @format A list with two objects:
#' \describe{
#'    \item{emat}{matrix of normalized expression values, where each row
#'      corresponds to a gene (rownames correspond to Entrez gene IDs) and each
#'      column corresponds to a sample}
#'    \item{groupVec}{vector of condition (tumor or non-tumor) for each sample}
#' }
#'
#' @source <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19188>
#'
#' @seealso [getRefCor()], [calcCCD()], [calcDeltaCCD()]
'GSE19188'
