#' @importFrom foreach foreach %do% %dopar%
#' @importFrom doRNG %dorng%
#' @importFrom data.table data.table :=
#' @importFrom rlang .data
NULL


#' Retrieve the reference correlation matrix for circadian gene co-expression.
#'
#' The pan-tissue reference matrix is based on a fixed-effects meta-analysis of
#' eight circadian transcriptome datasets from mice, as described in
#' Shilts et al. 2018(\doi{https://doi.org/10.7717/peerj.4327}). The human blood
#' reference matrix is based an analysis of three microarray datasets
#' (manuscript in preparation).
#'
#' @param species Currently either 'human' or 'mouse'. Only affects the row and
#'   column names of the correlation matrix, not the actual values.
#' @param tissue One of either 'pan' or 'blood'.
#' @param useEntrezGeneId If `FALSE`, row and column names of correlation matrix
#'   will correspond to gene symbols (e.g., PER2).
#'
#' @return A matrix of Spearman correlation values.
#'
#' @seealso [GSE19188], [plotRefHeatmap()], [calcCCD()], [calcDeltaCCD()]
#'
#' @export
getRefCor = function(
  species = c('human', 'mouse'), tissue = c('pan', 'blood'),
  useEntrezGeneId = TRUE) {

  species =  match.arg(species)
  tissue = match.arg(tissue)

  if (species == 'mouse' && tissue == 'blood') {
    stop('Blood reference is only available for species = \'human\'.')}

  prefix = if (isTRUE(useEntrezGeneId)) 'entrez' else 'symbol'
  suffix = if (species == 'human') 'hs' else 'mm'
  colname = paste(prefix, suffix, sep = '_')

  if (tissue == 'pan') {
    ref = refCorMouseEntrez
    rownames(ref) = clockGenes[[colname]]
  } else {
    ref = refCorHumanBlood
    rownames(ref) = bloodGenes[[colname]]}

  colnames(ref) = rownames(ref)
  return(ref)}


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
#' @param scale Logical indicating whether to scale CCD by the number of gene
#'   pairs.
#'
#' @return A data.table with columns for group name, CCD, and p-value.
#'
#' @examples
#' set.seed(35813)
#'
#' refCor = getRefCor()
#' ccdResult = calcCCD(refCor, GSE19188$emat, GSE19188$groupVec, nPerm = 100)
#'
#' @seealso [getRefCor()], [calcDeltaCCD()], [plotHeatmap()]
#'
#' @export
calcCCD = function(
  refCor, emat, groupVec = NULL, refEmat = NULL, nPerm = 1000, geneNames = NULL,
  dopar = FALSE, scale = FALSE) {

  groupNow = NULL
  method = 'spearman'
  doOp = if (isTRUE(dopar)) `%dorng%` else `%do%`

  refCor = checkRefCor(refCor, refEmat, geneNames, method)
  geneNames = checkGenes(emat, refCor, geneNames)

  if (is.null(groupVec)) {
    groupVec = rep('all', ncol(emat))
  } else if (length(groupVec) != ncol(emat)) {
    stop('Length of groupVec does not match the number of columns in emat.')
  } else if (min(table(groupVec)) < 3) {
    stop('Each unique group in groupVec must have at least three samples.')}

  if (nPerm > 1) {
    checkVar(emat, groupVec)
  } else {
    checkVar(emat[geneNames, ], groupVec)}

  nComb = choose(nrow(emat), length(geneNames))

  if (nPerm > 1) {
    result = foreach(groupNow = sort(unique(groupVec)), .combine = rbind) %do% {

      ccdObs = calcCCDSimple(refCor, emat[geneNames, groupVec == groupNow],
                             method = method, scale = scale)

      ccdRand = doOp(foreach(i = 1:nPerm, .combine = c), {
        genesNow = rownames(emat)
        genesNowRand = genesNow[sample.int(length(genesNow))]
        idxRand = genesNowRand %in% geneNames
        calcCCDSimple(refCor, emat[idxRand, groupVec == groupNow],
                      method = method, scale = scale)})

      pvalue = statmod::permp(sum(ccdRand <= ccdObs), nperm = nPerm,
        total.nperm = nComb, twosided = FALSE, method = 'approximate')

      data.table(group = groupNow, CCD = ccdObs, Pvalue = pvalue)}

  } else {
    result = foreach(groupNow = sort(unique(groupVec)), .combine = rbind) %do% {

      ccdObs = calcCCDSimple(refCor, emat[geneNames, groupVec == groupNow],
                             method = method, scale = scale)

      data.table(group = groupNow, CCD = ccdObs, Pvalue = NA)}}

  return(result)}


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
#' @param scale Logical indicating whether to use scaled CCDs to calculate
#'   difference.
#'
#' @return A data.table with columns for group 1, group 2, deltaCCD, and
#'   p-value. In each row, the deltaCCD is the CCD of group 2 minus the CCD of
#'   group 1, so group 1 corresponds to `groupNormal`.
#'
#' @examples
#' set.seed(35813)
#'
#' refCor = getRefCor()
#' deltaCcdResult = calcDeltaCCD(
#'   refCor, GSE19188$emat, GSE19188$groupVec, 'healthy', nPerm = 100)
#'
#' @seealso [getRefCor()], [calcCCD()], [plotHeatmap()]
#'
#' @export
calcDeltaCCD = function(
  refCor, emat, groupVec, groupNormal, refEmat = NULL, nPerm = 1000,
  geneNames = NULL, dopar = FALSE, scale = FALSE) {
  group2Now = i = NULL

  method = 'spearman'
  doOp = if (isTRUE(dopar)) `%dorng%` else `%do%`

  refCor = checkRefCor(refCor, refEmat, geneNames, method)
  geneNames = checkGenes(emat, refCor)

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

  checkVar(emat[geneNames, ], groupVec)

  result = data.table(
    group1 = groupNormal,
    group2 = setdiff(sort(unique(groupVec)), groupNormal))

  if (nPerm > 1) {
    resultTmp = foreach(group2Now = result$group2, .combine = rbind) %do% {
      idx1 = groupVec %in% c(groupNormal, group2Now)
      idx2 = groupVec[idx1] == group2Now
      ematNow = emat[geneNames, idx1]
      deltaCcdObs = calcDeltaCCDSimple(
        refCor, ematNow, idx2, method = method, scale = scale)

      idxPerm = makePerms(idx2, nPerm = nPerm, dopar = dopar)
      deltaCcdRand = doOp(foreach(i = 1:nrow(idxPerm), .combine = c), {
        calcDeltaCCDSimple(
          refCor, ematNow, idxPerm[i,], method = method, scale = scale)})

      nComb = choose(length(idx2), sum(idx2))
      pvalue = statmod::permp(
        sum(deltaCcdRand >= deltaCcdObs), nperm = nPerm, total.nperm = nComb,
        twosided = FALSE, method = 'approximate')
      data.table(DeltaCCD = deltaCcdObs, Pvalue = pvalue)}

  } else {
    resultTmp = foreach(group2Now = result$group2, .combine = rbind) %do% {
      idx1 = groupVec %in% c(groupNormal, group2Now)
      idx2 = groupVec[idx1] == group2Now
      ematNow = emat[geneNames, idx1]
      deltaCcdObs = calcDeltaCCDSimple(
        refCor, ematNow, idx2, method = method, scale = scale)
      data.table(DeltaCCD = deltaCcdObs, Pvalue = NA)}}

  result = cbind(result, resultTmp)
  return(result)}


#' Gene expression data for GSE19188.
#'
#' Data of gene expression measured by microarray for samples from human
#' non-small cell lung cancer.
#'
#' @format A list with two objects:
#' \describe{
#'    \item{emat}{Matrix of normalized expression values, where each row
#'      corresponds to a gene (rownames are Entrez Gene IDs) and each column
#'      corresponds to a sample. To save space, genes have been downsampled.}
#'    \item{groupVec}{Vector of condition (tumor or healthy) for each sample.}
#' }
#'
#' @source <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19188>
#'
#' @seealso [getRefCor()], [calcCCD()], [calcDeltaCCD()]
'GSE19188'
