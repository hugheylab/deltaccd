calcCorrSimple = function(dt, method = 'spearman') {
  
  geneNames = setdiff(colnames(dt), 'group')
  
  dt1 = data.table(stats::cor(as.matrix(dt[, geneNames]), method = method),
                   gene1 = geneNames)
  
  dt1 = data.table::melt(dt1, id.vars = 'gene1', measure.vars = 'gene2', 
                         value.name = 'rho')
  
  dt1 = dt1[gene1 != gene2]
  
  return(dt1)}


calcCorr = function(ematNow, groupVec, method = 'spearman') {
  
  dt = data.table(t(ematNow), group = groupVec)
  
  dtFinal = foreach(grp = groupVec, .combine = rbind) %do% {
    
    dtTmp = calcCorrSimple(dt[group == grp], method = method)
  
    return(dtTmp)}

  dtFinal[, gene1 := factor(gene1, rownames(ematNow))]
  dtFinal[, gene2 := factor(gene2, rev(rownames(ematNow)))]
  
  return(dtFinal)}
  

calcColorLimits = function(vals, vLow = -1, vMid = 0, vHigh = 1,
                           cLow = '#e66101', cMid = '#f7f7f7',
                           cHigh = '#5e3c99') {
  valRange = seq(0, 1, length.out = 201)
  colorScale = scales::div_gradient_pal(low = cLow, mid = cMid,
                                        high = cHigh)(valRange)

  minVal = (min(vals, na.rm = TRUE) - vLow) / (vHigh - vLow)
  idxLow = which.min(abs(minVal - valRange))

  maxVal = (max(vals, na.rm = TRUE) - vLow) / (vHigh - vLow)
  idxHigh = which.min(abs(maxVal - valRange))
  return(colorScale[c(idxLow, idxHigh)])}


plotHeatmapSimple = function(ggObj, cLims) {
  p = ggObj +
    ggplot2::geom_tile(ggplot2::aes(x = gene1, y = gene2, fill = rho)) +
    ggplot2::labs(x = 'Gene', y = 'Gene') +
    ggplot2::scale_fill_gradient2(low = cLims[1], mid = '#f7f7f7',
                                  high = cLims[2], na.value = 'grey80') +
    ggplot2::theme_light() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5,
                                                       hjust = 1),
                   strip.text = ggplot2::element_text(color = 'black'),
                   axis.text = ggplot2::element_text(color = 'black'),
                   legend.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0,
                                                   unit = 'cm'))}


#' Visualize gene co-expression.
#'
#' Make heatmaps of the co-expression (Spearman correlation) between pairs of
#' selected genes in a dataset.
#'
#' @param geneNames Vector indicating the subset of genes in the rownames of
#'   `emat` for which to calculate the correlations in expression.
#' @param emat Matrix of expression values, where each row corresponds to a
#'   gene and each column corresponds to a sample. The elements of `geneNames`
#'   should be present in the rownames of `emat`.
#' @param groupVec Optional vector indicating the group to which group each
#'   sample belongs. If not provided, the function assumes all samples belong
#'   to the same group.
#'
#' @return A `ggplot` object, which can be saved using [ggplot2::ggsave()].
#'   Heatmap colors will be directly comparable to any heatmaps created by this
#'   function or by [plotRefHeatmap()].
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
#'										'non-tumor', dopar = TRUE)
#'
#' pRef = plotRefHeatmap(refCor)
#' pTest = plotHeatmap(rownames(refCor), GSE19188$emat, GSE19188$groupVec)
#' }
#'
#' @seealso [calcCCD()], [calcDeltaCCD()], [plotRefHeatmap()]
#'
#' @export
plotHeatmap = function(geneNames, emat, groupVec = NULL) {
  method = 'spearman'

  ematNow = emat[geneNames,]
  if (nrow(ematNow)<2) {
    stop('Fewer than two genes in the supplied vector are in the expression matrix.')
  } else if (nrow(ematNow) < length(geneNames)) {
    warning(sprintf('%d gene(s) in the supplied vector is/are not in the expression matrix.',
                    length(geneNames) - nrow(ematNow)))}

  if (is.null(groupVec)) {
    groupVec = rep('all', ncol(emat))
  } else if (length(groupVec) != ncol(emat)) {
    stop('Length of groupVec does not match the number of columns in emat.')
  } else if (min(table(groupVec)) < 3) {
    stop('Each unique group in groupVec must have at least three samples.')}

  dt = calcCorr(ematNow, groupVec, method)
  cLims = calcColorLimits(dt$rho)
  p = plotHeatmapSimple(ggplot2::ggplot(dt) + 
                          ggplot2::facet_wrap(ggplot2::vars(group)), cLims)}


#' Visualize the reference pattern of gene co-expression.
#'
#' Make a heatmap of the reference correlation matrix for gene co-expression.
#'
#' @param refCor Correlation matrix, such as comes from [getRefCor()].
#'
#' @return A `ggplot` object, which can be saved using [ggplot2::ggsave()].
#'   Heatmap colors will be directly comparable to any heatmaps created by this
#'   function or by [plotHeatmap()].
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
#' @seealso [getRefCor()], [plotHeatmap()]
#'
#' @export
plotRefHeatmap = function(refCor) {
  if (any(rownames(refCor) != colnames(refCor)) || !isSymmetric(refCor)) {
    stop('refCor must be a correlation matrix, with identical rownames and colnames.')}
  
  dt = data.table(refCor, gene1 = geneNames)
  
  dt = data.table::melt(dt, id.vars = 'gene1', measure.vars = 'gene2'
    , value.name = 'rho')
  dt = dt[gene1 != gene2]
  dt[, gene1 := factor(gene1, rownames(ematNow))]
  dt[, gene2 := factor(gene2, rev(rownames(ematNow)))]
 
  cLims = calcColorLimits(dt$rho)
  p = plotHeatmapSimple(ggplot2::ggplot(dt), cLims)}
