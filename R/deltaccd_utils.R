calcDist = function(r1, r2) sqrt(sum((r1 - r2)^2, na.rm = TRUE))


calcCCDSimple = function(ref, emat, method = 'spearman', scale = FALSE) {
  corVecRef = ref[upper.tri(ref)]
  corMatTest = stats::cor(t(emat), method = method)
  corVecTest = corMatTest[upper.tri(corMatTest)]
  ccd = calcDist(corVecRef, corVecTest)
  
  if (isTRUE(scale)) {
    nPairs = choose(ncol(ref), 2)
    ccd = ccd / nPairs}
  
  return(ccd)}


checkVar = function(emat, groupVec) {
  
  varCheck = foreach (groupNow = sort(unique(groupVec)), .combine = rbind) %do% {
    
    varVec = apply(emat[, groupVec == groupNow], MARGIN = 1, 
                   FUN = stats::var, na.rm = TRUE)
    varDt = data.table::as.data.table(varVec, keep.rownames = 'gene')
    data.table::setnames(varDt, 'varVec', 'variance')
    varDt[, group := groupNow]
    
    zeroVar = varDt[variance == 0]}
  varCheck[, variance := NULL]  
  
  if (nrow(varCheck) > 0) {
    stop('Zero variance in the following gene-group pairs:\n', 
         paste(utils::capture.output(print(varCheck)), collapse = '\n'))}
  
  invisible()}


checkGenes = function(emat, refCor) {
  geneNames = rownames(refCor)[rownames(refCor) %in% rownames(emat)]
  if (length(geneNames) < nrow(refCor)) {
    missingGenes = setdiff(rownames(refCor), geneNames)
    stop(paste0('The following gene(s) is/are not in the expression matrix:\n',
                paste0(missingGenes, collapse = '\n')))} 
  
  return(geneNames)}


checkRefCor = function(refCor, refEmat = NULL, geneNames = NULL, method = 'spearman') {
  if (missing(refCor)) {
    if (is.null(refEmat)) {
      stop('Either refCor or refEmat must be supplied.')}
    refCor = stats::cor(t(refEmat[geneNames,]), method = method)
  } else if (any(rownames(refCor) != colnames(refCor)) || !isSymmetric(refCor)) {
    stop('refCor must be a correlation matrix, with identical rownames and colnames.')}
  
  return(refCor)}

calcDeltaCCDSimple = function(ref, emat, idx, method = 'spearman', scale = FALSE) {
  d0 = calcCCDSimple(ref, emat[, -idx], method = method, scale = scale)
  d1 = calcCCDSimple(ref, emat[, idx], method = method, scale = scale)
  d = d1 - d0
  return(d)}


makePerms = function(idx, nPerm = 1000, dopar = FALSE) {
  doOp = if (isTRUE(dopar)) `%dorng%` else `%do%`
  r = doOp(foreach(i = 1:nPerm, .combine = rbind), {
    sample(idx, length(idx))})
  return(r)}