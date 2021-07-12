test_that('calcDist', {
  d = calcDist(c(1, 0), c(1, 0))
  expect_equal(d, 0)
  
  d = calcDist(c(1, 0), c(0, 1))
  expect_equal(d, sqrt(2))
})

test_that('calcCCDSimple', {
  ccd = calcCCDSimple(diag(1, 2), diag(1, 2))
  expect_equal(ccd, 1)
  
  ref = diag(1, 2)
  emat = rbind(c(1, 0, 1), c(0, 1, 0))
  ccd = calcCCDSimple(ref, emat)
  expect_equal(ccd, 1)
  
  ref = diag(1, 3)
  emat = rbind(c(1, 0, 0), c(0, 1, -1), c(0, 1, 1))
  ccd = calcCCDSimple(ref, emat, scale = TRUE)
  expect_equal(ccd, 1/3)
})

test_that('checkVar', {
  groupVec = c('a', 'a', 'b', 'b')
  emat = cbind(diag(1, 2), diag(1, 2))
  expect_invisible(checkVar(emat, groupVec))
  
  emat = cbind(diag(1, 2), matrix(1, 2, 2))
  rownames(emat) = 1:2
  expect_error(checkVar(emat, c('a', 'a', 'b', 'b')), 
               paste0('Zero variance in the following gene-group pairs:\n', 
                     paste(utils::capture.output(
                       data.table(gene = c(1, 2), group = c('b', 'b'))), 
                       collapse = '\n')), fixed = TRUE)
})

test_that('checkGenes', {
  emat = diag(1, 2)
  rownames(emat) = 1:2
  refCor = diag(1, 2)
  rownames(refCor) = 1:2
  expect_visible(checkGenes(emat, refCor))
  
  refCor = diag(1, 4)
  rownames(refCor) = 1:4
  expect_error(checkGenes(emat, refCor), 
               paste0('The following gene(s) is/are not in the expression matrix:\n', 
                      paste0(3:4, collapse = '\n')), fixed = TRUE)
})

test_that('checkRefCor', {
  refCor = diag(1, 2)
  rownames(refCor) = 1:2
  colnames(refCor) = rownames(refCor)
  expect_visible(checkRefCor(refCor))
  
  expect_error(checkRefCor(), 'Either refCor or refEmat must be supplied.')
  
  rownames(refCor) = 3:4
  expect_error(checkRefCor(refCor),
               'refCor must be a correlation matrix, with identical rownames and colnames.')
  
  refCor = rbind(c(1, 0, 0), c(1, 0, 0))
  expect_error(checkRefCor(refCor),
               'refCor must be a correlation matrix, with identical rownames and colnames.')
})

test_that('calcDeltaCCDSimple', {
  ref = diag(1, 2)
  emat = cbind(diag(1, 2), diag(1, 2))
  expect_equal(calcDeltaCCDSimple(ref, emat, 1:2), 0)
})

test_that('makePerms', {
  perms = makePerms(1:10)
  expect_equal(nrow(perms), 1000)
  expect_equal(ncol(perms), 10)
})

# test_that('calcCorr', {
#   ematNow = rbind(c(1, -1, -1, 1), c(-1, 1, 1, -1))
#   rownames(ematNow) = 1:2
#   groupVec = c('a' ,'a', 'b', 'b')
#   
#   expect_
# })
